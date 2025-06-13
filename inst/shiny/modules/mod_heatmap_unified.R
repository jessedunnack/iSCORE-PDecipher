# Module: Unified Heatmap Visualization
# Integrates the unified heatmap system into Shiny

# UI function
mod_heatmap_unified_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        # Quick preset buttons
        div(class = "well well-sm", style = "text-align: center;",
          h4("Quick Comparisons", style = "margin-top: 0;"),
          actionButton(ns("preset_methods"), "Compare Methods", 
                       icon = icon("exchange-alt"), class = "btn-primary btn-sm"),
          actionButton(ns("preset_metrics"), "Compare Metrics", 
                       icon = icon("chart-bar"), class = "btn-primary btn-sm"),
          actionButton(ns("preset_directions"), "Compare Directions", 
                       icon = icon("arrows-alt-v"), class = "btn-primary btn-sm"),
          actionButton(ns("preset_comprehensive"), "Comprehensive View", 
                       icon = icon("th"), class = "btn-success btn-sm")
        )
      )
    ),
    
    fluidRow(
      # Configuration panel
      column(3,
        div(class = "box box-primary",
          div(class = "box-header with-border",
              h3("Heatmap Configuration", class = "box-title")
          ),
          div(class = "box-body",
          
          h5("Basic Settings", style = "font-weight: bold;"),
          
          selectInput(ns("heatmap_metric"),
                      "Primary Metric:",
                      choices = c("P-value (-log10)" = "pvalue",
                                "Fold Enrichment" = "foldenrichment",
                                "Z-score" = "zscore",
                                "GSEA NES" = "gsea"),
                      selected = "pvalue"),
          
          selectInput(ns("method_comparison"),
                      "Method Filter:",
                      choices = c("All Methods" = "all",
                                "MAST Only" = "MAST",
                                "MixScale Only" = "MixScale",
                                "Intersection (Both)" = "intersection",
                                "Union (Either)" = "union"),
                      selected = "all"),
          
          conditionalPanel(
            condition = sprintf("input['%s'] == 'gsea'", ns("heatmap_metric")),
            numericInput(ns("min_nes"),
                        "Min |NES| Magnitude:",
                        value = 1.0,
                        min = 0,
                        max = 3,
                        step = 0.1)
          ),
          
          hr(),
          h5("Display Options", style = "font-weight: bold;"),
          
          checkboxInput(ns("show_directions"),
                       "Show UP/DOWN/ALL Annotations",
                       value = TRUE),
          
          numericInput(ns("max_terms"),
                      "Maximum Terms to Display:",
                      value = 25,
                      min = 10,
                      max = 100,
                      step = 5),
          
          selectInput(ns("color_scheme"),
                      "Color Scheme:",
                      choices = c("Default" = "default_complexheatmap",
                                "Fixed Breaks" = "fixed_user",
                                "Adaptive" = "adaptive_quantile"),
                      selected = "default_complexheatmap"),
          
          hr(),
          h5("Advanced Options", style = "font-weight: bold;"),
          
          checkboxInput(ns("cluster_rows"),
                       "Cluster Rows (Terms)",
                       value = TRUE),
          
          checkboxInput(ns("cluster_columns"),
                       "Cluster Columns (Genes)",
                       value = TRUE),
          
          checkboxInput(ns("show_row_names"),
                       "Show Term Names",
                       value = TRUE),
          
          checkboxInput(ns("show_column_names"),
                       "Show Gene Names",
                       value = TRUE),
          
          hr(),
          
          # Custom title input
          textInput(ns("custom_title"),
                   "Custom Title (optional):",
                   value = "",
                   placeholder = "Leave empty for auto-generated title"),
          
          # Refresh button (in case reactive updates are paused)
          actionButton(ns("refresh_heatmap"),
                      "Refresh Heatmap",
                      icon = icon("sync"),
                      class = "btn-primary btn-block")
          )
        )
      ),
      
      # Main heatmap display
      column(9,
        div(class = "box box-primary",
          div(class = "box-header with-border",
              h3("Heatmap Visualization", class = "box-title")
          ),
          div(class = "box-body",
          
          # Heatmap output with dynamic height
          withSpinner(
            plotOutput(ns("heatmap_display"), height = "auto"),
            type = 6,
            color = "#3c8dbc"
          ),
          
          br(),
          
          # Method breakdown log (collapsible)
          tags$details(
            tags$summary(icon("info-circle"), " Method Details & Filtering Log"),
            verbatimTextOutput(ns("method_log"))
          ),
          
          br(),
          
          # Download options
          div(style = "text-align: center;",
            downloadButton(ns("download_pdf"), "Download PDF", class = "btn-sm"),
            downloadButton(ns("download_png"), "Download PNG", class = "btn-sm"),
            downloadButton(ns("download_svg"), "Download SVG", class = "btn-sm")
          )
          )
        )
      )
    )
  )
}

# Server function
mod_heatmap_unified_server <- function(id, app_data, global_selection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Load required packages and functions
    # Try ComplexHeatmap first, fallback to pheatmap
    use_complex_heatmap <- FALSE
    if (requireNamespace("ComplexHeatmap", quietly = TRUE) && 
        requireNamespace("circlize", quietly = TRUE)) {
      tryCatch({
        if (!exists("create_unified_enrichment_heatmap")) {
          source("unified_enrichment_heatmaps.R")
          cat("[HEATMAP] Sourced unified heatmap functions (ComplexHeatmap)\n")
        }
        use_complex_heatmap <- TRUE
      }, error = function(e) {
        cat("[HEATMAP] Error loading ComplexHeatmap functions:", e$message, "\n")
        use_complex_heatmap <<- FALSE
      })
    }
    
    if (!use_complex_heatmap) {
      cat("[HEATMAP] ComplexHeatmap not available, using pheatmap fallback\n")
      if (!requireNamespace("pheatmap", quietly = TRUE)) {
        install.packages("pheatmap")
      }
      library(pheatmap)
      library(RColorBrewer)
    }
    
    # Reactive values for this module
    module_data <- reactiveValues(
      current_heatmap = NULL,
      log_output = NULL
    )
    
    # Get filtered data based on global selection
    filtered_data <- reactive({
      req(app_data$data_loaded)
      req(app_data$consolidated_data)
      selection <- global_selection()
      
      # For heatmaps, we want data across multiple genes to create a meaningful visualization
      # The heatmap shows enrichment patterns across genes for the selected enrichment type
      if (input$method_comparison %in% c("intersection", "union", "all")) {
        # Get data across methods for all genes in the selected cluster
        data <- app_data$consolidated_data %>%
          filter(
            cluster == selection$cluster,
            enrichment_type == selection$enrichment_type,
            p.adjust <= selection$pval_threshold
          )
      } else {
        # Single method across multiple genes
        data <- app_data$consolidated_data %>%
          filter(
            method == input$method_comparison,
            cluster == selection$cluster,
            enrichment_type == selection$enrichment_type,
            p.adjust <= selection$pval_threshold
          )
      }
      
      # Debug logging
      cat("[HEATMAP] Filtered data rows:", nrow(data), "\n")
      cat("[HEATMAP] Method comparison:", input$method_comparison, "\n")
      cat("[HEATMAP] Selection cluster:", selection$cluster, "\n")
      cat("[HEATMAP] Selection enrichment type:", selection$enrichment_type, "\n")
      cat("[HEATMAP] Unique genes in filtered data:", n_distinct(data$mutation_perturbation), "\n")
      
      return(data)
    })
    
    # Main reactive heatmap generation
    current_heatmap <- reactive({
      data <- filtered_data()
      
      cat("[HEATMAP] Starting heatmap generation\n")
      cat("[HEATMAP] Data rows:", ifelse(is.null(data), "NULL", nrow(data)), "\n")
      
      if (is.null(data) || nrow(data) == 0) {
        cat("[HEATMAP] No data available, returning NULL\n")
        return(NULL)
      }
      
      selection <- global_selection()
      cat("[HEATMAP] Selected metric:", input$heatmap_metric, "\n")
      
      # Capture console output for logging
      log_output <- capture.output({
        
        if (use_complex_heatmap) {
          # Use ComplexHeatmap functions if available
          # Note: Since data is already filtered to one cluster, set cluster_val to NULL
          # to prevent double-filtering
          tryCatch({
            cat("[HEATMAP] Creating heatmap with:\n")
            cat("[HEATMAP] - Data rows:", nrow(data), "\n")
            cat("[HEATMAP] - Unique terms:", n_distinct(data$Description), "\n")
            cat("[HEATMAP] - Unique genes:", n_distinct(data$mutation_perturbation), "\n")
            cat("[HEATMAP] - Required columns present:\n")
            cat("[HEATMAP]   - p.adjust:", "p.adjust" %in% names(data), "\n")
            cat("[HEATMAP]   - Description:", "Description" %in% names(data), "\n")
            cat("[HEATMAP]   - mutation_perturbation:", "mutation_perturbation" %in% names(data), "\n")
            
            if (input$heatmap_metric == "pvalue") {
              cat("[HEATMAP] Calling create_pvalue_heatmap\n")
              ht <- create_pvalue_heatmap(
                data,
                cluster_val = selection$cluster,  # Pass cluster for proper title generation
                method_val = input$method_comparison,
                enrichment_types_to_include = selection$enrichment_type,
                max_terms_overall = input$max_terms,
                show_direction_annotation = input$show_directions,
                show_method_breakdown = TRUE,
                color_scaling_strategy = input$color_scheme,
                custom_title = if (input$custom_title != "") input$custom_title else NULL
              )
              cat("[HEATMAP] create_pvalue_heatmap returned:", class(ht)[1], "\n")
            } else if (input$heatmap_metric == "foldenrichment") {
              ht <- create_fold_enrichment_heatmap(
                data,
                cluster_val = selection$cluster,  # Pass cluster for proper title generation
                method_val = input$method_comparison,
                enrichment_types_to_include = selection$enrichment_type,
                max_terms_overall = input$max_terms,
                show_direction_annotation = input$show_directions,
                color_scaling_strategy = input$color_scheme,
                custom_title = if (input$custom_title != "") input$custom_title else NULL
              )
            } else if (input$heatmap_metric == "zscore") {
              ht <- create_zscore_heatmap(
                data,
                cluster_val = selection$cluster,  # Pass cluster for proper title generation
                method_val = input$method_comparison,
                enrichment_types_to_include = selection$enrichment_type,
                max_terms_overall = input$max_terms,
                show_direction_annotation = input$show_directions,
                color_scaling_strategy = input$color_scheme,
                custom_title = if (input$custom_title != "") input$custom_title else NULL
              )
            } else if (input$heatmap_metric == "gsea") {
              ht <- create_gsea_nes_heatmap(
                data,
                cluster_val = selection$cluster,  # Pass cluster for proper title generation
                method_val = input$method_comparison,
                max_terms_overall = input$max_terms,
                min_nes_magnitude = input$min_nes,
                custom_title = if (input$custom_title != "") input$custom_title else NULL
              )
            }
          }, error = function(e) {
            cat("[HEATMAP] Error creating heatmap:", e$message, "\n")
            ht <<- NULL
          })
        } else {
          # Fallback to pheatmap-based visualization
          tryCatch({
            ht <- create_simple_heatmap(
              data,
              metric = input$heatmap_metric,
              max_terms = input$max_terms,
              method_comparison = input$method_comparison,
              custom_title = if (input$custom_title != "") input$custom_title else NULL
            )
          }, error = function(e) {
            cat("[HEATMAP] Error creating simple heatmap:", e$message, "\n")
            ht <<- NULL
          })
        }
        
      })
      
      # Store in module data
      module_data$current_heatmap <- ht
      module_data$log_output <- log_output
      
      return(list(heatmap = ht, log = log_output))
    })
    
    # Simple heatmap creation function for pheatmap fallback
    create_simple_heatmap <- function(data, metric = "pvalue", max_terms = 25, 
                                    method_comparison = "all", custom_title = NULL) {
      
      cat("Creating simple heatmap with", nrow(data), "data points\n")
      
      if (nrow(data) == 0) return(NULL)
      
      # Prepare data matrix
      metric_col <- switch(metric,
        "pvalue" = "p.adjust",
        "foldenrichment" = "FoldEnrichment", 
        "zscore" = "zScore",
        "gsea" = "NES",
        "p.adjust"
      )
      
      if (!metric_col %in% names(data)) {
        cat("Metric column", metric_col, "not found in data\n")
        return(NULL)
      }
      
      # Transform values
      if (metric == "pvalue") {
        data$plot_value <- -log10(data[[metric_col]] + 1e-10)
      } else {
        data$plot_value <- data[[metric_col]]
      }
      
      # Select top terms
      top_terms <- data %>%
        group_by(Description) %>%
        summarise(max_val = max(plot_value, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_val)) %>%
        head(max_terms) %>%
        pull(Description)
      
      plot_data <- data %>%
        filter(Description %in% top_terms) %>%
        select(Description, mutation_perturbation, plot_value) %>%
        pivot_wider(names_from = mutation_perturbation, values_from = plot_value, 
                   values_fill = 0)
      
      # Create matrix
      mat <- as.matrix(plot_data[, -1])
      rownames(mat) <- make.unique(substr(plot_data$Description, 1, 50))
      
      if (nrow(mat) < 2 || ncol(mat) < 1) {
        cat("Matrix too small for heatmap\n")
        return(NULL)
      }
      
      # Color scale
      if (metric == "gsea") {
        colors <- RColorBrewer::brewer.pal(11, "RdBu")
        breaks <- seq(-max(abs(mat), na.rm = TRUE), max(abs(mat), na.rm = TRUE), length.out = 12)
      } else {
        colors <- RColorBrewer::brewer.pal(9, "YlOrRd")
        breaks <- seq(0, max(mat, na.rm = TRUE), length.out = 10)
      }
      
      # Create title
      title <- if (!is.null(custom_title)) {
        custom_title
      } else {
        paste(metric, "heatmap -", method_comparison)
      }
      
      return(list(
        type = "pheatmap",
        matrix = mat,
        colors = colors,
        breaks = breaks,
        title = title
      ))
    }
    
    # Render heatmap with dynamic height
    output$heatmap_display <- renderPlot({
      ht_data <- current_heatmap()
      
      # Check if we have valid heatmap data
      if (!is.null(ht_data) && !is.null(ht_data$heatmap)) {
        if (use_complex_heatmap) {
          # Additional safety check for ComplexHeatmap objects
          if (inherits(ht_data$heatmap, "Heatmap") || inherits(ht_data$heatmap, "HeatmapList")) {
            ComplexHeatmap::draw(ht_data$heatmap)
          } else {
            plot.new()
            text(0.5, 0.5, "Invalid heatmap object generated", 
                 cex = 1.5, col = "red")
          }
        } else {
          # Render pheatmap
          ht <- ht_data$heatmap
          if (!is.null(ht) && is.list(ht) && !is.null(ht$type) && ht$type == "pheatmap") {
            pheatmap(ht$matrix,
                    color = ht$colors,
                    breaks = ht$breaks,
                    main = ht$title,
                    cluster_rows = input$cluster_rows,
                    cluster_cols = input$cluster_columns,
                    show_rownames = input$show_row_names,
                    show_colnames = input$show_column_names,
                    fontsize_row = 8,
                    fontsize_col = 10)
          } else {
            plot.new()
            text(0.5, 0.5, "Invalid pheatmap object generated", 
                 cex = 1.5, col = "red")
          }
        }
      } else {
        plot.new()
        text(0.5, 0.5, "No data to display with current settings\n\nTry:\n• Different gene selection\n• Different enrichment type\n• Lower p-value threshold", 
             cex = 1.2, col = "gray50")
      }
    }, height = function() {
      # Dynamic height based on number of terms
      data <- filtered_data()
      if (is.null(data)) return(400)
      
      n_terms <- min(nrow(unique(data[, c("Description", "enrichment_type")])), 
                     input$max_terms)
      base_height <- 300
      per_term_height <- 15
      
      # Cap at reasonable maximum
      min(base_height + (n_terms * per_term_height), 1000)
    })
    
    # Method log output
    output$method_log <- renderPrint({
      ht_data <- current_heatmap()
      if (!is.null(ht_data$log)) {
        cat(paste(ht_data$log, collapse = "\n"))
      }
    })
    
    # Preset button handlers
    observeEvent(input$preset_methods, {
      selection <- global_selection()
      
      # Create comparison configs
      configs <- list(
        list(
          heatmap_type = "standard",
          metric_type = "p.adjust",
          method_val = "MAST",
          custom_title = "MAST Analysis"
        ),
        list(
          heatmap_type = "standard",
          metric_type = "p.adjust",
          method_val = "MixScale",
          custom_title = "MixScale Analysis"
        ),
        list(
          heatmap_type = "standard",
          metric_type = "p.adjust",
          method_val = "intersection",
          custom_title = "Shared Terms"
        )
      )
      
      # Create comparison heatmaps
      data <- filtered_data()
      if (!is.null(data) && nrow(data) > 0) {
        combined_ht <- tryCatch({
          create_comparison_heatmaps(
            data,
            configs,
            output_file = NULL  # Don't auto-save
          )
        }, error = function(e) {
          cat("[HEATMAP] Error creating method comparison:", e$message, "\n")
          NULL
        })
        
        # Update display
        output$heatmap_display <- renderPlot({
          if (!is.null(combined_ht)) {
            ComplexHeatmap::draw(combined_ht, ht_gap = unit(0.5, "cm"))
          } else {
            plot.new()
            text(0.5, 0.5, "Could not generate method comparison heatmaps", 
                 cex = 1.5, col = "red")
          }
        }, height = 600)
      } else {
        showNotification("No data available for method comparison with current selection", 
                        type = "warning", duration = 3)
      }
    })
    
    observeEvent(input$preset_metrics, {
      # Create metric comparison
      configs <- list(
        list(
          heatmap_type = "standard",
          metric_type = "p.adjust",
          custom_title = "Statistical Significance"
        ),
        list(
          heatmap_type = "standard",
          metric_type = "FoldEnrichment",
          custom_title = "Fold Enrichment"
        ),
        list(
          heatmap_type = "standard",
          metric_type = "zScore",
          custom_title = "Z-score"
        )
      )
      
      data <- filtered_data()
      if (!is.null(data) && nrow(data) > 0) {
        combined_ht <- tryCatch({
          create_comparison_heatmaps(
            data,
            configs,
            output_file = NULL
          )
        }, error = function(e) {
          cat("[HEATMAP] Error creating metric comparison:", e$message, "\n")
          NULL
        })
        
        output$heatmap_display <- renderPlot({
          if (!is.null(combined_ht)) {
            ComplexHeatmap::draw(combined_ht, ht_gap = unit(0.5, "cm"))
          } else {
            plot.new()
            text(0.5, 0.5, "Could not generate metric comparison heatmaps", 
                 cex = 1.5, col = "red")
          }
        }, height = 600)
      } else {
        showNotification("No data available for metric comparison with current selection", 
                        type = "warning", duration = 3)
      }
    })
    
    observeEvent(input$preset_directions, {
      # Update to show all directions
      showNotification("Generating direction comparison...", type = "message")
      # This would filter data by direction and create comparison
    })
    
    observeEvent(input$preset_comprehensive, {
      # Create comprehensive multi-panel view
      showNotification("Generating comprehensive view...", type = "message")
      # This would create a large multi-metric, multi-method comparison
    })
    
    # Manual refresh
    observeEvent(input$refresh_heatmap, {
      # Force re-calculation
      current_heatmap()
      showNotification("Heatmap refreshed", type = "message", duration = 2)
    })
    
    # Download handlers
    output$download_pdf <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.pdf")
      },
      content = function(file) {
        if (!is.null(module_data$current_heatmap)) {
          pdf(file, width = 12, height = 10)
          ComplexHeatmap::draw(module_data$current_heatmap)
          dev.off()
        } else {
          # Create empty PDF with error message
          pdf(file, width = 8, height = 6)
          plot.new()
          text(0.5, 0.5, "No heatmap available for download", cex = 1.5, col = "red")
          dev.off()
        }
      }
    )
    
    output$download_png <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.png")
      },
      content = function(file) {
        if (!is.null(module_data$current_heatmap)) {
          png(file, width = 1200, height = 1000, res = 150)
          ComplexHeatmap::draw(module_data$current_heatmap)
          dev.off()
        } else {
          # Create empty PNG with error message
          png(file, width = 800, height = 600, res = 150)
          plot.new()
          text(0.5, 0.5, "No heatmap available for download", cex = 1.5, col = "red")
          dev.off()
        }
      }
    )
    
    output$download_svg <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.svg")
      },
      content = function(file) {
        if (!is.null(module_data$current_heatmap)) {
          svg(file, width = 12, height = 10)
          ComplexHeatmap::draw(module_data$current_heatmap)
          dev.off()
        } else {
          # Create empty SVG with error message
          svg(file, width = 8, height = 6)
          plot.new()
          text(0.5, 0.5, "No heatmap available for download", cex = 1.5, col = "red")
          dev.off()
        }
      }
    )
    
    # Return current heatmap for use by other modules
    return(reactive({
      module_data$current_heatmap
    }))
  })
}