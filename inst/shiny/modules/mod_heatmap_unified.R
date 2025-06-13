# Module: Unified Heatmap Visualization
# Integrates the unified heatmap system into Shiny

# Required libraries
library(dplyr)
library(tidyr)

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
      
      # Create simplified heatmap using backup approach
      tryCatch({
        df <- data
        cat("[HEATMAP] Creating heatmap with ", nrow(df), " rows\n")
        cat("[HEATMAP] Available columns: ", paste(names(df), collapse = ", "), "\n")
        
        # Determine what to use as columns (x-axis) - use backup logic
        if ("mutation_perturbation" %in% names(df)) {
          x_var <- "mutation_perturbation"
        } else if ("gene" %in% names(df)) {
          x_var <- "gene"
        } else if ("cluster" %in% names(df)) {
          x_var <- "cluster"
        } else {
          # Create a composite identifier
          df$condition <- paste(df$method, df$cluster, sep = "_")
          x_var <- "condition"
        }
        
        # Use Description or term_name for y-axis
        if ("Description" %in% names(df)) {
          y_var <- "Description"
        } else if ("term_name" %in% names(df)) {
          y_var <- "term_name"
        } else {
          # Create generic term names
          df$term_id <- paste("Term", 1:nrow(df))
          y_var <- "term_id"
        }
        
        # Determine value column based on metric selection
        if (input$heatmap_metric == "pvalue" && "p.adjust" %in% names(df)) {
          df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
          value_var <- "heatmap_value"
          legend_title <- "-log10(p-value)"
        } else if (input$heatmap_metric == "foldenrichment" && "FoldEnrichment" %in% names(df)) {
          df$heatmap_value <- df$FoldEnrichment
          value_var <- "heatmap_value"
          legend_title <- "Fold Enrichment"
        } else if (input$heatmap_metric == "zscore" && "zScore" %in% names(df)) {
          df$heatmap_value <- df$zScore
          value_var <- "heatmap_value"
          legend_title <- "Z-Score"
        } else if (input$heatmap_metric == "gsea" && "NES" %in% names(df)) {
          df$heatmap_value <- df$NES
          value_var <- "heatmap_value"
          legend_title <- "NES"
        } else {
          # Default to a simple presence/absence
          df$heatmap_value <- 1
          value_var <- "heatmap_value"
          legend_title <- "Presence"
        }
        
        cat("[HEATMAP] Using x_var:", x_var, ", y_var:", y_var, ", value_var:", value_var, "\n")
        cat("[HEATMAP] Unique x values:", length(unique(df[[x_var]])), "\n")
        cat("[HEATMAP] Unique y values:", length(unique(df[[y_var]])), "\n")
        
        # Create the heatmap result
        ht_result <- list(
          type = "simple",
          data = df,
          x_var = x_var,
          y_var = y_var,
          value_var = value_var,
          legend_title = legend_title,
          custom_title = if (input$custom_title != "") input$custom_title else 
                        paste(legend_title, "- Cluster", selection$cluster)
        )
        
        cat("[HEATMAP] Successfully created simple heatmap data\n")
        return(list(heatmap = ht_result, log = "Simple heatmap created successfully"))
        
      }, error = function(e) {
        cat("[HEATMAP] Error creating heatmap:", e$message, "\n")
        cat("[HEATMAP] Error details:", toString(e), "\n")
        return(list(heatmap = NULL, log = paste("Error:", e$message)))
      })
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
        ht <- ht_data$heatmap
        
        # Handle simplified heatmap from backup approach
        if (is.list(ht) && !is.null(ht$type) && ht$type == "simple") {
          tryCatch({
            df <- ht$data
            x_var <- ht$x_var
            y_var <- ht$y_var
            value_var <- ht$value_var
            legend_title <- ht$legend_title
            
            cat("[RENDER] Creating heatmap plot with simple approach\n")
            cat("[RENDER] Matrix dimensions: ", length(unique(df[[y_var]])), " x ", length(unique(df[[x_var]])), "\n")
            
            # Check if we have enough data for a proper heatmap
            if (length(unique(df[[x_var]])) < 2) {
              # Create bar plot instead
              top_terms <- head(df[order(df[[value_var]], decreasing = TRUE), ], input$max_terms)
              
              barplot(top_terms[[value_var]], 
                      names.arg = substr(top_terms[[y_var]], 1, 30),
                      las = 2,
                      main = ht$custom_title,
                      ylab = legend_title,
                      cex.names = 0.7,
                      col = heat.colors(nrow(top_terms)))
              
            } else {
              # Create proper heatmap matrix using backup logic
              data_wide <- df %>%
                dplyr::select(all_of(c(y_var, x_var, value_var))) %>%
                tidyr::pivot_wider(names_from = all_of(x_var), 
                                  values_from = all_of(value_var), 
                                  values_fill = 0,
                                  values_fn = mean)  # Handle duplicates by averaging
              
              # Convert to matrix
              mat <- as.matrix(data_wide[,-1])
              rownames(mat) <- substr(data_wide[[y_var]], 1, 50)  # Truncate long names
              
              # Limit to max terms
              if (nrow(mat) > input$max_terms) {
                # Keep most significant terms (highest values)
                row_means <- rowMeans(mat, na.rm = TRUE)
                mat <- mat[order(row_means, decreasing = TRUE)[1:input$max_terms], ]
              }
              
              # Create color palette
              if (input$heatmap_metric == "pvalue") {
                col_fun <- colorRampPalette(c("white", "red"))(100)
              } else if (input$heatmap_metric == "gsea") {
                col_fun <- colorRampPalette(c("blue", "white", "red"))(100)
              } else {
                col_fun <- colorRampPalette(c("white", "blue"))(100)
              }
              
              # Basic heatmap using base R
              heatmap(mat,
                      col = col_fun,
                      scale = "none",
                      Rowv = if(input$cluster_rows) TRUE else NA,
                      Colv = if(input$cluster_columns) TRUE else NA,
                      margins = c(10, 15),
                      cexRow = 0.7,
                      cexCol = 0.8,
                      main = ht$custom_title)
            }
            
          }, error = function(e) {
            cat("[RENDER] Error in simple heatmap rendering:", e$message, "\n")
            plot.new()
            text(0.5, 0.5, paste("Error creating heatmap:\n", e$message), cex = 1.2, col = "red")
          })
          
        } else if (use_complex_heatmap && (inherits(ht, "Heatmap") || inherits(ht, "HeatmapList"))) {
          # ComplexHeatmap objects
          ComplexHeatmap::draw(ht)
          
        } else if (!use_complex_heatmap && is.list(ht) && !is.null(ht$type) && ht$type == "pheatmap") {
          # pheatmap objects
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
          text(0.5, 0.5, "Invalid heatmap object generated", 
               cex = 1.5, col = "red")
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