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
    
    # Source unified heatmap functions if not already loaded
    if (!exists("create_unified_enrichment_heatmap")) {
      source("unified_enrichment_heatmaps.R")
    }
    
    # Reactive values for this module
    module_data <- reactiveValues(
      current_heatmap = NULL,
      log_output = NULL
    )
    
    # Get filtered data based on global selection
    filtered_data <- reactive({
      req(app_data$data_loaded)
      selection <- global_selection()
      
      # For heatmaps, we might want broader data than just one selection
      # Depending on the comparison type
      if (input$method_comparison %in% c("intersection", "union", "all")) {
        # Get data across methods
        data <- app_data$consolidated_data %>%
          filter(
            gene == selection$gene,
            enrichment_type == selection$enrichment_type,
            p.adjust <= selection$pval_threshold
          )
      } else {
        # Single method
        data <- app_data$consolidated_data %>%
          filter(
            method == input$method_comparison,
            gene == selection$gene,
            enrichment_type == selection$enrichment_type,
            p.adjust <= selection$pval_threshold
          )
      }
      
      return(data)
    })
    
    # Main reactive heatmap generation
    current_heatmap <- reactive({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      selection <- global_selection()
      
      # Capture console output for logging
      log_output <- capture.output({
        
        # Determine which heatmap function to use
        if (input$heatmap_metric == "pvalue") {
          ht <- create_pvalue_heatmap(
            data,
            cluster_val = selection$cluster,
            method_val = input$method_comparison,
            enrichment_types_to_include = selection$enrichment_type,
            max_terms_overall = input$max_terms,
            show_direction_annotation = input$show_directions,
            show_method_breakdown = TRUE,
            color_scaling_strategy = input$color_scheme,
            custom_title = if (input$custom_title != "") input$custom_title else NULL
          )
        } else if (input$heatmap_metric == "foldenrichment") {
          ht <- create_fold_enrichment_heatmap(
            data,
            cluster_val = selection$cluster,
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
            cluster_val = selection$cluster,
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
            cluster_val = selection$cluster,
            method_val = input$method_comparison,
            max_terms_overall = input$max_terms,
            min_nes_magnitude = input$min_nes,
            custom_title = if (input$custom_title != "") input$custom_title else NULL
          )
        }
        
      })
      
      # Store in module data
      module_data$current_heatmap <- ht
      module_data$log_output <- log_output
      
      return(list(heatmap = ht, log = log_output))
    })
    
    # Render heatmap with dynamic height
    output$heatmap_display <- renderPlot({
      ht_data <- current_heatmap()
      if (!is.null(ht_data$heatmap)) {
        ComplexHeatmap::draw(ht_data$heatmap)
      } else {
        plot.new()
        text(0.5, 0.5, "No data to display with current settings", 
             cex = 1.5, col = "gray50")
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
      if (!is.null(data)) {
        combined_ht <- create_comparison_heatmaps(
          data,
          configs,
          output_file = NULL  # Don't auto-save
        )
        
        # Update display
        output$heatmap_display <- renderPlot({
          ComplexHeatmap::draw(combined_ht, ht_gap = unit(0.5, "cm"))
        }, height = 600)
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
      if (!is.null(data)) {
        combined_ht <- create_comparison_heatmaps(
          data,
          configs,
          output_file = NULL
        )
        
        output$heatmap_display <- renderPlot({
          ComplexHeatmap::draw(combined_ht, ht_gap = unit(0.5, "cm"))
        }, height = 600)
      }
    })
    
    observeEvent(input$preset_directions, {
      # Update to show all directions
      showNotification("Generating direction comparison...", type = "info")
      # This would filter data by direction and create comparison
    })
    
    observeEvent(input$preset_comprehensive, {
      # Create comprehensive multi-panel view
      showNotification("Generating comprehensive view...", type = "info")
      # This would create a large multi-metric, multi-method comparison
    })
    
    # Manual refresh
    observeEvent(input$refresh_heatmap, {
      # Force re-calculation
      current_heatmap()
      showNotification("Heatmap refreshed", type = "success", duration = 2)
    })
    
    # Download handlers
    output$download_pdf <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.pdf")
      },
      content = function(file) {
        pdf(file, width = 12, height = 10)
        ComplexHeatmap::draw(module_data$current_heatmap)
        dev.off()
      }
    )
    
    output$download_png <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.png")
      },
      content = function(file) {
        png(file, width = 1200, height = 1000, res = 150)
        ComplexHeatmap::draw(module_data$current_heatmap)
        dev.off()
      }
    )
    
    output$download_svg <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$enrichment_type, "_",
               input$heatmap_metric, "_heatmap.svg")
      },
      content = function(file) {
        svg(file, width = 12, height = 10)
        ComplexHeatmap::draw(module_data$current_heatmap)
        dev.off()
      }
    )
    
    # Return current heatmap for use by other modules
    return(reactive({
      module_data$current_heatmap
    }))
  })
}