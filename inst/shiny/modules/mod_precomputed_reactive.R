# Module: Pre-computed Results Browser (Reactive Version)
# Fully reactive implementation without manual load/update buttons

# UI function
mod_precomputed_reactive_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Status indicators
    fluidRow(
      column(12,
        uiOutput(ns("status_bar"))
      )
    ),
    
    fluidRow(
      column(12,
        div(
          class = "box box-primary",
          div(class = "box-header with-border",
            h3(class = "box-title", "Enrichment Analysis Results")
          ),
          div(class = "box-body",
            # Tab panels for different views
            tabsetPanel(
              id = ns("viz_tabs"),
              selected = "visualization",  # Make visualization the default tab
          
          # Statistical Analysis (Tier 2)
          tabPanel(
            "Statistical Analysis",
            value = "stats",
            icon = icon("chart-bar"),
            br(),
            fluidRow(
              column(3,
                wellPanel(
                  h4("Heatmap Settings"),
                  
                  selectInput(ns("heatmap_type"),
                              "Visualization Type:",
                              choices = c("P-value Significance" = "pvalue",
                                        "Fold Enrichment" = "foldenrichment",
                                        "Z-score" = "zscore"),
                              selected = "pvalue"),
                  
                  selectInput(ns("method_comparison"),
                              "Method Comparison:",
                              choices = c("All Methods" = "all",
                                        "MAST Only" = "MAST",
                                        "MixScale Only" = "MixScale",
                                        "Intersection" = "intersection",
                                        "Union" = "union"),
                              selected = "all"),
                  
                  checkboxInput(ns("show_directions"),
                               "Show UP/DOWN/ALL Annotations",
                               value = TRUE),
                  
                  numericInput(ns("max_terms_heatmap"),
                              "Maximum Terms:",
                              value = 25,
                              min = 10,
                              max = 50),
                  
                  # Preset buttons
                  hr(),
                  h5("Quick Comparisons"),
                  actionButton(ns("preset_methods"), "Compare Methods", 
                              icon = icon("exchange-alt"), class = "btn-sm btn-block"),
                  actionButton(ns("preset_directions"), "Compare Directions", 
                              icon = icon("arrows-alt-v"), class = "btn-sm btn-block")
                )
              ),
              column(9,
                # Live heatmap
                withSpinner(
                  plotOutput(ns("heatmap_display"), height = "auto"),
                  type = 6,
                  color = "#3c8dbc"
                ),
                br(),
                # Method breakdown log
                verbatimTextOutput(ns("heatmap_log"))
              )
            )
          ),
          
          # Results Visualization (formerly Detailed Analysis)
          tabPanel(
            "Results Visualization",
            value = "visualization",
            icon = icon("chart-bar"),
            br(),
            fluidRow(
              column(3,
                wellPanel(
                  h4("Plot Settings"),
                  
                  # Plot type selection
                  radioButtons(ns("plot_type"),
                              "Plot Type:",
                              choices = c("Dot Plot" = "dotplot",
                                        "Bar Plot" = "barplot",
                                        "Lollipop Plot" = "lollipop"),
                              selected = "dotplot"),
                  
                  hr(),
                  
                  # X-axis selection (for dot and lollipop plots)
                  conditionalPanel(
                    condition = "input.plot_type == 'dotplot' || input.plot_type == 'lollipop'",
                    ns = ns,
                    selectInput(ns("x_axis"),
                               "X-axis:",
                               choices = c("Gene Ratio" = "GeneRatio",
                                         "Count" = "Count",
                                         "Fold Enrichment" = "FoldEnrichment"),
                               selected = "GeneRatio")
                  ),
                  
                  # Dot size (for dotplot)
                  conditionalPanel(
                    condition = "input.plot_type == 'dotplot'",
                    ns = ns,
                    selectInput(ns("dot_size"),
                               "Dot Size:",
                               choices = c("Count" = "Count",
                                         "-log10(p-value)" = "neglog10p",
                                         "-log10(adjusted p-value)" = "neglog10padj"),
                               selected = "Count"),
                    
                    sliderInput(ns("size_range"),
                               "Size Range:",
                               min = 1, max = 10,
                               value = c(2, 6),
                               step = 0.5)
                  ),
                  
                  # Color mapping
                  selectInput(ns("color_by"),
                             "Color By:",
                             choices = c("Adjusted P-value" = "p.adjust",
                                       "P-value" = "pvalue",
                                       "Fold Enrichment" = "FoldEnrichment"),
                             selected = "p.adjust"),
                  
                  # Color scale
                  fluidRow(
                    column(6,
                      colourInput(ns("low_color"), "Low Color:", "#FF0000")
                    ),
                    column(6,
                      colourInput(ns("high_color"), "High Color:", "#0000FF")
                    )
                  ),
                  
                  hr(),
                  
                  # Number of terms
                  numericInput(ns("n_terms"),
                              "Number of Terms:",
                              value = 20,
                              min = 5,
                              max = 50,
                              step = 5),
                  
                  # Order by
                  selectInput(ns("order_by"),
                             "Order By:",
                             choices = c("P-value" = "pvalue",
                                       "Adjusted P-value" = "p.adjust",
                                       "Count" = "Count",
                                       "Fold Enrichment" = "FoldEnrichment"),
                             selected = "p.adjust"),
                  
                  checkboxInput(ns("show_id"),
                               "Show Term IDs",
                               value = FALSE)
                )
              ),
              column(9,
                withSpinner(
                  plotOutput(ns("enrichplot_display"), height = "600px"),
                  type = 4,
                  color = "#3c8dbc"
                ),
                br(),
                div(style = "text-align: center;",
                  downloadButton(ns("download_plot"), "Download Plot", class = "btn-primary"),
                  downloadButton(ns("download_data"), "Download Data", class = "btn-info")
                )
              )
            )
          ),
          
          # Data Table
          tabPanel(
            "Data Table",
            value = "table",
            icon = icon("table"),
            br(),
            fluidRow(
              column(12,
                h4("Significant Enrichment Terms"),
                br(),
                # Export options
                div(style = "margin-bottom: 10px;",
                  downloadButton(ns("download_csv"), "Download CSV", class = "btn-sm"),
                  downloadButton(ns("download_excel"), "Download Excel", class = "btn-sm")
                ),
                DT::dataTableOutput(ns("enrichment_table"))
              )
            )
          )
        ) # end tabsetPanel
        ) # end box-body div
      ) # end box div
      ) # end column
    ) # end fluidRow
  )
}

# Server function
mod_precomputed_reactive_server <- function(id, app_data, global_selection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Source required functions
    if (!exists("create_tier1_visualizations")) {
      source("R/visualization_tiers.R")
    }
    if (!exists("create_unified_enrichment_heatmap")) {
      source("unified_enrichment_heatmaps.R")
    }
    
    # Reactive filtered data
    filtered_data <- reactive({
      req(app_data$data_loaded)
      selection <- global_selection()
      
      get_significant_terms_from_consolidated(
        app_data$consolidated_data,
        analysis_type = selection$analysis_type,
        gene = selection$gene,
        cluster = selection$cluster,
        experiment = selection$experiment,
        enrichment_type = selection$enrichment_type,
        direction = selection$direction,
        pval_threshold = selection$pval_threshold
      )
    })
    
    # Status bar
    output$status_bar <- renderUI({
      data <- filtered_data()
      selection <- global_selection()
      
      if (!is.null(data) && nrow(data) > 0) {
        # Calculate unique genes if geneID column exists
        n_genes <- if ("geneID" %in% names(data)) {
          length(unique(unlist(strsplit(data$geneID, "/"))))
        } else NA
        
        tags$div(class = "well well-sm",
          tags$span(class = "text-primary", icon("database"), 
                   sprintf("%d significant terms", nrow(data))),
          tags$span(style = "margin-left: 20px;", icon("filter"),
                   sprintf("p < %g", selection$pval_threshold)),
          if (!is.na(n_genes)) {
            tags$span(style = "margin-left: 20px;", icon("dna"),
                     sprintf("%d unique genes", n_genes))
          },
          tags$span(style = "margin-left: 20px;", class = "text-success",
                   icon("bolt"), "Summary Data")
        )
      } else {
        # Check what enrichment types are available for this selection
        selection <- global_selection()
        available_msg <- NULL
        
        if (!is.null(app_data$consolidated_data)) {
          available_types <- unique(app_data$consolidated_data[
            app_data$consolidated_data$method == selection$analysis_type &
            app_data$consolidated_data$gene == selection$gene &
            app_data$consolidated_data$cluster == selection$cluster &
            app_data$consolidated_data$experiment == selection$experiment,
            "enrichment_type"
          ])
          
          if (length(available_types) > 0 && !(selection$enrichment_type %in% available_types)) {
            available_msg <- tags$div(
              tags$br(),
              tags$strong("Available enrichment types for this selection: "),
              paste(sort(available_types), collapse = ", "),
              tags$br(),
              tags$em("Try selecting one of these from the dropdown menu.")
            )
          }
        }
        
        tags$div(class = "alert alert-warning",
          icon("exclamation-triangle"),
          "No significant terms found with current settings",
          available_msg
        )
      }
    })
    
    # Tier 1: Quick View Visualizations
    tier1_plots <- reactive({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      create_tier1_visualizations(data)
    })
    
    output$summary_table <- DT::renderDataTable({
      plots <- tier1_plots()
      if (!is.null(plots$summary_table)) {
        DT::datatable(plots$summary_table,
                     options = list(pageLength = 10, dom = 'tip'),
                     rownames = FALSE)
      }
    })
    
    output$term_counts <- renderPlot({
      plots <- tier1_plots()
      plots$term_counts
    })
    
    output$pvalue_dist <- renderPlot({
      plots <- tier1_plots()
      plots$pvalue_distribution
    })
    
    output$top_terms_bar <- renderPlot({
      plots <- tier1_plots()
      plots$top_terms_bar
    })
    
    # Tier 2: Statistical Analysis (Heatmaps)
    current_heatmap <- reactive({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      selection <- global_selection()
      
      # Capture logging output
      log_output <- capture.output({
        if (input$heatmap_type == "pvalue") {
          ht <- create_pvalue_heatmap(
            data,
            cluster_val = selection$cluster,
            method_val = input$method_comparison,
            enrichment_types_to_include = selection$enrichment_type,
            max_terms_overall = input$max_terms_heatmap,
            show_direction_annotation = input$show_directions,
            show_method_breakdown = TRUE
          )
        } else if (input$heatmap_type == "foldenrichment") {
          ht <- create_fold_enrichment_heatmap(
            data,
            cluster_val = selection$cluster,
            method_val = input$method_comparison,
            enrichment_types_to_include = selection$enrichment_type,
            max_terms_overall = input$max_terms_heatmap,
            show_direction_annotation = input$show_directions
          )
        } else if (input$heatmap_type == "zscore") {
          ht <- create_zscore_heatmap(
            data,
            cluster_val = selection$cluster,
            method_val = input$method_comparison,
            enrichment_types_to_include = selection$enrichment_type,
            max_terms_overall = input$max_terms_heatmap,
            show_direction_annotation = input$show_directions
          )
        }
      })
      
      list(heatmap = ht, log = log_output)
    })
    
    # Dynamic height for heatmap
    output$heatmap_display <- renderPlot({
      ht_data <- current_heatmap()
      if (!is.null(ht_data$heatmap)) {
        ComplexHeatmap::draw(ht_data$heatmap)
      }
    }, height = function() {
      data <- filtered_data()
      if (is.null(data)) return(400)
      n_terms <- min(nrow(data), input$max_terms_heatmap)
      base_height <- 300
      per_term_height <- 15
      min(base_height + (n_terms * per_term_height), 800)
    })
    
    output$heatmap_log <- renderPrint({
      ht_data <- current_heatmap()
      if (!is.null(ht_data$log)) {
        cat(paste(ht_data$log, collapse = "\n"))
      }
    })
    
    # Preset comparisons
    observeEvent(input$preset_methods, {
      showNotification("Generating method comparison...", type = "message")
      updateSelectInput(session, "method_comparison", selected = "all")
      # Could implement side-by-side comparison here
    })
    
    observeEvent(input$preset_directions, {
      showNotification("Generating direction comparison...", type = "message")
      # Could implement UP vs DOWN comparison here
    })
    
    # Tier 3: Detailed Analysis (requires RDS loading)
    observeEvent(input$plot_type, {
      if (input$viz_tabs == "detailed") {
        selection <- global_selection()
        
        # Check if we need to load RDS data
        showNotification(
          "Loading full enrichment data for detailed visualization...",
          type = "info",
          duration = 2
        )
        
        # Construct file path
        file_path <- file.path(
          APP_CONFIG$enrichment_results_path,
          selection$analysis_type,
          selection$gene,
          selection$cluster,
          selection$experiment,
          selection$enrichment_type,
          paste0(selection$enrichment_type, "_", selection$direction, ".rds")
        )
        
        # Load with cache
        if (exists("load_enrichment_with_cache")) {
          rds_data <- load_enrichment_with_cache(
            selection$analysis_type,
            selection$gene,
            selection$cluster,
            selection$experiment,
            selection$enrichment_type,
            selection$direction
          )
        } else {
          rds_data <- load_enrichment_safe(file_path)
        }
        
        if (!is.null(rds_data)) {
          # Store the loaded data for plotting
          app_data$current_plot_data <- rds_data
        } else {
          showNotification(
            "Could not load full enrichment data",
            type = "error"
          )
        }
      }
    })
    
    # Enhanced visualization plot (ShinyGO-style)
    output$enrichplot_display <- renderPlot({
      data <- filtered_data()
      
      if (is.null(data) || nrow(data) == 0) {
        return(NULL)
      }
      
      # Order the data
      data <- data[order(data[[input$order_by]]), ]
      
      # Take top N terms
      if (nrow(data) > input$n_terms) {
        data <- data[1:input$n_terms, ]
      }
      
      # Reverse order for better display (best terms at top)
      data <- data[nrow(data):1, ]
      
      # Calculate fold enrichment if not present
      if (!"FoldEnrichment" %in% names(data) && "GeneRatio" %in% names(data)) {
        # Parse GeneRatio
        gr_parts <- strsplit(data$GeneRatio, "/")
        data$FoldEnrichment <- sapply(gr_parts, function(x) {
          if(length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) * 100
          else NA
        })
      }
      
      # Add -log10 p-values
      data$neglog10p <- -log10(data$pvalue)
      data$neglog10padj <- -log10(data$p.adjust)
      
      # Prepare term labels
      if (input$show_id && "ID" %in% names(data)) {
        data$term_label <- paste0(data$ID, ": ", data$Description)
      } else {
        data$term_label <- data$Description
      }
      
      # Truncate long labels
      data$term_label <- ifelse(nchar(data$term_label) > 50,
                               paste0(substr(data$term_label, 1, 47), "..."),
                               data$term_label)
      
      # Create the plot based on type
      if (input$plot_type == "dotplot") {
        # Dotplot
        p <- ggplot(data, aes_string(x = input$x_axis, y = "term_label")) +
          geom_point(aes_string(size = input$dot_size, color = input$color_by)) +
          scale_size_continuous(range = input$size_range) +
          scale_color_gradient(low = input$low_color, high = input$high_color) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.position = "right",
            panel.grid.major.y = element_blank()
          ) +
          labs(x = gsub("_", " ", input$x_axis),
               y = "",
               size = gsub("_", " ", input$dot_size),
               color = gsub("\\.", " ", input$color_by))
        
      } else if (input$plot_type == "barplot") {
        # Barplot
        p <- ggplot(data, aes_string(x = input$order_by, y = "term_label")) +
          geom_col(aes_string(fill = input$color_by)) +
          scale_fill_gradient(low = input$low_color, high = input$high_color) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.position = "right"
          ) +
          labs(x = gsub("\\.", " ", input$order_by),
               y = "",
               fill = gsub("\\.", " ", input$color_by))
        
      } else if (input$plot_type == "lollipop") {
        # Lollipop plot
        p <- ggplot(data, aes_string(x = input$x_axis, y = "term_label")) +
          geom_segment(aes_string(x = 0, xend = input$x_axis, 
                                 y = "term_label", yend = "term_label",
                                 color = input$color_by)) +
          geom_point(aes_string(color = input$color_by), size = 3) +
          scale_color_gradient(low = input$low_color, high = input$high_color) +
          theme_minimal() +
          theme(
            axis.text.y = element_text(size = 10),
            axis.text.x = element_text(size = 10),
            legend.position = "right",
            panel.grid.major.y = element_blank()
          ) +
          labs(x = gsub("_", " ", input$x_axis),
               y = "",
               color = gsub("\\.", " ", input$color_by))
      }
      
      # Set plot height based on number of terms
      plot_height <- max(400, input$n_terms * 20)
      
      p
    }, height = function() {
      # Dynamic height based on number of terms
      max(400, input$n_terms * 20)
    })
    
    # Download plot handler
    output$download_plot <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$cluster, "_",
               selection$enrichment_type, "_", input$plot_type, ".pdf")
      },
      content = function(file) {
        data <- filtered_data()
        if (!is.null(data) && nrow(data) > 0) {
          # Recreate the plot for download
          # (Same plotting code as above, just saved to file)
          pdf(file, width = 10, height = max(6, input$n_terms * 0.3))
          
          # ... (abbreviated - would include full plot recreation code)
          
          dev.off()
        }
      }
    )
    
    # Data table
    output$enrichment_table <- DT::renderDataTable({
      data <- filtered_data()
      if (is.null(data) || nrow(data) == 0) return(NULL)
      
      # Select relevant columns for display
      display_cols <- c("Description", "p.adjust", "pvalue", "Count", "GeneRatio", 
                       "enrichment_type", "direction", "geneID")
      display_cols <- intersect(display_cols, names(data))
      
      # Rename geneID column if present for clarity
      if ("geneID" %in% display_cols) {
        names(data)[names(data) == "geneID"] <- "Genes"
      }
      
      DT::datatable(
        data[, display_cols],
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          order = list(list(2, 'asc')),  # Order by p.adjust
          columnDefs = list(
            list(
              targets = which(display_cols == "Genes") - 1,  # Adjust for 0-based indexing
              render = JS(
                "function(data, type, row, meta) {
                  if(type === 'display' && data && data.length > 50) {
                    return data.substr(0, 50) + '...';
                  }
                  return data;
                }"
              )
            )
          )
        ),
        filter = 'top',
        rownames = FALSE,
        selection = 'single'
      ) %>%
        DT::formatSignif(columns = c("p.adjust", "pvalue"), digits = 3)
    })
    
    # Download data handler for visualization
    output$download_data <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$cluster, "_",
               selection$enrichment_type, "_data.csv")
      },
      content = function(file) {
        data <- filtered_data()
        if (!is.null(data) && nrow(data) > 0) {
          # Order and filter data same as plot
          data <- data[order(data[[input$order_by]]), ]
          if (nrow(data) > input$n_terms) {
            data <- data[1:input$n_terms, ]
          }
          write.csv(data, file, row.names = FALSE)
        }
      }
    )
    
    # Download handlers
    output$download_csv <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$cluster, "_", 
               selection$enrichment_type, "_enrichment.csv")
      },
      content = function(file) {
        data <- filtered_data()
        write.csv(data, file, row.names = FALSE)
      }
    )
    
    output$download_excel <- downloadHandler(
      filename = function() {
        selection <- global_selection()
        paste0(selection$gene, "_", selection$cluster, "_",
               selection$enrichment_type, "_enrichment.xlsx")
      },
      content = function(file) {
        data <- filtered_data()
        openxlsx::write.xlsx(data, file)
      }
    )
    
    # Return filtered data for other modules
    return(reactive({
      filtered_data()
    }))
  })
}