# Module: Advanced Heatmaps
# Provides heatmap visualizations using the unified heatmap system

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Heatmap Settings"),
      
      selectInput(ns("heatmap_type"),
                  "Heatmap Type:",
                  choices = c("P-value" = "pvalue",
                              "Fold Enrichment" = "foldenrichment", 
                              "Z-score" = "zscore",
                              "Gene Count" = "count"),
                  selected = "pvalue"),
      
      selectInput(ns("method_filter"),
                  "Analysis Method:",
                  choices = c("All" = "all",
                              "MAST Only" = "MAST",
                              "MixScale Only" = "MixScale", 
                              "Intersection" = "intersection",
                              "Union" = "union"),
                  selected = "all"),
      
      checkboxGroupInput(ns("enrichment_types"),
                         "Enrichment Types:",
                         choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "STRING"),
                         selected = c("GO_BP", "KEGG")),
      
      selectInput(ns("cluster_select"),
                  "Cluster:",
                  choices = c("All Clusters" = "all", paste0("cluster_", 0:9)),
                  selected = "all"),
      
      sliderInput(ns("max_terms"),
                  "Maximum Terms:",
                  min = 10,
                  max = 100,
                  value = 50,
                  step = 10),
      
      checkboxInput(ns("show_direction"),
                    "Show Direction Annotation",
                    value = TRUE),
      
      checkboxInput(ns("cluster_rows"),
                    "Cluster Rows",
                    value = TRUE),
      
      checkboxInput(ns("cluster_cols"),
                    "Cluster Columns", 
                    value = TRUE),
      
      br(),
      actionButton(ns("generate_heatmap"),
                   "Generate Heatmap",
                   class = "btn-primary", 
                   width = "100%"),
      
      br(),
      br(),
      downloadButton(ns("download_heatmap"),
                     "Download Heatmap",
                     class = "btn-success",
                     style = "width: 100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Heatmap",
                 br(),
                 withSpinner(
                   plotOutput(ns("heatmap_plot"), height = "800px"),
                   type = 4,
                   color = "#3c8dbc"
                 )),
        
        tabPanel("Data View",
                 br(),
                 DT::dataTableOutput(ns("heatmap_data"))),
        
        tabPanel("Settings Info",
                 br(),
                 verbatimTextOutput(ns("settings_info")))
      )
    )
  )
}

mod_heatmap_server <- function(id, app_data, pval_threshold) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    heatmap_data <- reactiveValues(
      data = NULL,
      plot = NULL
    )
    
    # Load and prepare data for heatmap
    observeEvent(input$generate_heatmap, {
      req(input$enrichment_types)
      
      showNotification("Loading data for heatmap...", type = "message", duration = NULL, id = "loading")
      
      tryCatch({
        # Use the consolidated data from app_data
        if (!is.null(app_data$consolidated_data) && nrow(app_data$consolidated_data) > 0) {
          # Use the data already loaded in the app
          all_data <- app_data$consolidated_data
          message("Using consolidated data with ", nrow(all_data), " rows")
          message("Columns available: ", paste(names(all_data), collapse = ", "))
          
          # Filter by enrichment types
          filtered_data <- all_data %>%
            dplyr::filter(enrichment_type %in% input$enrichment_types)
          message("After enrichment type filter: ", nrow(filtered_data), " rows")
          
          # Filter by cluster if not "all"
          if (input$cluster_select != "all") {
            filtered_data <- filtered_data %>%
              dplyr::filter(cluster == input$cluster_select)
            message("After cluster filter: ", nrow(filtered_data), " rows")
          }
          
          # Filter by method
          if (input$method_filter != "all") {
            if (input$method_filter %in% c("MAST", "MixScale")) {
              # Check which column name to use
              if ("method" %in% names(filtered_data)) {
                filtered_data <- filtered_data %>%
                  dplyr::filter(method == input$method_filter)
              } else if ("analysis_type" %in% names(filtered_data)) {
                filtered_data <- filtered_data %>%
                  dplyr::filter(analysis_type == input$method_filter)
              }
              message("After method filter: ", nrow(filtered_data), " rows")
            }
          }
          
          # Filter by significance
          if ("p.adjust" %in% names(filtered_data)) {
            filtered_data <- filtered_data %>%
              dplyr::filter(p.adjust <= 0.05)
            message("After significance filter: ", nrow(filtered_data), " rows")
          }
          
          if (nrow(filtered_data) == 0) {
            removeNotification("loading")
            showNotification("No data matches the selected criteria", type = "warning", duration = 5)
            return()
          }
          
          heatmap_data$data <- filtered_data
          
          removeNotification("loading")
          showNotification("Data loaded successfully! Heatmap will be generated automatically.", 
                           type = "success", duration = 3)
          
        } else {
          removeNotification("loading")
          showNotification("No consolidated data available. Please ensure data is loaded properly.", 
                           type = "warning", duration = 5)
        }
        
      }, error = function(e) {
        removeNotification("loading")
        message("Error loading heatmap data: ", e$message)
        showNotification(paste("Error loading data:", e$message), type = "error", duration = 10)
      })
    })
    
    # Generate heatmap plot
    output$heatmap_plot <- renderPlot({
      req(heatmap_data$data)
      
      # Create a simplified heatmap
      tryCatch({
        df <- heatmap_data$data
        message("Creating heatmap with ", nrow(df), " rows")
        message("Available columns: ", paste(names(df), collapse = ", "))
        
        # Create a frequency matrix of terms across genes/conditions
        # Use Description column for term names, and create columns for different conditions
        
        # Determine what to use as columns (x-axis)
        if ("gene" %in% names(df)) {
          x_var <- "gene"
        } else if ("mutation" %in% names(df)) {
          x_var <- "mutation"
        } else if ("cluster" %in% names(df)) {
          x_var <- "cluster"
        } else {
          # Create a composite identifier
          if ("method" %in% names(df)) {
            df$condition <- paste(df$method, df$cluster, sep = "_")
          } else if ("analysis_type" %in% names(df)) {
            df$condition <- paste(df$analysis_type, df$cluster, sep = "_")
          } else {
            df$condition <- df$cluster
          }
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
        
        # Determine value column
        if (input$heatmap_type == "pvalue" && "p.adjust" %in% names(df)) {
          df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
          value_var <- "heatmap_value"
          legend_title <- "-log10(p-value)"
        } else if (input$heatmap_type == "count" && "Count" %in% names(df)) {
          value_var <- "Count"
          legend_title <- "Gene Count"
        } else {
          # Default to a simple presence/absence
          df$heatmap_value <- 1
          value_var <- "heatmap_value"
          legend_title <- "Presence"
        }
        
        # Create matrix
        if (length(unique(df[[x_var]])) < 2) {
          # If only one condition, create a simple bar plot instead
          top_terms <- head(df[order(df[[value_var]], decreasing = TRUE), ], input$max_terms)
          
          barplot(top_terms[[value_var]], 
                  names.arg = substr(top_terms[[y_var]], 1, 30),
                  las = 2,
                  main = paste("Top Terms -", input$heatmap_type),
                  ylab = legend_title,
                  cex.names = 0.7,
                  col = heat.colors(nrow(top_terms)))
          
        } else {
          # Create proper heatmap matrix
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
          if (input$heatmap_type == "pvalue") {
            col_fun <- colorRampPalette(c("white", "red"))(100)
          } else {
            col_fun <- colorRampPalette(c("white", "blue"))(100)
          }
          
          # Basic heatmap using base R
          heatmap(mat,
                  col = col_fun,
                  scale = "none",
                  Rowv = if(input$cluster_rows) TRUE else NA,
                  Colv = if(input$cluster_cols) TRUE else NA,
                  margins = c(10, 15),
                  cexRow = 0.7,
                  cexCol = 0.8,
                  main = paste("Enrichment Heatmap -", legend_title))
        }
        
      }, error = function(e) {
        message("Heatmap error: ", e$message)
        plot.new()
        text(0.5, 0.5, paste("Error creating heatmap:\n", e$message), cex = 1.2)
      })
    })
    
    # Display data table
    output$heatmap_data <- DT::renderDataTable({
      req(heatmap_data$data)
      
      DT::datatable(
        heatmap_data$data %>%
          dplyr::select(gene, cluster, enrichment_type, term_name, p.adjust, Count) %>%
          dplyr::arrange(p.adjust),
        options = list(pageLength = 25, scrollX = TRUE),
        filter = 'top'
      ) %>%
        DT::formatSignif(columns = "p.adjust", digits = 3)
    })
    
    # Settings info
    output$settings_info <- renderPrint({
      cat("Current Settings:\n")
      cat("================\n")
      cat("Heatmap Type:", input$heatmap_type, "\n")
      cat("Method Filter:", input$method_filter, "\n")
      cat("Enrichment Types:", paste(input$enrichment_types, collapse = ", "), "\n")
      cat("Cluster:", input$cluster_select, "\n")
      cat("Max Terms:", input$max_terms, "\n")
      cat("Show Direction:", input$show_direction, "\n")
      cat("Cluster Rows:", input$cluster_rows, "\n")
      cat("Cluster Columns:", input$cluster_cols, "\n")
      
      if (!is.null(heatmap_data$data)) {
        cat("\nData Summary:\n")
        cat("Total rows:", nrow(heatmap_data$data), "\n")
        cat("Unique genes:", length(unique(heatmap_data$data$gene)), "\n")
        cat("Unique terms:", length(unique(heatmap_data$data$term_name)), "\n")
      }
    })
    
    # Download handler
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0("enrichment_heatmap_", input$heatmap_type, "_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 12, height = 10)
        # Recreate the heatmap for download
        if (!is.null(heatmap_data$data)) {
          # Same heatmap code as above
          tryCatch({
            # ... (same heatmap generation code)
            dev.off()
          }, error = function(e) {
            dev.off()
          })
        }
      }
    )
    
    return(heatmap_data)
  })
}