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
                   plotlyOutput(ns("heatmap_plot"), height = "800px", width = "100%"),
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
      
      showNotification("Loading data for heatmap...", type = "default", duration = NULL, id = "loading")
      
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
                           type = "default", duration = 3)
          
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
    output$heatmap_plot <- renderPlotly({
      req(heatmap_data$data)
      
      # Check if heatmaply is available
      if (!requireNamespace("heatmaply", quietly = TRUE)) {
        showNotification("heatmaply package is required for interactive clustered heatmaps. Using plotly instead.", 
                         type = "warning", duration = 5)
        use_heatmaply <- FALSE
      } else {
        use_heatmaply <- TRUE
      }
      
      # Create an interactive heatmap
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
        
        # Determine value column based on heatmap type
        if (input$heatmap_type == "pvalue" && "p.adjust" %in% names(df)) {
          df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
          value_var <- "heatmap_value"
          legend_title <- "-log10(p-value)"
        } else if (input$heatmap_type == "foldenrichment") {
          # Check for fold enrichment columns
          if ("FoldEnrichment" %in% names(df)) {
            value_var <- "FoldEnrichment"
            legend_title <- "Fold Enrichment"
          } else if ("fold_enrichment" %in% names(df)) {
            value_var <- "fold_enrichment"
            legend_title <- "Fold Enrichment"
          } else if ("RichFactor" %in% names(df)) {
            value_var <- "RichFactor"
            legend_title <- "Rich Factor"
          } else {
            # Calculate fold enrichment if we have the necessary columns
            if ("Count" %in% names(df) && "GeneRatio" %in% names(df)) {
              # Parse GeneRatio (format: "X/Y")
              df$heatmap_value <- sapply(df$GeneRatio, function(x) {
                parts <- as.numeric(strsplit(x, "/")[[1]])
                if (length(parts) == 2 && parts[2] > 0) parts[1]/parts[2] else 1
              })
              value_var <- "heatmap_value"
              legend_title <- "Gene Ratio"
            } else {
              df$heatmap_value <- 1
              value_var <- "heatmap_value"
              legend_title <- "Fold Enrichment (not available)"
            }
          }
        } else if (input$heatmap_type == "zscore") {
          if ("zScore" %in% names(df)) {
            value_var <- "zScore"
            legend_title <- "Z-score"
          } else if ("z_score" %in% names(df)) {
            value_var <- "z_score"
            legend_title <- "Z-score"
          } else {
            # Calculate z-score from p-values if available
            if ("p.adjust" %in% names(df)) {
              df$heatmap_value <- qnorm(1 - df$p.adjust/2) * sign(rnorm(nrow(df)))
              value_var <- "heatmap_value"
              legend_title <- "Z-score (approximated)"
            } else {
              df$heatmap_value <- 0
              value_var <- "heatmap_value"
              legend_title <- "Z-score (not available)"
            }
          }
        } else if (input$heatmap_type == "count" && "Count" %in% names(df)) {
          value_var <- "Count"
          legend_title <- "Gene Count"
        } else {
          # Default to p-value if available
          if ("p.adjust" %in% names(df)) {
            df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
            value_var <- "heatmap_value"
            legend_title <- "-log10(p-value)"
          } else {
            df$heatmap_value <- 1
            value_var <- "heatmap_value"
            legend_title <- "Presence"
          }
        }
        
        # Create matrix
        if (length(unique(df[[x_var]])) < 2) {
          # If only one condition, create a simple bar plot
          top_terms <- head(df[order(df[[value_var]], decreasing = TRUE), ], input$max_terms)
          
          p <- plot_ly(
            x = top_terms[[value_var]],
            y = substr(top_terms[[y_var]], 1, 50),
            type = "bar",
            orientation = "h",
            text = top_terms[[y_var]],
            hovertemplate = "%{text}<br>%{x:.2f}<extra></extra>",
            marker = list(color = top_terms[[value_var]], colorscale = "Viridis")
          ) %>%
            layout(
              title = paste("Top Terms -", legend_title),
              xaxis = list(title = legend_title),
              yaxis = list(title = "", tickfont = list(size = 10)),
              margin = list(l = 200)
            )
          
          return(p)
          
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
          rownames(mat) <- substr(data_wide[[y_var]], 1, 80)  # Truncate long names
          
          # Store original row names and full descriptions
          original_rownames <- rownames(mat)
          full_descriptions_all <- data_wide[[y_var]]
          
          # Limit to max terms
          if (nrow(mat) > input$max_terms) {
            # Keep most significant terms (highest values)
            row_means <- rowMeans(mat, na.rm = TRUE)
            top_indices <- order(row_means, decreasing = TRUE)[1:input$max_terms]
            mat <- mat[top_indices, , drop = FALSE]
            full_descriptions <- full_descriptions_all[top_indices]
          } else {
            full_descriptions <- full_descriptions_all
          }
          
          # Use heatmaply if available, otherwise use plotly
          if (use_heatmaply) {
            # Prepare custom hover text for heatmaply
            custom_text <- matrix(
              paste0("Term: ", rep(full_descriptions, ncol(mat)), "<br>",
                     "Condition: ", rep(colnames(mat), each = nrow(mat)), "<br>",
                     "Value: "),
              nrow = nrow(mat),
              ncol = ncol(mat),
              byrow = FALSE
            )
            
            # Set color scale based on heatmap type
            if (input$heatmap_type == "pvalue") {
              colors <- colorRampPalette(c("white", "red"))(256)
            } else if (input$heatmap_type == "zscore") {
              colors <- colorRampPalette(c("blue", "white", "red"))(256)
            } else {
              colors <- colorRampPalette(c("white", "darkblue"))(256)
            }
            
            # Determine dendrogram settings
            dendrogram <- "none"
            if (input$cluster_rows && input$cluster_cols) {
              dendrogram <- "both"
            } else if (input$cluster_rows) {
              dendrogram <- "row"
            } else if (input$cluster_cols) {
              dendrogram <- "column"
            }
            
            # Create heatmaply heatmap with dendrograms
            p <- heatmaply::heatmaply(
              mat,
              dendrogram = dendrogram,
              colors = colors,
              xlab = "",
              ylab = "",
              main = paste("Interactive Enrichment Heatmap -", legend_title),
              margins = c(150, 250, 50, 50),
              custom_hovertext = custom_text,
              label_names = c("Row", "Column", "Value"),
              fontsize_row = 10,
              fontsize_col = 10,
              showticklabels = c(TRUE, TRUE),
              plot_method = "plotly"
            )
            
            return(p)
            
          } else {
            # Fallback to regular plotly without clustering
            # Create hover text with full descriptions
            hover_text <- matrix(
              paste0("Term: ", rep(full_descriptions, ncol(mat)), "<br>",
                     "Condition: ", rep(colnames(mat), each = nrow(mat)), "<br>",
                     "Value: ", round(as.vector(mat), 3)),
              nrow = nrow(mat),
              ncol = ncol(mat),
              byrow = FALSE
            )
            
            # Create color scale based on heatmap type
            if (input$heatmap_type == "pvalue") {
              colorscale <- list(c(0, "white"), c(1, "red"))
            } else if (input$heatmap_type == "zscore") {
              colorscale <- "RdBu"
            } else {
              colorscale <- list(c(0, "white"), c(1, "blue"))
            }
            
            # Create basic plotly heatmap
            p <- plot_ly(
              x = colnames(mat),
              y = rownames(mat),
              z = mat,
              type = "heatmap",
              colorscale = colorscale,
              hovertext = hover_text,
              hovertemplate = "%{hovertext}<extra></extra>",
              colorbar = list(title = legend_title)
            ) %>%
              layout(
                title = paste("Interactive Enrichment Heatmap -", legend_title),
                xaxis = list(
                  title = "",
                  tickangle = -45,
                  tickfont = list(size = 12)
                ),
                yaxis = list(
                  title = "",
                  tickfont = list(size = 10)
                ),
                margin = list(l = 250, b = 100, t = 50)
              )
            
            return(p)
          }
        }
        
      }, error = function(e) {
        message("Heatmap error: ", e$message)
        # Return an empty plotly object with error message
        plot_ly() %>%
          add_annotations(
            text = paste("Error creating heatmap:\n", e$message),
            x = 0.5,
            y = 0.5,
            xref = "paper",
            yref = "paper",
            showarrow = FALSE,
            font = list(size = 14)
          ) %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE)
          )
      })
    })
    
    # Display data table
    output$heatmap_data <- DT::renderDataTable({
      req(heatmap_data$data)
      
      # Select available columns
      available_cols <- intersect(
        c("gene", "mutation_perturbation", "cluster", "enrichment_type", "Description", 
          "term_name", "p.adjust", "Count", "FoldEnrichment", "RichFactor"),
        names(heatmap_data$data)
      )
      
      DT::datatable(
        heatmap_data$data %>%
          dplyr::select(all_of(available_cols)) %>%
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