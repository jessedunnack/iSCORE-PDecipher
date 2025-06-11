# Module: Pre-computed Results Browser
# Allows users to browse and visualize existing enrichment analysis results

# UI function
mod_precomputed_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      # Visualization panel (now full width)
      div(class = "box box-primary",
        div(class = "box-header with-border",
            h3("Visualization", class = "box-title")
        ),
        div(class = "box-body",
        
        # Tab panels for different views
        tabsetPanel(
          id = ns("viz_tabs"),
          
          # Main plot tab
          tabPanel(
            "Enrichment Plot",
            value = "plot",
            br(),
            fluidRow(
              column(3,
                wellPanel(
                  h4("Data Summary"),
                  
                  # Data source indicator
                  uiOutput(ns("data_source_indicator")),
                  
                  br(),
                  verbatimTextOutput(ns("data_summary")),
                  
                  hr(),
                  
                  h4("Plot Settings"),
                  selectInput(ns("plot_type"),
                              "Plot Type:",
                              choices = c("Dot Plot" = "dotplot",
                                          "Bar Chart" = "bar",
                                          "Lollipop" = "lollipop",
                                          "Bubble Chart" = "bubble",
                                          "TreeMap" = "treemap",
                                          "Network Plot" = "network"),
                              selected = "dotplot"),
                  
                  numericInput(ns("n_terms"),
                               "Number of Terms:",
                               value = 20,
                               min = 5,
                               max = 50,
                               step = 5),
                  
                  selectInput(ns("x_axis"),
                              "X-axis:",
                              choices = c("Gene Ratio" = "GeneRatio",
                                          "Count" = "Count",
                                          "-log10(p-value)" = "neg_log_p",
                                          "Fold Enrichment" = "FoldEnrichment"),
                              selected = "GeneRatio"),
                  
                  selectInput(ns("color_by"),
                              "Color By:",
                              choices = c("Adjusted P-value" = "p.adjust",
                                          "P-value" = "pvalue",
                                          "Count" = "Count"),
                              selected = "p.adjust"),
                  
                  selectInput(ns("order_by"),
                              "Order By:",
                              choices = c("P-value" = "pvalue",
                                          "Adjusted P-value" = "p.adjust",
                                          "Count" = "Count",
                                          "Gene Ratio" = "GeneRatio"),
                              selected = "p.adjust"),
                  
                  checkboxInput(ns("show_labels"),
                                "Show Gene Labels",
                                value = FALSE)
                )
              ),
              column(9,
                withSpinner(
                  plotOutput(ns("main_plot"), height = "600px"),
                  type = 4,
                  color = "#3c8dbc"
                ),
                br(),
                downloadButton(ns("download_plot"), "Download Plot")
              )
            )
          ),
          
          # Interactive heatmap tab
          tabPanel(
            "Interactive Heatmap",
            value = "heatmap",
            withSpinner(
              plotOutput(ns("heatmap_plot"), height = "600px"),
              type = 4,
              color = "#3c8dbc"
            ),
            hr(),
            fluidRow(
              column(4,
                selectInput(
                  ns("heatmap_scope"),
                  "Heatmap Scope",
                  choices = APP_CONFIG$heatmap_scope,
                  selected = "across_clusters"
                )
              ),
              column(4,
                selectInput(
                  ns("heatmap_metric"),
                  "Heatmap Metric",
                  choices = APP_CONFIG$heatmap_metrics,
                  selected = "p.adjust"
                )
              ),
              column(4,
                numericInput(
                  ns("max_terms"),
                  "Max Terms to Display",
                  value = 30,
                  min = 10,
                  max = 100,
                  step = 5
                )
              )
            )
          ),
          
          # Data table tab
          tabPanel(
            "Data Table",
            value = "table",
            br(),
            DT::dataTableOutput(ns("enrichment_table"))
          ),
          
          # Summary statistics tab
          tabPanel(
            "Statistics",
            value = "stats",
            br(),
            fluidRow(
              valueBoxOutput(ns("total_terms_box")),
              valueBoxOutput(ns("sig_terms_box")),
              valueBoxOutput(ns("top_pval_box"))
            ),
            br(),
            plotOutput(ns("stats_plot"), height = "400px")
          )
        )
        )
      )
    ),
    
    # Additional options row
    fluidRow(
      div(class = "box box-primary",
        div(class = "box-header with-border",
            h3("Additional Options", class = "box-title")
        ),
        div(class = "box-body",
        
        fluidRow(
          column(3,
            checkboxInput(
              ns("show_labels"),
              "Show Gene Labels",
              value = FALSE
            )
          ),
          column(3,
            checkboxInput(
              ns("use_fdr"),
              "Use FDR Correction",
              value = TRUE
            )
          ),
          column(3,
            numericInput(
              ns("min_genes"),
              "Min Genes per Term",
              value = 3,
              min = 1,
              max = 20
            )
          ),
          column(3,
            selectInput(
              ns("color_scheme"),
              "Color Scheme",
              choices = c("Default", "Viridis", "RdBu", "YlOrRd"),
              selected = "Default"
            )
          )
        )
        )
      )
    )
  )
}

# Server function
mod_precomputed_server <- function(id, app_data, global_selection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for this module
    module_data <- reactiveValues(
      current_result = NULL,
      file_path = NULL,
      consolidated_terms = NULL,
      is_gsea = FALSE
    )
    
    # Reactive data loading when global selection changes
    filtered_terms <- reactive({
      req(global_selection())
      selection <- global_selection()
      req(selection$gene, selection$cluster, selection$experiment, 
          selection$enrichment_type, selection$direction)
      
      # First, get significant terms from consolidated data
      message("Checking consolidated data for significant terms...")
      significant_terms <- get_significant_terms_from_consolidated(
        app_data$consolidated_data,
        analysis_type = selection$analysis_type,
        gene = selection$gene,
        cluster = selection$cluster,
        experiment = selection$experiment,
        enrichment_type = selection$enrichment_type,
        direction = selection$direction,
        pval_threshold = selection$pval_threshold
      )
      
      # Store the filtered consolidated data
      module_data$consolidated_terms <- significant_terms
      
      if (!is.null(significant_terms) && nrow(significant_terms) > 0) {
        message("Found ", nrow(significant_terms), " significant terms in consolidated data")
        showNotification(paste("Found", nrow(significant_terms), "significant terms"), 
                        type = "success", duration = 3)
      }
      
      # Only load the full enrichment result file if needed for enrichPlot visualization
      # This will be done when actually creating the plot
      if (selection$enrichment_type == "GSEA") {
        # GSEA has special directory structure
        showNotification("GSEA results have nested structure - using consolidated data", 
                        type = "info")
        module_data$current_result <- NULL  # Will load on demand
        module_data$is_gsea <- TRUE
      } else {
        # Store file path for later loading if needed
        file_name <- paste0(selection$enrichment_type, "_", selection$direction, ".rds")
        file_path <- file.path(
          APP_CONFIG$enrichment_results_path,
          selection$analysis_type,
          selection$gene,
          selection$cluster,
          selection$experiment,
          selection$enrichment_type,
          file_name
        )
        
        module_data$file_path <- file_path
        module_data$is_gsea <- FALSE
        
        # For now, we'll only load the full result when needed for plotting
        module_data$current_result <- NULL
        
        showNotification("Data selection updated. Click 'Update Plot' to visualize.", 
                        type = "message", duration = 3)
      }
      
      # Update app_data with the consolidated terms
      app_data$current_enrichment <- significant_terms
    })
    
    # Data summary
    output$data_summary <- renderPrint({
      # Use consolidated terms for summary
      if (!is.null(module_data$consolidated_terms) && nrow(module_data$consolidated_terms) > 0) {
        selection <- global_selection()
        terms_data <- module_data$consolidated_terms
        
        cat("Data Source: Consolidated\n")
        cat("Total significant terms:", nrow(terms_data), "\n")
        
        # Find top term by p-value
        pval_col <- if ("p.adjust" %in% names(terms_data)) "p.adjust" 
                    else if ("fdr" %in% names(terms_data)) "fdr" 
                    else if ("padj" %in% names(terms_data)) "padj"
                    else NULL
        
        if (!is.null(pval_col)) {
          top_idx <- which.min(terms_data[[pval_col]])
          if (length(top_idx) > 0) {
            desc_col <- if ("Description" %in% names(terms_data)) "Description" 
                       else if ("description" %in% names(terms_data)) "description"
                       else if ("term" %in% names(terms_data)) "term"
                       else NULL
            
            if (!is.null(desc_col)) {
              cat("Top term:", substr(terms_data[[desc_col]][top_idx], 1, 30), "...\n")
              cat("Top p-value:", format_pvalue(terms_data[[pval_col]][top_idx]), "\n")
            }
          }
        }
        
        cat("\nP-value threshold:", selection$pval_threshold)
      } else if (!is.null(global_selection())) {
        selection <- global_selection()
        cat("No significant terms found\n")
        cat("\nCurrent Selection:")
        cat("\nAnalysis:", selection$analysis_type)
        cat("\nGene:", selection$gene)
        cat("\nCluster:", selection$cluster)
        cat("\nEnrichment:", selection$enrichment_type)
        cat("\nDirection:", selection$direction)
        cat("\nP-value threshold:", selection$pval_threshold)
      } else {
        cat("No data loaded")
      }
    })
    
    # Main enrichment plot
    # Reactive value for plot data
    plot_data <- reactiveVal(NULL)
    
    # Update plot when button is clicked
    observeEvent(input$update_plot, {
      req(global_selection())
      req(module_data$file_path)  # Ensure we have a file path
      
      # Wrap entire process in error handler
      tryCatch({
        selection <- global_selection()
        message("=== Starting plot data preparation ===")
        message("Enrichment type: ", selection$enrichment_type)
        
        # Load the full enrichment result for plotting if not already loaded
        if (is.null(module_data$current_result)) {
          message("Loading full enrichment result from: ", module_data$file_path)
          
          withProgress(message = 'Loading data for plotting...', value = 0, {
            incProgress(0.5)
            module_data$current_result <- load_enrichment_safe(module_data$file_path)
            incProgress(0.5)
          })
          
          if (is.null(module_data$current_result)) {
            showNotification("Failed to load enrichment data for plotting", type = "error")
            return()
          }
        }
        
        # Initialize df
        df <- NULL
        
        # Prepare data based on enrichment type
        if (selection$enrichment_type == "STRING") {
          df <- module_data$current_result$enrichment
          message("STRING data - Total terms before filtering: ", nrow(df))
          message("P-value threshold: ", selection$pval_threshold)
          message("FDR values range: ", min(df$fdr, na.rm = TRUE), " to ", max(df$fdr, na.rm = TRUE))
          df <- df[df$fdr <= selection$pval_threshold, ]
          message("Terms after filtering: ", nrow(df))
          # Add calculated columns
          df$neg_log_p <- -log10(df$p_value)
          df$GeneRatio <- df$number_of_genes / df$number_of_genes_in_background
          # Rename columns for consistency
          names(df)[names(df) == "fdr"] <- "p.adjust"
          names(df)[names(df) == "p_value"] <- "pvalue"
          names(df)[names(df) == "description"] <- "Description"
          names(df)[names(df) == "number_of_genes"] <- "Count"
        } else if (methods::is(module_data$current_result, "enrichResult")) {
          df <- module_data$current_result@result
          message("enrichResult data - Total terms before filtering: ", nrow(df))
          message("P-value threshold: ", selection$pval_threshold)
          
          # Ensure p.adjust column exists and is numeric
          if (!"p.adjust" %in% names(df)) {
            message("Warning: p.adjust column missing from enrichResult")
            df$p.adjust <- 1  # Default high p-value
          } else {
            df$p.adjust <- as.numeric(df$p.adjust)
            df$p.adjust[is.na(df$p.adjust)] <- 1  # Replace NA with high p-value
            message("p.adjust values range: ", min(df$p.adjust, na.rm = TRUE), " to ", max(df$p.adjust, na.rm = TRUE))
          }
          
          # Safe filtering
          threshold <- selection$pval_threshold
          if (is.null(threshold) || is.na(threshold) || !is.numeric(threshold)) {
            threshold <- 0.05  # Default
          }
          
          df <- df[df$p.adjust <= threshold, ]
          message("Terms after filtering: ", nrow(df))
          
          # Check if we have any data left after filtering
          if (nrow(df) > 0) {
            # Parse GeneRatio safely
            if ("GeneRatio" %in% names(df)) {
              df$GeneRatio <- sapply(strsplit(as.character(df$GeneRatio), "/"), function(x) {
                if (length(x) == 2 && !is.na(x[1]) && !is.na(x[2])) {
                  as.numeric(x[1]) / as.numeric(x[2])
                } else {
                  0.1  # Default small value
                }
              })
            } else {
              df$GeneRatio <- 0.1  # Default if column missing
            }
            
            # Create neg_log_p safely
            if ("pvalue" %in% names(df)) {
              pvals <- df$pvalue
              pvals[is.na(pvals)] <- 1  # Replace NA with 1
              pvals[pvals <= 0] <- 1e-300  # Replace 0 or negative with small value
              df$neg_log_p <- -log10(pvals)
            } else {
              df$neg_log_p <- 1  # Default
            }
            
            # Add Fold Enrichment if available
            if ("qvalue" %in% names(df) && "BgRatio" %in% names(df)) {
              bg_ratio <- sapply(strsplit(as.character(df$BgRatio), "/"), function(x) {
                if (length(x) == 2 && !is.na(x[1]) && !is.na(x[2])) {
                  result <- as.numeric(x[1]) / as.numeric(x[2])
                  if (is.na(result) || !is.finite(result) || result <= 0) {
                    return(0.1)  # Default small value
                  }
                  return(result)
                } else {
                  return(0.1)  # Default small value
                }
              })
              # Safe division
              fold_enrich <- df$GeneRatio / bg_ratio
              fold_enrich[is.na(fold_enrich) | !is.finite(fold_enrich)] <- 1
              df$FoldEnrichment <- fold_enrich
            } else {
              df$FoldEnrichment <- df$GeneRatio * 10  # Rough approximation
            }
          }
        }
        
        message("=== Data processing complete ===")
        message("df is null: ", is.null(df))
        message("df is data.frame: ", is.data.frame(df))
        if (!is.null(df)) {
          message("df nrow: ", nrow(df))
          message("df class: ", class(df))
          message("df columns: ", paste(names(df), collapse = ", "))
        }
      
      # Only process if we have data
      if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
        # Ensure all required columns exist with safe defaults
        required_cols <- c("Description", "Count", "p.adjust", "pvalue")
        for (col in required_cols) {
          if (!col %in% names(df)) {
            if (col == "Description") df$Description <- paste("Term", 1:nrow(df))
            if (col == "Count") df$Count <- 1
            if (col == "p.adjust") df$p.adjust <- 0.05
            if (col == "pvalue") df$pvalue <- 0.05
          }
        }
        
        # Safe column ordering
        order_col <- input$order_by
        if (is.null(order_col) || is.na(order_col) || !order_col %in% names(df)) {
          # Use p.adjust as default ordering column
          order_col <- "p.adjust"
          if (!order_col %in% names(df)) {
            order_col <- names(df)[1]  # Use first column as last resort
          }
        }
        
        # Safe ordering with error handling
        tryCatch({
          if (order_col %in% names(df)) {
            order_values <- df[[order_col]]
            # Handle non-numeric values
            if (!is.numeric(order_values)) {
              order_values <- as.numeric(as.factor(order_values))
            }
            # Handle NA values
            order_values[is.na(order_values)] <- max(order_values, na.rm = TRUE) + 1
            df <- df[order(order_values), ]
          }
        }, error = function(e) {
          message("Error ordering data: ", e$message)
          # Keep original order if sorting fails
        })
        
        # Take top N terms
        n_terms <- input$n_terms
        if (is.null(n_terms) || is.na(n_terms) || n_terms <= 0) {
          n_terms <- 20  # Default
        }
        df <- head(df, min(n_terms, nrow(df)))
        
        plot_data(df)
      } else {
        # No significant results
        plot_data(NULL)
        showNotification("No significant terms found with current p-value threshold", 
                        type = "warning", duration = 5)
      }
      
      }, error = function(e) {
        message("=== ERROR in plot data preparation ===")
        message("Error: ", e$message)
        message("Call: ", paste(deparse(e$call), collapse = " "))
        plot_data(NULL)
        showNotification(paste("Error preparing plot data:", e$message), 
                        type = "error", duration = 10)
      })
    }, ignoreNULL = FALSE)
    
    output$main_plot <- renderPlot({
      req(plot_data())
      df <- plot_data()
      
      # Check if we have data
      if (is.null(df) || nrow(df) == 0) {
        plot.new()
        text(0.5, 0.5, "No significant terms to display", cex = 1.5, col = "gray50")
        return()
      }
      
      # Create plot based on type
      if (input$plot_type == "dotplot") {
        p <- ggplot(df, aes_string(x = input$x_axis, 
                                   y = paste0("reorder(Description, ", input$x_axis, ")"))) +
          geom_point(aes_string(size = "Count", color = input$color_by)) +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_continuous(range = c(3, 10)) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10),
                axis.text.x = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold")) +
          labs(x = get_axis_label(input$x_axis),
               y = "",
               size = "Gene Count",
               color = get_axis_label(input$color_by),
               title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
        
      } else if (input$plot_type == "bar") {
        p <- ggplot(df, aes_string(x = paste0("reorder(Description, ", input$x_axis, ")"), 
                                   y = input$x_axis)) +
          geom_bar(stat = "identity", aes_string(fill = input$color_by)) +
          scale_fill_gradient(low = "lightblue", high = "darkred") +
          coord_flip() +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold")) +
          labs(x = "",
               y = get_axis_label(input$x_axis),
               fill = get_axis_label(input$color_by),
               title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
        
      } else if (input$plot_type == "lollipop") {
        p <- ggplot(df, aes_string(x = input$x_axis, 
                                   y = paste0("reorder(Description, ", input$x_axis, ")"))) +
          geom_segment(aes_string(x = 0, xend = input$x_axis, 
                                 y = paste0("reorder(Description, ", input$x_axis, ")"),
                                 yend = paste0("reorder(Description, ", input$x_axis, ")")),
                      color = "grey50") +
          geom_point(aes_string(size = "Count", color = input$color_by)) +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_continuous(range = c(3, 10)) +
          theme_minimal() +
          theme(axis.text.y = element_text(size = 10),
                plot.title = element_text(size = 14, face = "bold")) +
          labs(x = get_axis_label(input$x_axis),
               y = "",
               size = "Gene Count",
               color = get_axis_label(input$color_by),
               title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
               
      } else if (input$plot_type == "bubble") {
        # Enhanced bubble chart with size AND color mapping
        p <- ggplot(df, aes_string(x = input$x_axis, 
                                   y = "Count")) +
          geom_point(aes_string(size = "Count", color = input$color_by), alpha = 0.7) +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_continuous(range = c(5, 20)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 14, face = "bold")) +
          labs(x = get_axis_label(input$x_axis),
               y = "Gene Count",
               size = "Gene Count",
               color = get_axis_label(input$color_by),
               title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
               
      } else if (input$plot_type == "treemap") {
        # TreeMap visualization (simplified version)
        if (requireNamespace("treemap", quietly = TRUE)) {
          # Use treemap if available
          p <- ggplot(df, aes_string(area = "Count", fill = input$color_by)) +
            geom_rect(aes(xmin = 0, xmax = Count, ymin = 0, ymax = 1), 
                     color = "white", size = 0.5) +
            scale_fill_gradient(low = "lightblue", high = "darkred") +
            theme_void() +
            labs(title = paste("TreeMap -", input$enrichment_type, "-", input$direction, "regulated genes"))
        } else {
          # Fallback to tile plot
          df$tile_x <- rep(1:ceiling(sqrt(nrow(df))), length.out = nrow(df))
          df$tile_y <- rep(1:ceiling(nrow(df)/ceiling(sqrt(nrow(df)))), each = ceiling(sqrt(nrow(df))))[1:nrow(df)]
          
          p <- ggplot(df, aes(x = tile_x, y = tile_y)) +
            geom_tile(aes_string(fill = input$color_by), color = "white", size = 0.5) +
            scale_fill_gradient(low = "lightblue", high = "darkred") +
            theme_void() +
            labs(title = paste("Tile Map -", input$enrichment_type, "-", input$direction, "regulated genes"))
        }
        
      } else if (input$plot_type == "network") {
        # Simple network-style plot
        if (requireNamespace("igraph", quietly = TRUE) && requireNamespace("ggraph", quietly = TRUE)) {
          # Network plot if packages available
          showNotification("Network plot requires additional packages", type = "info")
          p <- ggplot(df, aes_string(x = input$x_axis, 
                                     y = paste0("reorder(Description, ", input$x_axis, ")"))) +
            geom_point(aes_string(size = "Count", color = input$color_by)) +
            scale_color_gradient(low = "blue", high = "red") +
            theme_minimal() +
            labs(title = "Network visualization coming soon...")
        } else {
          # Fallback to connected dot plot
          p <- ggplot(df, aes_string(x = input$x_axis, 
                                     y = paste0("reorder(Description, ", input$x_axis, ")"))) +
            geom_line(aes(group = 1), color = "gray70", alpha = 0.5) +
            geom_point(aes_string(size = "Count", color = input$color_by)) +
            scale_color_gradient(low = "blue", high = "red") +
            scale_size_continuous(range = c(3, 10)) +
            theme_minimal() +
            labs(x = get_axis_label(input$x_axis),
                 y = "",
                 title = paste("Connected Plot -", input$enrichment_type, "-", input$direction, "regulated genes"))
        }
      }
      
      # Add gene labels if requested
      if (input$show_labels && "geneID" %in% names(df)) {
        gene_labels <- sapply(strsplit(df$geneID, "/"), function(x) {
          if (length(x) > 3) {
            paste(c(x[1:3], "..."), collapse = "/")
          } else {
            paste(x, collapse = "/")
          }
        })
        p <- p + geom_text(aes(label = gene_labels), hjust = -0.1, size = 3)
      }
      
      print(p)
    })
    
    # Helper function for axis labels
    get_axis_label <- function(var) {
      labels <- c(
        "GeneRatio" = "Gene Ratio",
        "Count" = "Gene Count",
        "neg_log_p" = "-log10(p-value)",
        "FoldEnrichment" = "Fold Enrichment",
        "p.adjust" = "Adjusted P-value",
        "pvalue" = "P-value"
      )
      if (var %in% names(labels)) {
        return(labels[var])
      } else {
        return(var)
      }
    }
    
    # Download plot handler
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0(global_selection()$enrichment_type, "_", global_selection()$direction, "_enrichment_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        pdf(file, width = 10, height = 8)
        
        # Recreate the plot
        df <- plot_data()
        if (!is.null(df)) {
          # Create plot based on type
          if (input$plot_type == "dotplot") {
            p <- ggplot(df, aes_string(x = input$x_axis, 
                                       y = paste0("reorder(Description, ", input$x_axis, ")"))) +
              geom_point(aes_string(size = "Count", color = input$color_by)) +
              scale_color_gradient(low = "blue", high = "red") +
              scale_size_continuous(range = c(3, 10)) +
              theme_minimal() +
              theme(axis.text.y = element_text(size = 10),
                    axis.text.x = element_text(size = 10),
                    plot.title = element_text(size = 14, face = "bold")) +
              labs(x = get_axis_label(input$x_axis),
                   y = "",
                   size = "Gene Count",
                   color = get_axis_label(input$color_by),
                   title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
            
          } else if (input$plot_type == "bar") {
            p <- ggplot(df, aes_string(x = paste0("reorder(Description, ", input$x_axis, ")"), 
                                       y = input$x_axis)) +
              geom_bar(stat = "identity", aes_string(fill = input$color_by)) +
              scale_fill_gradient(low = "lightblue", high = "darkred") +
              coord_flip() +
              theme_minimal() +
              theme(axis.text.y = element_text(size = 10),
                    plot.title = element_text(size = 14, face = "bold")) +
              labs(x = "",
                   y = get_axis_label(input$x_axis),
                   fill = get_axis_label(input$color_by),
                   title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
            
          } else if (input$plot_type == "lollipop") {
            p <- ggplot(df, aes_string(x = input$x_axis, 
                                       y = paste0("reorder(Description, ", input$x_axis, ")"))) +
              geom_segment(aes_string(x = 0, xend = input$x_axis, 
                                     y = paste0("reorder(Description, ", input$x_axis, ")"),
                                     yend = paste0("reorder(Description, ", input$x_axis, ")")),
                          color = "grey50") +
              geom_point(aes_string(size = "Count", color = input$color_by)) +
              scale_color_gradient(low = "blue", high = "red") +
              scale_size_continuous(range = c(3, 10)) +
              theme_minimal() +
              theme(axis.text.y = element_text(size = 10),
                    plot.title = element_text(size = 14, face = "bold")) +
              labs(x = get_axis_label(input$x_axis),
                   y = "",
                   size = "Gene Count",
                   color = get_axis_label(input$color_by),
                   title = paste(input$enrichment_type, "-", input$direction, "regulated genes"))
          }
          
          # Add gene labels if requested
          if (input$show_labels && "geneID" %in% names(df)) {
            gene_labels <- sapply(strsplit(df$geneID, "/"), function(x) {
              if (length(x) > 3) {
                paste(c(x[1:3], "..."), collapse = "/")
              } else {
                paste(x, collapse = "/")
              }
            })
            p <- p + geom_text(aes(label = gene_labels), hjust = -0.1, size = 3)
          }
          
          print(p)
        }
        
        dev.off()
      }
    )
    
    # Heatmap plot
    output$heatmap_plot <- renderPlot({
      req(module_data$current_result)
      req(input$heatmap_scope)
      
      tryCatch({
        # Load data based on scope
        if (input$heatmap_scope == "across_clusters") {
          # Show terms across all clusters for current mutation
          showNotification("Loading data across all clusters...", type = "message", id = "heatmap_loading")
          
          # Get all clusters for current gene
          clusters <- get_available_clusters(input$gene)
          
          # Collect terms from all clusters
          term_data <- list()
          
          for (clust in clusters) {
            file_path <- file.path(
              APP_CONFIG$enrichment_results_path,
              input$analysis_type,
              input$gene,
              clust,
              input$experiment,
              input$enrichment_type,
              paste0(input$enrichment_type, "_", input$direction, ".rds")
            )
            
            if (file.exists(file_path)) {
              result <- readRDS(file_path)
              
              # Extract terms based on type
              if (input$enrichment_type == "STRING") {
                df <- result$enrichment
                df <- df[df$fdr <= pval_threshold(), ]
                if (nrow(df) > 0) {
                  term_data[[clust]] <- data.frame(
                    term = df$description,
                    pvalue = df$fdr,
                    count = df$number_of_genes,
                    stringsAsFactors = FALSE
                  )
                }
              } else if (methods::is(result, "enrichResult")) {
                df <- result@result
                df <- df[df$p.adjust <= pval_threshold(), ]
                if (nrow(df) > 0) {
                  term_data[[clust]] <- data.frame(
                    term = df$Description,
                    pvalue = df$p.adjust,
                    count = df$Count,
                    stringsAsFactors = FALSE
                  )
                }
              }
            }
          }
          
          removeNotification("heatmap_loading")
          
        } else {
          # Show terms across all mutations for current cluster
          showNotification("Loading data across all mutations...", type = "message", id = "heatmap_loading")
          
          # Get all genes
          genes <- get_available_genes()
          
          # Collect terms from all genes
          term_data <- list()
          
          for (g in genes) {
            # Check if this gene/cluster/experiment combination exists
            exp_list <- get_available_experiments(g, input$cluster, input$analysis_type)
            exp_to_use <- if (input$analysis_type == "MAST") "default" else if ("combined" %in% exp_list) "combined" else exp_list[1]
            
            if (!is.null(exp_to_use) && length(exp_to_use) > 0) {
              file_path <- file.path(
                APP_CONFIG$enrichment_results_path,
                input$analysis_type,
                g,
                input$cluster,
                exp_to_use,
                input$enrichment_type,
                paste0(input$enrichment_type, "_", input$direction, ".rds")
              )
              
              if (file.exists(file_path)) {
                result <- readRDS(file_path)
                
                # Extract terms based on type
                if (input$enrichment_type == "STRING") {
                  df <- result$enrichment
                  df <- df[df$fdr <= pval_threshold(), ]
                  if (nrow(df) > 0) {
                    term_data[[g]] <- data.frame(
                      term = df$description,
                      pvalue = df$fdr,
                      count = df$number_of_genes,
                      stringsAsFactors = FALSE
                    )
                  }
                } else if (methods::is(result, "enrichResult")) {
                  df <- result@result
                  df <- df[df$p.adjust <= pval_threshold(), ]
                  if (nrow(df) > 0) {
                    term_data[[g]] <- data.frame(
                      term = df$Description,
                      pvalue = df$p.adjust,
                      count = df$Count,
                      stringsAsFactors = FALSE
                    )
                  }
                }
              }
            }
          }
          
          removeNotification("heatmap_loading")
        }
        
        # Create matrix from collected data
        if (length(term_data) > 0) {
          # Get all unique terms
          all_terms <- unique(unlist(lapply(term_data, function(x) x$term)))
          
          # Limit to top terms by frequency
          term_freq <- table(unlist(lapply(term_data, function(x) x$term)))
          top_terms <- names(sort(term_freq, decreasing = TRUE))[1:min(input$max_terms, length(term_freq))]
          
          # Create matrix
          mat <- matrix(0, nrow = length(top_terms), ncol = length(term_data))
          rownames(mat) <- substr(top_terms, 1, 60)
          colnames(mat) <- names(term_data)
          
          # Fill matrix based on selected metric
          for (i in 1:length(top_terms)) {
            for (j in 1:length(term_data)) {
              term_idx <- which(term_data[[j]]$term == top_terms[i])
              if (length(term_idx) > 0) {
                if (input$heatmap_metric == "p.adjust") {
                  mat[i, j] <- -log10(term_data[[j]]$pvalue[term_idx[1]])
                } else if (input$heatmap_metric == "count") {
                  mat[i, j] <- term_data[[j]]$count[term_idx[1]]
                }
              }
            }
          }
          
          # Create heatmap
          heatmap(mat,
                  col = colorRampPalette(c("white", "red"))(100),
                  scale = "none",
                  margins = c(10, 15),
                  cexRow = 0.7,
                  cexCol = 0.8,
                  main = paste(input$enrichment_type, "-", 
                               ifelse(input$heatmap_scope == "across_clusters", 
                                      paste("Terms across clusters for", input$gene),
                                      paste("Terms across mutations in", input$cluster))))
          
        } else {
          plot.new()
          text(0.5, 0.5, "No significant terms found for heatmap", cex = 1.2)
        }
        
      }, error = function(e) {
        removeNotification("heatmap_loading")
        plot.new()
        text(0.5, 0.5, paste("Error creating heatmap:\n", e$message), cex = 1.2)
      })
    })
    
    # Data table - now uses consolidated data first
    output$enrichment_table <- DT::renderDataTable({
      # First try to use consolidated data
      if (!is.null(module_data$consolidated_terms) && nrow(module_data$consolidated_terms) > 0) {
        df <- module_data$consolidated_terms
        
        # Standardize column names based on enrichment type
        selection <- global_selection()
        
        # Identify columns and rename them consistently
        if (selection$enrichment_type == "STRING") {
          # STRING-specific columns
          col_mapping <- c(
            "term" = "ID",
            "description" = "Description",
            "number_of_genes" = "Gene Count",
            "Count" = "Gene Count",  # Alternative column name
            "geneID" = "Gene Symbols",
            "preferredNames" = "Gene Symbols",
            "fdr" = "Adj. P-value",
            "p_value" = "P-value",
            "pvalue" = "P-value"
          )
        } else {
          # Standard enrichment result columns
          col_mapping <- c(
            "ID" = "ID",
            "term_id" = "ID",
            "Description" = "Description",
            "term" = "Description",
            "Count" = "Gene Count",
            "number_of_genes" = "Gene Count",
            "geneID" = "Gene Symbols",
            "gene_list" = "Gene Symbols",
            "p.adjust" = "Adj. P-value",
            "padj" = "Adj. P-value",
            "fdr" = "Adj. P-value",
            "pvalue" = "P-value",
            "p_value" = "P-value"
          )
        }
        
        # Apply column mapping
        for (old_name in names(col_mapping)) {
          if (old_name %in% names(df) && !(col_mapping[old_name] %in% names(df))) {
            names(df)[names(df) == old_name] <- col_mapping[old_name]
          }
        }
        
        # Ensure required columns exist
        required_cols <- c("ID", "Description", "Gene Count", "Adj. P-value", "P-value")
        existing_cols <- intersect(required_cols, names(df))
        
        # Add Gene Symbols if available
        if ("Gene Symbols" %in% names(df)) {
          existing_cols <- c(existing_cols[1:3], "Gene Symbols", existing_cols[4:5])
          # Clean up gene symbols display
          df$`Gene Symbols` <- gsub("/", ", ", df$`Gene Symbols`)
        }
        
        # Select and order columns
        df <- df[, existing_cols, drop = FALSE]
        
        message("Displaying ", nrow(df), " significant terms from consolidated data")
        
      } else if (!is.null(module_data$current_result)) {
        # Fall back to using the full enrichment result if loaded
        selection <- global_selection()
        
        if (selection$enrichment_type == "STRING") {
          df <- module_data$current_result$enrichment
          df <- df[df$fdr <= selection$pval_threshold, ]
          # Add gene symbols column
          if ("preferredNames" %in% names(df)) {
            df$gene_symbols <- df$preferredNames
          } else {
            df$gene_symbols <- "N/A"
          }
          df <- df[, c("term", "description", "number_of_genes", "gene_symbols", "fdr", "p_value")]
          names(df) <- c("ID", "Description", "Gene Count", "Gene Symbols", "FDR", "P-value")
        } else if (methods::is(module_data$current_result, "enrichResult")) {
          df <- module_data$current_result@result
          df <- df[df$p.adjust <= selection$pval_threshold, ]
          
          # Add gene symbols handling (same as before)
          if ("geneID" %in% names(df)) {
            df$gene_symbols <- gsub("/", ", ", df$geneID)
          } else {
            df$gene_symbols <- "N/A"
          }
          
          df <- df[, c("ID", "Description", "Count", "gene_symbols", "p.adjust", "pvalue")]
          names(df) <- c("ID", "Description", "Gene Count", "Gene Symbols", "Adj. P-value", "P-value")
        }
      } else {
        # No data available
        return(NULL)
      }
      
      # Create the data table
      DT::datatable(
        df,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          order = list(list(ncol(df) - 1, 'asc')),  # Sort by adjusted p-value
          columnDefs = list(
            list(width = '300px', targets = which(names(df) == "Gene Symbols") - 1),  # Make gene symbols column wider
            list(className = 'dt-left', targets = '_all')
          )
        ),
        filter = 'top'
      ) %>%
        DT::formatSignif(columns = c("Adj. P-value", "P-value"), digits = 3)
    })
    
    # Value boxes for statistics - use consolidated data when available
    output$total_terms_box <- renderValueBox({
      total_terms <- 0
      
      # First check consolidated data
      if (!is.null(module_data$consolidated_terms)) {
        # For consolidated data, we only have significant terms
        # So we'll just show the count of significant terms as total
        total_terms <- nrow(module_data$consolidated_terms)
      } else if (!is.null(module_data$current_result)) {
        # Use full result if loaded
        total_terms <- summarize_enrichment(module_data$current_result, 
                                          global_selection()$enrichment_type)$total_terms
      }
      
      valueBox(
        value = total_terms,
        subtitle = "Total Terms",
        icon = icon("list"),
        color = "blue"
      )
    })
    
    output$sig_terms_box <- renderValueBox({
      sig_terms <- 0
      
      # First check consolidated data
      if (!is.null(module_data$consolidated_terms)) {
        sig_terms <- nrow(module_data$consolidated_terms)
      } else if (!is.null(module_data$current_result)) {
        # Use full result if loaded
        sig_terms <- summarize_enrichment(module_data$current_result, 
                                        global_selection()$enrichment_type,
                                        global_selection()$pval_threshold)$significant_terms
      }
      
      valueBox(
        value = sig_terms,
        subtitle = paste("Significant Terms (p <", global_selection()$pval_threshold, ")"),
        icon = icon("check"),
        color = "green"
      )
    })
    
    output$top_pval_box <- renderValueBox({
      top_pval <- 1
      
      # First check consolidated data
      if (!is.null(module_data$consolidated_terms) && nrow(module_data$consolidated_terms) > 0) {
        # Find the best p-value from consolidated data
        pval_cols <- c("p.adjust", "fdr", "padj")
        for (col in pval_cols) {
          if (col %in% names(module_data$consolidated_terms)) {
            top_pval <- min(module_data$consolidated_terms[[col]], na.rm = TRUE)
            break
          }
        }
      } else if (!is.null(module_data$current_result)) {
        # Use full result if loaded
        summary_stats <- summarize_enrichment(module_data$current_result, 
                                            global_selection()$enrichment_type)
        top_pval <- summary_stats$top_pvalue
      }
      
      valueBox(
        value = format_pvalue(top_pval),
        subtitle = "Best P-value",
        icon = icon("star"),
        color = "yellow"
      )
    })
    
    # Statistics plot
    output$stats_plot <- renderPlot({
      # Try to get p-values from either source
      pvals <- NULL
      
      if (!is.null(module_data$consolidated_terms) && nrow(module_data$consolidated_terms) > 0) {
        # Use consolidated data
        pval_cols <- c("p.adjust", "fdr", "padj")
        for (col in pval_cols) {
          if (col %in% names(module_data$consolidated_terms)) {
            pvals <- module_data$consolidated_terms[[col]]
            break
          }
        }
      } else if (!is.null(module_data$current_result)) {
        # Use full result if loaded
        selection <- global_selection()
        if (selection$enrichment_type == "STRING") {
          pvals <- module_data$current_result$enrichment$fdr
        } else if (methods::is(module_data$current_result, "enrichResult")) {
          pvals <- module_data$current_result@result$p.adjust
        }
      }
      
      if (!is.null(pvals)) {
        df <- data.frame(pvalue = pvals)
        
        ggplot(df, aes(x = pvalue)) +
          geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7) +
          geom_vline(xintercept = global_selection()$pval_threshold, color = "red", linetype = "dashed") +
          theme_minimal() +
          labs(x = "Adjusted P-value", y = "Count", 
               title = "P-value Distribution") +
          scale_x_continuous(limits = c(0, max(pvals, global_selection()$pval_threshold) * 1.1))
      } else {
        # No data available
        plot(1, type = "n", xlab = "", ylab = "", axes = FALSE)
        text(1, 1, "No data available for statistics plot", cex = 1.2)
      }
    })
    
    # Return the current data for use in other modules
    reactive({
      selection <- global_selection()
      list(
        data = module_data$current_result,
        consolidated_terms = module_data$consolidated_terms,
        file_path = module_data$file_path,
        metadata = list(
          analysis_type = selection$analysis_type,
          gene = selection$gene,
          cluster = selection$cluster,
          experiment = selection$experiment,
          enrichment_type = selection$enrichment_type,
          direction = selection$direction,
          pval_threshold = selection$pval_threshold
        )
      )
    })
  })
}