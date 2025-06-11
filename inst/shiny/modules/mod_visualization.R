# Module: Enhanced Visualization
# Advanced visualizations for enrichment results

mod_visualization_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Visualization Settings"),
      
      selectInput(ns("plot_type"),
                  "Visualization Type:",
                  choices = list(
                    "Dot Plot" = "dotplot"
                  ),
                  selected = "dotplot"),
      
      numericInput(ns("top_terms"),
                   "Number of Top Terms:",
                   value = 20,
                   min = 5,
                   max = 50,
                   step = 5),
      
      br(),
      
      selectInput(ns("x_axis"),
                  "X-axis:",
                  choices = c("-log10(p.adjust)" = "neg_log10_pval",
                            "Fold Enrichment" = "FoldEnrichment", 
                            "Gene Count" = "Count",
                            "Gene Ratio" = "GeneRatio",
                            "Rich Factor" = "RichFactor"),
                  selected = "neg_log10_pval"),
      
      selectInput(ns("color_by"),
                  "Color By:",
                  choices = c("-log10(adjusted p-value)" = "p-value", 
                             "Fold Enrichment" = "Fold Enrichment"),
                  selected = "p-value"),
      
      checkboxInput(ns("show_labels"),
                    "Show Labels",
                    value = TRUE),
      
      br(),
      helpText("Plot updates automatically when you change global settings"),
      
      br(),
      downloadButton(ns("download_plot"), "Download Plot",
                     class = "btn-success", style = "width: 100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = ns("main_tabs"),
        selected = "interactive",  # Make interactive the default tab
        
        tabPanel("Interactive Version",
                 value = "interactive",
                 br(),
                 plotlyOutput(ns("interactive_plot"), height = "700px") %>%
                   withSpinner(color = "#0073b7")),
        
        tabPanel("Static Plot",
                 value = "static",
                 br(),
                 plotOutput(ns("main_plot"), height = "700px") %>%
                   withSpinner(color = "#0073b7")),
        
        tabPanel("Plot Details",
                 value = "details",
                 br(),
                 verbatimTextOutput(ns("plot_info")),
                 br(),
                 DT::dataTableOutput(ns("plot_data")))
      )
    )
  )
}

mod_visualization_server <- function(id, global_selection, enrichment_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for plot data
    plot_data <- reactiveValues(
      data = NULL,
      plot = NULL,
      interactive = NULL
    )
    
    # Reactive data processing - automatically updates when global selection changes
    observe({
      req(global_selection(), enrichment_data())
      
      # Get filtered data from enrichment_data based on global selection
      selection <- global_selection()
      data <- enrichment_data()
      
      if (!is.null(data) && !is.null(selection)) {
        message("=== Visualization Data Filtering ===")
        message("Total data rows before filtering: ", nrow(data))
        message("Selection - Analysis type: ", selection$analysis_type)
        message("Selection - Gene: ", selection$gene)
        message("Selection - Cluster: ", selection$cluster)
        message("Selection - Experiment: ", selection$experiment)
        message("Selection - Enrichment type: ", selection$enrichment_type)
        message("Selection - Direction: ", selection$direction)
        
        # Check unique values in data
        if ("method" %in% names(data)) {
          message("Unique analysis types in data: ", paste(unique(data$method), collapse=", "))
        } else if ("analysis_type" %in% names(data)) {
          message("Unique analysis types in data: ", paste(unique(data$analysis_type), collapse=", "))
        }
        message("Unique genes in data: ", paste(head(unique(data$gene), 5), collapse=", "), "...")
        
        # Filter data based on selection
        # Handle both method and analysis_type column names
        if ("method" %in% names(data)) {
          analysis_filter <- data$method == selection$analysis_type
        } else {
          analysis_filter <- data$analysis_type == selection$analysis_type
        }
        
        # Apply all filters
        # Start with analysis type filter
        filtered_data <- data[analysis_filter, ]
        
        # Apply other filters only if they are not empty
        if (!is.null(selection$gene) && selection$gene != "") {
          filtered_data <- filtered_data[filtered_data$gene == selection$gene, ]
        }
        
        if (!is.null(selection$cluster) && selection$cluster != "") {
          filtered_data <- filtered_data[filtered_data$cluster == selection$cluster, ]
        }
        
        if (!is.null(selection$experiment) && selection$experiment != "") {
          filtered_data <- filtered_data[filtered_data$experiment == selection$experiment, ]
        }
        
        if (!is.null(selection$enrichment_type) && selection$enrichment_type != "") {
          filtered_data <- filtered_data[filtered_data$enrichment_type == selection$enrichment_type, ]
        }
        
        if (!is.null(selection$direction) && selection$direction != "ALL") {
          filtered_data <- filtered_data[filtered_data$direction == selection$direction, ]
        }
        
        data <- filtered_data
        
        message("Filtered data: ", nrow(data), " rows for visualization")
        message("=================================")
      } else {
        data <- NULL
        message("Data or selection is NULL")
      }
      
      # Process based on plot type - default to dotplot
      plot_type <- if(is.null(input$plot_type)) "dotplot" else input$plot_type
      top_terms <- if(is.null(input$top_terms)) 20 else input$top_terms
      
      if (plot_type == "dotplot") {
        plot_data$data <- prepare_dotplot_data(data, top_terms)
        x_axis_val <- if(is.null(input$x_axis)) "neg_log10_pval" else input$x_axis
        plot_data$plot <- create_dotplot(plot_data$data, input$color_by, input$show_labels, x_axis_val, selection)
      } else {
        # For other plot types, create a simple dotplot as fallback
        plot_data$data <- prepare_dotplot_data(data, top_terms)
        plot_data$plot <- create_dotplot(plot_data$data, input$color_by, input$show_labels, "neg_log10_pval", selection)
      }
      
      message("Plot updated for ", plot_type, " with ", if(!is.null(plot_data$data)) nrow(plot_data$data) else 0, " terms")
    })
    
    # Helper function to load gene lists from individual RDS files
    load_gene_lists <- function(data) {
      if (is.null(data) || nrow(data) == 0) return(data)
      
      # Add gene lists from individual RDS files
      data$gene_list <- vector("list", nrow(data))
      data$gene_symbols <- character(nrow(data))
      
      for (i in 1:nrow(data)) {
        tryCatch({
          # Construct file path for this specific result
          file_path <- file.path(
            "../enrichment_results",
            data$analysis_type[i],
            data$gene[i], 
            data$cluster[i],
            data$experiment[i],
            data$enrichment_type[i],
            paste0(data$enrichment_type[i], "_", data$direction[i], ".rds")
          )
          
          if (file.exists(file_path)) {
            enrichment_obj <- readRDS(file_path)
            
            # Extract genes for this specific term
            if (class(enrichment_obj)[1] == "enrichResult") {
              result_df <- enrichment_obj@result
              term_row <- result_df[result_df$ID == data$ID[i] | result_df$Description == data$Description[i], ]
              if (nrow(term_row) > 0 && !is.na(term_row$geneID[1])) {
                genes <- unlist(strsplit(term_row$geneID[1], "/"))
                data$gene_list[[i]] <- genes
                # Limit to first 10 genes for display
                display_genes <- if (length(genes) > 10) {
                  paste(c(genes[1:10], paste("... and", length(genes) - 10, "more")), collapse = ", ")
                } else {
                  paste(genes, collapse = ", ")
                }
                data$gene_symbols[i] <- display_genes
              }
            } else if (is.list(enrichment_obj) && "enrichment" %in% names(enrichment_obj)) {
              # STRING format
              result_df <- enrichment_obj$enrichment
              term_row <- result_df[result_df$term == data$ID[i] | result_df$description == data$Description[i], ]
              if (nrow(term_row) > 0 && !is.na(term_row$inputGenes[1])) {
                genes <- unlist(strsplit(term_row$inputGenes[1], ","))
                data$gene_list[[i]] <- genes
                display_genes <- if (length(genes) > 10) {
                  paste(c(genes[1:10], paste("... and", length(genes) - 10, "more")), collapse = ", ")
                } else {
                  paste(genes, collapse = ", ")
                }
                data$gene_symbols[i] <- display_genes
              }
            }
          }
          
          # Fallback if no genes found
          if (is.null(data$gene_list[[i]])) {
            data$gene_symbols[i] <- "Genes not available"
          }
          
        }, error = function(e) {
          data$gene_symbols[i] <- "Error loading genes"
        })
      }
      
      return(data)
    }

    # Helper functions for data preparation
    prepare_dotplot_data <- function(data, n_terms) {
      if (is.null(data) || nrow(data) == 0) {
        message("prepare_dotplot_data: No data provided")
        return(NULL)
      }
      
      message("prepare_dotplot_data: Input data has ", nrow(data), " rows")
      
      # Ensure p.adjust column exists and is numeric
      if (!"p.adjust" %in% names(data)) {
        message("Error: p.adjust column missing from data")
        return(NULL)
      }
      
      # Convert p.adjust to numeric if it's not
      data$p.adjust <- as.numeric(data$p.adjust)
      data <- data[!is.na(data$p.adjust), ]  # Remove any rows with NA p-values
      
      message("After cleaning p.adjust: ", nrow(data), " rows")
      
      # Sort by p-value and take top terms
      data <- data[order(data$p.adjust), ]
      if (nrow(data) > n_terms) {
        data <- data[1:n_terms, ]
      }
      
      message("After filtering to top ", n_terms, " terms: ", nrow(data), " rows")
      
      # Add -log10 p-value with safety checks
      pvals <- data$p.adjust
      pvals[pvals <= 0] <- 1e-300  # Replace zero or negative values
      data$neg_log10_pval <- -log10(pvals)
      
      # Ensure Count column is numeric
      if ("Count" %in% names(data)) {
        data$Count <- as.numeric(data$Count)
        data$Count[is.na(data$Count)] <- 1  # Replace NA with 1
      } else {
        data$Count <- 1  # Default count
      }
      
      # Ensure Description column exists
      if (!"Description" %in% names(data)) {
        if ("description" %in% names(data)) {
          data$Description <- data$description
        } else {
          data$Description <- paste("Term", 1:nrow(data))
        }
      }
      
      # Calculate additional metrics for visualization options
      if ("GeneRatio" %in% names(data) && is.character(data$GeneRatio)) {
        # Parse GeneRatio if it's in "x/y" format
        ratio_parts <- strsplit(data$GeneRatio, "/")
        data$GeneRatio_numeric <- sapply(ratio_parts, function(x) {
          if(length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA
        })
      } else if (!("GeneRatio" %in% names(data))) {
        data$GeneRatio_numeric <- data$Count / 100  # Rough estimate
      }
      
      # Calculate Fold Enrichment if not present
      if (!("FoldEnrichment" %in% names(data))) {
        if ("RichFactor" %in% names(data)) {
          data$FoldEnrichment <- data$RichFactor
        } else if ("GeneRatio_numeric" %in% names(data)) {
          data$FoldEnrichment <- data$GeneRatio_numeric * 10  # Rough approximation
        } else {
          data$FoldEnrichment <- 1
        }
      }
      
      # Ensure RichFactor exists
      if (!("RichFactor" %in% names(data))) {
        data$RichFactor <- data$FoldEnrichment
      }
      
      # Format DE genes if present (convert / to comma-separated)
      if ("geneID" %in% names(data)) {
        data$gene_list_formatted <- gsub("/", ", ", data$geneID)
      }
      
      # Create hover text here so it's available for both static and interactive plots
      data$hover_text <- paste0(
        "<b>", data$Description, "</b><br>",
        "ID: ", ifelse("ID" %in% names(data), data$ID, "N/A"), "<br>",
        "Gene Count: ", data$Count, "<br>",
        "Adj. p-value: ", sprintf("%.2e", data$p.adjust), "<br>",
        "-log10(adj. p-value): ", sprintf("%.2f", data$neg_log10_pval), "<br>",
        "Fold Enrichment: ", sprintf("%.2f", data$FoldEnrichment)
      )
      
      # Add DE genes if available
      if ("gene_list_formatted" %in% names(data)) {
        # Limit gene list display to prevent huge tooltips
        gene_display <- sapply(data$gene_list_formatted, function(genes) {
          gene_vec <- unlist(strsplit(genes, ", "))
          if (length(gene_vec) > 10) {
            paste(c(gene_vec[1:10], paste("... and", length(gene_vec) - 10, "more")), collapse = ", ")
          } else {
            genes
          }
        })
        data$hover_text <- paste0(data$hover_text, "<br><br><b>DE Genes:</b><br>", gene_display)
      }
      
      message("Final data ready with columns: ", paste(names(data), collapse = ", "))
      
      # Load gene lists for hover information (optional, may be slow)
      tryCatch({
        # For now, skip gene list loading to focus on basic plotting
        message("Skipping gene list loading for faster performance")
        # data <- load_gene_lists(data)
      }, error = function(e) {
        message("Error loading gene lists: ", e$message)
      })
      
      return(data)
    }
    
    # Removed other prepare functions - keeping only dotplot for simplicity
    
    # Plotting functions
    create_dotplot <- function(data, color_var, show_labels, x_axis = "neg_log10_pval", selection = NULL) {
      if (is.null(data) || nrow(data) == 0) {
        return(plot_placeholder("No data available for dot plot"))
      }
      
      # Debug: Check data structure
      message("=== Dotplot Debug ===")
      message("Data dimensions: ", nrow(data), " x ", ncol(data))
      message("X-axis selected: ", x_axis)
      message("Color by: ", color_var)
      
      # Map x_axis selection to actual column names
      x_col_map <- list(
        "neg_log10_pval" = "neg_log10_pval",
        "FoldEnrichment" = "FoldEnrichment",
        "Count" = "Count",
        "GeneRatio" = "GeneRatio_numeric",
        "RichFactor" = "RichFactor"
      )
      
      x_column <- if(x_axis %in% names(x_col_map)) x_col_map[[x_axis]] else "neg_log10_pval"
      
      # Ensure the selected x column exists
      if (!(x_column %in% names(data))) {
        message("Warning: ", x_column, " not found, using neg_log10_pval")
        x_column <- "neg_log10_pval"
      }
      
      # Map color_var to actual columns
      color_col_map <- list(
        "p-value" = "neg_log10_pval",
        "Fold Enrichment" = "FoldEnrichment"
      )
      
      color_column <- if(color_var %in% names(color_col_map)) color_col_map[[color_var]] else "neg_log10_pval"
      
      if (!(color_column %in% names(data))) {
        color_column <- "neg_log10_pval"
      }
      
      # Hover text is now created in prepare_dotplot_data, so no need to recreate it here
      # Just verify it exists
      if (!"hover_text" %in% names(data)) {
        message("Warning: hover_text column missing, creating basic version")
        data$hover_text <- paste0(
          "<b>", data$Description, "</b><br>",
          "Gene Count: ", data$Count, "<br>",
          "Adj. p-value: ", sprintf("%.2e", data$p.adjust)
        )
      }
      
      # Truncate descriptions for better display
      data$Description_short <- ifelse(nchar(data$Description) > 50, 
                                      paste0(substr(data$Description, 1, 47), "..."),
                                      data$Description)
      
      # Sort data by the selected x variable for y-axis ordering
      data <- data[order(data[[x_column]], decreasing = TRUE), ]
      
      # Get x-axis label
      x_labels <- list(
        "neg_log10_pval" = "-log10(adjusted p-value)",
        "FoldEnrichment" = "Fold Enrichment",
        "Count" = "Gene Count",
        "GeneRatio_numeric" = "Gene Ratio",
        "RichFactor" = "Rich Factor"
      )
      x_label <- if(x_column %in% names(x_labels)) x_labels[[x_column]] else x_column
      
      # Create appropriate color scale based on variable type
      if (color_column %in% c("direction", "enrichment_type")) {
        # Categorical color scale
        p <- ggplot(data, aes_string(x = x_column, y = paste0("reorder(Description_short, ", x_column, ")"))) +
          geom_point(aes_string(size = "Count", color = color_column)) +
          scale_size_continuous(name = "Gene Count", range = c(3, 10))
      } else {
        # Continuous color scale
        p <- ggplot(data, aes_string(x = x_column, y = paste0("reorder(Description_short, ", x_column, ")"))) +
          geom_point(aes_string(size = "Count", color = color_column)) +
          scale_color_gradient(low = "blue", high = "red", name = color_var) +
          scale_size_continuous(name = "Gene Count", range = c(3, 10))
      }
      
      # Create descriptive title
      plot_title <- if (!is.null(selection)) {
        paste0(
          selection$analysis_type, " Analysis: ",
          selection$gene, " (", selection$cluster, ")\n",
          selection$enrichment_type, " - ", selection$direction, " genes",
          if(selection$experiment != "default") paste0(" - ", selection$experiment) else ""
        )
      } else {
        "Enrichment Dot Plot"
      }
      
      p <- p +
        labs(x = x_label, y = "", title = plot_title) +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10),
              plot.title = element_text(size = 12, face = "bold", hjust = 0.5))
      
      if (!show_labels) {
        p <- p + theme(axis.text.y = element_blank())
      }
      
      # Store the hover text in the data for plotly
      attr(p, "hover_data") <- data$hover_text
      
      message("Plot created successfully")
      return(p)
    }
    
    create_barplot <- function(data, color_var) {
      if (is.null(data)) {
        return(plot_placeholder("No data available for bar plot"))
      }
      
      # Prepare hover text with gene information (same as dotplot)
      if ("gene_symbols" %in% colnames(data)) {
        data$hover_text <- paste0(
          "<b>", data$Description, "</b><br>",
          "ID: ", data$ID, "<br>",
          "Gene Count: ", data$Count, "<br>",
          "p-value: ", sprintf("%.2e", data$p.adjust), "<br>",
          "-log10(p-value): ", sprintf("%.2f", data$neg_log10_pval), "<br>",
          "<b>DE Genes:</b><br>", data$gene_symbols
        )
      } else {
        data$hover_text <- paste0(
          "<b>", data$Description, "</b><br>",
          "ID: ", data$ID, "<br>",
          "Gene Count: ", data$Count, "<br>",
          "p-value: ", sprintf("%.2e", data$p.adjust), "<br>",
          "-log10(p-value): ", sprintf("%.2f", data$neg_log10_pval)
        )
      }
      
      # Create bar plot with hover text
      p <- ggplot(data, aes(x = reorder(Description, neg_log10_pval), y = neg_log10_pval,
                           text = hover_text)) +
        geom_bar(stat = "identity", aes(fill = neg_log10_pval)) +
        scale_fill_gradient(low = "lightblue", high = "darkred", name = "-log10(p-value)") +
        coord_flip() +
        labs(x = "", y = "-log10(adjusted p-value)", 
             title = "Interactive Enrichment Bar Chart") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10),
              plot.title = element_text(size = 14, face = "bold"))
      
      return(p)
    }
    
    create_bubble_heatmap <- function(data, metric) {
      if (is.null(data)) {
        return(plot_placeholder("No data available for bubble heatmap"))
      }
      
      # Source bubble heatmap functions if not already loaded
      if (!exists("create_bubble_heatmap")) {
        source("R/bubble_heatmap_functions.R")
      }
      
      # Create the bubble heatmap
      tryCatch({
        ht <- create_bubble_heatmap(
          data = data,
          max_terms = 30,
          cluster_rows = TRUE,
          cluster_cols = TRUE,
          color_scale = "red",
          size_encoding = metric,  # "count" or "pvalue"
          title = paste("Enrichment Bubble Heatmap -", 
                       ifelse(metric == "count", "Size: Gene Count", "Size: Significance"))
        )
        
        # Return as a grob for plotting
        grid.newpage()
        draw(ht)
        return(recordPlot())
        
      }, error = function(e) {
        return(plot_placeholder(paste("Error creating bubble heatmap:\n", e$message)))
      })
    }
    
    create_network_plot <- function(data) {
      if (is.null(data)) {
        return(plot_placeholder("No data available for network plot"))
      }
      
      # Placeholder for network visualization
      return(plot_placeholder("Network Visualization\n(To be implemented)"))
    }
    
    create_enrichment_map <- function(data) {
      if (is.null(data)) {
        return(plot_placeholder("No data available for enrichment map"))
      }
      
      # Placeholder for enrichment map
      return(plot_placeholder("Enrichment Map\n(To be implemented)"))
    }
    
    create_upset_plot <- function(data) {
      if (is.null(data)) {
        return(plot_placeholder("No data available for upset plot"))
      }
      
      # Placeholder for upset plot
      return(plot_placeholder("Upset Plot\n(Requires UpSetR package)"))
    }
    
    plot_placeholder <- function(text) {
      ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = text, 
                 size = 8, color = "gray50") +
        theme_void() +
        xlim(0, 1) + ylim(0, 1)
    }
    
    # Render main plot
    output$main_plot <- renderPlot({
      if (!is.null(plot_data$plot)) {
        print(plot_data$plot)
      } else {
        plot(1, type = "n", xlab = "", ylab = "", main = "Loading enrichment visualization...")
      }
    })
    
    # Render interactive plot - always try to show something useful
    output$interactive_plot <- renderPlotly({
      # Default to dotplot if no plot type specified
      current_plot_type <- if(is.null(input$plot_type)) "dotplot" else input$plot_type
      
      if (!is.null(plot_data$data) && nrow(plot_data$data) > 0) {
        tryCatch({
          message("=== Creating interactive plot ===")
          message("Plot type: ", current_plot_type)
          message("Data rows: ", nrow(plot_data$data))
          
          # Use the data directly for plotly
          df <- plot_data$data
          
          # Create native plotly plot
          if (current_plot_type == "dotplot" || current_plot_type == "bar") {
            # The hover_text should already exist from prepare_dotplot_data
            # If not, create a basic one but this shouldn't happen
            if (!"hover_text" %in% names(df)) {
              message("Warning: hover_text missing, creating basic version")
              df$hover_text <- paste0("<b>", df$Description, "</b><br>",
                                    "Count: ", df$Count, "<br>",
                                    "-log10(adj. p): ", round(df$neg_log10_pval, 2))
            }
            
            # Get x-axis selection
            x_axis_selected <- if(!is.null(input$x_axis)) input$x_axis else "neg_log10_pval"
            x_col_map <- list(
              "neg_log10_pval" = "neg_log10_pval",
              "FoldEnrichment" = "FoldEnrichment",
              "Count" = "Count",
              "GeneRatio" = "GeneRatio_numeric",
              "RichFactor" = "RichFactor"
            )
            x_column <- if(x_axis_selected %in% names(x_col_map)) x_col_map[[x_axis_selected]] else "neg_log10_pval"
            
            # Ensure column exists
            if (!(x_column %in% names(df))) {
              x_column <- "neg_log10_pval"
            }
            
            # Get color selection
            color_selected <- if(!is.null(input$color_by)) input$color_by else "p-value"
            color_col_map <- list(
              "p-value" = "neg_log10_pval",
              "Fold Enrichment" = "FoldEnrichment"
            )
            color_column <- if(color_selected %in% names(color_col_map)) color_col_map[[color_selected]] else "neg_log10_pval"
            
            if (!(color_column %in% names(df))) {
              color_column <- "neg_log10_pval"
            }
            
            # Sort data for y-axis ordering based on x variable
            df <- df[order(df[[x_column]]), ]
            df$Description_ordered <- factor(df$Description, levels = df$Description)
            
            # Get x-axis label
            x_labels <- list(
              "neg_log10_pval" = "-log10(adjusted p-value)",
              "FoldEnrichment" = "Fold Enrichment",
              "Count" = "Gene Count",
              "GeneRatio_numeric" = "Gene Ratio",
              "RichFactor" = "Rich Factor"
            )
            x_label <- if(x_column %in% names(x_labels)) x_labels[[x_column]] else x_column
            
            # Calculate better size scaling (smaller dots)
            size_scale <- sqrt(df$Count)  # Square root for better scaling
            max_size <- max(size_scale, na.rm = TRUE)
            min_size <- min(size_scale, na.rm = TRUE)
            
            # Normalize sizes to reasonable range (4-20 pixels)
            df$plot_size <- 4 + (size_scale - min_size) / (max_size - min_size) * 16
            
            # Create descriptive title for interactive plot
            plot_title <- if (!is.null(global_selection())) {
              selection <- global_selection()
              paste0(
                "<b>", selection$analysis_type, " Analysis: ",
                selection$gene, " (", selection$cluster, ")</b><br>",
                "<span style='font-size:14px'>", selection$enrichment_type, " - ", selection$direction, " genes",
                if(selection$experiment != "default") paste0(" - ", selection$experiment) else "",
                "</span>"
              )
            } else {
              "<b>Interactive Enrichment Dot Plot</b>"
            }
            
            p <- plot_ly(df, 
                       x = as.formula(paste0("~", x_column)), 
                       y = ~Description_ordered,
                       type = 'scatter',
                       mode = 'markers',
                       marker = list(
                         size = ~plot_size,
                         color = as.formula(paste0("~", color_column)),
                         colorscale = if(color_column %in% c("direction", "enrichment_type")) 'Viridis' else list(c(0, 'blue'), c(1, 'red')),
                         colorbar = list(title = if(color_selected == "p-value") "-log10(adj. p-value)" else color_selected),
                         line = list(width = 1, color = 'rgba(0,0,0,0.3)')
                       ),
                       text = ~hover_text,
                       hovertemplate = '%{text}<extra></extra>') %>%
              layout(
                title = list(text = plot_title,
                            font = list(size = 16)),
                xaxis = list(title = x_label),
                yaxis = list(title = "", tickfont = list(size = 10)),
                hovermode = "closest",
                showlegend = FALSE,
                margin = list(l = 200, t = 80)  # More space for y-axis labels and title
              )
              
            message("Interactive plot created successfully")
            return(p)
          } else {
            # For other plot types, show message
            plotly_empty() %>%
              layout(title = list(text = "Interactive version coming soon for this plot type"))
          }
          
        }, error = function(e) {
          message("Error creating interactive plot: ", e$message)
          # Return error message
          if (!is.null(plot_data$data)) {
            df <- plot_data$data
            if (all(c("neg_log10_pval", "Description", "Count") %in% names(df))) {
              plot_ly(df, 
                     x = ~neg_log10_pval, 
                     y = ~reorder(Description, neg_log10_pval),
                     size = ~Count,
                     type = 'scatter',
                     mode = 'markers',
                     hovertemplate = ~paste0("<b>", Description, "</b><br>",
                                           "Gene Count: ", Count, "<br>",
                                           "-log10(p-value): ", round(neg_log10_pval, 2),
                                           "<extra></extra>")) %>%
                layout(title = "Enrichment Dot Plot",
                       xaxis = list(title = "-log10(adjusted p-value)"),
                       yaxis = list(title = ""))
            } else {
              plotly_empty() %>%
                layout(title = list(text = "No data available - check your selections"))
            }
          } else {
            plotly_empty() %>%
              layout(title = list(text = "Loading... Please wait"))
          }
        })
      } else {
        # Show loading message instead of error
        plotly_empty() %>%
          layout(
            title = list(text = "Loading enrichment data...", font = list(size = 14)),
            xaxis = list(title = ""),
            yaxis = list(title = "")
          )
      }
    })
    
    # Plot information
    output$plot_info <- renderPrint({
      if (!is.null(plot_data$data)) {
        cat("Plot Type:", input$plot_type, "\n")
        cat("Number of terms:", nrow(plot_data$data), "\n")
        cat("Color by:", input$color_by, "\n")
        if (input$plot_type == "bubble") {
          cat("Bubble metric:", input$bubble_metric, "\n")
        }
        if (input$plot_type == "network") {
          cat("Similarity cutoff:", input$similarity_cutoff, "\n")
        }
      } else {
        cat("No plot data available")
      }
    })
    
    # Data table
    output$plot_data <- DT::renderDataTable({
      if (!is.null(plot_data$data)) {
        df <- plot_data$data
        
        # Define columns to exclude
        exclude_cols <- c("method", "mutation_perturbation", "cluster", "enrichment_type", 
                         "direction", "gene", "geneID", "analysis_type", "GeneRatio_numeric", 
                         "neg_log10_pval", "hover_text", "Description_short", "plot_size",
                         "Description_ordered")
        
        # Remove excluded columns
        df <- df[, !(names(df) %in% exclude_cols), drop = FALSE]
        
        # Remove columns that are all NA or empty
        df <- df[, sapply(df, function(x) !all(is.na(x) | x == "")), drop = FALSE]
        
        # Reorder columns to put experiment and gene_list_formatted after qvalue
        col_order <- names(df)
        
        # Find positions
        qvalue_pos <- which(col_order == "qvalue")
        experiment_pos <- which(col_order == "experiment")
        gene_list_pos <- which(col_order == "gene_list_formatted")
        
        # Create new order
        if (length(qvalue_pos) > 0) {
          # Get columns before qvalue
          before_qvalue <- col_order[1:qvalue_pos]
          
          # Get columns to move
          cols_to_move <- c()
          if (length(experiment_pos) > 0) cols_to_move <- c(cols_to_move, "experiment")
          if (length(gene_list_pos) > 0) cols_to_move <- c(cols_to_move, "gene_list_formatted")
          
          # Get remaining columns
          remaining <- col_order[!(col_order %in% c(before_qvalue, cols_to_move))]
          
          # New order
          new_order <- c(before_qvalue, cols_to_move, remaining)
          df <- df[, new_order, drop = FALSE]
        }
        
        # Rename gene_list_formatted for display
        if ("gene_list_formatted" %in% names(df)) {
          names(df)[names(df) == "gene_list_formatted"] <- "DE_Genes"
        }
        
        DT::datatable(df,
                      options = list(pageLength = 10, scrollX = TRUE),
                      rownames = FALSE)
      }
    })
    
    # Download handler
    output$download_plot <- downloadHandler(
      filename = function() {
        if (!is.null(global_selection())) {
          selection <- global_selection()
          # Create descriptive filename
          paste0(
            selection$analysis_type, "_",
            selection$gene, "_",
            selection$cluster, "_",
            selection$enrichment_type, "_",
            selection$direction, "_",
            if(selection$experiment != "default") paste0(selection$experiment, "_") else "",
            Sys.Date(), ".pdf"
          )
        } else {
          paste0("enrichment_", input$plot_type, "_", Sys.Date(), ".pdf")
        }
      },
      content = function(file) {
        if (!is.null(plot_data$plot)) {
          ggsave(file, plot = plot_data$plot, width = 10, height = 8, device = "pdf")
        }
      }
    )
    
    return(plot_data)
  })
}