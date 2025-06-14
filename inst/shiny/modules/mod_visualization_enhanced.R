# Module: Enhanced Visualization with GSEA Support
# Advanced visualizations for enrichment results including GSEA plots

mod_visualization_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Visualization Settings"),
      
      # Conditional UI based on enrichment type
      conditionalPanel(
        condition = "output.is_gsea == false",
        ns = ns,
        
        selectInput(ns("plot_type"),
                    "Visualization Type:",
                    choices = list(
                      "Dot Plot" = "dotplot",
                      "Bar Plot" = "barplot",
                      "Lollipop Plot" = "lollipop"
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
                      value = TRUE)
      ),
      
      # GSEA-specific settings
      conditionalPanel(
        condition = "output.is_gsea == true",
        ns = ns,
        
        selectInput(ns("gsea_plot_type"),
                    "GSEA Plot Type:",
                    choices = list(
                      "Enrichment Plot" = "gseaplot",
                      "Dot Plot" = "dotplot",
                      "Ridge Plot" = "ridgeplot",
                      "GSEA Table" = "table"
                    ),
                    selected = "gseaplot"),
        
        # For enrichment plot
        conditionalPanel(
          condition = "input.gsea_plot_type == 'gseaplot'",
          ns = ns,
          
          selectInput(ns("gsea_term_select"),
                      "Select Gene Set:",
                      choices = NULL),
          
          checkboxInput(ns("gsea_show_pval"),
                        "Show p-value on plot",
                        value = TRUE),
          
          checkboxInput(ns("gsea_show_running"),
                        "Show running score",
                        value = TRUE),
          
          sliderInput(ns("gsea_base_size"),
                      "Base font size:",
                      min = 8,
                      max = 16,
                      value = 11,
                      step = 1)
        ),
        
        # For other GSEA plots
        conditionalPanel(
          condition = "input.gsea_plot_type != 'gseaplot'",
          ns = ns,
          
          numericInput(ns("gsea_top_terms"),
                       "Number of Top Terms:",
                       value = 20,
                       min = 5,
                       max = 50,
                       step = 5)
        )
      ),
      
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
        selected = "primary",
        
        tabPanel("Primary Plot",
                 value = "primary",
                 br(),
                 conditionalPanel(
                   condition = "output.is_gsea == false",
                   ns = ns,
                   plotlyOutput(ns("interactive_plot"), height = "700px") %>%
                     withSpinner(color = "#0073b7")
                 ),
                 conditionalPanel(
                   condition = "output.is_gsea == true",
                   ns = ns,
                   plotOutput(ns("gsea_plot"), height = "700px") %>%
                     withSpinner(color = "#0073b7")
                 )),
        
        
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
      interactive = NULL,
      is_gsea = FALSE
    )
    
    # Check if current selection is GSEA
    output$is_gsea <- reactive({
      selection <- global_selection()
      isTRUE(selection$enrichment_type == "GSEA")
    })
    outputOptions(output, "is_gsea", suspendWhenHidden = FALSE)
    
    # Data preparation function (from working backup version)
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
      
      message("Final data ready with columns: ", paste(names(data), collapse = ", "))
      return(data)
    }

    # Reactive data processing
    observe({
      req(global_selection(), enrichment_data())
      
      selection <- global_selection()
      data <- enrichment_data()
      
      
      # Check if this is GSEA data
      plot_data$is_gsea <- selection$enrichment_type == "GSEA"
      
      # Process the data
      if (nrow(data) > 0) {
        # Apply data preparation for better plotting
        if (!plot_data$is_gsea) {
          plot_data$data <- prepare_dotplot_data(data, input$top_terms %||% 20)
        } else {
          plot_data$data <- data
        }
        
        # Update GSEA term choices if GSEA
        if (plot_data$is_gsea) {
          gsea_terms <- unique(data$Description)
          updateSelectInput(session, "gsea_term_select",
                          choices = gsea_terms,
                          selected = gsea_terms[1])
        }
      } else {
        # Clear plot data when no data available
        plot_data$data <- NULL
      }
    })
    
    # Create standard enrichment plots
    create_standard_plot <- function(data, plot_type) {
      # Select top terms
      top_data <- data %>%
        arrange(p.adjust) %>%
        head(input$top_terms)
      
      # Check if we have any data
      if (nrow(top_data) == 0) {
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
               theme_void())
      }
      
      # Create appropriate plot based on type
      if (plot_type == "dotplot") {
        # Prepare data for dot plot - ensure all required columns exist
        plot_df <- top_data %>%
          mutate(
            # Handle duplicates by making descriptions unique
            Description = make.unique(Description, sep = " "),
            Description = factor(Description, levels = rev(unique(Description))),
            neg_log10_pval = -log10(pmax(p.adjust, 1e-100)),  # Prevent log(0)
            # Ensure Count column exists
            Count = if("Count" %in% names(.)) Count else if("count" %in% names(.)) count else 10,
            # Ensure FoldEnrichment exists 
            FoldEnrichment = if("FoldEnrichment" %in% names(.)) FoldEnrichment else 
                           if("fold_enrichment" %in% names(.)) fold_enrichment else 
                           if("enrichment_score" %in% names(.)) enrichment_score else 2
          )
        
        # CRITICAL FIX: Convert GeneRatio from "x/y" format to numeric
        if ("GeneRatio" %in% names(plot_df) && is.character(plot_df$GeneRatio)) {
          ratio_parts <- strsplit(plot_df$GeneRatio, "/")
          plot_df$GeneRatio_numeric <- sapply(ratio_parts, function(x) {
            if(length(x) == 2 && !is.na(x[1]) && !is.na(x[2])) {
              as.numeric(x[1]) / as.numeric(x[2])
            } else {
              NA
            }
          })
        } else {
          plot_df$GeneRatio_numeric <- plot_df$Count / 100  # Fallback estimate
        }
        
        # Determine x-axis variable with validation
        x_var <- switch(input$x_axis,
          "neg_log10_pval" = "neg_log10_pval",
          "FoldEnrichment" = if("FoldEnrichment" %in% names(plot_df)) "FoldEnrichment" else "neg_log10_pval",
          "Count" = if("Count" %in% names(plot_df)) "Count" else "neg_log10_pval",
          "GeneRatio" = if("GeneRatio_numeric" %in% names(plot_df)) "GeneRatio_numeric" else "neg_log10_pval",  # FIXED: Use numeric version
          "RichFactor" = if("RichFactor" %in% names(plot_df)) "RichFactor" else "neg_log10_pval",
          "neg_log10_pval"
        )
        
        # Validate that we have numeric data for plotting
        if (!is.numeric(plot_df[[x_var]]) || all(is.na(plot_df[[x_var]]))) {
          x_var <- "neg_log10_pval"  # Fallback to p-value
        }
        
        # Create base plot using tidy evaluation
        p <- ggplot(plot_df, aes(x = .data[[x_var]], y = Description))
        
        # Add points with validation
        if (input$color_by == "p-value") {
          p <- p + geom_point(aes(color = neg_log10_pval, size = Count), alpha = 0.8)
        } else {
          # Validate FoldEnrichment exists and is numeric
          if ("FoldEnrichment" %in% names(plot_df) && is.numeric(plot_df$FoldEnrichment)) {
            p <- p + geom_point(aes(color = FoldEnrichment, size = Count), alpha = 0.8)
          } else {
            p <- p + geom_point(aes(color = neg_log10_pval, size = Count), alpha = 0.8)
          }
        }
        
        # Customize appearance
        p <- p +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_continuous(range = c(3, 10), guide = guide_legend(title = "Count")) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 10),
            panel.grid.minor = element_blank()
          ) +
          labs(x = gsub("_", " ", tools::toTitleCase(gsub("_", " ", x_var))), y = "")
        
        return(p)
        
      } else if (plot_type == "barplot") {
        # Bar plot implementation
        plot_df <- top_data %>%
          mutate(
            # Handle duplicates by making descriptions unique
            Description = make.unique(Description, sep = " "),
            Description = factor(Description, levels = rev(unique(Description))),
            neg_log10_pval = -log10(p.adjust)
          )
        
        p <- ggplot(plot_df, aes(x = neg_log10_pval, y = reorder(Description, neg_log10_pval))) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_bw() +
          labs(x = "-log10(adjusted p-value)", y = "")
        
        return(p)
        
      } else if (plot_type == "lollipop") {
        # Lollipop plot implementation
        plot_df <- top_data %>%
          mutate(
            # Handle duplicates by making descriptions unique
            Description = make.unique(Description, sep = " "),
            Description = factor(Description, levels = rev(unique(Description))),
            neg_log10_pval = -log10(p.adjust)
          )
        
        p <- ggplot(plot_df, aes(x = neg_log10_pval, y = reorder(Description, neg_log10_pval))) +
          geom_segment(aes(x = 0, xend = neg_log10_pval, 
                          y = Description, yend = Description),
                      color = "grey50") +
          geom_point(size = 4, color = "steelblue") +
          theme_bw() +
          labs(x = "-log10(adjusted p-value)", y = "")
        
        return(p)
      }
    }
    
    # Create GSEA plots
    create_gsea_plot <- function(data, plot_type) {
      # Check if enrichplot is available
      if (!requireNamespace("enrichplot", quietly = TRUE)) {
        showNotification("enrichplot package is required for GSEA visualization. Install with: BiocManager::install('enrichplot')",
                        type = "error", duration = NULL)
        return(NULL)
      }
      
      library(enrichplot)
      
      if (plot_type == "gseaplot" && !is.null(input$gsea_term_select)) {
        # For gseaplot2, we need the original GSEA result object
        # Since we only have the data frame, we'll create a basic enrichment plot
        
        # Get the selected term data
        term_data <- data %>%
          filter(Description == input$gsea_term_select) %>%
          slice(1)
        
        if (nrow(term_data) == 0) return(NULL)
        
        # Create a simple enrichment plot visualization
        # Note: This is a simplified version since we don't have the full GSEA object
        p <- ggplot() +
          theme_bw() +
          labs(title = paste0("GSEA Enrichment Plot\n", input$gsea_term_select),
               subtitle = paste0("NES = ", round(term_data$NES, 3), 
                               ", p.adjust = ", format(term_data$p.adjust, scientific = TRUE, digits = 3))) +
          theme(plot.title = element_text(size = 14, face = "bold"),
                plot.subtitle = element_text(size = 12))
        
        # Add a note that full gseaplot2 requires original GSEA object
        p <- p + 
          annotate("text", x = 0.5, y = 0.5, 
                  label = "Note: Full enrichment plot requires original GSEA object.\nThis is a summary view.",
                  hjust = 0.5, vjust = 0.5, size = 5, color = "gray40")
        
        return(p)
        
      } else if (plot_type == "dotplot") {
        # GSEA dot plot
        top_data <- data %>%
          arrange(desc(abs(NES))) %>%
          head(input$gsea_top_terms)
        
        p <- ggplot(top_data, aes(x = NES, y = reorder(Description, NES))) +
          geom_point(aes(size = setSize, color = p.adjust)) +
          scale_color_gradient(low = "red", high = "blue") +
          scale_size_continuous(range = c(3, 10)) +
          theme_bw() +
          labs(x = "Normalized Enrichment Score", y = "",
               size = "Set Size", color = "Adjusted\np-value") +
          geom_vline(xintercept = 0, linetype = "dashed", color = "gray50")
        
        return(p)
        
      } else if (plot_type == "ridgeplot") {
        # Ridge plot for GSEA
        if (!requireNamespace("ggridges", quietly = TRUE)) {
          showNotification("ggridges package is required for ridge plots. Install with: install.packages('ggridges')",
                          type = "warning")
          return(NULL)
        }
        
        library(ggridges)
        
        top_data <- data %>%
          arrange(desc(abs(NES))) %>%
          head(input$gsea_top_terms)
        
        # Create a simplified ridge plot
        p <- ggplot(top_data, aes(x = NES, y = reorder(Description, NES))) +
          geom_density_ridges(aes(fill = p.adjust), alpha = 0.7) +
          scale_fill_gradient(low = "red", high = "blue") +
          theme_bw() +
          labs(x = "Normalized Enrichment Score", y = "",
               fill = "Adjusted\np-value")
        
        return(p)
        
      } else if (plot_type == "table") {
        # Return NULL - we'll handle table separately
        return(NULL)
      }
    }
    
    # Render plots based on data type
    output$interactive_plot <- renderPlotly({
      req(plot_data$data)
      req(!plot_data$is_gsea)
      
      p <- create_standard_plot(plot_data$data, input$plot_type)
      if (!is.null(p)) {
        # Convert to plotly with explicit parameters to avoid warnings
        plotly_obj <- ggplotly(p, tooltip = c("x", "y", "color", "size")) %>%
          config(displayModeBar = TRUE, displaylogo = FALSE) %>%
          layout(
            # Ensure proper axis labeling
            xaxis = list(title = p$labels$x),
            yaxis = list(title = p$labels$y),
            # Fix plotly warnings by ensuring proper trace configuration
            margin = list(l = 200, r = 50, t = 50, b = 50)
          )
        
        # Explicitly set trace mode for scatter plots to avoid warnings
        for (i in seq_along(plotly_obj$x$data)) {
          if (is.null(plotly_obj$x$data[[i]]$type)) {
            plotly_obj$x$data[[i]]$type <- "scatter"
          }
          if (is.null(plotly_obj$x$data[[i]]$mode)) {
            plotly_obj$x$data[[i]]$mode <- "markers"
          }
        }
        
        plotly_obj
      }
    })
    
    output$gsea_plot <- renderPlot({
      req(plot_data$data)
      req(plot_data$is_gsea)
      
      create_gsea_plot(plot_data$data, input$gsea_plot_type)
    })
    
    
    # Plot information
    output$plot_info <- renderPrint({
      req(plot_data$data)
      
      cat("Plot Information:\n")
      cat("================\n")
      cat("Data rows:", nrow(plot_data$data), "\n")
      cat("Enrichment type:", unique(plot_data$data$enrichment_type), "\n")
      
      if (plot_data$is_gsea) {
        cat("\nGSEA Statistics:\n")
        cat("NES range:", range(plot_data$data$NES), "\n")
        cat("Significant (|NES| > 1):", sum(abs(plot_data$data$NES) > 1), "\n")
      } else {
        cat("\nEnrichment Statistics:\n")
        cat("P-value range:", range(plot_data$data$p.adjust), "\n")
        cat("Significant (p < 0.05):", sum(plot_data$data$p.adjust < 0.05), "\n")
      }
    })
    
    # Data table
    output$plot_data <- DT::renderDataTable({
      req(plot_data$data)
      
      if (plot_data$is_gsea) {
        # GSEA-specific columns
        display_data <- plot_data$data %>%
          select(Description, NES, p.adjust, setSize, any_of(c("enrichmentScore", "rank"))) %>%
          arrange(desc(abs(NES)))
      } else {
        # Standard enrichment columns
        display_data <- plot_data$data %>%
          select(Description, p.adjust, Count, any_of(c("FoldEnrichment", "GeneRatio"))) %>%
          arrange(p.adjust)
      }
      
      DT::datatable(display_data,
                    options = list(pageLength = 20),
                    rownames = FALSE) %>%
        DT::formatSignif(columns = "p.adjust", digits = 3)
    })
    
    # Download handler
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("enrichment_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        if (plot_data$is_gsea) {
          p <- create_gsea_plot(plot_data$data, 
                               ifelse(input$gsea_plot_type == "table", "dotplot", input$gsea_plot_type))
        } else {
          p <- create_standard_plot(plot_data$data, input$plot_type)
        }
        
        if (!is.null(p)) {
          ggsave(file, p, width = 10, height = 8)
        }
      }
    )
    
    return(plot_data)
  })
}