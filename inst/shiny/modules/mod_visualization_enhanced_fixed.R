# Module: Enhanced Visualization with Gene Display and Hover Tooltips
# Complete fix for both table gene display and interactive hover tooltips

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
                      value = TRUE),
        
        checkboxInput(ns("show_genes_in_hover"),
                      "Show Genes in Hover",
                      value = TRUE)
      ),
      
      # GSEA-specific settings (keeping original)
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
      helpText("🧬 Gene lists are shown in Plot Details table and hover tooltips"),
      
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
                 div(
                   style = "background-color: #e8f4fd; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                   icon("dna", style = "color: #0073b7;"),
                   strong(" Gene Lists: ", style = "color: #0073b7;"),
                   "Each enrichment term shows its associated DE genes in the 'Associated_Genes' column below."
                 ),
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
    
    # Enhanced gene association loading
    gene_associations_loaded <- reactiveVal(FALSE)
    gene_data <- reactiveVal(NULL)
    
    observe({
      tryCatch({
        # Load gene association data
        gene_file <- "../../inst/extdata/gene_term_associations.rds"
        if (file.exists(gene_file)) {
          data <- readRDS(gene_file)
          gene_data(data)
          gene_associations_loaded(TRUE)
          message("Gene associations loaded: ", nrow(data), " associations")
        } else {
          message("Gene association file not found at: ", gene_file)
          gene_associations_loaded(FALSE)
        }
      }, error = function(e) {
        message("Failed to load gene associations: ", e$message)
        gene_associations_loaded(FALSE)
      })
    })
    
    # Enhanced gene lookup function with better matching
    get_genes_for_term_enhanced <- function(term_id, current_selection) {
      if (!gene_associations_loaded() || is.null(gene_data())) {
        return(list(genes = NULL, error = "Gene data unavailable"))
      }
      
      data <- gene_data()
      
      # Create composite key to match
      composite_key <- paste(
        current_selection$analysis_type,
        current_selection$gene,
        current_selection$cluster,
        current_selection$enrichment_type,
        current_selection$direction,
        "default",
        term_id,
        sep = "|"
      )
      
      # Try exact match
      match_idx <- which(data$composite_key == composite_key)
      
      if (length(match_idx) > 0) {
        genes <- unlist(strsplit(data$associated_genes[match_idx[1]], "/"))
        return(list(genes = genes, error = NULL))
      }
      
      # Fallback: find any match for this term
      term_matches <- data[data$term_id == term_id, ]
      if (nrow(term_matches) > 0) {
        genes <- unlist(strsplit(term_matches$associated_genes[1], "/"))
        return(list(genes = genes, error = NULL))
      }
      
      return(list(genes = NULL, error = "No genes found"))
    }
    
    # Format gene list for table display
    format_gene_list_table <- function(genes, max_genes = 20) {
      if (is.null(genes) || length(genes) == 0) return("No genes found")
      
      if (length(genes) <= max_genes) {
        return(paste(genes, collapse = ", "))
      } else {
        displayed_genes <- head(genes, max_genes)
        remaining <- length(genes) - max_genes
        return(paste0(paste(displayed_genes, collapse = ", "), ", ... ", remaining, " more"))
      }
    }
    
    # Format gene list for hover tooltips
    format_gene_list_hover <- function(genes, max_genes = 15, genes_per_line = 6) {
      if (is.null(genes) || length(genes) == 0) return("No genes")
      
      # Limit genes for readability
      display_genes <- head(genes, max_genes)
      
      # Split into lines
      gene_lines <- split(display_genes, ceiling(seq_along(display_genes) / genes_per_line))
      formatted_lines <- sapply(gene_lines, function(line) paste(line, collapse = ", "))
      
      result <- paste(formatted_lines, collapse = "<br>")
      
      # Add total count and remaining
      if (length(genes) > max_genes) {
        remaining <- length(genes) - max_genes
        result <- paste0("Genes (", length(genes), "):<br>", result, "<br>... ", remaining, " more")
      } else {
        result <- paste0("Genes (", length(genes), "):<br>", result)
      }
      
      return(result)
    }
    
    # Reactive values for plot data
    plot_data <- reactiveValues(
      data = NULL,
      plot = NULL,
      interactive = NULL,
      is_gsea = FALSE
    )
    
    # Check if current selection is GSEA
    output$is_gsea <- reactive({
      result <- processed_data()
      isTRUE(result$is_gsea)
    })
    outputOptions(output, "is_gsea", suspendWhenHidden = FALSE)
    
    # Data preparation function (keeping existing logic)
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
      data <- data[!is.na(data$p.adjust), ]
      
      message("After cleaning p.adjust: ", nrow(data), " rows")
      
      # Sort by p-value and take top terms
      data <- data[order(data$p.adjust), ]
      if (nrow(data) > n_terms) {
        data <- data[1:n_terms, ]
      }
      
      message("After filtering to top ", n_terms, " terms: ", nrow(data), " rows")
      
      # Add -log10 p-value with safety checks
      pvals <- data$p.adjust
      pvals[pvals <= 0] <- 1e-300
      data$neg_log10_pval <- -log10(pvals)
      
      # Ensure required columns
      if ("Count" %in% names(data)) {
        data$Count <- as.numeric(data$Count)
        data$Count[is.na(data$Count)] <- 1
      } else {
        data$Count <- 1
      }
      
      if (!"Description" %in% names(data)) {
        if ("description" %in% names(data)) {
          data$Description <- data$description
        } else {
          data$Description <- paste("Term", 1:nrow(data))
        }
      }
      
      # Calculate GeneRatio_numeric for plotting
      if ("GeneRatio" %in% names(data) && is.character(data$GeneRatio)) {
        ratio_parts <- strsplit(data$GeneRatio, "/")
        data$GeneRatio_numeric <- sapply(ratio_parts, function(x) {
          if(length(x) == 2) as.numeric(x[1]) / as.numeric(x[2]) else NA
        })
      } else if (!("GeneRatio" %in% names(data))) {
        data$GeneRatio_numeric <- data$Count / 100
      }
      
      # Calculate FoldEnrichment if not present
      if (!("FoldEnrichment" %in% names(data))) {
        if ("RichFactor" %in% names(data)) {
          data$FoldEnrichment <- data$RichFactor
        } else if ("GeneRatio_numeric" %in% names(data)) {
          data$FoldEnrichment <- data$GeneRatio_numeric * 10
        } else {
          data$FoldEnrichment <- 1
        }
      }
      
      if (!("RichFactor" %in% names(data))) {
        data$RichFactor <- data$FoldEnrichment
      }
      
      message("Final data ready with columns: ", paste(names(data), collapse = ", "))
      return(data)
    }

    # Reactive data processing
    processed_data <- reactive({
      req(global_selection(), enrichment_data())
      
      selection <- global_selection()
      data <- enrichment_data()
      
      # Check if this is GSEA data
      is_gsea <- selection$enrichment_type == "GSEA"
      
      # Process the data and return result
      if (nrow(data) > 0) {
        if (!is_gsea) {
          prepared_data <- prepare_dotplot_data(data, input$top_terms %||% 20)
        } else {
          prepared_data <- data
        }
        
        return(list(
          data = prepared_data,
          is_gsea = is_gsea
        ))
      } else {
        return(list(
          data = NULL,
          is_gsea = is_gsea
        ))
      }
    })
    
    # Enhanced plotting function with gene hover tooltips
    create_standard_plot_with_genes <- function(data, plot_type, show_genes_hover = TRUE) {
      top_data <- data %>%
        arrange(p.adjust) %>%
        head(input$top_terms)
      
      if (nrow(top_data) == 0) {
        return(ggplot() + 
               annotate("text", x = 0.5, y = 0.5, label = "No data available", size = 5) +
               theme_void())
      }
      
      # Add gene information for hover if enabled
      if (show_genes_hover && gene_associations_loaded()) {
        current_selection <- global_selection()
        
        # Get genes for each term
        top_data$hover_genes <- sapply(top_data$ID, function(term_id) {
          gene_result <- get_genes_for_term_enhanced(term_id, current_selection)
          if (!is.null(gene_result$genes)) {
            return(format_gene_list_hover(gene_result$genes))
          } else {
            return("No genes found")
          }
        })
      } else {
        top_data$hover_genes <- "Gene display disabled"
      }
      
      # Generate descriptive title
      selection <- global_selection()
      gene_part <- if (!is.null(selection$gene) && selection$gene != "" && selection$gene != "All") {
        selection$gene
      } else {
        NULL
      }
      
      if (!is.null(selection$analysis_type) && selection$analysis_type == "MAST") {
        if (!is.null(gene_part)) {
          comparison_part <- paste0(gene_part, " mutation vs isogenic eWT controls (MAST)")
        } else {
          comparison_part <- "MAST mutation analysis vs isogenic eWT controls"
        }
      } else if (!is.null(selection$analysis_type) && grepl("MixScale", selection$analysis_type)) {
        if (!is.null(gene_part)) {
          comparison_part <- paste0(gene_part, " CRISPRi knockdown vs Non-Targeting")
        } else {
          comparison_part <- "CRISPRi knockdown vs Non-Targeting"
        }
      } else {
        comparison_part <- if (!is.null(gene_part)) paste0(gene_part, " analysis") else "Enrichment analysis"
      }
      
      enrichment_part <- if (!is.null(selection$enrichment_type) && selection$enrichment_type != "") {
        selection$enrichment_type
      } else {
        "Enrichment"
      }
      
      direction_part <- if (!is.null(selection$direction) && selection$direction != "ALL") {
        paste0(selection$direction, "-regulated genes")
      } else {
        "All genes"
      }
      
      cluster_part <- if (!is.null(selection$cluster) && selection$cluster != "" && selection$cluster != "All") {
        paste0("Cluster ", gsub("cluster_", "", selection$cluster))
      } else {
        "All clusters"
      }
      
      plot_title <- paste0(comparison_part, "\n", 
                          enrichment_part, " (", direction_part, ")\n", 
                          cluster_part)
      
      # Create plot based on type
      if (plot_type == "dotplot") {
        plot_df <- top_data %>%
          mutate(
            Description = make.unique(Description, sep = " "),
            Description = factor(Description, levels = rev(unique(Description))),
            neg_log10_pval = -log10(pmax(p.adjust, 1e-100)),
            Count = if("Count" %in% names(.)) Count else 10,
            FoldEnrichment = if("FoldEnrichment" %in% names(.)) FoldEnrichment else 2
          )
        
        # Convert GeneRatio for plotting
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
          plot_df$GeneRatio_numeric <- plot_df$Count / 100
        }
        
        # Determine x-axis variable
        x_var <- switch(input$x_axis,
          "neg_log10_pval" = "neg_log10_pval",
          "FoldEnrichment" = if("FoldEnrichment" %in% names(plot_df)) "FoldEnrichment" else "neg_log10_pval",
          "Count" = if("Count" %in% names(plot_df)) "Count" else "neg_log10_pval",
          "GeneRatio" = if("GeneRatio_numeric" %in% names(plot_df)) "GeneRatio_numeric" else "neg_log10_pval",
          "RichFactor" = if("RichFactor" %in% names(plot_df)) "RichFactor" else "neg_log10_pval",
          "neg_log10_pval"
        )
        
        if (!is.numeric(plot_df[[x_var]]) || all(is.na(plot_df[[x_var]]))) {
          x_var <- "neg_log10_pval"
        }
        
        # Create custom hover text
        if (show_genes_hover) {
          plot_df$custom_hover <- paste0(
            "<b>", plot_df$Description, "</b><br>",
            "p.adjust: ", format(plot_df$p.adjust, scientific = TRUE, digits = 3), "<br>",
            "Count: ", plot_df$Count, "<br>",
            gsub("_", " ", tools::toTitleCase(x_var)), ": ", round(plot_df[[x_var]], 3), "<br><br>",
            plot_df$hover_genes
          )
        } else {
          plot_df$custom_hover <- paste0(
            "<b>", plot_df$Description, "</b><br>",
            "p.adjust: ", format(plot_df$p.adjust, scientific = TRUE, digits = 3), "<br>",
            "Count: ", plot_df$Count, "<br>",
            gsub("_", " ", tools::toTitleCase(x_var)), ": ", round(plot_df[[x_var]], 3)
          )
        }
        
        # Create base plot
        p <- ggplot(plot_df, aes(x = .data[[x_var]], y = Description, text = custom_hover))
        
        if (input$color_by == "p-value") {
          p <- p + geom_point(aes(color = neg_log10_pval, size = Count), alpha = 0.8)
        } else {
          if ("FoldEnrichment" %in% names(plot_df) && is.numeric(plot_df$FoldEnrichment)) {
            p <- p + geom_point(aes(color = FoldEnrichment, size = Count), alpha = 0.8)
          } else {
            p <- p + geom_point(aes(color = neg_log10_pval, size = Count), alpha = 0.8)
          }
        }
        
        p <- p +
          scale_color_gradient(low = "blue", high = "red") +
          scale_size_continuous(range = c(3, 10), guide = guide_legend(title = "Count")) +
          theme_bw() +
          theme(
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 12),
            legend.title = element_text(size = 10),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 11, margin = margin(b = 20))
          ) +
          labs(x = gsub("_", " ", tools::toTitleCase(gsub("_", " ", x_var))), y = "", title = plot_title)
        
        return(p)
        
      } else if (plot_type == "barplot") {
        plot_df <- top_data %>%
          mutate(
            Description = make.unique(Description, sep = " "),
            Description = factor(Description, levels = rev(unique(Description))),
            neg_log10_pval = -log10(p.adjust)
          )
        
        p <- ggplot(plot_df, aes(x = neg_log10_pval, y = reorder(Description, neg_log10_pval))) +
          geom_bar(stat = "identity", fill = "steelblue") +
          theme_bw() +
          theme(
            plot.title = element_text(hjust = 0.5, size = 11, margin = margin(b = 20))
          ) +
          labs(x = "-log10(adjusted p-value)", y = "", title = plot_title)
        
        return(p)
        
      } else if (plot_type == "lollipop") {
        plot_df <- top_data %>%
          mutate(
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
          theme(
            plot.title = element_text(hjust = 0.5, size = 11, margin = margin(b = 20))
          ) +
          labs(x = "-log10(adjusted p-value)", y = "", title = plot_title)
        
        return(p)
      }
    }
    
    # Render interactive plot with enhanced hover
    output$interactive_plot <- renderPlotly({
      result <- processed_data()
      req(result$data)
      req(!result$is_gsea)
      
      p <- create_standard_plot_with_genes(result$data, input$plot_type, input$show_genes_in_hover)
      if (!is.null(p)) {
        if (input$show_genes_in_hover) {
          # Use custom hover text
          plotly_obj <- ggplotly(p, tooltip = "text") %>%
            config(displayModeBar = TRUE, displaylogo = FALSE) %>%
            layout(
              xaxis = list(title = p$labels$x),
              yaxis = list(title = p$labels$y),
              margin = list(l = 200, r = 50, t = 50, b = 50)
            )
        } else {
          # Use default hover
          plotly_obj <- ggplotly(p, tooltip = c("x", "y", "color", "size")) %>%
            config(displayModeBar = TRUE, displaylogo = FALSE) %>%
            layout(
              xaxis = list(title = p$labels$x),
              yaxis = list(title = p$labels$y),
              margin = list(l = 200, r = 50, t = 50, b = 50)
            )
        }
        
        plotly_obj
      }
    })
    
    # Enhanced data table with prominent gene column
    output$plot_data <- DT::renderDataTable({
      result <- processed_data()
      req(result$data)
      
      current_selection <- global_selection()
      
      if (result$is_gsea) {
        display_data <- result$data %>%
          select(ID, Description, NES, p.adjust, setSize, any_of(c("enrichmentScore", "rank"))) %>%
          arrange(desc(abs(NES)))
      } else {
        display_data <- result$data %>%
          select(ID, Description, p.adjust, Count, any_of(c("FoldEnrichment", "GeneRatio"))) %>%
          arrange(p.adjust)
      }
      
      # Add prominent gene column
      if (gene_associations_loaded() && !is.null(current_selection)) {
        message("Adding gene lists to table...")
        
        display_data$Associated_Genes <- sapply(display_data$ID, function(term_id) {
          gene_result <- get_genes_for_term_enhanced(term_id, current_selection)
          format_gene_list_table(gene_result$genes)
        })
        
        message("Gene lists added for ", nrow(display_data), " terms")
      } else {
        display_data$Associated_Genes <- "Gene data unavailable"
      }
      
      # Create enhanced datatable
      dt <- DT::datatable(
        display_data,
        options = list(
          pageLength = 20,
          scrollX = TRUE,
          columnDefs = list(
            list(width = "200px", targets = which(names(display_data) == "Description") - 1),
            list(width = "350px", targets = which(names(display_data) == "Associated_Genes") - 1),
            list(className = "dt-center", targets = 0:(ncol(display_data)-2))
          ),
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        extensions = 'Buttons',
        rownames = FALSE,
        selection = "single"
      ) %>%
        DT::formatSignif(columns = "p.adjust", digits = 3)
      
      # Highlight gene column
      if ("Associated_Genes" %in% names(display_data)) {
        dt <- dt %>%
          DT::formatStyle("Associated_Genes",
                         fontSize = "90%",
                         color = "#0066cc",
                         backgroundColor = "#f8f9fa",
                         fontWeight = "500")
      }
      
      return(dt)
    })
    
    # Plot information
    output$plot_info <- renderPrint({
      result <- processed_data()
      req(result$data)
      
      cat("Plot Information:\n")
      cat("================\n")
      cat("Data rows:", nrow(result$data), "\n")
      cat("Gene associations loaded:", gene_associations_loaded(), "\n")
      
      if (gene_associations_loaded()) {
        cat("Gene data entries:", nrow(gene_data()), "\n")
      }
      
      if (result$is_gsea) {
        cat("\nGSEA Statistics:\n")
        cat("NES range:", range(result$data$NES), "\n")
        cat("Significant (|NES| > 1):", sum(abs(result$data$NES) > 1), "\n")
      } else {
        cat("\nEnrichment Statistics:\n")
        cat("P-value range:", range(result$data$p.adjust), "\n")
        cat("Significant (p < 0.05):", sum(result$data$p.adjust < 0.05), "\n")
      }
    })
    
    # Download handler (keeping existing)
    output$download_plot <- downloadHandler(
      filename = function() {
        paste0("enrichment_plot_", Sys.Date(), ".pdf")
      },
      content = function(file) {
        result <- processed_data()
        if (result$is_gsea) {
          # Handle GSEA plots here if needed
          p <- NULL
        } else {
          p <- create_standard_plot_with_genes(result$data, input$plot_type, FALSE)
        }
        
        if (!is.null(p)) {
          ggsave(file, p, width = 10, height = 8)
        }
      }
    )
    
    return(processed_data)
  })
}