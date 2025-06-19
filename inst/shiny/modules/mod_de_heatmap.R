# Module: DE Gene Heatmaps
# Interactive heatmap builder with advanced controls for DE genes

library(dplyr)
library(ggplot2)
library(plotly)
if (requireNamespace("heatmaply", quietly = TRUE)) {
  library(heatmaply)
}
if (requireNamespace("pheatmap", quietly = TRUE)) {
  library(pheatmap)
}

# Load shared functions with error handling
tryCatch({
  source("modules/shared/de_processing.R")
}, error = function(e) {
  message("Could not load de_processing.R: ", e$message)
  
  # Define minimal fallback function
  process_de_data_with_metadata <<- function(de_results, metadata = NULL, max_genes_per_condition = 300) {
    message("Using fallback DE processing")
    return(data.frame(
      gene = character(0),
      cluster = character(0),
      gene_name = character(0),
      log2FC = numeric(0),
      pvalue = numeric(0),
      method = character(0),
      source = character(0),
      experiment = character(0)
    ))
  }
})

tryCatch({
  source("modules/shared/source_labeling.R")
}, error = function(e) {
  message("Could not load source_labeling.R: ", e$message)
  
  # Define minimal fallback functions
  add_source_labels <<- function(data, ...) return(data)
  create_source_label <<- function(gene, method, modality = NULL) return(as.character(gene))
})

mod_de_heatmap_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("DE Gene Heatmap Builder"),
      
      # Data scope
      wellPanel(
        h5("Data Scope"),
        selectInput(ns("cluster_select"),
                    "Cluster:",
                    choices = paste0("cluster_", 0:13),
                    selected = "cluster_0"),
        
        helpText("Cluster selection is synced with global settings"),
        
        # Add gene filtering option
        checkboxInput(ns("filter_by_genes"),
                      "Filter by specific genes",
                      value = FALSE),
        
        conditionalPanel(
          condition = paste0("input['", ns("filter_by_genes"), "'] == true"),
          selectizeInput(ns("gene_filter"),
                         "Genes to include:",
                         choices = NULL,
                         multiple = TRUE,
                         options = list(placeholder = "Select genes...",
                                       maxItems = 50))
        )
      ),
      
      # Source filter
      checkboxGroupInput(ns("source_filter"),
                         "Include Sources:",
                         choices = c("iSCORE-PD" = "iSCORE-PD",
                                   "CRISPRi" = "CRISPRi",
                                   "CRISPRa" = "CRISPRa"),
                         selected = c("iSCORE-PD", "CRISPRi")),
      
      # Gene selection
      wellPanel(
        h5("Gene Selection"),
        selectInput(ns("gene_selection_mode"),
                    "Gene Selection:",
                    choices = c("All Significant" = "all",
                                "Top N by Significance" = "top_n",
                                "Custom List" = "custom"),
                    selected = "top_n"),
        
        conditionalPanel(
          condition = paste0("input['", ns("gene_selection_mode"), "'] == 'top_n'"),
          numericInput(ns("top_n_genes"),
                       "Number of Top Genes:",
                       value = 50,
                       min = 10,
                       max = 500,
                       step = 10)
        ),
        
        conditionalPanel(
          condition = paste0("input['", ns("gene_selection_mode"), "'] == 'custom'"),
          textAreaInput(ns("custom_genes"),
                        "Gene Names (one per line):",
                        placeholder = "GENE1\nGENE2\nGENE3",
                        height = "100px")
        )
      ),
      
      # Cutoff controls
      wellPanel(
        h5("Significance Cutoffs"),
        sliderInput(ns("log2fc_cutoff"),
                    "Min |log2FC|:",
                    min = 0,
                    max = 3,
                    value = 0.5,
                    step = 0.1),
        
        sliderInput(ns("pvalue_cutoff"),
                    "Max p-value:",
                    min = 0.0001,
                    max = 0.1,
                    value = 0.01,
                    step = 0.001)
      ),
      
      # Display options
      wellPanel(
        h5("Display Options"),
        selectInput(ns("heatmap_metric"),
                    "Heatmap Values:",
                    choices = c("log2FC" = "log2FC",
                                "-log10(p-value)" = "neg_log10_pval",
                                "Z-score" = "zscore"),
                    selected = "log2FC"),
        
        selectInput(ns("color_scheme"),
                    "Color Scheme:",
                    choices = c("Red-Blue (diverging)" = "RdBu",
                                "Blue-White-Red" = "BWR",
                                "Viridis" = "viridis",
                                "Plasma" = "plasma"),
                    selected = "RdBu"),
        
        checkboxInput(ns("cluster_genes"),
                      "Cluster Genes (rows)",
                      value = TRUE),
        
        checkboxInput(ns("cluster_conditions"),
                      "Cluster Conditions (columns)",
                      value = TRUE),
        
        checkboxInput(ns("show_gene_annotations"),
                      "Show Gene Annotations",
                      value = FALSE)
      ),
      
      # Performance options
      wellPanel(
        h5("Performance"),
        checkboxInput(ns("preview_mode"),
                      "Preview Mode (faster)",
                      value = TRUE),
        
        conditionalPanel(
          condition = paste0("input['", ns("preview_mode"), "'] == true"),
          helpText("Preview limits to 200 genes for faster rendering")
        ),
        
        sliderInput(ns("max_display_genes"),
                    "Max Genes to Display:",
                    min = 50,
                    max = 1000,
                    value = 200,
                    step = 50)
      ),
      
      # Data loading controls
      wellPanel(
        h5("Data Management"),
        actionButton(ns("load_data_btn"),
                     "Load DE Data",
                     class = "btn-info",
                     width = "100%"),
        br(), br(),
        actionButton(ns("refresh_data_btn"),
                     "Refresh Data",
                     class = "btn-warning",
                     width = "100%"),
        br(), br(),
        helpText("Click 'Load DE Data' first to load differential expression results. This may take 1-2 minutes for large datasets.")
      ),
      
      br(),
      actionButton(ns("generate_heatmap"),
                   "Generate Heatmap",
                   class = "btn-primary",
                   width = "100%"),
      
      br(), br(),
      downloadButton(ns("download_heatmap_png"),
                     "Download PNG",
                     class = "btn-success",
                     width = "100%"),
      
      br(),
      downloadButton(ns("download_heatmap_pdf"),
                     "Download PDF",
                     class = "btn-success",
                     width = "100%")
    ),
    
    mainPanel(
      width = 9,
      
      # Status panel
      fluidRow(
        column(12,
          wellPanel(
            style = "background-color: #f8f9fa; margin-bottom: 20px;",
            h4("Heatmap Status", style = "margin-top: 0;"),
            verbatimTextOutput(ns("heatmap_status"))
          )
        )
      ),
      
      # Main heatmap display
      tabsetPanel(
        id = ns("heatmap_tabs"),
        
        tabPanel("Interactive Heatmap",
                 value = "interactive",
                 br(),
                 conditionalPanel(
                   condition = "output.has_heatmaply == true",
                   ns = ns,
                   plotlyOutput(ns("interactive_heatmap"), height = "700px")
                 ),
                 conditionalPanel(
                   condition = "output.has_heatmaply == false",
                   ns = ns,
                   plotOutput(ns("static_heatmap"), height = "700px")
                 )),
        
        tabPanel("Data Table",
                 value = "table",
                 br(),
                 h4("Heatmap Data"),
                 DT::dataTableOutput(ns("heatmap_data_table")),
                 br(),
                 h4("Gene Statistics"),
                 DT::dataTableOutput(ns("gene_stats_table"))),
        
        tabPanel("Settings Summary",
                 value = "settings",
                 br(),
                 h4("Current Settings"),
                 verbatimTextOutput(ns("settings_summary")),
                 br(),
                 h4("Data Processing Log"),
                 verbatimTextOutput(ns("processing_log")))
      )
    )
  )
}

mod_de_heatmap_server <- function(id, app_data, global_selection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values
    heatmap_data <- reactiveValues(
      de_data = NULL,
      processed_data = NULL,
      heatmap_matrix = NULL,
      plot = NULL,
      processing_log = character(),
      de_loaded = FALSE,  # Flag to prevent repeated loading
      local_updating = FALSE  # Flag to prevent update loops
    )
    
    # Sync cluster selection with global settings
    observe({
      req(global_selection()$cluster)
      if (!heatmap_data$local_updating && input$cluster_select != global_selection()$cluster) {
        heatmap_data$local_updating <- TRUE
        updateSelectInput(session, "cluster_select", selected = global_selection()$cluster)
        heatmap_data$local_updating <- FALSE
      }
    })
    
    # Update global settings when local cluster changes
    observeEvent(input$cluster_select, {
      if (!heatmap_data$local_updating) {
        heatmap_data$local_updating <- TRUE
        # Send message to update global cluster
        session$sendInputMessage("update_cluster_from_module", 
                               list(value = input$cluster_select))
        heatmap_data$local_updating <- FALSE
      }
    })
    
    # Load DE data ONLY when user clicks "Load Data" button
    observeEvent(input$load_data_btn, {
      # Prevent multiple loading
      if (heatmap_data$de_loaded) {
        showNotification("Data already loaded. Use 'Refresh Data' if needed.", type = "message", duration = 3)
        return()
      }
      
      req(app_data$data_loaded)
      
      # Use the same DE file path logic as the existing DE Results module
      data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
      if (data_dir == "") {
        heatmap_data$processing_log <- c(heatmap_data$processing_log,
                                        "WARNING: No data directory configured")
        showNotification(
          "No data directory configured. DE heatmaps not available.",
          type = "warning"
        )
        return()
      }
      
      de_file_path <- file.path(data_dir, "full_DE_results.rds")
      
      if (!file.exists(de_file_path)) {
        heatmap_data$processing_log <- c(heatmap_data$processing_log,
                                        "WARNING: DE results file not found")
        showNotification(
          "DE results file not found in the configured data directory. DE heatmaps not available.",
          type = "warning"
        )
        return()
      }
      
      # Load and process actual DE results with performance optimization
      heatmap_data$processing_log <- c(heatmap_data$processing_log, "Starting DE data loading...")
      showNotification("Loading DE results for heatmaps (this may take 1-2 minutes)...", type = "message", duration = 10)
      
      tryCatch({
        de_results <- readRDS(de_file_path)
        heatmap_data$processing_log <- c(heatmap_data$processing_log, "DE file loaded, starting processing...")
        
        # Process only the selected cluster to improve performance
        # Add timeout to prevent hanging
        setTimeLimit(cpu = 30, elapsed = 60, transient = TRUE)  # 1 minute timeout
        
        on.exit(setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE))  # Reset timeout
        
        # Filter to only the selected cluster before processing
        selected_cluster <- input$cluster_select
        heatmap_data$processing_log <- c(heatmap_data$processing_log, 
                                        paste("Processing DE data for", selected_cluster, "only"))
        
        processed_de <- process_de_data_with_metadata(de_results, 
                                                     max_genes_per_condition = 50,
                                                     cluster_filter = selected_cluster)
        
        if (nrow(processed_de) == 0) {
          heatmap_data$processing_log <- c(heatmap_data$processing_log,
                                          "WARNING: No significant DE results found")
          showNotification("No significant DE results found in the data file", type = "warning")
          return()
        }
        
        heatmap_data$de_data <- processed_de
        heatmap_data$processing_log <- c(heatmap_data$processing_log, 
                                        paste("Loaded", nrow(processed_de), "significant DE gene records"))
        
        # Update gene filter choices
        available_genes <- unique(processed_de$gene)
        updateSelectizeInput(session, "gene_filter",
                            choices = available_genes,
                            selected = head(available_genes, 5))
        
        showNotification(
          paste("Loaded", nrow(processed_de), "significant DE genes from", length(available_genes), "genes/conditions"),
          type = "message",
          duration = 3
        )
        
        # Set flag to prevent repeated loading
        heatmap_data$de_loaded <- TRUE
        
      }, error = function(e) {
        heatmap_data$processing_log <- c(heatmap_data$processing_log,
                                        paste("ERROR loading DE results:", e$message))
        showNotification(
          paste("Error loading DE results:", e$message),
          type = "error",
          duration = 10
        )
      })
    })
    
    # Refresh data button - allows reloading
    observeEvent(input$refresh_data_btn, {
      # Reset the loaded flag and clear data
      heatmap_data$de_loaded <- FALSE
      heatmap_data$de_data <- NULL
      heatmap_data$processed_data <- NULL
      heatmap_data$heatmap_matrix <- NULL
      heatmap_data$plot <- NULL
      heatmap_data$processing_log <- character()
      
      showNotification("Data cleared. Click 'Load DE Data' to reload.", type = "message", duration = 5)
      
      # Clear gene filter choices
      updateSelectizeInput(session, "gene_filter", choices = NULL, selected = NULL)
    })
    
    # Process data for heatmap
    processed_heatmap_data <- reactive({
      req(heatmap_data$de_loaded)  # Only process when data is loaded
      req(heatmap_data$de_data)
      req(input$source_filter)
      
      data <- heatmap_data$de_data
      
      # Apply source filter
      filtered_data <- data %>%
        dplyr::filter(source %in% input$source_filter)
      
      # Apply cluster filter (always use the selected cluster)
      filtered_data <- filtered_data %>%
        dplyr::filter(cluster == input$cluster_select)
      
      # Apply gene filter if enabled
      if (input$filter_by_genes && length(input$gene_filter) > 0) {
        filtered_data <- filtered_data %>%
          dplyr::filter(gene_name %in% input$gene_filter)
      }
      
      # Apply significance cutoffs
      filtered_data <- filtered_data %>%
        dplyr::filter(
          abs(log2FC) >= input$log2fc_cutoff,
          pvalue <= input$pvalue_cutoff
        )
      
      # Gene selection
      if (input$gene_selection_mode == "top_n") {
        filtered_data <- filtered_data %>%
          dplyr::arrange(pvalue) %>%
          dplyr::group_by(gene, cluster, source) %>%
          dplyr::slice_head(n = input$top_n_genes) %>%
          dplyr::ungroup()
      } else if (input$gene_selection_mode == "custom" && input$custom_genes != "") {
        custom_gene_list <- trimws(strsplit(input$custom_genes, "\n")[[1]])
        custom_gene_list <- custom_gene_list[custom_gene_list != ""]
        filtered_data <- filtered_data %>%
          dplyr::filter(gene_name %in% custom_gene_list)
      }
      
      # Apply performance limits
      if (input$preview_mode && nrow(filtered_data) > input$max_display_genes) {
        filtered_data <- sample_de_for_preview(filtered_data, input$max_display_genes)
      }
      
      heatmap_data$processing_log <- c(heatmap_data$processing_log,
                                      paste("Processed to", nrow(filtered_data), "DE records"))
      
      return(filtered_data)
    })
    
    # Generate heatmap matrix
    observeEvent(input$generate_heatmap, {
      req(processed_heatmap_data())
      
      data <- processed_heatmap_data()
      
      if (nrow(data) == 0) {
        showNotification("No data available with current filters", type = "warning")
        return()
      }
      
      # Create condition labels with source distinction
      data$condition_label <- paste0(data$gene, " (", data$source, ") - ", data$cluster)
      
      # Prepare heatmap value
      if (input$heatmap_metric == "log2FC") {
        data$heatmap_value <- data$log2FC
      } else if (input$heatmap_metric == "neg_log10_pval") {
        data$heatmap_value <- -log10(pmax(data$pvalue, 1e-300))
      } else if (input$heatmap_metric == "zscore") {
        # Calculate z-score within each condition
        data <- data %>%
          dplyr::group_by(condition_label) %>%
          dplyr::mutate(heatmap_value = scale(log2FC)[,1]) %>%
          dplyr::ungroup()
      }
      
      # Create matrix for heatmap
      heatmap_matrix <- data %>%
        dplyr::select(gene_name, condition_label, heatmap_value) %>%
        tidyr::pivot_wider(names_from = condition_label, 
                          values_from = heatmap_value,
                          values_fill = 0,
                          values_fn = mean)
      
      # Convert to matrix
      mat <- as.matrix(heatmap_matrix[,-1])
      rownames(mat) <- heatmap_matrix$gene_name
      
      # Limit matrix size if needed
      if (nrow(mat) > 500) {
        # Keep top genes by variance
        gene_vars <- apply(mat, 1, var, na.rm = TRUE)
        top_genes <- order(gene_vars, decreasing = TRUE)[1:500]
        mat <- mat[top_genes, , drop = FALSE]
      }
      
      heatmap_data$heatmap_matrix <- mat
      heatmap_data$processed_data <- data
      
      showNotification("Heatmap generated successfully!", type = "message")
    })
    
    # Check if heatmaply is available
    output$has_heatmaply <- reactive({
      requireNamespace("heatmaply", quietly = TRUE) && !is.null(heatmap_data$heatmap_matrix)
    })
    outputOptions(output, "has_heatmaply", suspendWhenHidden = FALSE)
    
    # Interactive heatmap
    output$interactive_heatmap <- renderPlotly({
      req(heatmap_data$heatmap_matrix)
      
      if (requireNamespace("heatmaply", quietly = TRUE)) {
        mat <- heatmap_data$heatmap_matrix
        
        # Set color scale
        colors <- switch(input$color_scheme,
          "RdBu" = colorRampPalette(c("blue", "white", "red"))(256),
          "BWR" = colorRampPalette(c("blue", "white", "red"))(256),
          "viridis" = viridis::viridis(256),
          "plasma" = viridis::plasma(256),
          colorRampPalette(c("blue", "white", "red"))(256)
        )
        
        # Determine clustering
        dendrogram <- "none"
        if (input$cluster_genes && input$cluster_conditions) {
          dendrogram <- "both"
        } else if (input$cluster_genes) {
          dendrogram <- "row"
        } else if (input$cluster_conditions) {
          dendrogram <- "column"
        }
        
        # Create heatmap
        tryCatch({
          heatmaply::heatmaply(
            mat,
            dendrogram = dendrogram,
            colors = colors,
            main = paste("DE Gene Heatmap -", input$heatmap_metric),
            xlab = "Conditions (Gene Source - Cluster)",
            ylab = "Gene Names",
            showticklabels = c(TRUE, TRUE),
            margins = c(100, 200, 40, 20)
          )
        }, error = function(e) {
          # Fallback to simpler version
          heatmaply::heatmaply(
            mat,
            dendrogram = "none",
            colors = colors,
            main = "DE Gene Heatmap",
            showticklabels = c(FALSE, FALSE)
          )
        })
      }
    })
    
    # Static heatmap fallback
    output$static_heatmap <- renderPlot({
      req(heatmap_data$heatmap_matrix)
      
      mat <- heatmap_data$heatmap_matrix
      
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(
          mat,
          cluster_rows = input$cluster_genes,
          cluster_cols = input$cluster_conditions,
          main = paste("DE Gene Heatmap -", input$heatmap_metric),
          show_rownames = nrow(mat) <= 50,
          show_colnames = ncol(mat) <= 20,
          fontsize = 8
        )
      } else {
        # Base R heatmap
        heatmap(mat,
                main = paste("DE Gene Heatmap -", input$heatmap_metric),
                cexRow = 0.7,
                cexCol = 0.7)
      }
    })
    
    # Status output
    output$heatmap_status <- renderPrint({
      if (!heatmap_data$de_loaded) {
        cat("Status: Ready to load data\n")
        cat("=========================\n")
        cat("Click 'Load DE Data' to begin.\n")
        cat("This will load differential expression results from your dataset.\n")
        cat("Expected loading time: 1-2 minutes for large datasets.\n")
        return()
      }
      
      if (is.null(heatmap_data$de_data)) {
        cat("Status: Data loading in progress...\n")
        cat("===================================\n")
        cat("Please wait while DE results are processed.\n")
        return()
      }
      
      tryCatch({
        data <- processed_heatmap_data()
      }, error = function(e) {
        cat("Status: Waiting for data initialization\n")
        cat("=====================================\n")
        cat("Click 'Load DE Data' when ready.\n")
        return()
      })
      
      if (is.null(data)) {
        cat("Status: Data loaded, waiting for settings\n")
        cat("=========================================\n")
        cat("DE data loaded successfully!\n")
        cat("Configure your settings and click 'Generate Heatmap'.\n")
        return()
      }
      
      cat("Heatmap Status:\n")
      cat("==============\n")
      cat("DE records:", nrow(data), "\n")
      cat("Unique genes:", length(unique(data$gene_name)), "\n")
      cat("Conditions:", length(unique(paste(data$gene, data$source, data$cluster))), "\n")
      cat("Sources:", paste(unique(data$source), collapse = ", "), "\n")
      
      if (!is.null(heatmap_data$heatmap_matrix)) {
        mat <- heatmap_data$heatmap_matrix
        cat("\nHeatmap Matrix:\n")
        cat("Genes (rows):", nrow(mat), "\n")
        cat("Conditions (cols):", ncol(mat), "\n")
        cat("Preview mode:", input$preview_mode, "\n")
      }
    })
    
    # Data table
    output$heatmap_data_table <- DT::renderDataTable({
      req(processed_heatmap_data())
      
      display_data <- processed_heatmap_data() %>%
        dplyr::select(gene_name, gene, source, cluster, log2FC, pvalue) %>%
        dplyr::arrange(pvalue)
      
      DT::datatable(display_data,
                    options = list(pageLength = 20, scrollX = TRUE),
                    rownames = FALSE) %>%
        DT::formatSignif(columns = c("log2FC", "pvalue"), digits = 3)
    })
    
    # Settings summary
    output$settings_summary <- renderPrint({
      cat("Current Heatmap Settings:\n")
      cat("========================\n")
      cat("Cluster:", input$cluster_select, "\n")
      cat("Sources:", paste(input$source_filter, collapse = ", "), "\n")
      cat("Gene selection:", input$gene_selection_mode, "\n")
      cat("log2FC cutoff:", input$log2fc_cutoff, "\n")
      cat("P-value cutoff:", input$pvalue_cutoff, "\n")
      cat("Heatmap metric:", input$heatmap_metric, "\n")
      cat("Color scheme:", input$color_scheme, "\n")
      cat("Preview mode:", input$preview_mode, "\n")
      cat("Max display genes:", input$max_display_genes, "\n")
    })
    
    # Processing log
    output$processing_log <- renderPrint({
      if (length(heatmap_data$processing_log) > 0) {
        cat("Processing Log:\n")
        cat("===============\n")
        for (i in seq_along(heatmap_data$processing_log)) {
          cat(paste0("[", i, "] ", heatmap_data$processing_log[i], "\n"))
        }
      } else {
        cat("No processing steps logged yet.")
      }
    })
    
    return(reactive({ processed_heatmap_data() }))
  })
}