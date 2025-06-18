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

# Load shared functions
source("modules/shared/de_processing.R")
source("modules/shared/source_labeling.R")

mod_de_heatmap_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("DE Gene Heatmap Builder"),
      
      # Data scope
      wellPanel(
        h5("Data Scope"),
        selectInput(ns("data_scope"),
                    "Dataset Scope:",
                    choices = c("Full Dataset" = "full",
                                "Selected Clusters" = "clusters",
                                "Selected Genes" = "genes"),
                    selected = "clusters"),
        
        conditionalPanel(
          condition = paste0("input['", ns("data_scope"), "'] == 'clusters'"),
          checkboxGroupInput(ns("cluster_filter"),
                             "Clusters:",
                             choices = paste0("cluster_", 0:9),
                             selected = paste0("cluster_", 0:2))
        ),
        
        conditionalPanel(
          condition = paste0("input['", ns("data_scope"), "'] == 'genes'"),
          selectizeInput(ns("gene_filter"),
                         "Genes:",
                         choices = NULL,
                         multiple = TRUE,
                         options = list(placeholder = "Select genes...",
                                       maxItems = 20))
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

mod_de_heatmap_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values
    heatmap_data <- reactiveValues(
      de_data = NULL,
      processed_data = NULL,
      heatmap_matrix = NULL,
      plot = NULL,
      processing_log = character()
    )
    
    # Load DE data (mock for now)
    observe({
      req(app_data$data)
      
      # Create more comprehensive mock DE data
      genes <- c("LRRK2", "PINK1", "SNCA", "GBA", "PARK7", "ATP13A2", "DNAJC6", "FBXO7", "SYNJ1", "VPS13C")
      clusters <- paste0("cluster_", 0:4)
      sources <- c("iSCORE-PD", "CRISPRi", "CRISPRa")
      
      mock_de_data <- expand.grid(
        gene = genes,
        cluster = clusters, 
        source = sources,
        stringsAsFactors = FALSE
      )
      
      # Add random gene names and statistics
      n_genes_per_condition <- 100
      full_mock_data <- data.frame()
      
      for (i in 1:nrow(mock_de_data)) {
        condition_data <- data.frame(
          gene = mock_de_data$gene[i],
          cluster = mock_de_data$cluster[i],
          source = mock_de_data$source[i],
          gene_name = paste0("ENSG", sprintf("%08d", 1:n_genes_per_condition)),
          log2FC = rnorm(n_genes_per_condition, 0, 1.5),
          pvalue = runif(n_genes_per_condition, 0.0001, 0.1),
          method = ifelse(mock_de_data$source[i] == "iSCORE-PD", "MAST", "MixScale"),
          experiment = "default",
          stringsAsFactors = FALSE
        )
        full_mock_data <- rbind(full_mock_data, condition_data)
      }
      
      # Add source labels
      full_mock_data$source_label <- paste0(full_mock_data$gene, " (", full_mock_data$source, ")")
      
      heatmap_data$de_data <- full_mock_data
      heatmap_data$processing_log <- c(heatmap_data$processing_log, 
                                      paste("Loaded", nrow(full_mock_data), "DE gene records"))
      
      # Update gene filter choices
      available_genes <- unique(full_mock_data$gene)
      updateSelectizeInput(session, "gene_filter",
                          choices = available_genes,
                          selected = available_genes[1:5])
    })
    
    # Process data for heatmap
    processed_heatmap_data <- reactive({
      req(heatmap_data$de_data)
      req(input$source_filter)
      
      data <- heatmap_data$de_data
      
      # Apply source filter
      filtered_data <- data %>%
        dplyr::filter(source %in% input$source_filter)
      
      # Apply data scope filter
      if (input$data_scope == "clusters") {
        req(input$cluster_filter)
        filtered_data <- filtered_data %>%
          dplyr::filter(cluster %in% input$cluster_filter)
      } else if (input$data_scope == "genes") {
        req(input$gene_filter)
        filtered_data <- filtered_data %>%
          dplyr::filter(gene %in% input$gene_filter)
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
      data <- processed_heatmap_data()
      
      if (is.null(data)) {
        cat("No data loaded yet.")
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
      cat("Data scope:", input$data_scope, "\n")
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