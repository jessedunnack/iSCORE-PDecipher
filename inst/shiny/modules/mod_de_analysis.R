# Module: DE Gene Analysis
# Cross-condition comparisons, overlaps, and correlations

# Load required libraries
library(dplyr)
library(ggplot2)
library(plotly)

# Load optional libraries
venn_available <- requireNamespace("VennDiagram", quietly = TRUE)
if (venn_available) {
  library(VennDiagram)
}

upset_available <- requireNamespace("UpSetR", quietly = TRUE)
if (upset_available) {
  library(UpSetR)
}

# Load shared functions
source("modules/shared/de_processing.R")
source("modules/shared/source_labeling.R")

mod_de_analysis_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("DE Gene Analysis Settings"),
      
      # Dataset selection
      selectInput(ns("dataset_filter"),
                  "Dataset Filter:",
                  choices = c("All Datasets" = "all",
                              "iSCORE-PD Only" = "iSCORE-PD",
                              "PerturbSeq Only" = "PerturbSeq"),
                  selected = "all"),
      
      # Gene/mutation selection
      selectizeInput(ns("gene_selection"),
                     "Genes/Mutations:",
                     choices = NULL,
                     multiple = TRUE,
                     options = list(placeholder = "Select genes...",
                                   maxItems = 10)),
      
      # Cluster selection
      checkboxGroupInput(ns("cluster_selection"),
                         "Clusters:",
                         choices = paste0("cluster_", 0:13),
                         selected = paste0("cluster_", 0:2)),
      
      # Significance cutoffs
      wellPanel(
        h5("Significance Cutoffs"),
        sliderInput(ns("log2fc_threshold"),
                    "Min |log2FC|:",
                    min = 0,
                    max = 2,
                    value = 0.25,
                    step = 0.05),
        
        sliderInput(ns("pvalue_threshold"),
                    "Max p-value:",
                    min = 0.001,
                    max = 0.1,
                    value = 0.05,
                    step = 0.005)
      ),
      
      # Analysis type
      selectInput(ns("analysis_type"),
                  "Analysis Type:",
                  choices = c("Gene Overlaps" = "overlaps",
                              "Correlations" = "correlations",
                              "Scatter Plots" = "scatter",
                              "Top Genes" = "top_genes"),
                  selected = "overlaps"),
      
      br(),
      actionButton(ns("run_analysis"),
                   "Run Analysis",
                   class = "btn-primary",
                   width = "100%"),
      
      br(), br(),
      downloadButton(ns("download_results"),
                     "Download Results",
                     class = "btn-success",
                     width = "100%")
    ),
    
    mainPanel(
      width = 9,
      
      # Status and summary
      fluidRow(
        column(12,
          wellPanel(
            style = "background-color: #f8f9fa; margin-bottom: 20px;",
            h4("Analysis Summary", style = "margin-top: 0;"),
            verbatimTextOutput(ns("analysis_summary"))
          )
        )
      ),
      
      # Main visualization panel
      tabsetPanel(
        id = ns("viz_tabs"),
        
        tabPanel("Gene Overlaps",
                 value = "overlaps",
                 br(),
                 conditionalPanel(
                   condition = "output.has_upset == true",
                   ns = ns,
                   h4("UpSet Plot"),
                   plotOutput(ns("upset_plot"), height = "500px")
                 ),
                 h4("Venn Diagram"),
                 plotOutput(ns("venn_plot"), height = "400px"),
                 h4("Overlap Statistics"),
                 DT::dataTableOutput(ns("overlap_table"))),
        
        tabPanel("Correlations", 
                 value = "correlations",
                 br(),
                 h4("Cross-Condition log2FC Correlations"),
                 plotlyOutput(ns("correlation_plot"), height = "600px"),
                 br(),
                 h4("Correlation Matrix"),
                 DT::dataTableOutput(ns("correlation_table"))),
        
        tabPanel("Scatter Plots",
                 value = "scatter", 
                 br(),
                 fluidRow(
                   column(6,
                     selectInput(ns("scatter_x"),
                                 "X-axis Condition:",
                                 choices = NULL)
                   ),
                   column(6,
                     selectInput(ns("scatter_y"),
                                 "Y-axis Condition:",
                                 choices = NULL)
                   )
                 ),
                 plotlyOutput(ns("scatter_plot"), height = "600px")),
        
        tabPanel("Top Genes",
                 value = "top_genes",
                 br(),
                 fluidRow(
                   column(4,
                     numericInput(ns("top_n"),
                                  "Number of Top Genes:",
                                  value = 20,
                                  min = 5,
                                  max = 100,
                                  step = 5)
                   ),
                   column(4,
                     selectInput(ns("top_metric"),
                                 "Ranking Metric:",
                                 choices = c("P-value" = "pvalue",
                                           "|log2FC|" = "abs_log2fc",
                                           "Combined Score" = "combined"),
                                 selected = "pvalue")
                   ),
                   column(4,
                     selectInput(ns("top_direction"),
                                 "Direction:",
                                 choices = c("All" = "all",
                                           "Up-regulated" = "up",
                                           "Down-regulated" = "down"),
                                 selected = "all")
                   )
                 ),
                 plotlyOutput(ns("top_genes_plot"), height = "600px"))
      )
    )
  )
}

mod_de_analysis_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values
    analysis_data <- reactiveValues(
      de_data = NULL,
      processed_data = NULL,
      overlaps = NULL,
      correlations = NULL
    )
    
    # Load and process real DE data using same path as existing DE module
    observe({
      req(app_data$data)
      
      # Use the same DE file path logic as the existing DE Results module
      data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
      if (data_dir == "") {
        showNotification(
          "No data directory configured. DE analysis not available.",
          type = "warning"
        )
        return()
      }
      
      de_file_path <- file.path(data_dir, "full_DE_results.rds")
      
      if (!file.exists(de_file_path)) {
        showNotification(
          "DE results file not found in the configured data directory. DE analysis not available.",
          type = "warning"
        )
        return()
      }
      
      # Load and process actual DE results
      tryCatch({
        de_results <- readRDS(de_file_path)
        processed_de <- process_de_data_with_metadata(de_results)
        
        if (nrow(processed_de) == 0) {
          showNotification("No DE results found in the data file", type = "warning")
          return()
        }
        
        analysis_data$de_data <- processed_de
        
        # Update gene choices
        available_genes <- unique(processed_de$gene)
        updateSelectizeInput(session, "gene_selection",
                            choices = available_genes,
                            selected = head(available_genes, 3))
        
        showNotification(
          paste("Loaded", nrow(processed_de), "DE gene records from", length(available_genes), "genes"),
          type = "message"
        )
        
      }, error = function(e) {
        showNotification(
          paste("Error loading DE results:", e$message),
          type = "error",
          duration = NULL
        )
      })
    })
    
    # Process data based on user selections
    processed_data <- reactive({
      req(analysis_data$de_data)
      req(input$gene_selection)
      req(input$cluster_selection)
      
      data <- analysis_data$de_data
      
      # Apply filters
      filtered_data <- data %>%
        dplyr::filter(
          gene %in% input$gene_selection,
          cluster %in% input$cluster_selection,
          abs(log2FC) >= input$log2fc_threshold,
          pvalue <= input$pvalue_threshold
        )
      
      # Dataset filter
      if (input$dataset_filter != "all") {
        filtered_data <- filtered_data %>%
          dplyr::filter(source == input$dataset_filter)
      }
      
      return(filtered_data)
    })
    
    # Analysis summary
    output$analysis_summary <- renderPrint({
      data <- processed_data()
      if (is.null(data) || nrow(data) == 0) {
        cat("No data available with current filters.")
        return()
      }
      
      cat("Analysis Summary:\n")
      cat("================\n")
      cat("Total DE genes:", nrow(data), "\n")
      cat("Unique genes:", length(unique(data$gene_name)), "\n")
      cat("Conditions:", length(unique(paste(data$gene, data$cluster, data$source))), "\n")
      cat("Sources:", paste(unique(data$source), collapse = ", "), "\n")
      cat("Clusters:", paste(unique(data$cluster), collapse = ", "), "\n")
      cat("\nTop significant genes:\n")
      
      top_genes <- data %>%
        dplyr::arrange(pvalue) %>%
        dplyr::slice_head(n = 5) %>%
        dplyr::select(gene_name, gene, source, log2FC, pvalue)
      
      print(top_genes)
    })
    
    # Gene overlap analysis
    observe({
      req(processed_data())
      
      if (input$analysis_type == "overlaps") {
        data <- processed_data()
        
        # Create comparison groups
        comparison_groups <- list()
        conditions <- unique(paste(data$gene, data$cluster, data$source, sep = "_"))
        
        for (condition in conditions[1:min(5, length(conditions))]) {  # Limit to 5 for visualization
          parts <- strsplit(condition, "_")[[1]]
          if (length(parts) >= 3) {
            group_name <- paste(parts[1], parts[3], sep = " ")
            comparison_groups[[group_name]] <- list(
              gene = parts[1],
              cluster = paste(parts[2:3], collapse = "_"),
              source = parts[length(parts)]
            )
          }
        }
        
        overlaps <- calculate_de_overlaps(data, comparison_groups, 
                                         list(log2FC = input$log2fc_threshold, 
                                              pvalue = input$pvalue_threshold))
        analysis_data$overlaps <- overlaps
      }
    })
    
    # Check if UpSetR is available
    output$has_upset <- reactive({
      upset_available && !is.null(analysis_data$overlaps)
    })
    outputOptions(output, "has_upset", suspendWhenHidden = FALSE)
    
    # UpSet plot
    output$upset_plot <- renderPlot({
      req(analysis_data$overlaps)
      
      if (upset_available) {
        gene_lists <- analysis_data$overlaps$gene_lists
        if (length(gene_lists) >= 2) {
          tryCatch({
            UpSetR::upset(UpSetR::fromList(gene_lists), nsets = length(gene_lists))
          }, error = function(e) {
            plot.new()
            text(0.5, 0.5, "Error generating UpSet plot", cex = 1.5)
          })
        }
      } else {
        plot.new()
        text(0.5, 0.5, "UpSetR package not available.\nInstall with: install.packages('UpSetR')", cex = 1.2)
      }
    })
    
    # Venn diagram
    output$venn_plot <- renderPlot({
      req(analysis_data$overlaps)
      
      gene_lists <- analysis_data$overlaps$gene_lists
      if (length(gene_lists) >= 2) {
        if (!venn_available) {
          plot.new()
          text(0.5, 0.5, "VennDiagram package not available.\nInstall with: install.packages('VennDiagram')", cex = 1.2)
          return()
        }
        
        tryCatch({
          if (length(gene_lists) == 2) {
            venn.plot <- VennDiagram::venn.diagram(
              x = gene_lists,
              category.names = names(gene_lists),
              filename = NULL,
              output = TRUE,
              col = "transparent",
              fill = c("lightblue", "lightcoral"),
              alpha = 0.50,
              cex = 1.5,
              fontfamily = "serif",
              cat.cex = 1.5,
              cat.fontfamily = "serif"
            )
            grid::grid.draw(venn.plot)
          } else if (length(gene_lists) == 3) {
            venn.plot <- VennDiagram::venn.diagram(
              x = gene_lists,
              category.names = names(gene_lists),
              filename = NULL,
              output = TRUE,
              col = "transparent",
              fill = c("lightblue", "lightcoral", "lightgreen"),
              alpha = 0.50,
              cex = 1.5,
              fontfamily = "serif",
              cat.cex = 1.5,
              cat.fontfamily = "serif"
            )
            grid::grid.draw(venn.plot)
          } else {
            plot.new()
            text(0.5, 0.5, "Venn diagrams support 2-3 sets only.\nSee UpSet plot for more comparisons.", cex = 1.2)
          }
        }, error = function(e) {
          plot.new()
          text(0.5, 0.5, paste("Error generating Venn diagram:", e$message), cex = 1)
        })
      } else {
        plot.new()
        text(0.5, 0.5, "Need at least 2 conditions for overlap analysis", cex = 1.5)
      }
    })
    
    # Overlap table
    output$overlap_table <- DT::renderDataTable({
      req(analysis_data$overlaps)
      
      overlaps <- analysis_data$overlaps
      summary_data <- overlaps$summary
      
      # Add intersection counts
      if (length(overlaps$intersections) > 0) {
        intersection_df <- data.frame(
          group = names(overlaps$intersections),
          gene_count = sapply(overlaps$intersections, length),
          type = "Intersection",
          stringsAsFactors = FALSE
        )
        
        summary_data$type <- "Individual"
        display_data <- rbind(summary_data, intersection_df)
      } else {
        summary_data$type <- "Individual"
        display_data <- summary_data
      }
      
      DT::datatable(display_data,
                    options = list(pageLength = 15),
                    rownames = FALSE) %>%
        DT::formatStyle("type",
                       backgroundColor = DT::styleEqual("Intersection", "#e8f4fd"))
    })
    
    # Return reactive data for potential use by other modules
    return(reactive({ processed_data() }))
  })
}