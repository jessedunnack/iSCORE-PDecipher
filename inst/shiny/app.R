# iSCORE-PDecipher Shiny Visualization App
# A ShinyGO-inspired interface for Parkinson's Disease functional enrichment analysis
# Combines reference app design with flexible data loading

# Load required libraries
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(plotly)
library(dplyr)
library(ggplot2)
library(glue)
library(tibble)
library(colourpicker)

# Check if we have data file provided or need upload interface
has_data <- Sys.getenv("ISCORE_HAS_DATA", unset = "FALSE") == "TRUE"
data_file <- Sys.getenv("ISCORE_DATA_FILE", unset = "")

# Source global functions and configurations
if (file.exists("global_minimal.R")) {
  source("global_minimal.R")
} else {
  # Create minimal config if file doesn't exist
  APP_CONFIG <- list(
    analysis_types = c("MAST", "MixScale", "MixScale_CRISPRa"),
    enrichment_types = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA")
  )
}

# Source modules if they exist
safe_source <- function(file) {
  if (file.exists(file)) {
    tryCatch(source(file), error = function(e) {
      cat("Warning: Could not source", file, ":", e$message, "\n")
    })
  }
}

safe_source("R/startup_manager.R")
safe_source("R/cache_manager.R")
safe_source("modules/mod_landing_page.R")
safe_source("modules/mod_precomputed_reactive.R")
safe_source("modules/mod_visualization.R")
safe_source("modules/mod_comparison.R")
safe_source("modules/mod_heatmap.R")
safe_source("modules/mod_pathview.R")
safe_source("modules/mod_export.R")

# Define UI with persistent sidebar
ui <- fluidPage(
  title = "iSCORE-PDecipher",
  
  # Header content (shinyjs and CSS)
  useShinyjs(),
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "custom.css"),
    tags$link(rel = "stylesheet", type = "text/css", href = "custom_enhanced.css"),
    tags$script(src = "app.js"),
    tags$script(HTML("
      // Show/hide loading overlay
      Shiny.addCustomMessageHandler('showLoading', function(message) {
        document.getElementById('loadingOverlay').classList.add('active');
      });
      
      Shiny.addCustomMessageHandler('hideLoading', function(message) {
        document.getElementById('loadingOverlay').classList.remove('active');
      });
      
      // Monitor Shiny busy state
      $(document).on('shiny:busy', function(event) {
        setTimeout(function() {
          if ($('.shiny-busy').length > 0) {
            $('#loadingOverlay').addClass('active');
          }
        }, 200); // Small delay to avoid flashing on quick operations
      });
      
      $(document).on('shiny:idle', function(event) {
        $('#loadingOverlay').removeClass('active');
      });
    ")),
    tags$style(HTML("
      /* Sidebar styling */
      .sidebar-fixed {
        position: fixed;
        top: 60px;
        bottom: 0;
        left: 0;
        width: 300px;
        padding: 20px;
        background-color: #f8f9fa;
        border-right: 1px solid #dee2e6;
        overflow-y: auto;
        z-index: 1000;
      }
      
      /* Main content offset */
      .main-content {
        margin-left: 320px;
        padding: 20px;
        position: relative;
        min-height: calc(100vh - 60px);
      }
      
      /* Loading overlay */
      .loading-overlay {
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background-color: rgba(255, 255, 255, 0.9);
        z-index: 999;
        display: none;
        align-items: center;
        justify-content: center;
      }
      
      .loading-overlay.active {
        display: flex;
      }
      
      .loading-content {
        text-align: center;
        padding: 40px;
        background: white;
        border-radius: 8px;
        box-shadow: 0 2px 10px rgba(0,0,0,0.1);
      }
      
      .loading-spinner {
        border: 4px solid #f3f3f3;
        border-top: 4px solid #3c8dbc;
        border-radius: 50%;
        width: 50px;
        height: 50px;
        animation: spin 1s linear infinite;
        margin: 0 auto 20px;
      }
      
      @keyframes spin {
        0% { transform: rotate(0deg); }
        100% { transform: rotate(360deg); }
      }
      
      /* Header styling */
      .app-header {
        position: fixed;
        top: 0;
        left: 0;
        right: 0;
        height: 50px;
        background-color: #3c8dbc;
        color: white;
        padding: 10px 20px;
        z-index: 1001;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      .app-header h2 {
        margin: 0;
        font-size: 24px;
        font-weight: 500;
      }
      
      /* Content area adjustment */
      .content-area {
        margin-top: 60px;
      }
      
      /* Tab styling */
      .nav-tabs {
        border-bottom: 2px solid #dee2e6;
        margin-bottom: 20px;
      }
      
      .nav-tabs .nav-link {
        color: #495057;
        border: none;
        border-bottom: 3px solid transparent;
        padding: 10px 20px;
      }
      
      .nav-tabs .nav-link:hover {
        border-color: transparent;
        border-bottom-color: #dee2e6;
      }
      
      .nav-tabs .nav-link.active {
        color: #3c8dbc;
        background-color: transparent;
        border-color: transparent;
        border-bottom-color: #3c8dbc;
        font-weight: 500;
      }
      
      /* Upload area styling */
      .upload-area {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 40px;
        border-radius: 15px;
        text-align: center;
        margin: 20px auto;
        max-width: 600px;
        box-shadow: 0 10px 30px rgba(0,0,0,0.2);
      }
      
      .upload-box {
        border: 3px dashed rgba(255,255,255,0.5);
        border-radius: 10px;
        padding: 30px;
        transition: all 0.3s ease;
      }
      
      .upload-box:hover {
        border-color: rgba(255,255,255,0.8);
        background-color: rgba(255,255,255,0.1);
      }
      
      /* Responsive */
      @media (max-width: 768px) {
        .sidebar-fixed {
          position: relative;
          width: 100%;
          top: auto;
          border-right: none;
          border-bottom: 1px solid #dee2e6;
        }
        
        .main-content {
          margin-left: 0;
        }
      }
    "))
  ),
  
  # Loading overlay
  div(id = "loadingOverlay", class = "loading-overlay",
    div(class = "loading-content",
      div(class = "loading-spinner"),
      h4("Loading data..."),
      p("Please wait while we process your enrichment data.")
    )
  ),
  
  # App header
  div(class = "app-header",
    h2(icon("dna"), "iSCORE-PDecipher: PD Enrichment Explorer")
  ),
  
  div(class = "content-area",
    # Conditional UI based on whether data is loaded
    conditionalPanel(
      condition = "!output.dataLoaded",
      
      div(class = "upload-area",
        h2(icon("upload"), "Load Your Data"),
        p("Upload your consolidated enrichment results to begin analysis"),
        
        div(class = "upload-box",
          fileInput("dataFile", 
                    label = NULL,
                    accept = c(".rds", ".RDS"),
                    placeholder = "Choose RDS file...",
                    buttonLabel = "Browse Files",
                    multiple = FALSE),
          
          br(),
          p(icon("info-circle"), "Accepts RDS files with enrichment analysis results", 
            style = "font-size: 14px; opacity: 0.8;")
        ),
        
        # Loading indicator
        conditionalPanel(
          condition = "output.uploading",
          div(
            style = "margin-top: 20px;",
            div(class = "loading-spinner", style = "width: 30px; height: 30px;"),
            "Processing your data..."
          )
        )
      )
    ),
    
    # Main interface (shown when data is loaded)
    conditionalPanel(
      condition = "output.dataLoaded",
      
      # Sidebar with global settings
      div(class = "sidebar-fixed",
        h3(icon("cogs"), "Global Settings"),
        hr(),
        
        # Data overview
        div(id = "dataOverview",
          h5("Dataset Overview", style = "margin-bottom: 10px; font-weight: bold; color: #337ab7;"),
          verbatimTextOutput("dataOverviewText", placeholder = TRUE),
          hr()
        ),
        
        # Data Selection
        h5("Data Selection", style = "margin-bottom: 15px; font-weight: bold; color: #337ab7;"),
        
        selectInput("global_analysis_type",
                    "Analysis Type",
                    choices = character(0),
                    width = "100%"),
        
        selectInput("global_gene",
                    "Gene/Mutation",
                    choices = character(0),
                    width = "100%"),
        
        selectInput("global_cluster",
                    "Cell Cluster", 
                    choices = character(0),
                    width = "100%"),
        
        selectInput("global_experiment",
                    "Experiment",
                    choices = character(0),
                    width = "100%"),
        
        selectInput("global_enrichment_type",
                    "Enrichment Database",
                    choices = character(0),
                    width = "100%"),
        
        # P-value threshold
        hr(),
        numericInput("global_pval",
                     "P-value Threshold",
                     value = 0.05,
                     min = 0.001,
                     max = 0.5,
                     step = 0.001,
                     width = "100%"),
        
        br(),
        actionButton("resetToDefaults",
                     "Reset to Defaults",
                     class = "btn-warning",
                     style = "width: 100%;")
      ),
      
      # Main content area
      div(class = "main-content",
        tabsetPanel(
          id = "main_tabs",
          
          # Overview/Landing Page
          tabPanel(
            "Overview",
            icon = icon("home"),
            value = "overview",
            br(),
            div(
              h2("Welcome to iSCORE-PDecipher"),
              p("Explore Parkinson's disease functional enrichment results from MAST and MixScale analyses."),
              
              fluidRow(
                column(6,
                  wellPanel(
                    h4(icon("chart-bar"), "Basic Visualization"),
                    p("Create dotplots and other standard visualizations"),
                    actionButton("goToVisualization", "Start Exploring", class = "btn-primary")
                  )
                ),
                column(6,
                  wellPanel(
                    h4(icon("th"), "Heatmap Analysis"),
                    p("Generate clustered heatmaps for pathway comparison"),
                    actionButton("goToHeatmap", "Create Heatmaps", class = "btn-info")
                  )
                )
              ),
              
              fluidRow(
                column(6,
                  wellPanel(
                    h4(icon("chart-line"), "Method Comparison"),
                    p("Compare MAST vs MixScale enrichment results"),
                    actionButton("goToComparison", "Compare Methods", class = "btn-success")
                  )
                ),
                column(6,
                  wellPanel(
                    h4(icon("download"), "Export Results"),
                    p("Download figures and filtered datasets"),
                    actionButton("goToExport", "Export Data", class = "btn-warning")
                  )
                )
              )
            )
          ),
          
          # Basic Visualization (main plotting interface)
          tabPanel(
            "Dotplot Visualization",
            icon = icon("chart-bar"),
            value = "visualization",
            br(),
            
            # Simple version of visualization since modules might not load
            sidebarLayout(
              sidebarPanel(
                width = 3,
                h4("Plot Settings"),
                
                selectInput("plot_type",
                            "Visualization Type:",
                            choices = list("Dot Plot" = "dotplot"),
                            selected = "dotplot"),
                
                numericInput("top_terms",
                             "Number of Terms:",
                             value = 20,
                             min = 5,
                             max = 50,
                             step = 5),
                
                selectInput("direction_filter",
                            "Direction:",
                            choices = c("All" = "ALL", "Up" = "UP", "Down" = "DOWN"),
                            selected = "ALL"),
                
                br(),
                actionButton("updatePlot", "Update Plot", class = "btn-primary", width = "100%")
              ),
              
              mainPanel(
                width = 9,
                tabsetPanel(
                  tabPanel("Plot",
                           withSpinner(plotlyOutput("mainPlot", height = "700px"))),
                  tabPanel("Data Table",
                           withSpinner(DT::dataTableOutput("plotData"))),
                  tabPanel("Summary",
                           verbatimTextOutput("plotSummary"))
                )
              )
            )
          ),
          
          # Heatmap Visualization  
          tabPanel(
            "Clustered Heatmaps",
            icon = icon("th"),
            value = "heatmap",
            br(),
            
            sidebarLayout(
              sidebarPanel(
                width = 3,
                h4("Heatmap Settings"),
                
                selectInput("heatmap_metric",
                            "Color Metric:",
                            choices = c("P-value" = "pvalue", "Fold Enrichment" = "foldenrich"),
                            selected = "pvalue"),
                
                selectInput("heatmap_genes",
                            "Genes to Include:",
                            choices = character(0),
                            multiple = TRUE),
                
                numericInput("heatmap_max_terms",
                             "Max Terms per Gene:",
                             value = 20,
                             min = 5,
                             max = 50),
                
                checkboxInput("cluster_rows", "Cluster Rows", value = TRUE),
                checkboxInput("cluster_cols", "Cluster Columns", value = TRUE),
                
                br(),
                actionButton("generateHeatmap", "Generate Heatmap", class = "btn-primary", width = "100%")
              ),
              
              mainPanel(
                width = 9,
                withSpinner(plotOutput("heatmapPlot", height = "800px"))
              )
            )
          ),
          
          # Export Results
          tabPanel(
            "Export & Download",
            icon = icon("download"),
            value = "export",
            br(),
            
            fluidRow(
              column(6,
                wellPanel(
                  h4("Download Filtered Data"),
                  p("Export the currently filtered dataset"),
                  downloadButton("downloadFilteredData", "Download CSV", class = "btn-success", width = "100%")
                )
              ),
              column(6,
                wellPanel(
                  h4("Download Plots"),
                  p("Save current visualizations"),
                  downloadButton("downloadPlot", "Download Plot (PNG)", class = "btn-info", width = "100%")
                )
              )
            )
          )
        )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    full_data = NULL,
    filtered_data = NULL,
    data_loaded = FALSE,
    current_plot = NULL
  )
  
  # Check if data was provided via launch function
  observe({
    if (has_data && file.exists(data_file)) {
      cat("Loading provided data file:", data_file, "\n")
      
      tryCatch({
        values$full_data <- readRDS(data_file)
        values$data_loaded <- TRUE
        
        # Initialize UI choices
        initialize_ui_choices()
        
        showNotification("Data loaded successfully!", type = "success")
        
      }, error = function(e) {
        showNotification(paste("Error loading data:", e$message), type = "error")
      })
    }
  })
  
  # Handle file upload
  observeEvent(input$dataFile, {
    req(input$dataFile)
    
    output$uploading <- reactive({ TRUE })
    outputOptions(output, "uploading", suspendWhenHidden = FALSE)
    
    tryCatch({
      values$full_data <- readRDS(input$dataFile$datapath)
      values$data_loaded <- TRUE
      
      # Initialize UI choices
      initialize_ui_choices()
      
      showNotification("Data uploaded successfully!", type = "success")
      output$uploading <- reactive({ FALSE })
      
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
      output$uploading <- reactive({ FALSE })
    })
  })
  
  # Function to initialize UI choices
  initialize_ui_choices <- function() {
    req(values$full_data)
    
    # Fix column names if needed
    if ("mutation_perturbation" %in% names(values$full_data) && !"gene" %in% names(values$full_data)) {
      values$full_data$gene <- values$full_data$mutation_perturbation
    }
    
    # Update analysis type choices
    if ("method" %in% names(values$full_data)) {
      analysis_types <- sort(unique(values$full_data$method))
      updateSelectInput(session, "global_analysis_type", choices = analysis_types, selected = analysis_types[1])
    }
    
    # Update enrichment type choices
    if ("enrichment_type" %in% names(values$full_data)) {
      enrich_types <- sort(unique(values$full_data$enrichment_type))
      updateSelectInput(session, "global_enrichment_type", choices = enrich_types, selected = enrich_types[1])
    }
  }
  
  # Update gene choices based on analysis type
  observe({
    req(values$full_data, input$global_analysis_type)
    
    genes <- unique(values$full_data[values$full_data$method == input$global_analysis_type, "gene"])
    genes <- sort(genes)
    
    updateSelectInput(session, "global_gene", choices = genes, selected = genes[1])
    updateSelectInput(session, "heatmap_genes", choices = genes, selected = head(genes, 3))
  })
  
  # Update cluster choices
  observe({
    req(values$full_data, input$global_analysis_type, input$global_gene)
    
    clusters <- unique(values$full_data[
      values$full_data$method == input$global_analysis_type &
      values$full_data$gene == input$global_gene,
      "cluster"
    ])
    clusters <- sort(clusters)
    
    updateSelectInput(session, "global_cluster", choices = clusters, selected = clusters[1])
  })
  
  # Output for conditional UI
  output$dataLoaded <- reactive({
    return(values$data_loaded)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  output$uploading <- reactive({ FALSE })
  outputOptions(output, "uploading", suspendWhenHidden = FALSE)
  
  # Data overview
  output$dataOverviewText <- renderText({
    req(values$full_data)
    
    paste(
      paste("Total terms:", format(nrow(values$full_data), big.mark = ",")),
      paste("Genes:", length(unique(values$full_data$gene))),
      paste("Methods:", length(unique(values$full_data$method))),
      paste("Enrichment types:", length(unique(values$full_data$enrichment_type))),
      sep = "\n"
    )
  })
  
  # Filter data for current plot
  filtered_plot_data <- reactive({
    req(values$full_data)
    req(input$global_analysis_type, input$global_gene, input$global_enrichment_type)
    
    data <- values$full_data
    
    # Apply filters
    data <- data[data$method == input$global_analysis_type, ]
    data <- data[data$gene == input$global_gene, ]
    data <- data[data$enrichment_type == input$global_enrichment_type, ]
    
    if (input$direction_filter != "ALL") {
      data <- data[data$direction == input$direction_filter, ]
    }
    
    # Filter by p-value
    if ("p.adjust" %in% names(data)) {
      data <- data[data$p.adjust <= input$global_pval, ]
    }
    
    # Take top terms
    if (nrow(data) > input$top_terms) {
      data <- data[order(data$p.adjust), ][1:input$top_terms, ]
    }
    
    return(data)
  })
  
  # Main dotplot
  output$mainPlot <- renderPlotly({
    req(filtered_plot_data())
    
    data <- filtered_plot_data()
    
    if (nrow(data) == 0) {
      return(plotly_empty("No data available for current selection"))
    }
    
    # Create dotplot
    p <- ggplot(data, aes(x = -log10(p.adjust), y = reorder(Description, -log10(p.adjust)))) +
      geom_point(aes(size = ifelse("Count" %in% names(data), Count, 10),
                     color = -log10(p.adjust)), alpha = 0.7) +
      scale_color_viridis_c(name = "-log10(p.adj)") +
      scale_size_continuous(name = "Gene Count", range = c(3, 10)) +
      theme_minimal() +
      theme(
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "bottom"
      ) +
      labs(
        title = paste(input$global_gene, "-", input$global_enrichment_type, "Enrichment"),
        x = "-log10(adjusted p-value)",
        y = "Pathway"
      )
    
    ggplotly(p, tooltip = c("x", "y", "size", "color"))
  })
  
  # Heatmap
  output$heatmapPlot <- renderPlot({
    req(input$generateHeatmap)
    req(values$full_data, input$heatmap_genes)
    
    isolate({
      data <- values$full_data
      
      # Filter data
      data <- data[data$gene %in% input$heatmap_genes, ]
      data <- data[data$enrichment_type == input$global_enrichment_type, ]
      data <- data[data$p.adjust <= input$global_pval, ]
      
      if (nrow(data) == 0) {
        plot.new()
        text(0.5, 0.5, "No data available for current selection", cex = 1.5)
        return()
      }
      
      # Get top terms per gene
      top_terms <- data %>%
        group_by(gene) %>%
        slice_min(p.adjust, n = input$heatmap_max_terms) %>%
        ungroup()
      
      # Create matrix for heatmap
      if (input$heatmap_metric == "pvalue") {
        plot_data <- top_terms %>%
          select(gene, Description, p.adjust) %>%
          mutate(value = -log10(p.adjust)) %>%
          select(-p.adjust)
      } else {
        # Use fold enrichment if available
        if ("FoldEnrichment" %in% names(top_terms)) {
          plot_data <- top_terms %>%
            select(gene, Description, FoldEnrichment) %>%
            rename(value = FoldEnrichment)
        } else {
          plot_data <- top_terms %>%
            select(gene, Description, p.adjust) %>%
            mutate(value = -log10(p.adjust)) %>%
            select(-p.adjust)
        }
      }
      
      # Pivot to matrix
      mat <- plot_data %>%
        pivot_wider(names_from = gene, values_from = value, values_fill = 0) %>%
        column_to_rownames("Description") %>%
        as.matrix()
      
      if (ncol(mat) > 1 && nrow(mat) > 1) {
        pheatmap::pheatmap(
          mat,
          clustering_distance_rows = "euclidean",
          clustering_distance_cols = "euclidean",
          cluster_rows = input$cluster_rows,
          cluster_cols = input$cluster_cols,
          color = viridis::viridis(100),
          fontsize_row = 8,
          fontsize_col = 10,
          main = paste("Clustered Heatmap -", input$global_enrichment_type)
        )
      } else {
        plot.new()
        text(0.5, 0.5, "Need multiple genes and terms for clustering", cex = 1.2)
      }
    })
  })
  
  # Data table
  output$plotData <- DT::renderDataTable({
    req(filtered_plot_data())
    
    display_cols <- c("gene", "enrichment_type", "Description", "p.adjust")
    available_cols <- intersect(display_cols, names(filtered_plot_data()))
    
    DT::datatable(
      filtered_plot_data()[, available_cols, drop = FALSE],
      options = list(pageLength = 15, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      DT::formatSignif(columns = "p.adjust", digits = 3)
  })
  
  # Navigation buttons
  observeEvent(input$goToVisualization, {
    updateTabsetPanel(session, "main_tabs", selected = "visualization")
  })
  
  observeEvent(input$goToHeatmap, {
    updateTabsetPanel(session, "main_tabs", selected = "heatmap")
  })
  
  # Download handlers
  output$downloadFilteredData <- downloadHandler(
    filename = function() {
      paste0("iscore_enrichment_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_plot_data(), file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)