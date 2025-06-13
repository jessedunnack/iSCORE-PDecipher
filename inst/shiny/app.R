# iSCORE-PDecipher Shiny Visualization App
# A ShinyGO-inspired interface for Parkinson's Disease functional enrichment analysis

# Note: global.R is automatically loaded by Shiny before app.R
# This includes all required library loading

# Ensure critical libraries are loaded (fallback)
critical_packages <- c("shinyjs", "shinycssloaders", "shinyWidgets", "DT", "plotly")

for (pkg in critical_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(pkg, "package is required but not installed"))
  }
  if (!paste0("package:", pkg) %in% search()) {
    library(pkg, character.only = TRUE)
  }
}

# Ensure APP_CONFIG is available from global.R
if (!exists("APP_CONFIG")) {
  # If global.R wasn't loaded properly, create minimal config
  APP_CONFIG <- list(
    analysis_types = c("MAST", "MixScale"),
    enrichment_types = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
    directions = c("ALL", "UP", "DOWN")
  )
}

# Source core data management (NEW - centralized data loading)
source("R/data_manager.R")

# Source startup manager (NO immediate execution now)
source("R/startup_manager.R")

# Source cache manager
source("R/cache_manager.R")

# Source modules
# source("modules/mod_landing_page.R")  # Original landing page
# source("modules/mod_landing_page_with_umap.R")  # Enhanced with UMAP viewer v1
source("modules/mod_landing_page_with_umap_v2.R")  # Compact layout with auto-dataset selection
source("modules/mod_precomputed_reactive.R")
# source("modules/mod_visualization.R")  # Original version
source("modules/mod_visualization_enhanced.R")  # Enhanced with GSEA support
source("modules/mod_comparison.R")
source("modules/mod_heatmap_unified.R")
source("modules/mod_pathview.R")
source("modules/mod_export.R")

# Ensure UI functions are available (fallback assignments)
ui_functions <- list(
  useShinyjs = "shinyjs",
  withSpinner = "shinycssloaders",
  pickerInput = "shinyWidgets",
  dataTableOutput = "DT",
  plotlyOutput = "plotly"
)

for (func_name in names(ui_functions)) {
  if (!exists(func_name)) {
    pkg_name <- ui_functions[[func_name]]
    if (requireNamespace(pkg_name, quietly = TRUE)) {
      assign(func_name, get(func_name, envir = asNamespace(pkg_name)))
    } else {
      stop(paste("Function", func_name, "not available from package", pkg_name))
    }
  }
}

# Define UI with persistent sidebar
ui <- fluidPage(
  title = "PD Enrichment Explorer",
  
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
  
  # App header
  div(class = "app-header",
    h2(icon("dna"), "PD Enrichment Explorer")
  ),
  
  div(class = "content-area",
    # Sidebar with global settings
    div(class = "sidebar-fixed",
      h3(icon("cogs"), "Global Settings"),
      hr(),
      
      # Data Selection
      h5("Data Selection", style = "margin-bottom: 15px; font-weight: bold; color: #337ab7;"),
      
      selectInput("global_analysis_type",
                  "Analysis Type",
                  choices = c("MAST", "MixScale"),
                  selected = "MAST",
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
                  choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
                  selected = "GO_BP",
                  width = "100%"),
      
      # GSEA-specific options
      conditionalPanel(
        condition = "input.global_enrichment_type == 'GSEA'",
        wellPanel(
          style = "background-color: #fff; border: 1px solid #ddd; margin-top: 10px;",
          h5("GSEA Options", style = "margin-top: 0;"),
          checkboxGroupInput(
            "global_gsea_databases",
            "Gene Set Collections:",
            choices = c("C1: Positional" = "c1.all",
                       "C2: Curated (CP)" = "c2.cp",
                       "C5: GO" = "c5.all",
                       "H: Hallmark" = "h.all"),
            selected = c("c2.cp", "h.all")
          ),
          numericInput(
            "global_min_nes",
            "Min |NES|:",
            value = 1.0,
            min = 0,
            max = 3,
            step = 0.1,
            width = "100%"
          )
        )
      ),
      
      radioButtons("global_direction",
                   "Gene Regulation",
                   choices = c("ALL", "UP", "DOWN"),
                   selected = "ALL"),
      
      hr(),
      
      h5("Analysis Settings", style = "margin-bottom: 15px; font-weight: bold; color: #337ab7;"),
      
      sliderInput("global_pval",
                  "P-value Threshold",
                  min = 0.001,
                  max = 0.05,
                  value = 0.05,
                  step = 0.001,
                  width = "100%"),
      
      hr(),
      
      # Data status
      uiOutput("data_status")
    ),
      
    # Main content area
    div(class = "main-content",
      # Loading overlay
      div(id = "loadingOverlay", class = "loading-overlay",
        div(class = "loading-content",
          div(class = "loading-spinner"),
          h4("Loading data..."),
          p("Please wait while we fetch your results")
        )
      ),
      
      tabsetPanel(
        id = "main_tabs",
        
        # Landing Page
        tabPanel(
          "Overview",
          icon = icon("home"),
          value = "overview",
          br(),
          landingPageWithUmapUI("landing")
        ),
        
        # Basic Visualization (main plotting interface)
        tabPanel(
          "Basic Visualization",
          icon = icon("chart-bar"),
          value = "visualization",
          br(),
          mod_visualization_ui("visualization_module")
        ),
        
        # Method Comparison
        tabPanel(
          "Method Comparison",
          icon = icon("chart-line"),
          value = "comparison",
          br(),
          h2("Compare Analysis Methods"),
          p("Compare MAST vs MixScale results to identify convergent pathways"),
          mod_comparison_ui("comparison_module")
        ),
        
        # Heatmap Visualization
        tabPanel(
          "Heatmap Visualization",
          icon = icon("th"),
          value = "heatmap",
          br(),
          h2("Interactive Heatmaps"),
          p("Create customizable heatmaps of enrichment results"),
          mod_heatmap_unified_ui("heatmap_module")
        ),
        
        # KEGG Pathview
        tabPanel(
          "KEGG Pathview",
          icon = icon("project-diagram"),
          value = "pathview",
          br(),
          h2("KEGG Pathway Visualization"),
          p("Visualize differential expression on KEGG pathway diagrams"),
          mod_pathview_ui("pathview_module")
        ),
        
        # Export Results
        tabPanel(
          "Export Results",
          icon = icon("download"),
          value = "export",
          br(),
          h2("Export Results"),
          p("Download figures and data in various formats"),
          mod_export_ui("export_module")
        )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  
  # Initialize reactive values
  app_data <- reactiveValues(
    consolidated_data = NULL,
    data_loaded = FALSE,
    available_genes = character(0),
    default_method = "MAST",
    default_gene = "LRRK2",
    default_cluster = "cluster_0",
    default_experiment = "default",
    default_enrichment = "GO_BP",
    default_direction = "ALL"
  )
  
  # Global p-value threshold
  global_pval <- reactive({
    input$global_pval
  })
  
  # Initialize app with data - run once on startup
  observe({
    # Check environment variables first
    has_data <- Sys.getenv("ISCORE_HAS_DATA", unset = "FALSE") == "TRUE"
    data_file <- Sys.getenv("ISCORE_DATA_FILE", unset = "")
    
    # Only initialize once
    if (!app_data$data_loaded) {
      if (has_data && file.exists(data_file)) {
        cat("Loading provided data file:", data_file, "\n")
        initialize_app_with_data(app_data, data_file)
      } else {
        initialize_app_with_data(app_data)
      }
    }
  })
  
  # Update gene choices based on consolidated data
  observe({
    req(app_data$data_loaded)
    req(input$global_analysis_type)
    
    if (!is.null(app_data$consolidated_data)) {
      genes <- unique(app_data$consolidated_data[
        app_data$consolidated_data$method == input$global_analysis_type, 
        "gene"
      ])
      genes <- sort(genes)
      
      # Set default if available
      selected <- NULL
      if (length(genes) > 0) {
        if (!is.null(app_data$default_gene) && app_data$default_gene %in% genes) {
          selected <- app_data$default_gene
        } else {
          selected <- genes[1]
        }
      }
      
      updateSelectInput(session, "global_gene", 
                       choices = genes, 
                       selected = selected)
    }
  })
  
  # Populate cluster choices based on consolidated data
  observe({
    req(app_data$data_loaded)
    req(input$global_analysis_type, input$global_gene)
    
    if (!is.null(app_data$consolidated_data)) {
      clusters <- unique(app_data$consolidated_data[
        app_data$consolidated_data$method == input$global_analysis_type &
        app_data$consolidated_data$gene == input$global_gene,
        "cluster"
      ])
      clusters <- sort(clusters)
      
      # Set default if available
      selected <- NULL
      if (length(clusters) > 0) {
        if (!is.null(app_data$default_cluster) && app_data$default_cluster %in% clusters) {
          selected <- app_data$default_cluster
        } else {
          selected <- clusters[1]
        }
      }
      
      updateSelectInput(session, "global_cluster", 
                       choices = clusters, 
                       selected = selected)
    }
  })
  
  # Populate experiment choices based on consolidated data
  observe({
    req(app_data$data_loaded)
    req(input$global_analysis_type, input$global_gene, input$global_cluster)
    
    if (!is.null(app_data$consolidated_data)) {
      experiments <- unique(app_data$consolidated_data[
        app_data$consolidated_data$method == input$global_analysis_type &
        app_data$consolidated_data$gene == input$global_gene &
        app_data$consolidated_data$cluster == input$global_cluster,
        "experiment"
      ])
      experiments <- sort(experiments)
      
      # For MAST, default should typically be "default"
      selected <- NULL
      if (length(experiments) > 0) {
        if (input$global_analysis_type == "MAST" && "default" %in% experiments) {
          selected <- "default"
        } else {
          selected <- experiments[1]
        }
      }
      
      updateSelectInput(session, "global_experiment", 
                       choices = experiments, 
                       selected = selected)
    }
  })
  
  # Update enrichment types based on what's available for current selection
  observe({
    req(app_data$data_loaded)
    req(input$global_analysis_type, input$global_gene, input$global_cluster, input$global_experiment)
    
    if (!is.null(app_data$consolidated_data)) {
      # Get available enrichment types for current selection
      available_types <- unique(app_data$consolidated_data[
        app_data$consolidated_data$method == input$global_analysis_type &
        app_data$consolidated_data$gene == input$global_gene &
        app_data$consolidated_data$cluster == input$global_cluster &
        app_data$consolidated_data$experiment == input$global_experiment,
        "enrichment_type"
      ])
      available_types <- sort(available_types)
      
      if (length(available_types) > 0) {
        # Keep current selection if it's still available, otherwise pick first
        current <- input$global_enrichment_type
        selected <- if (current %in% available_types) {
          current
        } else {
          # Try common defaults in order
          defaults <- c("GO_BP", "GO_GOALL", "STRING", "GO_CC", "GO_MF")
          default_match <- defaults[defaults %in% available_types][1]
          if (!is.na(default_match)) default_match else available_types[1]
        }
        
        updateSelectInput(session, "global_enrichment_type",
                         choices = available_types,
                         selected = selected)
        
        # Show notification if enrichment type was changed
        if (!is.null(current) && !(current %in% available_types)) {
          showNotification(
            paste0("Note: ", current, " not available for this selection. Switched to ", selected),
            type = "warning",
            duration = 5
          )
        }
      }
    }
  })
  
  # Data status display
  output$data_status <- renderUI({
    if (app_data$data_loaded) {
      tags$div(
        class = "alert alert-success",
        style = "margin-top: 20px;",
        icon("check-circle"),
        " Data loaded successfully",
        br(),
        tags$small(
          sprintf("%d enrichment results", 
                  ifelse(is.null(app_data$consolidated_data), 0, nrow(app_data$consolidated_data)))
        )
      )
    } else {
      tags$div(
        class = "alert alert-warning",
        style = "margin-top: 20px;",
        icon("exclamation-triangle"),
        " Loading data..."
      )
    }
  })
  
  # Create reactive for global data selection
  global_data_selection <- reactive({
    list(
      analysis_type = input$global_analysis_type,
      gene = input$global_gene,
      cluster = input$global_cluster,
      experiment = input$global_experiment,
      enrichment_type = input$global_enrichment_type,
      direction = input$global_direction,
      pval_threshold = input$global_pval
    )
  })
  
  # Module servers
  # Landing page module
  landingPageWithUmapServer("landing", data = app_data)
  
  # Create filtered data reactive for other modules to use
  filtered_data <- reactive({
    req(app_data$data_loaded)
    selection <- global_data_selection()
    
    get_significant_terms_from_consolidated(
      app_data$consolidated_data,
      analysis_type = selection$analysis_type,
      gene = selection$gene,
      cluster = selection$cluster,
      experiment = selection$experiment,
      enrichment_type = selection$enrichment_type,
      direction = selection$direction,
      pval_threshold = selection$pval_threshold
    )
  })
  
  # Main visualization module (enhanced with interactivity)
  visualization_results <- mod_visualization_server(
    "visualization_module",
    global_selection = global_data_selection,
    enrichment_data = reactive({ app_data$consolidated_data })
  )
  
  comparison_results <- mod_comparison_server(
    "comparison_module",
    app_data = app_data,
    pval_threshold = global_pval
  )
  
  heatmap_results <- mod_heatmap_unified_server(
    "heatmap_module",
    app_data = app_data,
    global_selection = global_data_selection
  )
  
  pathview_results <- mod_pathview_server(
    "pathview_module",
    app_data = app_data,
    selected_enrichment_data = filtered_data,
    global_selection = global_data_selection
  )
  
  export_status <- mod_export_server(
    "export_module",
    app_data = app_data,
    precomputed_data = filtered_data,
    comparison_data = comparison_results,
    visualization_data = visualization_results
  )
}

# Run the app
shinyApp(ui = ui, server = server)