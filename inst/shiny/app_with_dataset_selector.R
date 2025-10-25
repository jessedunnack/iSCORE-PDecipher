# iSCORE-PDecipher Shiny App with Dataset Selection
# Enhanced version supporting multiple datasets

# Load required files
source("global.R")
source("R/dataset_selector.R")
source("R/dataset_modal_functions.R")
source("R/startup_config.R")
source("R/data_manager.R")

# Load Shiny modules
source("modules/mod_landing_page_with_umap_v2.R")  # Landing page with UMAP
source("modules/mod_precomputed_reactive.R")       # Precomputed data handling
source("modules/mod_enrichment_gene_display_v2.R") # Enrichment gene display
source("modules/mod_comparison.R")                 # Method comparison
source("modules/mod_de_results.R")                 # DE Results with volcano plots
source("modules/mod_heatmap.R")                    # Heatmap module
source("modules/mod_pathview.R")                   # Pathview module
source("modules/mod_export.R")                     # Export module
source("modules/mod_de_analysis.R")                # Cross-condition DE analysis
source("modules/mod_de_heatmap.R")                 # Interactive DE heatmaps
source("modules/mod_signature_nomination.R")       # Signature nomination
source("modules/mod_signature_trends.R")           # Signature trends
source("modules/mod_perturbseq_only.R")            # NEW v0.5.0: Perturb-seq only

# Create function aliases for landingPageWithUmap module
mod_landing_page_ui <- landingPageWithUmapUI
mod_landing_page_server <- landingPageWithUmapServer

# Create function aliases for enrichment module
mod_enrichment_analysis_ui <- enrichmentGeneDisplayUI
mod_enrichment_analysis_server <- enrichmentGeneDisplayServer

# Initialize startup configuration
startup_config <- create_startup_config()

# Define UI
ui <- fluidPage(
  title = "iSCORE-PDecipher: PD Analysis Platform",
  
  # Include shinyjs for dynamic UI
  useShinyjs(),
  
  # Add custom CSS
  tags$head(
    tags$style(HTML("
      .dataset-banner {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 15px;
        margin-bottom: 20px;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      
      .dataset-banner h3 {
        margin: 0 0 10px 0;
        color: white;
      }
      
      .dataset-info {
        display: flex;
        justify-content: space-between;
        align-items: center;
      }
      
      .dataset-stats {
        display: flex;
        gap: 20px;
      }
      
      .dataset-stat {
        display: flex;
        align-items: center;
        gap: 5px;
      }
      
      .change-dataset-btn {
        background: white;
        color: #667eea;
        border: none;
        padding: 8px 15px;
        border-radius: 3px;
        cursor: pointer;
        font-weight: bold;
      }
      
      .change-dataset-btn:hover {
        background: #f0f0f0;
      }
    "))
  ),
  
  # Dataset selection banner
  div(class = "dataset-banner",
    div(class = "dataset-info",
      div(
        h3(textOutput("current_dataset_name", inline = TRUE)),
        div(class = "dataset-stats",
          div(class = "dataset-stat",
            icon("cell"),
            textOutput("dataset_cells", inline = TRUE)
          ),
          div(class = "dataset-stat",
            icon("dna"),
            textOutput("dataset_genes", inline = TRUE)
          ),
          div(class = "dataset-stat",
            icon("project-diagram"),
            textOutput("dataset_clusters", inline = TRUE)
          ),
          conditionalPanel(
            condition = "output.has_crispr",
            div(class = "dataset-stat",
              icon("cut"),
              "CRISPRi data"
            )
          )
        )
      ),
      actionButton("change_dataset", "Change Dataset", class = "change-dataset-btn")
    )
  ),
  
  # Main app content (conditional on data loaded)
  conditionalPanel(
    condition = "output.data_loaded",
    
    # Navigation tabs
    navbarPage(
      "Analysis",
      id = "main_tabs",
      
      # Landing page tab
      tabPanel(
        "Overview",
        icon = icon("home"),
        mod_landing_page_ui("landing_page")
      ),
      
      # Differential Expression tab
      tabPanel(
        "DE Results",
        icon = icon("chart-line"),
        mod_de_results_ui("de_results")
      ),
      
      # Enrichment Analysis tab
      tabPanel(
        "Enrichment",
        icon = icon("sitemap"),
        mod_enrichment_analysis_ui("enrichment")
      ),
      
      # Heatmap tab
      tabPanel(
        "Heatmaps",
        icon = icon("th"),
        mod_heatmap_ui("heatmap")
      ),
      
      # Settings tab
      tabPanel(
        "Settings",
        icon = icon("cog"),
        wellPanel(
          h4("App Settings"),
          hr(),
          
          h5("Performance"),
          checkboxInput("use_cache", "Enable plot caching", value = TRUE),
          checkboxInput("use_preview", "Use preview mode for large datasets", value = TRUE),
          
          hr(),
          h5("Display"),
          sliderInput("max_terms", "Maximum enrichment terms to display",
                     min = 10, max = 100, value = 30, step = 10),
          
          hr(),
          h5("Data"),
          p("Current data directory:"),
          verbatimTextOutput("current_data_dir"),
          
          hr(),
          actionButton("clear_cache", "Clear Cache", icon = icon("trash")),
          actionButton("reload_data", "Reload Data", icon = icon("refresh"))
        )
      )
    )
  ),
  
  # Loading screen
  conditionalPanel(
    condition = "!output.data_loaded",
    div(style = "text-align: center; padding: 50px;",
      h3("Loading data..."),
      br(),
      div(class = "spinner-border", role = "status")
    )
  )
)

# Define Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    current_dataset = NULL,
    data_loaded = FALSE,
    enrichment_data = NULL,
    de_data = NULL,
    seurat_obj = NULL
  )
  
  # Initialize with default dataset
  observe({
    if (!values$data_loaded && !is.null(startup_config$default_dataset)) {
      load_dataset(startup_config$default_dataset, startup_config)
      
      # Load data
      withProgress(message = "Loading dataset...", value = 0, {
        incProgress(0.3, detail = "Loading enrichment data...")
        values$enrichment_data <- get_enrichment_data(force_reload = TRUE)
        
        incProgress(0.3, detail = "Loading DE results...")
        de_file <- Sys.getenv("ISCORE_DE_FILE")
        if (file.exists(de_file)) {
          values$de_data <- readRDS(de_file)
        }
        
        incProgress(0.3, detail = "Preparing UI...")
        values$current_dataset <- startup_config$default_dataset
        values$data_loaded <- TRUE
      })
    }
  })
  
  # Dataset information outputs
  output$current_dataset_name <- renderText({
    if (!is.null(values$current_dataset)) {
      dataset <- startup_config$available_datasets[[values$current_dataset]]
      dataset$name
    } else {
      "No dataset loaded"
    }
  })
  
  output$dataset_cells <- renderText({
    if (!is.null(values$current_dataset)) {
      dataset <- startup_config$available_datasets[[values$current_dataset]]
      format(dataset$stats$n_cells, big.mark = ",")
    }
  })
  
  output$dataset_genes <- renderText({
    if (!is.null(values$current_dataset)) {
      dataset <- startup_config$available_datasets[[values$current_dataset]]
      format(dataset$stats$n_genes, big.mark = ",")
    }
  })
  
  output$dataset_clusters <- renderText({
    if (!is.null(values$current_dataset)) {
      dataset <- startup_config$available_datasets[[values$current_dataset]]
      dataset$stats$n_clusters
    }
  })
  
  output$has_crispr <- reactive({
    if (!is.null(values$current_dataset)) {
      dataset <- startup_config$available_datasets[[values$current_dataset]]
      dataset$stats$has_crispr
    } else {
      FALSE
    }
  })
  outputOptions(output, "has_crispr", suspendWhenHidden = FALSE)
  
  output$data_loaded <- reactive({
    values$data_loaded
  })
  outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)
  
  # Change dataset button
  observeEvent(input$change_dataset, {
    showModal(
      create_dataset_selector_modal(
        startup_config$available_datasets,
        values$current_dataset
      )
    )
  })
  
  # Handle dataset selection
  observeEvent(input$load_dataset, {
    req(input$selected_dataset)
    
    if (input$selected_dataset != values$current_dataset) {
      
      # Load new dataset
      withProgress(message = "Switching dataset...", value = 0, {
        incProgress(0.2, detail = "Loading configuration...")
        load_dataset(input$selected_dataset, startup_config)
        
        incProgress(0.3, detail = "Loading enrichment data...")
        values$enrichment_data <- get_enrichment_data(force_reload = TRUE)
        
        incProgress(0.3, detail = "Loading DE results...")
        de_file <- Sys.getenv("ISCORE_DE_FILE")
        if (file.exists(de_file)) {
          values$de_data <- readRDS(de_file)
        }
        
        incProgress(0.2, detail = "Updating UI...")
        values$current_dataset <- input$selected_dataset
      })
      
      showNotification(
        paste("Switched to:", startup_config$available_datasets[[input$selected_dataset]]$name),
        type = "success"
      )
    }
    
    removeModal()
  })
  
  # Settings outputs
  output$current_data_dir <- renderText({
    Sys.getenv("ISCORE_DATA_DIR", "Not set")
  })
  
  # Clear cache button
  observeEvent(input$clear_cache, {
    if (exists("GLOBAL_UMAP_CACHE")) {
      GLOBAL_UMAP_CACHE$clear()
    }
    init_data_manager()
    showNotification("Cache cleared", type = "info")
  })
  
  # Reload data button
  observeEvent(input$reload_data, {
    values$enrichment_data <- get_enrichment_data(force_reload = TRUE)
    showNotification("Data reloaded", type = "success")
  })
  
  # Call modules (only if they exist)
  if (exists("mod_landing_page_server")) {
    callModule(mod_landing_page_server, "landing_page", 
              reactive(values$enrichment_data))
  }
  
  if (exists("mod_de_results_server")) {
    callModule(mod_de_results_server, "de_results",
              reactive(values$de_data))
  }
  
  if (exists("mod_enrichment_analysis_server")) {
    callModule(mod_enrichment_analysis_server, "enrichment",
              reactive(values$enrichment_data))
  }
  
  if (exists("mod_heatmap_server")) {
    callModule(mod_heatmap_server, "heatmap",
              reactive(values$enrichment_data))
  }
}

# Run the app
shinyApp(ui = ui, server = server)