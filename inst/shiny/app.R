# iSCORE-PDecipher Shiny Application
# Handles both direct data loading and file upload interface

library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

# Source modules and functions
source("R/heatmap_functions.R")

# Check if we have data file provided or need upload interface
has_data <- Sys.getenv("ISCORE_HAS_DATA", unset = "FALSE") == "TRUE"
data_file <- Sys.getenv("ISCORE_DATA_FILE", unset = "")

# UI
ui <- fluidPage(
  useShinyjs(),
  
  titlePanel("iSCORE-PDecipher: Parkinson's Disease Enrichment Analysis"),
  
  # Conditional UI based on whether data is provided
  conditionalPanel(
    condition = "!output.dataLoaded",
    
    wellPanel(
      h3("Load Enrichment Data"),
      p("Please upload your consolidated enrichment results file (RDS format):"),
      
      fileInput("dataFile", 
                label = "Choose RDS File",
                accept = c(".rds", ".RDS"),
                placeholder = "No file selected"),
      
      # Drag and drop area
      div(
        id = "dropArea",
        style = "border: 2px dashed #ccc; border-radius: 10px; padding: 20px; text-align: center; margin: 10px 0;",
        p("Or drag and drop your RDS file here", style = "color: #666; margin: 0;")
      ),
      
      # Loading indicator
      conditionalPanel(
        condition = "output.uploading",
        div(
          style = "text-align: center; margin: 20px;",
          tags$i(class = "fa fa-spinner fa-spin", style = "font-size: 24px;"),
          br(),
          "Loading data..."
        )
      )
    )
  ),
  
  # Main analysis interface (shown when data is loaded)
  conditionalPanel(
    condition = "output.dataLoaded",
    
    sidebarLayout(
      sidebarPanel(
        width = 3,
        
        h4("Data Overview"),
        verbatimTextOutput("dataOverview"),
        
        hr(),
        
        h4("Filter Options"),
        
        selectInput("selectedGenes", 
                   "Select Genes/Mutations:",
                   choices = NULL,
                   multiple = TRUE),
        
        selectInput("selectedClusters",
                   "Select Clusters:",
                   choices = NULL,
                   multiple = TRUE),
        
        selectInput("enrichmentType",
                   "Enrichment Type:",
                   choices = NULL),
        
        selectInput("direction",
                   "Direction:",
                   choices = c("All" = "ALL", "Up-regulated" = "UP", "Down-regulated" = "DOWN")),
        
        hr(),
        
        downloadButton("downloadData", "Download Filtered Data", class = "btn-primary")
      ),
      
      mainPanel(
        width = 9,
        
        tabsetPanel(
          type = "tabs",
          
          tabPanel("Data Table",
                  withSpinner(DT::dataTableOutput("enrichmentTable"))
          ),
          
          tabPanel("Heatmap",
                  withSpinner(plotlyOutput("enrichmentHeatmap", height = "600px"))
          ),
          
          tabPanel("Bar Plot",
                  withSpinner(plotlyOutput("enrichmentBarplot", height = "600px"))
          ),
          
          tabPanel("Summary",
                  withSpinner(verbatimTextOutput("summaryStats"))
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  values <- reactiveValues(
    full_data = NULL,
    filtered_data = NULL,
    data_loaded = FALSE
  )
  
  # Check if data was provided via launch function
  observe({
    if (has_data && file.exists(data_file)) {
      cat("Loading provided data file:", data_file, "\n")
      
      tryCatch({
        values$full_data <- readRDS(data_file)
        values$data_loaded <- TRUE
        cat("Loaded", nrow(values$full_data), "enrichment terms\n")
        
        # Initialize UI elements
        updateSelectInput(session, "selectedGenes",
                         choices = sort(unique(values$full_data$mutation_perturbation)))
        
        updateSelectInput(session, "selectedClusters", 
                         choices = sort(unique(values$full_data$cluster)))
        
        updateSelectInput(session, "enrichmentType",
                         choices = sort(unique(values$full_data$enrichment_type)))
        
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
      
      showNotification("Data loaded successfully!", type = "success")
      
      # Initialize UI elements
      updateSelectInput(session, "selectedGenes",
                       choices = sort(unique(values$full_data$mutation_perturbation)))
      
      updateSelectInput(session, "selectedClusters", 
                       choices = sort(unique(values$full_data$cluster)))
      
      updateSelectInput(session, "enrichmentType",
                       choices = sort(unique(values$full_data$enrichment_type)))
      
      output$uploading <- reactive({ FALSE })
      
    }, error = function(e) {
      showNotification(paste("Error loading file:", e$message), type = "error")
      output$uploading <- reactive({ FALSE })
    })
  })
  
  # Output for conditional UI
  output$dataLoaded <- reactive({
    return(values$data_loaded)
  })
  outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)
  
  output$uploading <- reactive({ FALSE })
  outputOptions(output, "uploading", suspendWhenHidden = FALSE)
  
  # Data overview
  output$dataOverview <- renderText({
    req(values$full_data)
    
    paste(
      paste("Total terms:", nrow(values$full_data)),
      paste("Genes:", length(unique(values$full_data$mutation_perturbation))),
      paste("Clusters:", length(unique(values$full_data$cluster))),
      paste("Enrichment types:", length(unique(values$full_data$enrichment_type))),
      paste("Methods:", length(unique(values$full_data$method))),
      sep = "\n"
    )
  })
  
  # Filter data based on selections
  observe({
    req(values$full_data)
    
    filtered <- values$full_data
    
    # Filter by genes
    if (!is.null(input$selectedGenes) && length(input$selectedGenes) > 0) {
      filtered <- filtered[filtered$mutation_perturbation %in% input$selectedGenes, ]
    }
    
    # Filter by clusters  
    if (!is.null(input$selectedClusters) && length(input$selectedClusters) > 0) {
      filtered <- filtered[filtered$cluster %in% input$selectedClusters, ]
    }
    
    # Filter by enrichment type
    if (!is.null(input$enrichmentType) && input$enrichmentType != "") {
      filtered <- filtered[filtered$enrichment_type == input$enrichmentType, ]
    }
    
    # Filter by direction
    if (input$direction != "ALL") {
      filtered <- filtered[filtered$direction == input$direction, ]
    }
    
    values$filtered_data <- filtered
  })
  
  # Data table
  output$enrichmentTable <- DT::renderDataTable({
    req(values$filtered_data)
    
    display_data <- values$filtered_data[, c("mutation_perturbation", "cluster", "enrichment_type", 
                                            "direction", "Description", "p.adjust")]
    
    DT::datatable(display_data,
                  options = list(
                    pageLength = 25,
                    scrollX = TRUE,
                    dom = 'Bfrtip'
                  )) %>%
      DT::formatSignif(columns = "p.adjust", digits = 3)
  })
  
  # Heatmap
  output$enrichmentHeatmap <- renderPlotly({
    req(values$filtered_data)
    req(nrow(values$filtered_data) > 0)
    
    # Simple heatmap showing p-values
    if (nrow(values$filtered_data) > 100) {
      plot_data <- values$filtered_data[1:100, ]  # Limit for performance
    } else {
      plot_data <- values$filtered_data
    }
    
    p <- ggplot(plot_data, aes(x = cluster, y = Description, fill = -log10(p.adjust))) +
      geom_tile() +
      scale_fill_viridis_c(name = "-log10(p.adj)") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8),
            axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Enrichment Heatmap", x = "Cluster", y = "Pathway")
    
    ggplotly(p)
  })
  
  # Bar plot
  output$enrichmentBarplot <- renderPlotly({
    req(values$filtered_data)
    req(nrow(values$filtered_data) > 0)
    
    # Top 20 most significant terms
    top_terms <- values$filtered_data[order(values$filtered_data$p.adjust), ][1:min(20, nrow(values$filtered_data)), ]
    
    p <- ggplot(top_terms, aes(x = reorder(Description, -log10(p.adjust)), y = -log10(p.adjust))) +
      geom_col(aes(fill = direction)) +
      coord_flip() +
      theme_minimal() +
      labs(title = "Top Enriched Terms", x = "Pathway", y = "-log10(p.adjust)")
    
    ggplotly(p)
  })
  
  # Summary statistics
  output$summaryStats <- renderText({
    req(values$filtered_data)
    
    paste(
      paste("Filtered results:", nrow(values$filtered_data), "terms"),
      paste("Significant (p.adj < 0.05):", sum(values$filtered_data$p.adjust < 0.05, na.rm = TRUE)),
      paste("Highly significant (p.adj < 0.01):", sum(values$filtered_data$p.adjust < 0.01, na.rm = TRUE)),
      paste("Mean -log10(p.adj):", round(mean(-log10(values$filtered_data$p.adjust), na.rm = TRUE), 2)),
      sep = "\n"
    )
  })
  
  # Download handler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("iscore_enrichment_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(values$filtered_data, file, row.names = FALSE)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)