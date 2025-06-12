# Simple UMAP Viewer Module
# Displays UMAP plot without user controls - dataset determined by app data

mod_umap_viewer_simple_ui <- function(id) {
  ns <- NS(id)
  
  withSpinner(
    plotOutput(ns("umap_plot"), height = "500px"),
    type = 4,
    color = "#3c8dbc"
  )
}

mod_umap_viewer_simple_server <- function(id, app_data, dataset_name = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for UMAP data
    umap_data <- reactiveValues(
      sce = NULL,
      loaded = FALSE
    )
    
    # Load appropriate UMAP data based on app context
    observe({
      req(app_data$data_loaded)
      
      # Determine which dataset to load
      if (!is.null(dataset_name)) {
        dataset_to_load <- dataset_name
      } else {
        # Auto-detect based on data
        has_crispri <- any(grepl("MixScale", app_data$consolidated_data$method))
        has_mutations <- any(grepl("MAST", app_data$consolidated_data$method))
        
        if (has_crispri && has_mutations) {
          dataset_to_load <- "Full_Dataset"
        } else if (has_crispri) {
          dataset_to_load <- "iSCORE_PD_CRISPRi"
        } else {
          dataset_to_load <- "iSCORE_PD"
        }
      }
      
      # Try to load UMAP data
      possible_paths <- c(
        system.file("extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds"), 
                    package = "iSCORE.PDecipher"),
        file.path("inst", "extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds")),
        paste0("../../inst/extdata/umap_data/", dataset_to_load, "_umap_data.rds")
      )
      
      for (path in possible_paths) {
        if (file.exists(path)) {
          tryCatch({
            umap_data$sce <- readRDS(path)
            umap_data$loaded <- TRUE
            message("Loaded UMAP data from: ", path)
            break
          }, error = function(e) {
            message("Failed to load from ", path, ": ", e$message)
          })
        }
      }
    })
    
    # Render UMAP plot
    output$umap_plot <- renderPlot({
      if (!umap_data$loaded || is.null(umap_data$sce)) {
        plot.new()
        text(0.5, 0.5, "UMAP visualization not available", 
             cex = 1.2, col = "gray60")
        return()
      }
      
      # Check if dittoSeq is available
      if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        plot.new()
        text(0.5, 0.5, "dittoSeq package required", 
             cex = 1.2, col = "red")
        return()
      }
      
      library(dittoSeq)
      
      # Create UMAP plot
      tryCatch({
        dittoDimPlot(
          umap_data$sce,
          var = "seurat_clusters",
          reduction.use = "UMAP",
          size = 0.3,
          do.label = TRUE,
          labels.size = 3,
          legend.show = TRUE,
          main = ""
        )
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "Error creating UMAP", 
             cex = 1.2, col = "red")
      })
    })
    
    # Return UMAP data for other uses
    return(umap_data)
  })
}