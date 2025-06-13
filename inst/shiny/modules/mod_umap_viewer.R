# Module: UMAP Viewer
# Displays UMAP plots using dittoseq from pre-extracted minimal data

mod_umap_viewer_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(4,
        wellPanel(
          h4("UMAP Visualization Settings"),
          selectInput(ns("dataset"),
                      "Dataset:",
                      choices = c("iSCORE-PD (Mutations)" = "iSCORE_PD",
                                  "iSCORE-PD + CRISPRi" = "iSCORE_PD_CRISPRi", 
                                  "Full Dataset (CRISPRi+a)" = "Full_Dataset"),
                      selected = "iSCORE_PD"),
          
          selectInput(ns("color_by"),
                      "Color by:",
                      choices = c("Coarse Clusters" = "seurat_clusters",
                                  "Fine Clusters" = "seurat_clusters_fine",
                                  "Perturbation" = "scMAGeCK_gene_assignment"),
                      selected = "seurat_clusters"),
          
          checkboxInput(ns("show_labels"),
                        "Show cluster labels",
                        value = TRUE),
          
          sliderInput(ns("point_size"),
                      "Point size:",
                      min = 0.1,
                      max = 2,
                      value = 0.5,
                      step = 0.1),
          
          sliderInput(ns("label_size"),
                      "Label size:",
                      min = 2,
                      max = 8,
                      value = 4,
                      step = 0.5),
          
          checkboxInput(ns("legend_show"),
                        "Show legend",
                        value = TRUE),
          
          actionButton(ns("refresh_plot"),
                       "Refresh Plot",
                       icon = icon("refresh"),
                       class = "btn-primary")
        )
      ),
      column(8,
        withSpinner(
          plotOutput(ns("umap_plot"), height = "600px"),
          type = 4,
          color = "#3c8dbc"
        ),
        br(),
        verbatimTextOutput(ns("plot_info"))
      )
    )
  )
}

mod_umap_viewer_server <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store data
    umap_data <- reactiveValues(
      sce_list = NULL,
      current_sce = NULL,
      data_loaded = FALSE
    )
    
    # Load UMAP data on module initialization
    observe({
      # Try to load from multiple possible locations
      possible_paths <- c(
        system.file("extdata", "umap_data", "all_umap_data_combined.rds", 
                    package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", "all_umap_data_combined.rds"),
        file.path(dirname(getwd()), "extdata", "umap_data", "all_umap_data_combined.rds"),
        "../../inst/extdata/umap_data/all_umap_data_combined.rds"
      )
      
      data_loaded <- FALSE
      for (path in possible_paths) {
        if (file.exists(path)) {
          tryCatch({
            umap_data$sce_list <- readRDS(path)
            umap_data$data_loaded <- TRUE
            data_loaded <- TRUE
            showNotification("UMAP data loaded successfully!", 
                             type = "message", 
                             duration = 3)
            break
          }, error = function(e) {
            message("Failed to load from ", path, ": ", e$message)
          })
        }
      }
      
      if (!data_loaded) {
        showNotification("UMAP data not found. Please run extract_umap_data.R first.", 
                         type = "error",
                         duration = NULL)
      }
    })
    
    # Update color_by choices based on selected dataset
    observeEvent(input$dataset, {
      req(umap_data$sce_list)
      
      if (input$dataset %in% names(umap_data$sce_list)) {
        sce <- umap_data$sce_list[[input$dataset]]
        available_cols <- colnames(colData(sce))
        
        # Define choices based on available columns
        color_choices <- c()
        if ("seurat_clusters" %in% available_cols) {
          color_choices <- c(color_choices, "Coarse Clusters" = "seurat_clusters")
        }
        if ("seurat_clusters_fine" %in% available_cols) {
          color_choices <- c(color_choices, "Fine Clusters" = "seurat_clusters_fine")
        }
        if ("scMAGeCK_gene_assignment" %in% available_cols) {
          color_choices <- c(color_choices, "Perturbation" = "scMAGeCK_gene_assignment")
        }
        
        updateSelectInput(session, "color_by", choices = color_choices)
        umap_data$current_sce <- sce
      }
    })
    
    # Generate UMAP plot
    output$umap_plot <- renderPlot({
      req(umap_data$current_sce)
      input$refresh_plot  # Trigger on refresh button
      
      # Check if dittoSeq is available
      if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        plot.new()
        text(0.5, 0.5, "dittoSeq package is required for UMAP visualization.\nPlease install with: BiocManager::install('dittoSeq')", 
             cex = 1.2, col = "red")
        return()
      }
      
      library(dittoSeq)
      
      # Check if selected color_by variable exists
      if (!input$color_by %in% colnames(colData(umap_data$current_sce))) {
        plot.new()
        text(0.5, 0.5, paste("Variable", input$color_by, "not found in this dataset"), 
             cex = 1.2, col = "red")
        return()
      }
      
      # Create the plot
      tryCatch({
        p <- dittoDimPlot(
          umap_data$current_sce,
          var = input$color_by,
          reduction.use = "UMAP",
          size = input$point_size,
          do.label = input$show_labels,
          labels.size = input$label_size,
          legend.show = input$legend_show,
          main = paste(input$dataset, "Dataset -", 
                       switch(input$color_by,
                              "seurat_clusters" = "Coarse Clusters",
                              "seurat_clusters_fine" = "Fine Clusters",
                              "scMAGeCK_gene_assignment" = "Perturbations",
                              input$color_by))
        )
        
        # For perturbation coloring, we might want to adjust the plot
        if (input$color_by == "scMAGeCK_gene_assignment") {
          # Get unique perturbations
          perturbs <- unique(colData(umap_data$current_sce)[[input$color_by]])
          n_perturbs <- length(perturbs)
          
          # If too many perturbations, hide legend and labels
          if (n_perturbs > 50 && input$legend_show) {
            p <- p + theme(legend.position = "none")
            showNotification("Legend hidden due to large number of perturbations", 
                             type = "message", duration = 3)
          }
        }
        
        return(p)
        
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating plot:\n", e$message), 
             cex = 1.2, col = "red")
      })
    })
    
    # Display plot information
    output$plot_info <- renderPrint({
      req(umap_data$current_sce)
      
      sce <- umap_data$current_sce
      meta <- metadata(sce)
      
      cat("Dataset Information:\n")
      cat("===================\n")
      cat("Dataset:", input$dataset, "\n")
      cat("Total cells:", meta$n_cells, "\n")
      cat("Coarse clusters:", meta$n_clusters_coarse, "\n")
      
      if (!is.na(meta$n_clusters_fine)) {
        cat("Fine clusters:", meta$n_clusters_fine, "\n")
      }
      
      if (input$color_by %in% colnames(colData(sce))) {
        var_data <- colData(sce)[[input$color_by]]
        n_unique <- length(unique(var_data))
        cat("\nCurrent visualization:\n")
        cat("Variable:", input$color_by, "\n")
        cat("Unique values:", n_unique, "\n")
        
        if (n_unique <= 20) {
          cat("Values:", paste(sort(unique(var_data)), collapse = ", "), "\n")
        }
      }
    })
    
    return(umap_data)
  })
}