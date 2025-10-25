#' Data Loader Module for Shiny App
#'
#' Handles loading of large datasets with preview mode option
#' Updated: October 24, 2025 - Adding FDR-corrected pooled data support
#'
#' @author Claude

#' Data Loader UI
#' @export
mod_data_loader_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      column(12,
        wellPanel(
          h4("Dataset Status"),
          
          # Dataset info display
          uiOutput(ns("dataset_info")),
          
          # Preview mode controls
          conditionalPanel(
            condition = "output.has_preview_available",
            ns = ns,
            hr(),
            h5("Performance Mode"),
            radioButtons(
              ns("data_mode"),
              label = NULL,
              choices = list(
                "Preview Mode (50K cells) - Fast" = "preview",
                "Full Dataset - Complete" = "full"
              ),
              selected = "preview",
              inline = TRUE
            ),
            
            # Load full data button
            conditionalPanel(
              condition = "input.data_mode == 'full'",
              ns = ns,
              actionButton(
                ns("load_full"),
                "Load Full Dataset",
                class = "btn-primary",
                icon = icon("database")
              ),
              br(),
              tags$small(
                class = "text-warning",
                "Loading full dataset may take 30-60 seconds and requires ~16GB RAM"
              )
            )
          ),
          
          # Memory usage indicator
          uiOutput(ns("memory_info"))
        )
      )
    )
  )
}

#' Data Loader Server
#' @export
mod_data_loader_server <- function(id, data_path, preview_cells = 50000) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values to store data
    values <- reactiveValues(
      seurat_obj = NULL,
      preview_obj = NULL,
      current_obj = NULL,
      is_preview = TRUE,
      loading = FALSE,
      memory_stats = NULL
    )
    
    # Check if preview is available
    output$has_preview_available <- reactive({
      if (is.null(values$seurat_obj)) return(FALSE)
      ncol(values$seurat_obj) > preview_cells
    })
    outputOptions(output, "has_preview_available", suspendWhenHidden = FALSE)
    
    # Initial data load
    observe({
      if (is.null(values$seurat_obj)) {
        values$loading <- TRUE
        
        withProgress(message = "Loading dataset...", value = 0, {
          
          # Load the full dataset
          incProgress(0.3, detail = "Reading file...")
          values$seurat_obj <- readRDS(data_path)
          
          # Check if we need preview mode
          if (ncol(values$seurat_obj) > preview_cells) {
            incProgress(0.3, detail = "Creating preview dataset...")
            
            # Try to load cached preview first
            cache_dir <- "cache/"
            dataset_hash <- digest::digest(list(
              n_cells = ncol(values$seurat_obj),
              n_genes = nrow(values$seurat_obj),
              preview_cells = preview_cells
            ))
            cache_file <- file.path(cache_dir, sprintf("preview_%s.rds", dataset_hash))
            
            if (file.exists(cache_file)) {
              incProgress(0.3, detail = "Loading cached preview...")
              values$preview_obj <- readRDS(cache_file)
            } else {
              # Create preview
              values$preview_obj <- sample_seurat_cells(
                values$seurat_obj, 
                n_cells = preview_cells,
                preserve_proportions = TRUE
              )
              
              # Cache it
              if (!dir.exists(cache_dir)) {
                dir.create(cache_dir, recursive = TRUE)
              }
              saveRDS(values$preview_obj, cache_file)
            }
            
            # Start with preview
            values$current_obj <- values$preview_obj
            values$is_preview <- TRUE
            
          } else {
            # Dataset is small enough, use full
            values$current_obj <- values$seurat_obj
            values$is_preview <- FALSE
          }
          
          # Calculate memory stats
          incProgress(0.1, detail = "Calculating memory usage...")
          values$memory_stats <- estimate_memory_usage(values$current_obj)
        })
        
        values$loading <- FALSE
      }
    })
    
    # Handle mode switching
    observeEvent(input$load_full, {
      if (!is.null(values$seurat_obj) && values$is_preview) {
        values$loading <- TRUE
        
        withProgress(message = "Loading full dataset...", value = 0.5, {
          values$current_obj <- values$seurat_obj
          values$is_preview <- FALSE
          values$memory_stats <- estimate_memory_usage(values$current_obj)
        })
        
        values$loading <- FALSE
        
        showNotification(
          "Full dataset loaded successfully!",
          type = "success",
          duration = 3
        )
      }
    })
    
    # Dataset info display
    output$dataset_info <- renderUI({
      if (is.null(values$current_obj)) {
        return(tags$p(class = "text-muted", "No dataset loaded"))
      }
      
      n_cells <- ncol(values$current_obj)
      n_genes <- nrow(values$current_obj)
      n_clusters <- length(unique(values$current_obj$seurat_clusters))
      
      mode_text <- if (values$is_preview) {
        tags$span(
          class = "badge badge-warning",
          sprintf("PREVIEW MODE (%s cells)", format(n_cells, big.mark = ","))
        )
      } else {
        tags$span(
          class = "badge badge-success", 
          "FULL DATASET"
        )
      }
      
      tagList(
        tags$p(
          tags$strong("Current Dataset: "), mode_text
        ),
        tags$ul(
          tags$li(sprintf("Cells: %s", format(n_cells, big.mark = ","))),
          tags$li(sprintf("Genes: %s", format(n_genes, big.mark = ","))),
          tags$li(sprintf("Clusters: %d", n_clusters))
        ),
        
        if (values$is_preview && !is.null(values$seurat_obj)) {
          tags$p(
            class = "text-info",
            sprintf("Full dataset available: %s cells",
                   format(ncol(values$seurat_obj), big.mark = ","))
          )
        }
      )
    })
    
    # Memory info display
    output$memory_info <- renderUI({
      if (is.null(values$memory_stats)) {
        return(NULL)
      }
      
      tagList(
        hr(),
        tags$p(
          tags$strong("Memory Usage:"),
          sprintf(" %.1f MB", values$memory_stats$total_mb)
        ),
        if (values$memory_stats$recommended_ram_gb > 8) {
          tags$p(
            class = "text-warning",
            icon("exclamation-triangle"),
            sprintf("Recommended RAM: %d GB", values$memory_stats$recommended_ram_gb)
          )
        }
      )
    })
    
    # Return current dataset for other modules
    return(reactive({
      list(
        data = values$current_obj,
        is_preview = values$is_preview,
        loading = values$loading
      )
    }))
  })
}

#' Quick Dataset Switcher UI
#' @export
quick_data_switcher_ui <- function(id) {
  ns <- NS(id)
  
  conditionalPanel(
    condition = "true",  # Always show in header
    tags$div(
      class = "navbar-text",
      style = "margin-right: 15px;",
      
      # Inline mode switcher
      tags$span(
        id = ns("mode_indicator"),
        style = "margin-right: 10px;"
      ),
      
      # Quick toggle
      actionLink(
        ns("toggle_mode"),
        label = NULL,
        icon = icon("exchange-alt"),
        style = "color: white;"
      )
    )
  )
}

#' Quick Dataset Switcher Server
#' @export
quick_data_switcher_server <- function(id, data_loader) {
  moduleServer(id, function(input, output, session) {
    
    # Mode indicator
    output$mode_indicator <- renderUI({
      data_info <- data_loader()
      
      if (is.null(data_info$data)) {
        return(NULL)
      }
      
      if (data_info$is_preview) {
        tags$span(
          class = "label label-warning",
          "Preview Mode"
        )
      } else {
        tags$span(
          class = "label label-success",
          "Full Dataset"
        )
      }
    })
    
    # Toggle mode
    observeEvent(input$toggle_mode, {
      # This would trigger mode change in parent module
      showModal(modalDialog(
        title = "Switch Dataset Mode",
        "Use the Dataset Status panel to switch between preview and full dataset modes.",
        footer = modalButton("OK")
      ))
    })
  })
}