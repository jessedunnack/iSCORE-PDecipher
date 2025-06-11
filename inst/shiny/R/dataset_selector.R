# Dataset Selector Module for PerturbSeq Analysis App
# Handles dataset selection UI and logic

#' Dataset Selector UI
#' @param id Module ID
#' @return Shiny UI elements
datasetSelectorUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Current dataset indicator
    div(id = ns("dataset_indicator"),
      style = "background: #f8f9fa; border: 1px solid #dee2e6; border-radius: 5px; padding: 10px; margin-bottom: 15px;",
      
      div(style = "display: flex; justify-content: space-between; align-items: center;",
        div(
          strong("Current Dataset: "),
          span(id = ns("current_dataset_name"), 
               style = "color: #0066cc; font-weight: bold;", 
               "Loading...")
        ),
        actionButton(ns("change_dataset"), 
                    "Change Dataset", 
                    class = "btn-outline-secondary btn-sm",
                    icon = icon("exchange-alt"))
      ),
      
      # Dataset info
      div(style = "margin-top: 8px; font-size: 0.9em; color: #666;",
        textOutput(ns("dataset_info_text"))
      )
    )
  )
}

#' Dataset Selector Server
#' @param id Module ID
#' @param startup_config Reactive startup configuration
#' @return List of reactive values
datasetSelectorServer <- function(id, startup_config) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values
    values <- reactiveValues(
      current_dataset = NULL,
      dataset_info = NULL,
      show_selector = FALSE
    )
    
    # Initialize with default dataset
    observe({
      req(startup_config())
      config <- startup_config()
      
      if (config$app_ready && !is.null(config$default_dataset)) {
        values$current_dataset <- config$default_dataset
        
        # Update display
        dataset_list <- config$available_datasets
        if (values$current_dataset %in% names(dataset_list)) {
          dataset_info <- dataset_list[[values$current_dataset]]
          updateTextOutput(session, "current_dataset_name", dataset_info$name)
        }
      }
    })
    
    # Show change dataset modal
    observeEvent(input$change_dataset, {
      config <- startup_config()
      if (!is.null(config) && config$app_ready) {
        showModal(
          create_dataset_selector_modal(
            config$available_datasets,
            values$current_dataset
          )
        )
      }
    })
    
    # Handle dataset selection from modal
    observeEvent(input$selected_dataset, {
      if (!is.null(input$selected_dataset)) {
        values$selected_for_loading <- input$selected_dataset
      }
    })
    
    # Load selected dataset
    observeEvent(input$load_dataset, {
      req(input$selected_dataset)
      
      # Update current dataset
      values$current_dataset <- input$selected_dataset
      
      # Update UI display
      config <- startup_config()
      dataset_list <- config$available_datasets
      if (values$current_dataset %in% names(dataset_list)) {
        dataset_info <- dataset_list[[values$current_dataset]]
        
        # Update dataset name display
        session$sendCustomMessage("updateDatasetName", dataset_info$name)
        
        # Update info text
        values$dataset_info <- dataset_info$description
      }
      
      # Close modal
      removeModal()
      
      # Signal that dataset changed
      values$dataset_changed <- Sys.time()
    })
    
    # Show upload option
    observeEvent(input$show_upload_option, {
      removeModal()
      showModal(create_custom_upload_modal())
    })
    
    # Handle custom file upload validation
    observeEvent(input$custom_dataset_file, {
      req(input$custom_dataset_file)
      
      validation_result <- validate_dataset(input$custom_dataset_file$datapath)
      
      if (validation_result$valid) {
        output$upload_validation_message <- renderUI({
          div(class = "alert alert-success",
            icon("check"),
            "Dataset validated successfully! Ready to load."
          )
        })
        
        # Enable load button
        shinyjs::enable("load_custom_dataset")
        
      } else {
        output$upload_validation_message <- renderUI({
          div(class = "alert alert-danger",
            icon("exclamation-triangle"),
            validation_result$message
          )
        })
        
        # Disable load button
        shinyjs::disable("load_custom_dataset")
      }
    })
    
    # Load custom dataset
    observeEvent(input$load_custom_dataset, {
      req(input$custom_dataset_file)
      
      # Set current dataset to custom file path
      values$current_dataset <- input$custom_dataset_file$datapath
      values$dataset_info <- paste("Custom upload:", input$custom_dataset_file$name)
      
      # Update UI
      session$sendCustomMessage("updateDatasetName", "Custom Dataset")
      
      # Close modal
      removeModal()
      
      # Signal dataset changed
      values$dataset_changed <- Sys.time()
    })
    
    # Back to dataset selector
    observeEvent(input$back_to_dataset_selector, {
      removeModal()
      config <- startup_config()
      if (!is.null(config) && config$app_ready) {
        showModal(
          create_dataset_selector_modal(
            config$available_datasets,
            values$current_dataset
          )
        )
      }
    })
    
    # Dataset info text output
    output$dataset_info_text <- renderText({
      if (!is.null(values$dataset_info)) {
        return(values$dataset_info)
      } else if (!is.null(values$current_dataset)) {
        return(paste("Dataset:", values$current_dataset))
      } else {
        return("No dataset loaded")
      }
    })
    
    # Current dataset name output  
    output$current_dataset_name <- renderText({
      if (!is.null(values$current_dataset)) {
        config <- startup_config()
        if (!is.null(config$available_datasets) && 
            values$current_dataset %in% names(config$available_datasets)) {
          return(config$available_datasets[[values$current_dataset]]$name)
        } else {
          return("Custom Dataset")
        }
      } else {
        return("Loading...")
      }
    })
    
    # Return reactive values for parent app
    return(reactive({
      list(
        current_dataset = values$current_dataset,
        dataset_changed = values$dataset_changed,
        dataset_info = values$dataset_info
      )
    }))
  })
}

cat("Dataset selector module loaded successfully\n")