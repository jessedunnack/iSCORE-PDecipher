# startup_manager.R
# Handles initial data loading and app initialization for the Shiny app

initialize_app_with_data <- function(app_data, data_file = NULL) {
  cat("Initializing app with data...\n")
  
  # Check for environment variable first
  if (is.null(data_file)) {
    has_data <- Sys.getenv("ISCORE_HAS_DATA", unset = "FALSE") == "TRUE"
    if (has_data) {
      data_file <- Sys.getenv("ISCORE_DATA_FILE", unset = "")
    }
  }
  
  # Try to load data file if provided
  if (!is.null(data_file) && data_file != "" && file.exists(data_file)) {
    cat("Loading provided data file:", data_file, "\n")
    
    tryCatch({
      # Load the data
      data <- readRDS(data_file)
      cat("Data loaded successfully. Rows:", nrow(data), "Columns:", ncol(data), "\n")
      
      # Fix column names if needed
      if ("mutation_perturbation" %in% names(data) && !"gene" %in% names(data)) {
        data$gene <- data$mutation_perturbation
      }
      
      # Map method column if needed
      if ("method" %in% names(data) && !"analysis_type" %in% names(data)) {
        data$analysis_type <- data$method
      }
      
      # Store in app_data
      app_data$consolidated_data <- data
      app_data$data_loaded <- TRUE
      app_data$available_genes <- sort(unique(data$gene))
      
      cat("App data initialized successfully.\n")
      
    }, error = function(e) {
      cat("Error loading data file:", e$message, "\n")
      app_data$data_loaded <- FALSE
    })
  } else {
    # Try to find data in common locations
    possible_paths <- c(
      "/Users/hockemeyer/Desktop/Functional Enrichment/all_enrichment_padj005_complete_with_direction.rds",
      file.path(getwd(), "data", "all_enrichment_padj005_complete_with_direction.rds"),
      file.path(getwd(), "all_enrichment_padj005_complete_with_direction.rds"),
      file.path("..", "all_enrichment_padj005_complete_with_direction.rds")
    )
    
    data_found <- FALSE
    for (path in possible_paths) {
      if (file.exists(path)) {
        cat("Found data file at:", path, "\n")
        initialize_app_with_data(app_data, path)
        data_found <- TRUE
        break
      }
    }
    
    if (!data_found) {
      cat("Consolidated data not found at expected location\n")
      app_data$data_loaded <- FALSE
    }
  }
}

#' Process uploaded enrichment file
#' 
#' @param file_info File info from fileInput
#' @return Processed data frame or NULL if error
process_uploaded_file <- function(file_info) {
  
  if (is.null(file_info)) {
    return(NULL)
  }
  
  # Check file extension
  ext <- tools::file_ext(file_info$datapath)
  if (ext != "rds") {
    showNotification(
      "Please upload an RDS file containing consolidated enrichment results",
      type = "error"
    )
    return(NULL)
  }
  
  cat("User selected file:", file_info$name, "\n")
  
  tryCatch({
    # Load the data
    data <- readRDS(file_info$datapath)
    
    # Process the data
    if ("source_file" %in% names(data)) {
      # Extract metadata from filename
      data$gene <- sapply(strsplit(data$source_file, "_"), function(x) x[1])
      data$cluster <- sapply(strsplit(data$source_file, "_"), function(x) x[2])
      data$comparison <- sapply(strsplit(data$source_file, "_"), function(x) x[3])
      data$experiment <- sapply(strsplit(data$source_file, "_"), function(x) x[4])
    }
    
    # Verify it's enrichment data
    required_cols <- c("term", "p.adjust")
    if (!all(required_cols %in% names(data))) {
      stop("This doesn't appear to be a consolidated enrichment file")
    }
    
    # Filter for significant results
    if ("p.adjust" %in% names(data)) {
      data <- data[data$p.adjust <= 0.05, ]
    }
    
    showNotification(
      paste("Successfully loaded", nrow(data), "enrichment results"),
      type = "success"
    )
    
    return(data)
    
  }, error = function(e) {
    showNotification(
      paste("Error loading file:", e$message),
      type = "error"
    )
    return(NULL)
  })
}

#' Check if data is loaded
#' 
#' @return TRUE if data is loaded, FALSE otherwise
is_data_loaded <- function() {
  !is.null(app_data$data) && nrow(app_data$data) > 0
}

#' Create file picker modal dialog
#' 
#' @return Modal dialog UI
create_file_picker_modal <- function() {
  modalDialog(
    title = "Select Enrichment Results File",
    
    div(
      class = "alert alert-info",
      icon("info-circle"),
      "The consolidated enrichment data file was not found at the expected location.",
      br(),
      "Please select the file 'consolidated_enrichment_results.rds'"
    ),
    
    fileInput(
      "enrichment_file",
      "Choose RDS file",
      accept = ".rds"
    ),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("load_file", "Load Data", class = "btn-primary")
    ),
    
    size = "m",
    easyClose = FALSE
  )
}