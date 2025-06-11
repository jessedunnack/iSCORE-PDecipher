# startup_manager.R
# Handles initial data loading and file selection for the Shiny app

# Initialize app_data in global environment if not exists
if (!exists("app_data")) {
  app_data <- list(
    data = NULL,
    available_genes = character(),
    available_clusters = character(),
    startup_message = NULL,
    file_picker_shown = FALSE
  )
}

# Try to load consolidated data on startup
if (!app_data$file_picker_shown) {
  default_path <- file.path(getwd(), "data", "consolidated_enrichment_results.rds")
  if (file.exists(default_path)) {
    cat("Loading consolidated data from", default_path, "\n")
    
    tryCatch({
      # Load the data
      data <- readRDS(default_path)
      
      # Process the data
      if ("source_file" %in% names(data)) {
        # Extract metadata from filename
        # Assuming filename format: gene_cluster_comparison_experiment_...
        data$gene <- sapply(strsplit(data$source_file, "_"), function(x) x[1])
        data$cluster <- sapply(strsplit(data$source_file, "_"), function(x) x[2])
        data$comparison <- sapply(strsplit(data$source_file, "_"), function(x) x[3])
        data$experiment <- sapply(strsplit(data$source_file, "_"), function(x) x[4])
      }
      
      # Verify required columns
      required_cols <- c("term", "p.adjust", "gene", "cluster", "comparison", "experiment")
      missing_cols <- setdiff(required_cols, names(data))
      
      if (length(missing_cols) > 0) {
        stop("Missing required columns after processing: ", paste(missing_cols, collapse = ", "))
      }
      
      # Filter for significant results
      if ("p.adjust" %in% names(data)) {
        data <- data[data$p.adjust <= 0.05, ]
        cat("Filtered to", nrow(data), "significant results (p.adjust ≤ 0.05)\n")
      }
      
      # Store in app_data
      app_data$data <- data
      app_data$available_genes <- unique(data$gene)
      app_data$available_clusters <- unique(data$cluster)
      app_data$startup_message <- paste("✓ Loaded", nrow(data), "significant enrichment results")
      
      # Create gene sets
      app_data$gene_sets <- split(data$gene_id, data$term)
      
      # Success message
      cat("Successfully loaded", nrow(data), "enrichment terms\n")
      cat("Available genes:", length(app_data$available_genes), "\n")
      cat("Available clusters:", length(app_data$available_clusters), "\n")
      
    }, error = function(e) {
      cat("Error loading consolidated data:", e$message, "\n")
      app_data$data <- NULL
      app_data$file_picker_shown <- TRUE
      app_data$startup_message <- paste("⚠ Error loading data:", e$message)
    })
    
  } else {
    cat("Consolidated data not found at expected location\n")
    app_data$file_picker_shown <- TRUE
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
    
    # Update app_data
    app_data$data <- data
    app_data$available_genes <- unique(data$gene)
    app_data$available_clusters <- unique(data$cluster)
    
    # Create gene sets
    if ("gene_id" %in% names(data) && "term" %in% names(data)) {
      app_data$gene_sets <- split(data$gene_id, data$term)
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