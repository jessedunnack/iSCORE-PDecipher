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

#' Initialize App Data (Call when ready)
#' Sets up app_data with centralized data management
initialize_app_data <- function() {
  
  # Only initialize once
  if (!is.null(app_data$data)) {
    return(TRUE)
  }
  
  tryCatch({
    # Use centralized data manager instead of loading directly
    data <- get_enrichment_data()
    
    if (is.null(data)) {
      app_data$startup_message <- "⚠ Failed to load enrichment data"
      return(FALSE)
    }
    
    # Add gene column for compatibility
    if (!"gene" %in% names(data)) {
      data$gene <- data$mutation_perturbation
    }
    
    # Store in app_data
    app_data$data <- data
    app_data$consolidated_data <- data
    app_data$available_genes <- get_available_genes()
    app_data$available_clusters <- get_available_clusters()
    app_data$startup_message <- paste("✓ Loaded", nrow(data), "significant enrichment results")
    
    cat("App data initialized with", nrow(data), "enrichment terms\n")
    return(TRUE)
    
  }, error = function(e) {
    cat("Error initializing app data:", e$message, "\n")
    app_data$startup_message <- paste("⚠ Error:", e$message)
    return(FALSE)
  })
}

# NOTE: No immediate execution - data loading is now on-demand

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

#' Initialize app with data
#' 
#' @param app_data Reactive values to update
#' @param data_file Optional path to data file
initialize_app_with_data <- function(app_data, data_file = NULL) {
  # Check if data already loaded from startup
  if (!is.null(app_data$data) && nrow(app_data$data) > 0) {
    # Data already loaded in startup_manager.R
    app_data$consolidated_data <- app_data$data
    app_data$data_loaded <- TRUE
    return(invisible(TRUE))
  }
  
  # Otherwise try to load data
  if (!is.null(data_file) && file.exists(data_file)) {
    tryCatch({
      data <- readRDS(data_file)
      
      # Process column names
      if ("mutation_perturbation" %in% names(data)) {
        data$gene <- data$mutation_perturbation
      }
      
      # Add method column if not present
      if (!"method" %in% names(data)) {
        if ("log2FC" %in% names(data)) {
          data$method <- "MixScale"
        } else {
          data$method <- "MAST"
        }
      }
      
      # Add experiment column if not present
      if (!"experiment" %in% names(data)) {
        data$experiment <- "Default"
      }
      
      # Store data
      app_data$data <- data
      app_data$consolidated_data <- data
      app_data$available_genes <- unique(data$gene)
      app_data$available_clusters <- unique(data$cluster)
      app_data$data_loaded <- TRUE
      
      cat("Successfully initialized app with", nrow(data), "enrichment results\n")
      
    }, error = function(e) {
      cat("Error loading data file:", e$message, "\n")
      app_data$data_loaded <- FALSE
    })
  } else {
    # Try default loading through load_enrichment_data
    tryCatch({
      data <- load_enrichment_data("all_modalities")
      
      # Process column names
      if ("mutation_perturbation" %in% names(data)) {
        data$gene <- data$mutation_perturbation
      }
      
      # Add method column if not present
      if (!"method" %in% names(data)) {
        if ("log2FC" %in% names(data)) {
          data$method <- "MixScale"
        } else {
          data$method <- "MAST"
        }
      }
      
      # Add experiment column if not present
      if (!"experiment" %in% names(data)) {
        data$experiment <- "Default"
      }
      
      # Store data
      app_data$data <- data
      app_data$consolidated_data <- data
      app_data$available_genes <- unique(data$gene)
      app_data$available_clusters <- unique(data$cluster)
      app_data$data_loaded <- TRUE
      
      cat("Successfully initialized app with", nrow(data), "enrichment results\n")
      
    }, error = function(e) {
      cat("Error in initialize_app_with_data:", e$message, "\n")
      app_data$data_loaded <- FALSE
    })
  }
}