#' Launch iSCORE-PDecipher Shiny Application
#'
#' This function launches the interactive Shiny application for exploring
#' Parkinson's disease mutation and perturbation enrichment results.
#' It now supports dataset selection and automatic file generation.
#'
#' @param data_dir Path to the dataset directory containing analysis results.
#'   If NULL, the app will show an interactive dataset selector.
#' @param port Port number for the Shiny app (default: auto-select)
#' @param launch.browser Whether to launch the app in browser (default: TRUE)
#' @param ... Additional arguments passed to shiny::runApp
#'
#' @return Launches the Shiny application
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch with specific dataset
#' launch_iscore_app(data_dir = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/")
#' 
#' # Launch with interactive dataset selection
#' launch_iscore_app()
#' }
launch_iscore_app <- function(data_dir = NULL, port = getOption("shiny.port"), 
                              launch.browser = getOption("shiny.launch.browser", TRUE), ...) {
  
  # Load required functions - handle both installed and development scenarios
  validator_path <- system.file("R", "dataset_validator.R", package = "iSCORE.PDecipher")
  if (validator_path == "") {
    # Development mode - use relative paths
    pkg_root <- if (basename(getwd()) == "iSCORE-PDecipher") {
      getwd()
    } else if (file.exists("iSCORE-PDecipher/R/dataset_validator.R")) {
      "iSCORE-PDecipher"
    } else {
      dirname(getwd())
    }
    validator_path <- file.path(pkg_root, "R", "dataset_validator.R")
    generator_path <- file.path(pkg_root, "R", "data_generator.R")
  } else {
    generator_path <- system.file("R", "data_generator.R", package = "iSCORE.PDecipher")
  }
  
  source(validator_path)
  source(generator_path)
  
  # If no data_dir provided, show dataset selector
  if (is.null(data_dir)) {
    data_dir <- select_dataset_directory()
    if (is.null(data_dir)) {
      stop("No dataset selected. Launch cancelled.", call. = FALSE)
    }
  }
  
  # Validate the dataset directory
  message("\n=== Validating dataset directory ===")
  validation <- validate_dataset_directory(data_dir)
  
  # Show validation messages
  for (msg in validation$messages) {
    message(msg)
  }
  
  # If not valid due to missing files, offer to generate them
  if (!validation$valid && length(validation$missing) > 0) {
    message("\n=== Missing required files ===")
    
    # Check if packages are installed
    if (!check_required_packages()) {
      stop("Please install missing packages before continuing.", call. = FALSE)
    }
    
    # Offer to generate missing files
    success <- generate_missing_files(data_dir, validation$missing)
    
    if (!success) {
      stop("Failed to generate required files. Please check error messages above.", call. = FALSE)
    }
    
    # Re-validate after generation
    message("\n=== Re-validating dataset ===")
    validation <- validate_dataset_directory(data_dir)
    
    if (!validation$valid) {
      stop("Dataset still not valid after generation. Please check the data directory.", call. = FALSE)
    }
  }
  
  message("\n✓ All required files present. Launching app...")
  
  # Get pkg_root for development mode (same as earlier in the function)
  pkg_root <- if (basename(getwd()) == "iSCORE-PDecipher") {
    getwd()
  } else if (file.exists("iSCORE-PDecipher/R/dataset_validator.R")) {
    "iSCORE-PDecipher"
  } else {
    dirname(getwd())
  }
  
  # Get the path to the Shiny app
  app_dir <- system.file("shiny", package = "iSCORE.PDecipher")
  
  if (app_dir == "") {
    # Development mode - try multiple locations
    if (file.exists("inst/shiny/app.R")) {
      app_dir <- "inst/shiny"
    } else if (file.exists("shiny/app.R")) {
      app_dir <- "shiny"
    } else if (file.exists(file.path(pkg_root, "inst/shiny/app.R"))) {
      app_dir <- file.path(pkg_root, "inst/shiny")
    } else {
      # Try the separate shiny_app directory
      shiny_app_path <- file.path(dirname(pkg_root), "shiny_app")
      if (dir.exists(shiny_app_path) && file.exists(file.path(shiny_app_path, "app.R"))) {
        app_dir <- shiny_app_path
      } else {
        stop("Could not find Shiny app in any expected location.", call. = FALSE)
      }
    }
  }
  
  # Set environment variables for the app
  Sys.setenv(ISCORE_DATA_DIR = normalizePath(data_dir))
  Sys.setenv(ISCORE_HAS_DATA = "TRUE")
  
  # Set paths to specific files
  Sys.setenv(ISCORE_DE_FILE = file.path(data_dir, "full_DE_results.rds"))
  Sys.setenv(ISCORE_ENRICHMENT_FILE = file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds"))
  Sys.setenv(ISCORE_ENRICHMENT_DIR = file.path(data_dir, "enrichment_results"))
  
  # Launch the app
  shiny::runApp(appDir = app_dir, 
                port = port, 
                launch.browser = launch.browser, 
                ...)
}

#' Select dataset directory interactively
#'
#' @return Selected directory path or NULL if cancelled
select_dataset_directory <- function() {
  # Get pre-configured options
  options <- get_dataset_options()
  
  message("\n=== Select Dataset ===")
  message("Available datasets:")
  
  for (i in seq_along(options)) {
    message(sprintf("[%d] %s", i, names(options)[i]))
    message(sprintf("    %s", options[[i]]))
  }
  
  message("[C] Choose custom directory")
  message("[Q] Quit")
  
  choice <- readline("Enter your choice: ")
  
  if (tolower(choice) == "q") {
    return(NULL)
  }
  
  if (tolower(choice) == "c") {
    # Use directory chooser if available
    if (requireNamespace("rstudioapi", quietly = TRUE) && 
        rstudioapi::isAvailable()) {
      data_dir <- rstudioapi::selectDirectory(
        caption = "Select dataset directory",
        label = "Select",
        path = dirname(options[[1]])
      )
    } else {
      data_dir <- readline("Enter the full path to your dataset directory: ")
    }
    
    if (data_dir == "" || is.null(data_dir)) {
      return(NULL)
    }
    
    return(normalizePath(data_dir, mustWork = FALSE))
  }
  
  # Check if numeric choice
  choice_num <- suppressWarnings(as.integer(choice))
  if (!is.na(choice_num) && choice_num >= 1 && choice_num <= length(options)) {
    return(options[[choice_num]])
  }
  
  message("Invalid choice. Please try again.")
  return(select_dataset_directory())
}

# Create alias for backward compatibility (not exported)
launch_app <- launch_iscore_app