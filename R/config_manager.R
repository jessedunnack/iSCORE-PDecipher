#' Configuration Management for iSCORE-PDecipher
#'
#' Functions to manage user configuration settings including dataset paths
#' in a cross-platform compatible way.
#'

#' Get the configuration file path
#'
#' @return Path to the configuration file
#' @export
get_config_path <- function() {
  config_dir <- rappdirs::user_config_dir("iSCORE.PDecipher", "jessedunnack")
  
  # Create config directory if it doesn't exist
  if (!dir.exists(config_dir)) {
    dir.create(config_dir, recursive = TRUE)
  }
  
  file.path(config_dir, "config.json")
}

#' Load configuration settings
#'
#' @return List containing configuration settings
#' @export
load_config <- function() {
  config_file <- get_config_path()
  
  if (file.exists(config_file)) {
    tryCatch({
      config <- jsonlite::fromJSON(config_file)
      return(config)
    }, error = function(e) {
      warning("Failed to load config file, using defaults: ", e$message)
      return(list())
    })
  }
  
  # Return empty list if no config file exists
  list()
}

#' Save configuration settings
#'
#' @param config List containing configuration settings
#' @export
save_config <- function(config) {
  config_file <- get_config_path()
  
  tryCatch({
    jsonlite::write_json(config, config_file, pretty = TRUE, auto_unbox = TRUE)
    return(TRUE)
  }, error = function(e) {
    warning("Failed to save config file: ", e$message)
    return(FALSE)
  })
}

#' Get the parent data directory path
#'
#' @return Path to the parent directory containing dataset folders, or NULL if not set
#' @export
get_parent_data_dir <- function() {
  config <- load_config()
  config$parent_data_dir %||% NULL
}

#' Set the parent data directory path
#'
#' @param path Path to the parent directory containing dataset folders
#' @export
set_parent_data_dir <- function(path) {
  config <- load_config()
  config$parent_data_dir <- normalizePath(path, mustWork = FALSE)
  save_config(config)
}

#' Check if this is the first launch (no config exists)
#'
#' @return Logical indicating if this is the first launch
#' @export
is_first_launch <- function() {
  config <- load_config()
  is.null(config$parent_data_dir) || !dir.exists(config$parent_data_dir)
}

#' Prompt user to select parent data directory
#'
#' @return Path to selected directory, or NULL if cancelled
#' @export
prompt_for_parent_dir <- function() {
  cat("\n")
  cat("=== iSCORE-PDecipher First-Time Setup ===\n")
  cat("\n")
  cat("Please specify the parent directory that contains your dataset folders.\n")
  cat("This directory should contain the following subdirectories:\n")
  cat("  - iSCORE-PD/\n")
  cat("  - iSCORE-PD_plus_CRISPRi/\n")
  cat("\n")
  
  if (interactive()) {
    # Try to use system file chooser if available
    if (requireNamespace("tcltk", quietly = TRUE)) {
      tryCatch({
        selected_dir <- tcltk::tk_choose.dir(
          caption = "Select parent directory containing iSCORE dataset folders"
        )
        if (!is.na(selected_dir) && nzchar(selected_dir)) {
          return(selected_dir)
        }
      }, error = function(e) {
        # Fall back to manual entry if GUI fails
      })
    }
  }
  
  # Manual path entry fallback
  while (TRUE) {
    cat("Enter the full path to your parent data directory: ")
    path <- readline()
    
    if (nzchar(path)) {
      # Expand ~ for home directory on Unix-like systems
      path <- path.expand(path)
      
      if (dir.exists(path)) {
        return(normalizePath(path))
      } else {
        cat("Directory does not exist:", path, "\n")
        cat("Please try again or press Ctrl+C to cancel.\n\n")
      }
    } else {
      cat("Please enter a valid path or press Ctrl+C to cancel.\n\n")
    }
  }
}

#' Validate that required dataset folders exist in parent directory
#'
#' @param parent_dir Path to parent directory
#' @return List with validation results
#' @export
validate_parent_dir <- function(parent_dir) {
  if (!dir.exists(parent_dir)) {
    return(list(
      valid = FALSE,
      message = paste("Parent directory does not exist:", parent_dir)
    ))
  }
  
  required_folders <- c(
    "iSCORE-PD",
    "iSCORE-PD_plus_CRISPRi"
  )
  
  existing_folders <- character(0)
  missing_folders <- character(0)
  
  for (folder in required_folders) {
    folder_path <- file.path(parent_dir, folder)
    if (dir.exists(folder_path)) {
      existing_folders <- c(existing_folders, folder)
    } else {
      missing_folders <- c(missing_folders, folder)
    }
  }
  
  if (length(existing_folders) == 0) {
    return(list(
      valid = FALSE,
      message = paste(
        "No dataset folders found in:", parent_dir,
        "\nExpected folders:", paste(required_folders, collapse = ", ")
      )
    ))
  }
  
  if (length(missing_folders) > 0) {
    message <- paste(
      "Found", length(existing_folders), "of", length(required_folders), "dataset folders.",
      "\nAvailable:", paste(existing_folders, collapse = ", "),
      "\nMissing:", paste(missing_folders, collapse = ", ")
    )
  } else {
    message <- paste("All", length(required_folders), "dataset folders found.")
  }
  
  return(list(
    valid = TRUE,
    message = message,
    existing_folders = existing_folders,
    missing_folders = missing_folders
  ))
}

#' Setup parent data directory with validation
#'
#' @param prompt_if_missing Logical, whether to prompt user if no valid config exists
#' @return Path to parent directory, or NULL if setup failed
#' @export
setup_parent_dir <- function(prompt_if_missing = TRUE) {
  # Check if we already have a valid parent directory
  parent_dir <- get_parent_data_dir()
  
  if (!is.null(parent_dir)) {
    validation <- validate_parent_dir(parent_dir)
    if (validation$valid) {
      return(parent_dir)
    } else {
      cat("Current parent directory is invalid:", validation$message, "\n")
      if (!prompt_if_missing) {
        return(NULL)
      }
    }
  }
  
  # Need to set up parent directory
  if (!prompt_if_missing) {
    return(NULL)
  }
  
  while (TRUE) {
    parent_dir <- prompt_for_parent_dir()
    
    if (is.null(parent_dir)) {
      cat("Setup cancelled.\n")
      return(NULL)
    }
    
    validation <- validate_parent_dir(parent_dir)
    
    if (validation$valid) {
      cat("\n", validation$message, "\n")
      set_parent_data_dir(parent_dir)
      cat("Configuration saved.\n\n")
      return(parent_dir)
    } else {
      cat("\nError:", validation$message, "\n")
      cat("Please select a different directory.\n\n")
    }
  }
}