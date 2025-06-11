# Dataset validation functions for iSCORE-PDecipher
# These functions check for required files and validate dataset structure

#' Check if a directory contains required source data
#'
#' @param data_dir Path to dataset directory
#' @return List with status and messages
check_source_data <- function(data_dir) {
  message("Checking source data in: ", data_dir)
  
  results <- list(
    valid = TRUE,
    messages = character(),
    has_mast = FALSE,
    has_mixscale = FALSE
  )
  
  # Check for MAST data
  mast_dir <- file.path(data_dir, "iSCORE-PD_MAST_analysis")
  if (dir.exists(mast_dir)) {
    mast_files <- list.files(mast_dir, pattern = "mutation.*\\.rds$", recursive = TRUE)
    if (length(mast_files) > 0) {
      results$has_mast <- TRUE
      results$messages <- c(results$messages, 
                           paste("✓ Found MAST data:", length(mast_files), "files"))
    } else {
      results$messages <- c(results$messages, "✗ MAST directory exists but contains no RDS files")
      results$valid <- FALSE
    }
  } else {
    results$messages <- c(results$messages, "ℹ No MAST data directory found")
  }
  
  # Check for MixScale data (multiple possible directory names)
  mixscale_patterns <- c("PerturbSeq_MixScale_analysis", "CRISPRi_PerturbSeq_Reports", "CRISPRa_PerturbSeq_Reports")
  mixscale_found <- FALSE
  
  for (pattern in mixscale_patterns) {
    mixscale_dirs <- list.dirs(data_dir, recursive = FALSE)
    matching_dirs <- mixscale_dirs[grep(pattern, basename(mixscale_dirs))]
    
    if (length(matching_dirs) > 0) {
      for (mdir in matching_dirs) {
        mixscale_files <- list.files(mdir, pattern = "DEGs\\.rds$", recursive = TRUE)
        if (length(mixscale_files) > 0) {
          mixscale_found <- TRUE
          results$has_mixscale <- TRUE
          results$messages <- c(results$messages, 
                               paste("✓ Found MixScale data in", basename(mdir), ":", 
                                     length(mixscale_files), "files"))
        }
      }
    }
  }
  
  if (!mixscale_found) {
    results$messages <- c(results$messages, "ℹ No MixScale data found")
  }
  
  # Check if we have at least one data source
  if (!results$has_mast && !results$has_mixscale) {
    results$valid <- FALSE
    results$messages <- c(results$messages, 
                         "✗ No valid source data found (need MAST and/or MixScale results)")
  }
  
  return(results)
}

#' Check which required files are missing
#'
#' @param data_dir Path to dataset directory
#' @return List of missing file types
check_missing_files <- function(data_dir) {
  missing <- character()
  
  # Check for full_DE_results.rds
  if (!file.exists(file.path(data_dir, "full_DE_results.rds"))) {
    missing <- c(missing, "full_DE_results")
  }
  
  # Check for enrichment_results directory
  enrichment_dir <- file.path(data_dir, "enrichment_results")
  if (!dir.exists(enrichment_dir) || length(list.files(enrichment_dir, recursive = TRUE)) == 0) {
    missing <- c(missing, "enrichment_results")
  }
  
  # Check for consolidated enrichment file
  if (!file.exists(file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds"))) {
    missing <- c(missing, "consolidated_enrichment")
  }
  
  return(missing)
}

#' Validate that a dataset directory is ready for the app
#'
#' @param data_dir Path to dataset directory
#' @return List with validation results
validate_dataset_directory <- function(data_dir) {
  if (!dir.exists(data_dir)) {
    return(list(
      valid = FALSE,
      messages = paste("Directory does not exist:", data_dir),
      missing = character()
    ))
  }
  
  # Check source data
  source_check <- check_source_data(data_dir)
  
  if (!source_check$valid) {
    return(list(
      valid = FALSE,
      messages = source_check$messages,
      missing = character()
    ))
  }
  
  # Check for missing files
  missing <- check_missing_files(data_dir)
  
  messages <- source_check$messages
  
  if (length(missing) == 0) {
    messages <- c(messages, "✓ All required files present")
    valid <- TRUE
  } else {
    valid <- FALSE
    if ("full_DE_results" %in% missing) {
      messages <- c(messages, "✗ Missing: full_DE_results.rds (differential expression compilation)")
    }
    if ("enrichment_results" %in% missing) {
      messages <- c(messages, "✗ Missing: enrichment_results/ directory (functional enrichment)")
    }
    if ("consolidated_enrichment" %in% missing) {
      messages <- c(messages, "✗ Missing: consolidated enrichment file")
    }
  }
  
  return(list(
    valid = valid,
    messages = messages,
    missing = missing,
    has_mast = source_check$has_mast,
    has_mixscale = source_check$has_mixscale
  ))
}

#' Get pre-configured dataset options
#'
#' @return Named list of dataset paths
get_dataset_options <- function() {
  # Check which OS we're on and adjust paths accordingly
  if (.Platform$OS.type == "windows") {
    base_path <- "E:/ASAP/scRNASeq/PerturbSeq/final"
  } else {
    base_path <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final"
  }
  
  list(
    "iSCORE-PD only" = file.path(base_path, "iSCORE-PD"),
    "iSCORE-PD + CRISPRi" = file.path(base_path, "iSCORE-PD_plus_CRISPRi"),
    "iSCORE-PD + CRISPRi + CRISPRa" = file.path(base_path, "iSCORE-PD_plus_CRISPRi_and_CRISPRa")
  )
}