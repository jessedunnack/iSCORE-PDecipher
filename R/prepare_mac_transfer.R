#' Prepare Files for Mac Transfer
#'
#' This function helps identify and prepare all required files for transferring
#' the iSCORE-PDecipher datasets to a Mac system.
#'

#' Check which files exist in a dataset directory
#'
#' @param dataset_dir Path to dataset directory
#' @return List with file existence status and sizes
#' @export
check_dataset_files <- function(dataset_dir) {
  if (!dir.exists(dataset_dir)) {
    return(list(
      valid = FALSE,
      message = paste("Dataset directory does not exist:", dataset_dir)
    ))
  }
  
  # Essential files
  essential_files <- c(
    "all_enrichment_padj005_complete_with_direction.rds"
  )
  
  # Optional but recommended files
  recommended_files <- c(
    "full_DE_results.rds"
  )
  
  # UMAP data files
  umap_dir <- file.path(dataset_dir, "inst", "extdata", "umap_data")
  umap_files <- c(
    "all_umap_data_combined.rds"
  )
  
  # Check for dataset-specific UMAP files
  dataset_name <- basename(dataset_dir)
  umap_pattern <- switch(dataset_name,
    "iSCORE-PD" = "iSCORE_PD_umap_data",
    "iSCORE-PD_plus_CRISPRi" = "iSCORE_PD_CRISPRi_umap_data",
    "iSCORE-PD_plus_CRISPRi_plus_CRISPRa" = "Full_Dataset_umap_data",
    "unknown_pattern"
  )
  
  if (umap_pattern != "unknown_pattern") {
    specific_umap_files <- paste0(umap_pattern, c("_30pc.rds", "_50pc.rds", "_100pc.rds"))
    cluster_markers_file <- paste0(gsub("umap_data", "cluster_markers", umap_pattern), ".rds")
    umap_files <- c(umap_files, specific_umap_files, cluster_markers_file)
  }
  
  # Check file existence and sizes
  file_status <- list()
  total_size <- 0
  
  # Check essential files
  for (file in essential_files) {
    file_path <- file.path(dataset_dir, file)
    if (file.exists(file_path)) {
      size_mb <- round(file.info(file_path)$size / 1024^2, 1)
      file_status[[file]] <- list(exists = TRUE, size_mb = size_mb, type = "essential")
      total_size <- total_size + size_mb
    } else {
      file_status[[file]] <- list(exists = FALSE, size_mb = 0, type = "essential")
    }
  }
  
  # Check recommended files
  for (file in recommended_files) {
    file_path <- file.path(dataset_dir, file)
    if (file.exists(file_path)) {
      size_mb <- round(file.info(file_path)$size / 1024^2, 1)
      file_status[[file]] <- list(exists = TRUE, size_mb = size_mb, type = "recommended")
      total_size <- total_size + size_mb
    } else {
      file_status[[file]] <- list(exists = FALSE, size_mb = 0, type = "recommended")
    }
  }
  
  # Check UMAP files
  for (file in umap_files) {
    file_path <- file.path(umap_dir, file)
    if (file.exists(file_path)) {
      size_mb <- round(file.info(file_path)$size / 1024^2, 1)
      file_status[[paste0("umap/", file)]] <- list(exists = TRUE, size_mb = size_mb, type = "umap")
      total_size <- total_size + size_mb
    } else {
      file_status[[paste0("umap/", file)]] <- list(exists = FALSE, size_mb = 0, type = "umap")
    }
  }
  
  return(list(
    valid = TRUE,
    dataset_name = dataset_name,
    file_status = file_status,
    total_size_mb = total_size,
    umap_dir = umap_dir
  ))
}

#' Generate transfer report for all datasets
#'
#' @param parent_dir Parent directory containing dataset folders
#' @return Comprehensive report of all files and sizes
#' @export
generate_transfer_report <- function(parent_dir = NULL) {
  if (is.null(parent_dir)) {
    parent_dir <- get_parent_data_dir()
    if (is.null(parent_dir)) {
      stop("No parent directory configured. Please run setup_parent_dir() first.")
    }
  }
  
  if (!dir.exists(parent_dir)) {
    stop("Parent directory does not exist: ", parent_dir)
  }
  
  # Check for dataset directories
  dataset_dirs <- c(
    "iSCORE-PD",
    "iSCORE-PD_plus_CRISPRi", 
    "iSCORE-PD_plus_CRISPRi_plus_CRISPRa"
  )
  
  report <- list()
  total_essential <- 0
  total_recommended <- 0
  total_umap <- 0
  
  cat("=== iSCORE-PDecipher Mac Transfer Report ===\n\n")
  cat("Parent Directory:", parent_dir, "\n\n")
  
  for (dataset_name in dataset_dirs) {
    dataset_path <- file.path(parent_dir, dataset_name)
    
    if (dir.exists(dataset_path)) {
      cat("Dataset:", dataset_name, "\n")
      cat("Path:", dataset_path, "\n")
      
      status <- check_dataset_files(dataset_path)
      report[[dataset_name]] <- status
      
      if (status$valid) {
        # Categorize files by type
        essential_size <- sum(sapply(status$file_status, function(x) {
          if (x$type == "essential" && x$exists) x$size_mb else 0
        }))
        
        recommended_size <- sum(sapply(status$file_status, function(x) {
          if (x$type == "recommended" && x$exists) x$size_mb else 0
        }))
        
        umap_size <- sum(sapply(status$file_status, function(x) {
          if (x$type == "umap" && x$exists) x$size_mb else 0
        }))
        
        total_essential <- total_essential + essential_size
        total_recommended <- total_recommended + recommended_size
        total_umap <- total_umap + umap_size
        
        cat("  Essential files:", essential_size, "MB\n")
        cat("  Recommended files:", recommended_size, "MB\n")
        cat("  UMAP files:", umap_size, "MB\n")
        cat("  Total:", status$total_size_mb, "MB\n")
        
        # List missing essential files
        missing_essential <- names(status$file_status)[sapply(status$file_status, function(x) {
          x$type == "essential" && !x$exists
        })]
        
        if (length(missing_essential) > 0) {
          cat("  ⚠️  Missing essential files:", paste(missing_essential, collapse = ", "), "\n")
        }
        
        # Check if UMAP directory exists
        if (!dir.exists(status$umap_dir)) {
          cat("  ⚠️  UMAP directory missing:", status$umap_dir, "\n")
        }
        
      } else {
        cat("  ❌ Error:", status$message, "\n")
      }
      
    } else {
      cat("Dataset:", dataset_name, "❌ Not found\n")
    }
    
    cat("\n")
  }
  
  # Summary
  cat("=== Transfer Summary ===\n")
  cat("Total essential files:", total_essential, "MB\n")
  cat("Total recommended files:", total_recommended, "MB\n") 
  cat("Total UMAP files:", total_umap, "MB\n")
  cat("Grand total:", total_essential + total_recommended + total_umap, "MB\n\n")
  
  # Recommendations
  cat("=== Recommendations ===\n")
  cat("Minimum transfer (essential + UMAP):", total_essential + total_umap, "MB\n")
  cat("Complete transfer (all files):", total_essential + total_recommended + total_umap, "MB\n\n")
  
  cat("For Mac compatibility:\n")
  cat("1. Copy the entire parent directory structure to Mac\n")
  cat("2. Install iSCORE.PDecipher package on Mac\n")
  cat("3. Run launch_iscore_app() and specify the Mac path when prompted\n")
  cat("4. Missing files can be auto-generated if source data is available\n\n")
  
  invisible(report)
}

#' Copy essential files for Mac transfer
#'
#' @param source_parent_dir Source parent directory
#' @param dest_parent_dir Destination parent directory  
#' @param copy_mode One of "minimal", "recommended", or "complete"
#' @export
prepare_mac_transfer_copy <- function(source_parent_dir, dest_parent_dir, copy_mode = "recommended") {
  
  if (!dir.exists(source_parent_dir)) {
    stop("Source directory does not exist: ", source_parent_dir)
  }
  
  # Create destination directory
  if (!dir.exists(dest_parent_dir)) {
    dir.create(dest_parent_dir, recursive = TRUE)
  }
  
  dataset_dirs <- c(
    "iSCORE-PD",
    "iSCORE-PD_plus_CRISPRi", 
    "iSCORE-PD_plus_CRISPRi_plus_CRISPRa"
  )
  
  for (dataset_name in dataset_dirs) {
    source_dataset <- file.path(source_parent_dir, dataset_name)
    dest_dataset <- file.path(dest_parent_dir, dataset_name)
    
    if (dir.exists(source_dataset)) {
      cat("Preparing", dataset_name, "...\n")
      
      # Create destination dataset directory
      if (!dir.exists(dest_dataset)) {
        dir.create(dest_dataset, recursive = TRUE)
      }
      
      status <- check_dataset_files(source_dataset)
      
      if (status$valid) {
        files_to_copy <- names(status$file_status)[sapply(status$file_status, function(x) {
          x$exists && (
            x$type == "essential" || 
            (copy_mode %in% c("recommended", "complete") && x$type == "recommended") ||
            (copy_mode == "complete" && x$type == "umap") ||
            (x$type == "umap" && grepl("combined", names(status$file_status)[1]))
          )
        })]
        
        for (file_key in files_to_copy) {
          if (startsWith(file_key, "umap/")) {
            # UMAP file
            file_name <- sub("umap/", "", file_key)
            source_file <- file.path(source_dataset, "inst", "extdata", "umap_data", file_name)
            dest_umap_dir <- file.path(dest_dataset, "inst", "extdata", "umap_data")
            
            if (!dir.exists(dest_umap_dir)) {
              dir.create(dest_umap_dir, recursive = TRUE)
            }
            
            dest_file <- file.path(dest_umap_dir, file_name)
          } else {
            # Root level file
            source_file <- file.path(source_dataset, file_key)
            dest_file <- file.path(dest_dataset, file_key)
          }
          
          if (file.exists(source_file)) {
            cat("  Copying", file_key, "...\n")
            file.copy(source_file, dest_file, overwrite = TRUE)
          }
        }
      }
    }
  }
  
  cat("\nTransfer preparation complete!\n")
  cat("Destination:", dest_parent_dir, "\n")
  
  # Generate report for destination
  cat("\nValidating copied files...\n")
  generate_transfer_report(dest_parent_dir)
}