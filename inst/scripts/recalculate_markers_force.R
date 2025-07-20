#!/usr/bin/env Rscript

# Script to force recalculation of cluster markers with new parameters
# This ensures old marker files are overwritten

cat("=====================================\n")
cat("Cluster Marker Recalculation Script\n")
cat("=====================================\n")
cat("\nThis script will:\n")
cat("1. Delete existing marker files\n")
cat("2. Recalculate markers with new parameters:\n")
cat("   - LFC threshold: 0.5 (was 0.25)\n")
cat("   - Min.pct: 0.25 (was 0.1)\n")
cat("   - Min.diff.pct: 0.2 (new)\n")
cat("   - Including both positive and negative markers\n")
cat("\n")

# Check if user wants to proceed
if (interactive()) {
  response <- readline("Do you want to proceed? (y/n): ")
  if (tolower(response) != "y") {
    cat("Cancelled.\n")
    quit("no")
  }
}

# Load required libraries
library(Seurat)
library(SingleCellExperiment)
library(dplyr)

# Set paths
SCRIPT_DIR <- dirname(sys.frame(1)$ofile)
if (is.null(SCRIPT_DIR)) {
  SCRIPT_DIR <- getwd()
}

# Source the main extraction script
source(file.path(SCRIPT_DIR, "extract_umap_data.R"))

# Path to marker files
MARKER_DIR <- file.path(dirname(dirname(SCRIPT_DIR)), "extdata", "umap_data")

# List existing marker files
marker_files <- list.files(MARKER_DIR, pattern = "_cluster_markers\\.rds$", full.names = TRUE)

if (length(marker_files) > 0) {
  cat("\nFound existing marker files:\n")
  for (f in marker_files) {
    file_info <- file.info(f)
    cat(sprintf("  - %s (modified: %s)\n", basename(f), file_info$mtime))
  }
  
  # Backup old files
  backup_dir <- file.path(MARKER_DIR, paste0("backup_", format(Sys.Date(), "%Y%m%d")))
  dir.create(backup_dir, showWarnings = FALSE)
  
  cat("\nBacking up old files to:", backup_dir, "\n")
  for (f in marker_files) {
    file.copy(f, file.path(backup_dir, basename(f)))
  }
  
  # Delete old files
  cat("\nDeleting old marker files...\n")
  file.remove(marker_files)
}

# Modify the main function to force recalculation
main_force_recalc <- function() {
  message("\nStarting FORCED UMAP data extraction with marker recalculation...")
  message("=======================================")
  
  sce_list <- list()
  
  for (dataset_name in names(DATASETS)) {
    dataset_info <- DATASETS[[dataset_name]]
    seurat_file <- dataset_info$file
    
    if (!file.exists(seurat_file)) {
      warning(sprintf("Seurat file not found: %s", seurat_file))
      next
    }
    
    message(sprintf("\nProcessing %s dataset...", dataset_name))
    message(sprintf("Loading from: %s", seurat_file))
    
    # Load Seurat object
    seurat_obj <- tryCatch({
      readRDS(seurat_file)
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", seurat_file, e$message))
      return(NULL)
    })
    
    if (is.null(seurat_obj)) next
    
    # Extract minimal data with FORCED MARKER RECALCULATION
    sce <- tryCatch({
      extract_minimal_umap_data(seurat_obj, dataset_name, 
                               expected_clusters = dataset_info$expected_clusters,
                               force_recalculate = TRUE)  # FORCE RECALCULATION
    }, error = function(e) {
      warning(sprintf("Failed to extract data from %s: %s", dataset_name, e$message))
      return(NULL)
    })
    
    if (!is.null(sce)) {
      sce_list[[dataset_name]] <- sce
      
      # Save individual dataset files
      for (pc_version in c("pc1", "pc2", "pc3")) {
        output_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data_%s.rds", dataset_name, pc_version))
        saveRDS(sce, output_file)
        message(sprintf("  - Saved %s data to %s", pc_version, output_file))
      }
      
      # Legacy filename for compatibility
      legacy_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data.rds", dataset_name))
      saveRDS(sce, legacy_file)
      message(sprintf("  - Saved legacy format to %s", legacy_file))
    }
  }
  
  if (length(sce_list) > 0) {
    # Create summary statistics
    summary_stats <- create_summary_stats(sce_list)
    summary_file <- file.path(OUTPUT_DIR, "umap_data_summary.csv")
    write.csv(summary_stats, summary_file, row.names = FALSE)
    message(sprintf("Saved summary statistics to: %s", summary_file))
    
    # Print summary
    message("\nExtraction Summary:")
    message("===================")
    print(summary_stats)
  }
  
  message("\nForced recalculation complete!")
  
  # List new marker files
  new_marker_files <- list.files(MARKER_DIR, pattern = "_cluster_markers\\.rds$", full.names = TRUE)
  if (length(new_marker_files) > 0) {
    cat("\nNew marker files created:\n")
    for (f in new_marker_files) {
      file_info <- file.info(f)
      cat(sprintf("  - %s (size: %d bytes)\n", basename(f), file_info$size))
    }
  }
}

# Run the forced recalculation
cat("\nStarting marker recalculation...\n")
cat("This may take 10-30 minutes depending on dataset size.\n\n")

main_force_recalc()

cat("\n✓ Marker recalculation complete!\n")
cat("\nNext steps:\n")
cat("1. Reinstall the package: remotes::install_github('jessedunnack/iSCORE-PDecipher', force = TRUE)\n")
cat("2. Launch the app: launch_iscore_app()\n")
cat("3. Check the Overview page to see new markers\n")