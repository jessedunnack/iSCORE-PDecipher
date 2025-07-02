#!/usr/bin/env Rscript

# Archive existing Dataset 2 UMAP files before regenerating

cat("=== ARCHIVING DATASET 2 UMAP FILES ===\n\n")

# Auto-detect platform and set paths accordingly
if (.Platform$OS.type == "windows") {
  BASE_PATH <- "E:/ASAP/scRNASeq/PerturbSeq/final"
} else {
  BASE_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final"
}

UMAP_OUTPUT_DIR <- file.path(BASE_PATH, "update_analysis_scripts", "iSCORE-PDecipher", "inst", "extdata", "umap_data")
DATASET_NAME <- "iSCORE_PD_CRISPRi"

# Create backup directory with timestamp
backup_dir <- file.path(UMAP_OUTPUT_DIR, paste0("backup_", format(Sys.Date(), "%Y%m%d")))
if (!dir.exists(backup_dir)) {
  dir.create(backup_dir, recursive = TRUE)
  cat("Created backup directory:", backup_dir, "\n")
}

# Find all iSCORE_PD_CRISPRi files
files_to_archive <- list.files(UMAP_OUTPUT_DIR, 
                               pattern = paste0("^", DATASET_NAME, ".*\\.rds$"),
                               full.names = TRUE)

if (length(files_to_archive) > 0) {
  cat("\nFiles to archive:\n")
  for (file in files_to_archive) {
    cat("  -", basename(file), "\n")
  }
  
  cat("\nArchiving files...\n")
  for (file in files_to_archive) {
    backup_file <- file.path(backup_dir, basename(file))
    file.copy(file, backup_file)
    cat("  - Copied", basename(file), "to backup\n")
  }
  
  cat("\nRemoving original files...\n")
  for (file in files_to_archive) {
    file.remove(file)
    cat("  - Removed", basename(file), "\n")
  }
  
  cat("\n✓ Successfully archived", length(files_to_archive), "files to:", backup_dir, "\n")
} else {
  cat("No", DATASET_NAME, "files found to archive.\n")
}

cat("\nDone. Ready for fresh file generation.\n")