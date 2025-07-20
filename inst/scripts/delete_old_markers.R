#!/usr/bin/env Rscript

# Simple script to delete old marker files
# After running this, the regular extract_umap_data.R will recalculate markers

cat("=====================================\n")
cat("Delete Old Marker Files\n")
cat("=====================================\n")

# Path to marker files
MARKER_DIR <- file.path(dirname(dirname(getwd())), "inst", "extdata", "umap_data")

# If running from scripts directory
if (basename(getwd()) == "scripts") {
  MARKER_DIR <- file.path(dirname(getwd()), "extdata", "umap_data")
}

cat("\nLooking for marker files in:", MARKER_DIR, "\n")

# List existing marker files
marker_files <- list.files(MARKER_DIR, pattern = "_cluster_markers\\.rds$", full.names = TRUE)

if (length(marker_files) == 0) {
  cat("\nNo marker files found.\n")
  quit("no")
}

cat("\nFound marker files:\n")
for (f in marker_files) {
  file_info <- file.info(f)
  cat(sprintf("  - %s (modified: %s, size: %.1f KB)\n", 
              basename(f), 
              format(file_info$mtime, "%Y-%m-%d %H:%M"), 
              file_info$size/1024))
}

# Check with user
if (interactive()) {
  cat("\nThese files contain markers calculated with OLD parameters:\n")
  cat("  - LFC threshold: 0.25\n")
  cat("  - Min.pct: 0.1\n") 
  cat("  - Only positive markers\n")
  cat("\nDo you want to delete them? (y/n): ")
  
  response <- readline()
  if (tolower(response) != "y") {
    cat("Cancelled.\n")
    quit("no")
  }
}

# Create backup
backup_dir <- file.path(MARKER_DIR, paste0("backup_markers_", format(Sys.time(), "%Y%m%d_%H%M%S")))
dir.create(backup_dir, showWarnings = FALSE)

cat("\nBacking up to:", backup_dir, "\n")
for (f in marker_files) {
  file.copy(f, file.path(backup_dir, basename(f)))
  cat("  - Backed up", basename(f), "\n")
}

# Delete files
cat("\nDeleting old marker files...\n")
for (f in marker_files) {
  file.remove(f)
  cat("  - Deleted", basename(f), "\n")
}

cat("\n✓ Old marker files deleted!\n")
cat("\nNext steps:\n")
cat("1. Run: source('extract_umap_data.R'); main()\n")
cat("2. This will recalculate markers with NEW parameters:\n")
cat("   - LFC threshold: 0.5\n")
cat("   - Min.pct: 0.25\n")
cat("   - Min.diff.pct: 0.2\n")
cat("   - Both positive and negative markers\n")