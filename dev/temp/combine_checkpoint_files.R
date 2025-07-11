# Combine Checkpoint Files from Gene Association Extraction
# This script combines the existing checkpoint files into the final dataset

library(dplyr)

cat("=== Combining Gene Association Checkpoint Files ===\n")

# Check what checkpoint files exist
extdata_dir <- "inst/extdata"
checkpoint_files <- list.files(extdata_dir, pattern = "gene_associations.*\\.rds$", full.names = TRUE)

cat("Found checkpoint files:\n")
for (file in checkpoint_files) {
  file_info <- file.info(file)
  size_mb <- round(file_info$size / (1024^2), 1)
  cat("  ", basename(file), ":", size_mb, "MB\n")
}

# Load and check the main datasets
cat("\n=== Loading Main Dataset Files ===\n")

# Load iSCORE-PD data (complete)
iscore_file <- file.path(extdata_dir, "gene_associations_iscore.rds")
if (file.exists(iscore_file)) {
  cat("Loading iSCORE-PD data...\n")
  iscore_data <- readRDS(iscore_file)
  cat("  iSCORE-PD associations:", nrow(iscore_data), "\n")
  cat("  Unique terms:", length(unique(iscore_data$term_id)), "\n")
  cat("  Unique genes:", length(unique(iscore_data$gene)), "\n")
} else {
  cat("Warning: iSCORE-PD file not found!\n")
  iscore_data <- NULL
}

# Find the latest CRISPRi checkpoint
crispi_temp_files <- list.files(extdata_dir, pattern = "gene_associations_iscore_crispi\\.rds_temp_.*\\.rds$", full.names = TRUE)
if (length(crispi_temp_files) > 0) {
  # Get the highest batch number
  batch_nums <- as.numeric(gsub(".*temp_(\\d+)\\.rds$", "\\1", crispi_temp_files))
  latest_batch <- max(batch_nums)
  latest_crispi_file <- crispi_temp_files[which.max(batch_nums)]
  
  cat("Loading iSCORE-PD_plus_CRISPRi data from batch", latest_batch, "...\n")
  crispi_data <- readRDS(latest_crispi_file)
  cat("  iSCORE-PD_plus_CRISPRi associations:", nrow(crispi_data), "\n")
  cat("  Unique terms:", length(unique(crispi_data$term_id)), "\n")
  cat("  Unique genes:", length(unique(crispi_data$gene)), "\n")
} else {
  cat("Warning: No CRISPRi checkpoint files found!\n")
  crispi_data <- NULL
}

# Combine datasets
cat("\n=== Combining Datasets ===\n")
if (!is.null(iscore_data) && !is.null(crispi_data)) {
  # Check column compatibility
  iscore_cols <- colnames(iscore_data)
  crispi_cols <- colnames(crispi_data)
  
  cat("iSCORE-PD columns:", paste(iscore_cols, collapse = ", "), "\n")
  cat("CRISPRi columns:", paste(crispi_cols, collapse = ", "), "\n")
  
  if (identical(iscore_cols, crispi_cols)) {
    cat("Columns are compatible, combining...\n")
    combined_data <- bind_rows(iscore_data, crispi_data)
  } else {
    cat("Columns differ, using intersection...\n")
    common_cols <- intersect(iscore_cols, crispi_cols)
    cat("Common columns:", paste(common_cols, collapse = ", "), "\n")
    combined_data <- bind_rows(
      iscore_data[, common_cols],
      crispi_data[, common_cols]
    )
  }
} else if (!is.null(iscore_data)) {
  cat("Using only iSCORE-PD data...\n")
  combined_data <- iscore_data
} else if (!is.null(crispi_data)) {
  cat("Using only CRISPRi data...\n")
  combined_data <- crispi_data
} else {
  stop("No data found in checkpoint files!")
}

cat("\n=== Combined Dataset Summary ===\n")
cat("Total associations:", nrow(combined_data), "\n")
cat("Unique terms:", length(unique(combined_data$term_id)), "\n")
cat("Unique genes:", length(unique(combined_data$gene)), "\n")
cat("Analysis types:", paste(sort(unique(combined_data$analysis_type)), collapse = ", "), "\n")
cat("Enrichment types:", paste(sort(unique(combined_data$enrichment_type)), collapse = ", "), "\n")
cat("Directions:", paste(sort(unique(combined_data$direction)), collapse = ", "), "\n")

# Save final combined dataset
output_file <- file.path(extdata_dir, "gene_term_associations.rds")
cat("\nSaving combined dataset to:", output_file, "\n")
saveRDS(combined_data, output_file, compress = "xz")

# Check file size
file_size_mb <- round(file.size(output_file) / (1024^2), 1)
cat("Final file size:", file_size_mb, "MB\n")

if (file_size_mb > 100) {
  cat("WARNING: File size may be too large for GitHub deployment!\n")
  cat("Consider further optimization or data reduction.\n")
} else {
  cat("SUCCESS: File size acceptable for GitHub deployment!\n")
}

cat("\n=== Sample Data ===\n")
if (nrow(combined_data) > 0) {
  sample_data <- head(combined_data, 3)
  print(sample_data)
}

cat("\n=== Combination Complete ===\n")