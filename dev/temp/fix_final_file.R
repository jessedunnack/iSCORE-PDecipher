# Fix Final Gene Associations File
# Recreate with standard compression

library(dplyr)

cat("=== Fixing Gene Associations File ===\n")

# Load from the working checkpoint files
extdata_dir <- "inst/extdata"

# Load iSCORE-PD data
cat("Loading iSCORE-PD data...\n")
iscore_data <- readRDS(file.path(extdata_dir, "gene_associations_iscore.rds"))
cat("  Loaded:", nrow(iscore_data), "associations\n")

# Load latest CRISPRi checkpoint
crispi_files <- list.files(extdata_dir, pattern = "gene_associations_iscore_crispi\\.rds_temp_\\d+\\.rds$", full.names = TRUE)
batch_nums <- as.numeric(gsub(".*temp_(\\d+)\\.rds$", "\\1", basename(crispi_files)))
latest_file <- crispi_files[which.max(batch_nums)]
latest_batch <- max(batch_nums)

cat("Loading CRISPRi data from batch", latest_batch, "...\n")
crispi_data <- readRDS(latest_file)
cat("  Loaded:", nrow(crispi_data), "associations\n")

# Add composite key if missing
if (!"composite_key" %in% colnames(iscore_data)) {
  cat("Adding composite keys to iSCORE data...\n")
  iscore_data$composite_key <- paste(
    iscore_data$analysis_type, 
    iscore_data$gene, 
    iscore_data$cluster,
    iscore_data$enrichment_type, 
    iscore_data$direction, 
    "default",  # experiment
    iscore_data$term_id, 
    sep = "|"
  )
}

if (!"composite_key" %in% colnames(crispi_data)) {
  cat("Adding composite keys to CRISPRi data...\n")
  crispi_data$composite_key <- paste(
    crispi_data$analysis_type, 
    crispi_data$gene, 
    crispi_data$cluster,
    crispi_data$enrichment_type, 
    crispi_data$direction, 
    "default",  # experiment
    crispi_data$term_id, 
    sep = "|"
  )
}

# Combine efficiently
cat("Combining datasets...\n")
combined_data <- rbind(iscore_data, crispi_data)

# Clean up memory
rm(iscore_data, crispi_data)
gc()

cat("Combined:", nrow(combined_data), "associations\n")

# Summary
cat("\n=== Dataset Summary ===\n")
cat("Total associations:", nrow(combined_data), "\n")
cat("Unique terms:", length(unique(combined_data$term_id)), "\n")
cat("Unique genes:", length(unique(combined_data$gene)), "\n")
cat("Analysis types:", paste(sort(unique(combined_data$analysis_type)), collapse = ", "), "\n")

# Save with standard compression (faster and more reliable)
output_file <- "inst/extdata/gene_term_associations.rds"
cat("\nSaving with standard compression...\n")
saveRDS(combined_data, output_file, compress = TRUE)

# Check file size
file_size_mb <- round(file.size(output_file) / (1024^2), 1)
cat("Final file size:", file_size_mb, "MB\n")

if (file_size_mb > 100) {
  cat("WARNING: File size may be too large for GitHub!\n")
} else {
  cat("SUCCESS: File size acceptable for GitHub!\n")
}

# Quick test
cat("\n=== Quick Test ===\n")
test_data <- readRDS(output_file)
cat("Successfully reloaded:", nrow(test_data), "associations\n")

# Show sample
sample_row <- test_data[1, ]
cat("\nSample data:\n")
cat("  Term:", sample_row$term_id, "\n")
cat("  Gene:", sample_row$gene, "| Cluster:", sample_row$cluster, "\n")
cat("  Associated genes:", substr(sample_row$associated_genes, 1, 50), "...\n")

cat("\n=== File Fixed Successfully ===\n")