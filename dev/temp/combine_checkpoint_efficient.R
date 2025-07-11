# Efficient Combination of Gene Association Checkpoint Files
# Optimized for large datasets with memory management

library(dplyr)

cat("=== Efficient Gene Association Combination ===\n")

# Check existing final file
final_file <- "inst/extdata/gene_term_associations.rds"
if (file.exists(final_file)) {
  cat("Final file already exists. Checking contents...\n")
  existing_data <- readRDS(final_file)
  cat("Existing data:", nrow(existing_data), "associations\n")
  cat("File size:", round(file.size(final_file) / (1024^2), 1), "MB\n")
  
  # Quick validation
  cat("Analysis types:", paste(unique(existing_data$analysis_type), collapse = ", "), "\n")
  cat("Unique terms:", length(unique(existing_data$term_id)), "\n")
  cat("Unique genes:", length(unique(existing_data$gene)), "\n")
  
  # Test lookup functionality
  cat("\n=== Testing Lookup Functionality ===\n")
  sample_rows <- head(existing_data, 3)
  cat("Sample entries:\n")
  for (i in 1:3) {
    row <- sample_rows[i, ]
    cat("  Term:", row$term_id, "| Gene:", row$gene, "| Cluster:", row$cluster, "\n")
    cat("    Genes:", substr(row$associated_genes, 1, 50), "...\n")
  }
  
  cat("\nFinal dataset already available!\n")
  quit(save = "no")
}

cat("Creating final dataset from checkpoint files...\n")

# Load datasets efficiently
extdata_dir <- "inst/extdata"

# Load iSCORE-PD (smaller dataset first)
iscore_file <- file.path(extdata_dir, "gene_associations_iscore.rds")
cat("Loading iSCORE-PD data (22.5 MB)...\n")
iscore_data <- readRDS(iscore_file)
cat("  Loaded:", nrow(iscore_data), "associations\n")

# Find and load latest CRISPRi checkpoint  
crispi_files <- list.files(extdata_dir, pattern = "gene_associations_iscore_crispi\\.rds_temp_\\d+\\.rds$", full.names = TRUE)
batch_nums <- as.numeric(gsub(".*temp_(\\d+)\\.rds$", "\\1", basename(crispi_files)))
latest_file <- crispi_files[which.max(batch_nums)]
latest_batch <- max(batch_nums)

cat("Loading CRISPRi data from batch", latest_batch, "(243.3 MB)...\n")
cat("This may take a moment...\n")
crispi_data <- readRDS(latest_file)
cat("  Loaded:", nrow(crispi_data), "associations\n")

# Efficient combination using data.table approach
cat("\nCombining datasets efficiently...\n")
combined_data <- rbind(iscore_data, crispi_data)

# Memory cleanup
rm(iscore_data, crispi_data)
gc()

cat("Combined dataset:", nrow(combined_data), "associations\n")

# Create summary stats before saving
cat("\n=== Final Dataset Summary ===\n")
stats <- list(
  total_associations = nrow(combined_data),
  unique_terms = length(unique(combined_data$term_id)),
  unique_genes = length(unique(combined_data$gene)),
  analysis_types = sort(unique(combined_data$analysis_type)),
  enrichment_types = sort(unique(combined_data$enrichment_type)),
  directions = sort(unique(combined_data$direction))
)

cat("Total associations:", stats$total_associations, "\n")
cat("Unique terms:", stats$unique_terms, "\n") 
cat("Unique genes:", stats$unique_genes, "\n")
cat("Analysis types:", paste(stats$analysis_types, collapse = ", "), "\n")
cat("Enrichment types:", paste(stats$enrichment_types, collapse = ", "), "\n")
cat("Directions:", paste(stats$directions, collapse = ", "), "\n")

# Save with maximum compression
cat("\nSaving final dataset with XZ compression...\n")
cat("This may take several minutes...\n")
saveRDS(combined_data, final_file, compress = "xz")

# Check final file size
file_size_mb <- round(file.size(final_file) / (1024^2), 1)
cat("Final file size:", file_size_mb, "MB\n")

if (file_size_mb > 100) {
  cat("WARNING: File size (", file_size_mb, "MB) exceeds GitHub 100MB limit!\n")
  cat("Consider data reduction strategies.\n")
} else {
  cat("SUCCESS: File size acceptable for GitHub deployment!\n")
}

# Test data access
cat("\n=== Sample Data Preview ===\n")
sample_data <- head(combined_data, 2)
for (i in 1:2) {
  row <- sample_data[i, ]
  cat("Example", i, ":\n")
  cat("  Term:", row$term_id, "\n")
  cat("  Description:", substr(row$description, 1, 60), "...\n")
  cat("  Gene:", row$gene, "| Cluster:", row$cluster, "\n")
  cat("  Associated genes:", substr(row$associated_genes, 1, 80), "...\n")
  cat("  Gene count:", row$gene_count, "\n\n")
}

cat("=== Gene Association Extraction Complete ===\n")
cat("Ready for deployment and testing!\n")