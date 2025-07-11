# Create Working Gene Associations Dataset
# Use the largest checkpoint file as our working dataset

library(dplyr)

cat("=== Creating Working Gene Associations Dataset ===\n")

# Use the largest and most complete checkpoint file
# This represents batch 190 of CRISPRi processing (8.6M associations)
extdata_dir <- "inst/extdata"
source_file <- file.path(extdata_dir, "gene_associations_iscore_crispi.rds_temp_190.rds")
target_file <- file.path(extdata_dir, "gene_term_associations.rds")

cat("Source file:", source_file, "\n")
cat("Target file:", target_file, "\n")

# Check source file
if (!file.exists(source_file)) {
  stop("Source file not found: ", source_file)
}

source_size <- round(file.size(source_file) / (1024^2), 1)
cat("Source file size:", source_size, "MB\n")

# Load and validate
cat("Loading source data...\n")
data <- readRDS(source_file)
cat("Loaded:", nrow(data), "associations\n")

# Ensure composite key exists
if (!"composite_key" %in% colnames(data)) {
  cat("Adding composite keys...\n")
  data$composite_key <- paste(
    data$analysis_type, 
    data$gene, 
    data$cluster,
    data$enrichment_type, 
    data$direction, 
    "default",  # experiment
    data$term_id, 
    sep = "|"
  )
}

# Summary stats
cat("\n=== Dataset Summary ===\n")
cat("Total associations:", nrow(data), "\n")
cat("Unique terms:", length(unique(data$term_id)), "\n")
cat("Unique genes:", length(unique(data$gene)), "\n")
cat("Analysis types:", paste(sort(unique(data$analysis_type)), collapse = ", "), "\n")
cat("Enrichment types:", paste(sort(unique(data$enrichment_type)), collapse = ", "), "\n")
cat("Directions:", paste(sort(unique(data$direction)), collapse = ", "), "\n")

# Sample genes by analysis type
cat("\nGenes by analysis type:\n")
for (analysis in unique(data$analysis_type)) {
  genes <- sort(unique(data$gene[data$analysis_type == analysis]))
  cat("  ", analysis, ":", paste(genes, collapse = ", "), "\n")
}

# Save as working dataset
cat("\nSaving as working dataset...\n")
saveRDS(data, target_file, compress = TRUE)

# Verify
final_size <- round(file.size(target_file) / (1024^2), 1)
cat("Final file size:", final_size, "MB\n")

# GitHub compatibility
if (final_size > 100) {
  cat("WARNING: File size exceeds GitHub 100MB limit\n")
} else {
  cat("SUCCESS: File size acceptable for GitHub (< 100MB)\n")
}

# Test reload
cat("\nTesting file integrity...\n")
test_data <- readRDS(target_file)
cat("Successfully reloaded:", nrow(test_data), "associations\n")

# Show sample data
cat("\n=== Sample Data ===\n")
sample_rows <- head(test_data, 2)
for (i in 1:2) {
  row <- sample_rows[i, ]
  cat("Example", i, ":\n")
  cat("  Term:", row$term_id, "\n")
  cat("  Description:", substr(row$description, 1, 50), "...\n")
  cat("  Gene:", row$gene, "| Cluster:", row$cluster, "| Type:", row$enrichment_type, "\n")
  cat("  Associated genes (", row$gene_count, "):", substr(row$associated_genes, 1, 60), "...\n")
  cat("  Composite key:", substr(row$composite_key, 1, 60), "...\n\n")
}

cat("=== Working Dataset Ready ===\n")
cat("Gene association lookup system is ready for testing!\n")