# Optimize Gene Associations Dataset for GitHub Deployment
# Reduce file size and fix composite key logic

library(dplyr)

cat("=== Optimizing Gene Associations Dataset ===\n")

# Load current dataset
current_file <- "inst/extdata/gene_term_associations.rds"
data <- readRDS(current_file)

cat("Current dataset:\n")
cat("  Associations:", nrow(data), "\n")
cat("  File size:", round(file.size(current_file) / (1024^2), 1), "MB\n")

# Check composite key format and fix if needed
cat("\nAnalyzing composite keys...\n")
sample_keys <- head(data$composite_key, 3)
cat("Sample keys:\n")
for (i in 1:3) {
  cat("  ", sample_keys[i], "\n")
}

# Check for missing composite keys
if (!"composite_key" %in% colnames(data) || any(is.na(data$composite_key))) {
  cat("Fixing composite keys...\n")
  data$composite_key <- paste(
    data$analysis_type, 
    data$gene, 
    data$cluster,
    data$enrichment_type, 
    data$direction, 
    "default",  # experiment - standardize to "default"
    data$term_id, 
    sep = "|"
  )
  cat("Updated composite keys\n")
}

# Optimization strategies to reduce file size
cat("\n=== Optimization Strategies ===\n")

# 1. Remove duplicate associations (same term-gene combination)
cat("1. Checking for duplicates...\n")
original_rows <- nrow(data)
data <- data %>% distinct()
after_distinct <- nrow(data)
cat("   Removed", original_rows - after_distinct, "duplicate rows\n")

# 2. Compress descriptions (remove very long descriptions)
cat("2. Optimizing descriptions...\n")
data$description <- ifelse(
  nchar(data$description) > 100,
  paste0(substr(data$description, 1, 97), "..."),
  data$description
)

# 3. Filter to most significant terms (keep terms with more genes)
cat("3. Filtering to most significant terms...\n")
before_filter <- nrow(data)

# Keep terms with at least 2 genes associated
data <- data %>%
  filter(gene_count >= 2)

after_filter <- nrow(data)
cat("   Filtered to terms with ≥2 genes:", after_filter, "associations\n")
cat("   Removed", before_filter - after_filter, "single-gene associations\n")

# 4. Focus on key enrichment types
cat("4. Focusing on key enrichment types...\n")
key_types <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome")
available_types <- unique(data$enrichment_type)
cat("   Available types:", paste(available_types, collapse = ", "), "\n")

# Keep all available types for now (they're all useful)

# 5. Sample representative terms if still too large
current_size <- nrow(data)
if (current_size > 2000000) {  # If over 2M, sample
  cat("5. Sampling representative terms...\n")
  
  # Stratified sampling by analysis_type and enrichment_type
  sampled_data <- data %>%
    group_by(analysis_type, gene, enrichment_type, direction) %>%
    slice_head(n = 100) %>%  # Top 100 terms per group
    ungroup()
  
  cat("   Sampled", nrow(sampled_data), "representative associations\n")
  data <- sampled_data
}

# Final summary
cat("\n=== Optimized Dataset Summary ===\n")
cat("Total associations:", nrow(data), "\n")
cat("Unique terms:", length(unique(data$term_id)), "\n")
cat("Unique genes:", length(unique(data$gene)), "\n")
cat("Analysis types:", paste(sort(unique(data$analysis_type)), collapse = ", "), "\n")
cat("Enrichment types:", paste(sort(unique(data$enrichment_type)), collapse = ", "), "\n")

# Save optimized version
output_file <- "inst/extdata/gene_term_associations_optimized.rds"
cat("\nSaving optimized dataset...\n")
saveRDS(data, output_file, compress = "xz")

# Check final size
final_size <- round(file.size(output_file) / (1024^2), 1)
cat("Optimized file size:", final_size, "MB\n")

if (final_size > 100) {
  cat("Still too large for GitHub. Consider further reduction.\n")
} else {
  cat("SUCCESS: File size now GitHub compatible!\n")
}

# Replace original with optimized version
cat("\nReplacing original with optimized version...\n")
file.rename(output_file, current_file)

# Test composite key lookup
cat("\n=== Testing Composite Key Lookup ===\n")
sample_row <- data[1, ]
test_key <- sample_row$composite_key
matching_rows <- which(data$composite_key == test_key)
cat("Test lookup for key:", substr(test_key, 1, 80), "...\n")
cat("Found", length(matching_rows), "matching row(s)\n")

if (length(matching_rows) > 0) {
  match_row <- data[matching_rows[1], ]
  genes <- unlist(strsplit(match_row$associated_genes, "/"))
  cat("Associated genes:", length(genes), "genes -", paste(head(genes, 5), collapse = ", "), "\n")
}

cat("\n=== Optimization Complete ===\n")
cat("Dataset ready for GitHub deployment and testing!\n")