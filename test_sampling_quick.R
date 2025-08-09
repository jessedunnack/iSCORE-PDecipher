#!/usr/bin/env Rscript

# Quick test of data sampling functions
cat("Quick Data Sampling Test\n")
cat("========================\n\n")

# Source functions
source("R/data_sampling.R")

# Test memory estimation function
cat("Testing memory estimation...\n")
test_stats <- estimate_memory_usage(list(
  meta.data = data.frame(a = 1:1000),
  reductions = list()
), include_assays = FALSE)

cat(sprintf("Memory stats calculated: %.2f MB\n", test_stats$metadata_mb))

# Test UMAP extraction
cat("\nTesting UMAP extraction function...\n")
# Create mock Seurat-like object
mock_obj <- list(
  reductions = list(
    umap = list(
      cell.embeddings = matrix(rnorm(200), ncol = 2, 
                               dimnames = list(paste0("Cell", 1:100), c("UMAP_1", "UMAP_2")))
    )
  ),
  meta.data = data.frame(
    seurat_clusters = factor(sample(0:3, 100, replace = TRUE)),
    orig.ident = rep("Sample1", 100)
  )
)
class(mock_obj) <- "Seurat"

# Extract UMAP
umap_data <- extract_umap_data(mock_obj, sample_n = 50)
cat(sprintf("Extracted %d cells with %d columns\n", nrow(umap_data), ncol(umap_data)))

cat("\n✅ Quick test passed!\n")