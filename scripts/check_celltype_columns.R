#!/usr/bin/env Rscript

# Simple script to check celltype columns and their alignment

library(Seurat)
library(dplyr)

cat("Loading Seurat object...\n")
x <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")

cat("\nTotal cells:", ncol(x), "\n")
cat("Unique clusters in seurat_clusters_fine:", paste(sort(unique(x$seurat_clusters_fine)), collapse=", "), "\n\n")

# Check all celltype-related columns
celltype_cols <- grep("celltype|validated", colnames(x@meta.data), value = TRUE)
cat("Celltype-related columns found:\n")
for (col in celltype_cols) {
  cat("\n", col, ":\n")
  # Check if column has any non-NA values
  non_na <- sum(!is.na(x@meta.data[[col]]) & x@meta.data[[col]] != "")
  cat("  Non-empty values:", non_na, "/", ncol(x), "\n")
  
  # Show unique values (up to 10)
  unique_vals <- unique(x@meta.data[[col]])
  unique_vals <- unique_vals[!is.na(unique_vals) & unique_vals != ""]
  if (length(unique_vals) > 0) {
    cat("  Unique values:", paste(head(unique_vals, 10), collapse=", "))
    if (length(unique_vals) > 10) cat("...")
    cat("\n")
  }
}

# Check specific clusters
cat("\n\nChecking specific clusters:\n")
cat("============================\n")

# Cluster 5 (should be Dopaminergic)
cat("\nCluster 5 analysis:\n")
cluster5_cells <- which(x$seurat_clusters_fine == 5)
cat("Number of cells:", length(cluster5_cells), "\n")

for (col in celltype_cols) {
  vals <- x@meta.data[cluster5_cells, col]
  non_empty <- vals[!is.na(vals) & vals != ""]
  if (length(non_empty) > 0) {
    cat("\n", col, ":\n")
    print(table(non_empty))
  }
}

# Cluster 4 (should be Choroid Plexus)
cat("\n\nCluster 4 analysis:\n")
cluster4_cells <- which(x$seurat_clusters_fine == 4)
cat("Number of cells:", length(cluster4_cells), "\n")

for (col in celltype_cols) {
  vals <- x@meta.data[cluster4_cells, col]
  non_empty <- vals[!is.na(vals) & vals != ""]
  if (length(non_empty) > 0) {
    cat("\n", col, ":\n")
    print(table(non_empty))
  }
}

# Check marker genes availability
cat("\n\nChecking key marker genes:\n")
markers <- c("TH", "TTR", "TUBB3", "SOX2", "FOXA2")
for (marker in markers) {
  if (marker %in% rownames(x)) {
    cat(marker, "- present\n")
  } else {
    cat(marker, "- NOT FOUND\n")
  }
}

# Check available assays
cat("\nAvailable assays:", paste(names(x@assays), collapse=", "), "\n")