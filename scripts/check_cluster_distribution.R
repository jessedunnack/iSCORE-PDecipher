#!/usr/bin/env Rscript

# Quick script to check cluster distribution
library(Seurat)
library(dplyr)

# Load Seurat object
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds")

# Check seurat_clusters_fine
cat("=== seurat_clusters_fine distribution ===\n")
cluster_table <- table(seurat_obj$seurat_clusters_fine)
cat("Number of unique clusters:", length(cluster_table), "\n")
cat("Cluster range:", min(as.numeric(names(cluster_table))), "-", max(as.numeric(names(cluster_table))), "\n\n")

# Show first 20 clusters with cell counts
cat("Cluster sizes (first 20):\n")
print(head(sort(cluster_table, decreasing = TRUE), 20))

# Check other potential fine cluster columns
cat("\n\n=== Checking other high-resolution clustering columns ===\n")
high_res_cols <- grep("res\\.(0\\.[5-9]|[1-9])", colnames(seurat_obj@meta.data), value = TRUE)

for (col in high_res_cols) {
  n_clusters <- length(unique(seurat_obj@meta.data[[col]]))
  cat(sprintf("%-20s: %d clusters\n", col, n_clusters))
}

# Find the column with exactly 36 clusters
cat("\n\n=== Columns with exactly 36 clusters ===\n")
all_cluster_cols <- grep("cluster|res\\.", colnames(seurat_obj@meta.data), value = TRUE)
for (col in all_cluster_cols) {
  n_clusters <- length(unique(seurat_obj@meta.data[[col]]))
  if (n_clusters == 36) {
    cat(sprintf("%-20s: %d clusters\n", col, n_clusters))
  }
}

# Check SCT_snn_res.0.8 specifically
if ("SCT_snn_res.0.8" %in% colnames(seurat_obj@meta.data)) {
  cat("\n\n=== SCT_snn_res.0.8 (36 clusters expected) ===\n")
  res08_table <- table(seurat_obj$SCT_snn_res.0.8)
  cat("Number of clusters:", length(res08_table), "\n")
  cat("Cluster range:", min(as.numeric(names(res08_table))), "-", max(as.numeric(names(res08_table))), "\n")
}