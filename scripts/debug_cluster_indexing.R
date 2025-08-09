#!/usr/bin/env Rscript

# Debug script to understand cluster indexing

library(Seurat)
library(dplyr)

# Load data
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds")
annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv")

cat("=== Debugging Cluster Indexing ===\n\n")

# Check Seurat clusters
cat("Seurat cluster info:\n")
cat("- Class of seurat_clusters_fine:", class(seurat_obj$seurat_clusters_fine), "\n")
cat("- Unique values:", sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_fine)))), "\n")
cat("- First 10 values:", head(seurat_obj$seurat_clusters_fine, 10), "\n")

# Check annotations
cat("\nAnnotation cluster info:\n")
cat("- Unique values in annotations:", sort(unique(annotations$fine_cluster)), "\n")

# Check if they match
seurat_clusters_numeric <- as.numeric(as.character(seurat_obj$seurat_clusters_fine))
cat("\nMatching check:\n")
cat("- Min Seurat cluster:", min(seurat_clusters_numeric), "\n")
cat("- Max Seurat cluster:", max(seurat_clusters_numeric), "\n")
cat("- Min annotation cluster:", min(annotations$fine_cluster), "\n")
cat("- Max annotation cluster:", max(annotations$fine_cluster), "\n")

# Test lookup
test_lookup <- setNames(annotations$cell_type, as.character(annotations$fine_cluster))
cat("\nTest lookup creation:\n")
cat("- Lookup names:", names(test_lookup)[1:5], "\n")
cat("- Seurat cluster values as character:", as.character(seurat_obj$seurat_clusters_fine)[1:5], "\n")

# Try the lookup
test_result <- test_lookup[as.character(seurat_obj$seurat_clusters_fine)]
cat("\nLookup test:\n")
cat("- NA count:", sum(is.na(test_result)), "\n")
cat("- Non-NA count:", sum(!is.na(test_result)), "\n")