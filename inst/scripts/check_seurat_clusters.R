#!/usr/bin/env Rscript

# Quick script to check the Seurat object's cluster status

library(Seurat)

cat("=== CHECKING SEURAT OBJECT CLUSTERS ===\n\n")

# Load the Seurat object
seurat_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"

cat("Loading Seurat object...\n")
seurat_obj <- readRDS(seurat_file)

cat("\n1. DEFAULT IDENTITIES:\n")
cat("   - Number of unique Idents:", length(unique(Idents(seurat_obj))), "\n")
cat("   - First 10 Idents:", paste(head(sort(unique(as.character(Idents(seurat_obj)))), 10), collapse=", "), "\n")

cat("\n2. SEURAT_CLUSTERS IN METADATA:\n")
cat("   - Number of unique clusters:", length(unique(seurat_obj@meta.data$seurat_clusters)), "\n")
cat("   - Cluster labels:", paste(sort(unique(as.character(seurat_obj@meta.data$seurat_clusters))), collapse=", "), "\n")

cat("\n3. AVAILABLE REDUCTIONS:\n")
cat("   - Reductions:", paste(names(seurat_obj@reductions), collapse=", "), "\n")

cat("\n4. RESOLUTION PARAMETERS:\n")
# Check for resolution parameters in misc slot
if ("FindClusters" %in% names(seurat_obj@misc)) {
  cat("   - FindClusters parameters found\n")
}

# Check commands history
if (length(seurat_obj@commands) > 0) {
  cluster_commands <- grep("FindClusters", names(seurat_obj@commands), value = TRUE)
  if (length(cluster_commands) > 0) {
    cat("   - FindClusters commands found:", length(cluster_commands), "\n")
    # Show the most recent one
    latest_cmd <- cluster_commands[length(cluster_commands)]
    cat("   - Latest clustering command:", latest_cmd, "\n")
  }
}

cat("\nDone.\n")