#!/usr/bin/env Rscript

# Inspect the structure of our UMAP file without SingleCellExperiment

cat("=== INSPECTING UMAP FILE STRUCTURE ===\n\n")

umap_file <- "inst/extdata/umap_data/iSCORE_PD_CRISPRi_umap_data_30pc.rds"

if (file.exists(umap_file)) {
  cat("Loading file:", umap_file, "\n")
  
  # Read the raw RDS structure
  data <- readRDS(umap_file)
  
  cat("\n1. OBJECT CLASS:\n")
  cat("   - Class:", class(data), "\n")
  
  if ("colData" %in% slotNames(data)) {
    cat("\n2. COLDATA STRUCTURE:\n")
    coldata <- data@colData
    cat("   - Number of rows:", nrow(coldata), "\n")
    cat("   - Column names:", paste(colnames(coldata), collapse=", "), "\n")
    
    if ("seurat_clusters" %in% colnames(coldata)) {
      clusters <- unique(coldata$seurat_clusters)
      cat("   - Number of unique clusters:", length(clusters), "\n")
      cat("   - First 10 clusters:", paste(head(sort(as.character(clusters)), 10), collapse=", "), "\n")
      cat("   - Cluster class:", class(clusters), "\n")
    } else {
      cat("   - seurat_clusters column NOT FOUND!\n")
    }
  }
  
  if ("reducedDims" %in% slotNames(data)) {
    cat("\n3. REDUCED DIMS:\n")
    rdims <- data@reducedDims
    cat("   - Available reductions:", paste(names(rdims), collapse=", "), "\n")
    if ("UMAP" %in% names(rdims)) {
      umap_coords <- rdims$UMAP
      cat("   - UMAP dimensions:", dim(umap_coords), "\n")
    }
  }
  
} else {
  cat("File not found:", umap_file, "\n")
}

cat("\nDone.\n")