#!/usr/bin/env Rscript

# Ultra-quick clustering test for exactly 3 resolutions
# Tests only: 0.125, 0.15, 0.175

library(Seurat)

test_three_resolutions <- function(seurat_file, target_clusters) {
  
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_file)
  
  # Current state
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    current_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
    cat("Current clusters:", current_clusters, "\n")
  }
  cat("Target clusters:", target_clusters, "\n\n")
  
  # Test exactly these 3 resolutions
  resolutions <- c(0.125, 0.15, 0.175)
  results <- list()
  match_found <- FALSE
  
  for (res in resolutions) {
    cat("Testing resolution", res, "... ")
    
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    n_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
    
    cat(n_clusters, "clusters")
    
    if (n_clusters == target_clusters) {
      cat(" ✓ MATCH!")
      match_found <- TRUE
      
      # Save immediately
      cat("\n\nFound match! Save updated object? (y/n): ")
      if (tolower(readline()) == "y") {
        saveRDS(seurat_obj, seurat_file)
        cat("✓ Saved!\n")
      }
      return(res)
    }
    cat("\n")
    
    results[[as.character(res)]] <- n_clusters
  }
  
  # No match found
  cat("\n✗ None of the resolutions (0.125, 0.15, 0.175) gave", target_clusters, "clusters\n")
  cat("\nResults:\n")
  for (res in names(results)) {
    cat("  Resolution", res, ":", results[[res]], "clusters\n")
  }
  
  return(NULL)
}

# Quick check what DE data expects
check_de_clusters <- function(dataset_name) {
  de_files <- list(
    iSCORE_PD = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/full_DE_results.rds",
    iSCORE_PD_CRISPRi = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds",
    Full_Dataset = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_DE_results.rds"
  )
  
  de_file <- de_files[[dataset_name]]
  if (file.exists(de_file)) {
    de_data <- readRDS(de_file)
    
    cat("Checking", dataset_name, "DE data...\n")
    
    # Get unique clusters from first gene
    if ("iSCORE_PD_MAST" %in% names(de_data) && length(de_data$iSCORE_PD_MAST) > 0) {
      first_gene <- names(de_data$iSCORE_PD_MAST)[1]
      clusters <- names(de_data$iSCORE_PD_MAST[[first_gene]])
      cat("Expected clusters:", length(clusters), "\n")
      cat("Clusters:", paste(clusters, collapse=", "), "\n")
      return(length(clusters))
    }
    
    if ("CRISPRi_Mixscale" %in% names(de_data) && length(de_data$CRISPRi_Mixscale) > 0) {
      first_gene <- names(de_data$CRISPRi_Mixscale)[1]
      clusters <- names(de_data$CRISPRi_Mixscale[[first_gene]])
      cat("Expected clusters:", length(clusters), "\n")
      cat("Clusters:", paste(clusters, collapse=", "), "\n")
      return(length(clusters))
    }
  }
  return(NULL)
}

# Main
if (interactive()) {
  cat("=== Test 3 Resolutions (0.125, 0.15, 0.175) ===\n\n")
  cat("Usage:\n")
  cat("  test_three_resolutions('path/to/seurat.rds', target_clusters)\n")
  cat("  check_de_clusters('dataset_name')\n")
  cat("\nExample:\n")
  cat("  # First check what DE expects\n")
  cat("  n <- check_de_clusters('iSCORE_PD')\n")
  cat("  # Then test the 3 resolutions\n")
  cat("  test_three_resolutions('E:/ASAP/.../final_iSCORE-PD.rds', n)\n")
}