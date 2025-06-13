#!/usr/bin/env Rscript

# Quick clustering fix for large datasets (200K+ cells)
# Focuses on resolution range 1.5-1.8 for efficiency

library(Seurat)

# Fast resolution finder
find_resolution_fast <- function(seurat_file, target_clusters, resolution_range = c(1.5, 1.8), step = 0.05) {
  
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_file)
  
  # Current state
  current_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
  cat("Current clusters:", current_clusters, "\n")
  cat("Target clusters:", target_clusters, "\n")
  
  if (current_clusters == target_clusters) {
    cat("✓ Already has correct number of clusters!\n")
    return(NULL)
  }
  
  cat("\nTesting resolutions from", resolution_range[1], "to", resolution_range[2], "...\n\n")
  
  # Test resolutions
  results <- data.frame()
  for (res in seq(resolution_range[1], resolution_range[2], by = step)) {
    cat("Resolution", sprintf("%.2f", res), ": ")
    
    # Run clustering
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    n_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
    
    cat(n_clusters, "clusters")
    if (n_clusters == target_clusters) cat(" ✓ MATCH!")
    cat("\n")
    
    results <- rbind(results, data.frame(resolution = res, n_clusters = n_clusters))
    
    # Stop if found
    if (n_clusters == target_clusters) {
      cat("\n✓ Found matching resolution:", res, "\n")
      
      # Save option
      cat("\nSave with this clustering? (y/n): ")
      if (tolower(readline()) == "y") {
        saveRDS(seurat_obj, seurat_file)
        cat("✓ Saved!\n")
      }
      return(res)
    }
  }
  
  # No exact match - show closest
  cat("\n✗ No exact match in range", resolution_range[1], "-", resolution_range[2], "\n")
  print(results)
  
  # Suggest next step
  if (all(results$n_clusters > target_clusters)) {
    cat("\nTry lower resolutions (< ", resolution_range[1], ")\n", sep="")
  } else if (all(results$n_clusters < target_clusters)) {
    cat("\nTry higher resolutions (> ", resolution_range[2], ")\n", sep="")
  }
  
  return(NULL)
}

# Quick fix for each dataset
fix_datasets <- function() {
  datasets <- list(
    iSCORE_PD = list(
      file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/final_iSCORE-PD.rds",
      expected = 15  # From your DE analysis
    ),
    iSCORE_PD_CRISPRi = list(
      file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
      expected = NULL  # Will check DE data
    ),
    Full_Dataset = list(
      file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_dataset.rds",
      expected = NULL  # Will check DE data
    )
  )
  
  for (name in names(datasets)) {
    cat("\n━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    cat("Dataset:", name, "\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    
    if (!is.null(datasets[[name]]$expected)) {
      find_resolution_fast(
        datasets[[name]]$file, 
        datasets[[name]]$expected,
        resolution_range = c(1.5, 1.8),
        step = 0.05
      )
    } else {
      cat("Need to check DE/enrichment data for expected clusters\n")
      cat("Run: check_expected_clusters('", name, "')\n", sep="")
    }
  }
}

# Check expected clusters from DE data
check_expected_clusters <- function(dataset_name) {
  de_files <- list(
    iSCORE_PD = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/full_DE_results.rds",
    iSCORE_PD_CRISPRi = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds",
    Full_Dataset = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_DE_results.rds"
  )
  
  de_file <- de_files[[dataset_name]]
  if (file.exists(de_file)) {
    de_data <- readRDS(de_file)
    
    # Check MAST clusters
    if ("iSCORE_PD_MAST" %in% names(de_data)) {
      clusters <- c()
      for (gene in names(de_data$iSCORE_PD_MAST)) {
        clusters <- union(clusters, names(de_data$iSCORE_PD_MAST[[gene]]))
      }
      cat("MAST clusters (", length(clusters), "):", paste(sort(clusters), collapse=", "), "\n")
    }
    
    # Check MixScale clusters
    if ("CRISPRi_Mixscale" %in% names(de_data)) {
      clusters <- c()
      for (gene in names(de_data$CRISPRi_Mixscale)) {
        clusters <- union(clusters, names(de_data$CRISPRi_Mixscale[[gene]]))
      }
      cat("MixScale clusters (", length(clusters), "):", paste(sort(clusters), collapse=", "), "\n")
    }
  }
}

# Main
if (interactive()) {
  cat("=== Quick Clustering Fix (200K+ cells) ===\n\n")
  cat("Functions:\n")
  cat("  find_resolution_fast(seurat_file, target_clusters) - Fix single dataset\n")
  cat("  fix_datasets() - Check all datasets\n")
  cat("  check_expected_clusters('dataset_name') - See what DE data expects\n")
  cat("\nExample:\n")
  cat("  find_resolution_fast('path/to/seurat.rds', 15)\n")
  cat("\nThis focuses on resolution range 1.5-1.8 for speed.\n")
}