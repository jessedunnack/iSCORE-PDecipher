#!/usr/bin/env Rscript

# Interactive script to fix clustering mismatches
# Run this when clusters in Seurat object don't match your analysis data

library(Seurat)

# Function to interactively fix clustering
fix_clustering_interactive <- function(seurat_file, target_clusters = NULL) {
  
  # Load Seurat object
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_file)
  
  # Check current state
  cat("\n=== Current Clustering State ===\n")
  cat("Total cells:", ncol(seurat_obj), "\n")
  
  # Check default identities
  default_idents <- levels(Idents(seurat_obj))
  cat("Default identities:", length(default_idents), "levels\n")
  
  # Check seurat_clusters
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    current_clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
    cat("Current seurat_clusters:", length(current_clusters), "\n")
    cat("Clusters:", paste(current_clusters, collapse = ", "), "\n")
  } else {
    cat("No seurat_clusters found in metadata!\n")
    current_clusters <- NULL
  }
  
  # Ask for target if not provided
  if (is.null(target_clusters)) {
    cat("\nHow many clusters should this dataset have?\n")
    cat("(Check your DE results or enrichment files for the expected number)\n")
    target_clusters <- as.numeric(readline("Enter target number of clusters: "))
  }
  
  cat("\nTarget clusters:", target_clusters, "\n")
  
  # If already matches, we're done
  if (!is.null(current_clusters) && length(current_clusters) == target_clusters) {
    cat("\n✓ Clustering already matches target!\n")
    return(seurat_obj)
  }
  
  # Find best resolution
  cat("\n=== Finding Best Resolution ===\n")
  cat("Testing resolutions from 0.1 to 2.0...\n\n")
  
  best_res <- NULL
  best_diff <- Inf
  results <- list()
  
  for (res in seq(0.1, 2.0, by = 0.1)) {
    cat(sprintf("Resolution %.1f: ", res))
    
    seurat_obj_test <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    n_clusters <- length(unique(seurat_obj_test@meta.data$seurat_clusters))
    
    cat(n_clusters, "clusters")
    
    diff <- abs(n_clusters - target_clusters)
    if (diff < best_diff) {
      best_diff <- diff
      best_res <- res
    }
    
    if (n_clusters == target_clusters) {
      cat(" ✓ MATCH!")
    }
    cat("\n")
    
    results[[as.character(res)]] <- n_clusters
    
    # Stop if we found exact match
    if (n_clusters == target_clusters) {
      seurat_obj <- seurat_obj_test
      break
    }
  }
  
  # If no exact match, offer options
  if (best_diff > 0) {
    cat("\n✗ No exact match found.\n")
    cat("Best resolution:", best_res, "gives", results[[as.character(best_res)]], "clusters\n")
    
    cat("\nWould you like to:\n")
    cat("1. Use best resolution (", best_res, ")\n")
    cat("2. Try a custom resolution\n")
    cat("3. Cancel\n")
    
    choice <- readline("Enter choice (1-3): ")
    
    if (choice == "1") {
      cat("\nApplying resolution", best_res, "...\n")
      seurat_obj <- FindClusters(seurat_obj, resolution = best_res, verbose = FALSE)
    } else if (choice == "2") {
      custom_res <- as.numeric(readline("Enter custom resolution: "))
      cat("\nApplying resolution", custom_res, "...\n")
      seurat_obj <- FindClusters(seurat_obj, resolution = custom_res, verbose = FALSE)
    } else {
      cat("\nCancelled. No changes made.\n")
      return(NULL)
    }
  }
  
  # Verify final clustering
  final_clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
  cat("\n=== Final Clustering ===\n")
  cat("Number of clusters:", length(final_clusters), "\n")
  cat("Clusters:", paste(final_clusters, collapse = ", "), "\n")
  
  # Ask to save
  cat("\nSave updated Seurat object? (y/n): ")
  if (tolower(readline()) == "y") {
    # Create backup
    backup_file <- sub("\\.rds$", "_backup.rds", seurat_file)
    cat("Creating backup at:", backup_file, "\n")
    file.copy(seurat_file, backup_file)
    
    # Save updated object
    cat("Saving updated object...\n")
    saveRDS(seurat_obj, seurat_file)
    cat("✓ Saved!\n")
  }
  
  return(seurat_obj)
}

# Quick check function
quick_check <- function(seurat_file) {
  seurat_obj <- readRDS(seurat_file)
  
  cat("File:", seurat_file, "\n")
  cat("Cells:", ncol(seurat_obj), "\n")
  
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
    cat("Clusters:", length(clusters), "-", paste(clusters[1:min(5, length(clusters))], collapse = ", "))
    if (length(clusters) > 5) cat("...")
    cat("\n")
  } else {
    cat("No seurat_clusters found!\n")
  }
}

# Main interactive function
if (interactive()) {
  cat("=== Clustering Fix Tool ===\n\n")
  cat("Functions available:\n")
  cat("  quick_check('path/to/seurat.rds') - Check current clustering\n")
  cat("  fix_clustering_interactive('path/to/seurat.rds', target_clusters = 15) - Fix clustering\n")
  cat("\nExample:\n")
  cat("  fix_clustering_interactive('E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/final_iSCORE-PD.rds', 15)\n")
} else {
  # If run as script, prompt for file
  cat("Enter path to Seurat RDS file: ")
  file_path <- readline()
  
  if (file.exists(file_path)) {
    fix_clustering_interactive(file_path)
  } else {
    cat("File not found:", file_path, "\n")
  }
}