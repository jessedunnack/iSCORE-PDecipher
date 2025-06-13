#!/usr/bin/env Rscript

# Script to check cluster consistency across Seurat objects and saved data
# Helps identify mismatches and find appropriate clustering resolution

library(Seurat)
library(SingleCellExperiment)
library(dplyr)

# Configuration - same as extract_umap_data.R
DATASETS <- list(
  iSCORE_PD = list(
    seurat_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/final_iSCORE-PD.rds",
    de_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/full_DE_results.rds",
    enrichment_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/all_enrichment_padj005_complete_with_direction.rds"
  ),
  iSCORE_PD_CRISPRi = list(
    seurat_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
    de_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds",
    enrichment_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
  ),
  Full_Dataset = list(
    seurat_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_dataset.rds",
    de_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_DE_results.rds",
    enrichment_file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/all_enrichment_padj005_complete_with_direction.rds"
  )
)

# Function to get cluster info from Seurat object
check_seurat_clusters <- function(seurat_file) {
  if (!file.exists(seurat_file)) {
    return(list(error = paste("File not found:", seurat_file)))
  }
  
  cat("Loading Seurat object from:", seurat_file, "\n")
  seurat_obj <- readRDS(seurat_file)
  
  # Check default identities
  default_idents <- levels(Idents(seurat_obj))
  n_default_idents <- length(default_idents)
  
  # Check seurat_clusters if available
  seurat_clusters <- NULL
  n_seurat_clusters <- NA
  if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
    seurat_clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
    n_seurat_clusters <- length(seurat_clusters)
  }
  
  # Get cell count
  n_cells <- ncol(seurat_obj)
  
  # Check available reductions
  reductions <- names(seurat_obj@reductions)
  
  return(list(
    n_cells = n_cells,
    default_idents = default_idents,
    n_default_idents = n_default_idents,
    seurat_clusters = seurat_clusters,
    n_seurat_clusters = n_seurat_clusters,
    reductions = reductions
  ))
}

# Function to check DE results clusters
check_de_clusters <- function(de_file) {
  if (!file.exists(de_file)) {
    return(list(error = paste("File not found:", de_file)))
  }
  
  cat("Loading DE results from:", de_file, "\n")
  de_results <- readRDS(de_file)
  
  clusters <- list()
  
  # Check MAST data
  if ("iSCORE_PD_MAST" %in% names(de_results)) {
    mast_clusters <- c()
    for (gene in names(de_results$iSCORE_PD_MAST)) {
      mast_clusters <- union(mast_clusters, names(de_results$iSCORE_PD_MAST[[gene]]))
    }
    clusters$MAST <- sort(mast_clusters)
  }
  
  # Check MixScale CRISPRi data
  if ("CRISPRi_Mixscale" %in% names(de_results)) {
    mixscale_clusters <- c()
    for (gene in names(de_results$CRISPRi_Mixscale)) {
      mixscale_clusters <- union(mixscale_clusters, names(de_results$CRISPRi_Mixscale[[gene]]))
    }
    clusters$MixScale_CRISPRi <- sort(mixscale_clusters)
  }
  
  # Check MixScale CRISPRa data
  if ("CRISPRa_Mixscale" %in% names(de_results)) {
    crispa_clusters <- c()
    for (gene in names(de_results$CRISPRa_Mixscale)) {
      crispa_clusters <- union(crispa_clusters, names(de_results$CRISPRa_Mixscale[[gene]]))
    }
    clusters$MixScale_CRISPRa <- sort(crispa_clusters)
  }
  
  return(clusters)
}

# Function to check enrichment results clusters
check_enrichment_clusters <- function(enrichment_file) {
  if (!file.exists(enrichment_file)) {
    return(list(error = paste("File not found:", enrichment_file)))
  }
  
  cat("Loading enrichment results from:", enrichment_file, "\n")
  enrichment_data <- readRDS(enrichment_file)
  
  clusters <- sort(unique(enrichment_data$cluster))
  methods <- unique(enrichment_data$method)
  
  return(list(
    all_clusters = clusters,
    n_clusters = length(clusters),
    methods = methods
  ))
}

# Function to try different resolutions
find_matching_resolution <- function(seurat_obj, target_n_clusters, resolution_range = c(0.1, 2.0)) {
  cat("\nSearching for resolution that gives", target_n_clusters, "clusters...\n")
  
  # Try a range of resolutions
  resolutions <- seq(resolution_range[1], resolution_range[2], by = 0.1)
  results <- data.frame(resolution = numeric(), n_clusters = numeric())
  
  for (res in resolutions) {
    cat("  Trying resolution", res, "... ")
    seurat_obj <- FindClusters(seurat_obj, resolution = res, verbose = FALSE)
    n_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
    cat(n_clusters, "clusters\n")
    
    results <- rbind(results, data.frame(resolution = res, n_clusters = n_clusters))
    
    if (n_clusters == target_n_clusters) {
      cat("\n✓ Found matching resolution:", res, "\n")
      return(list(resolution = res, n_clusters = n_clusters, success = TRUE))
    }
  }
  
  # If no exact match, suggest closest
  closest <- results[which.min(abs(results$n_clusters - target_n_clusters)), ]
  cat("\n✗ No exact match found. Closest:", closest$resolution, "gives", closest$n_clusters, "clusters\n")
  
  return(list(
    resolution = closest$resolution, 
    n_clusters = closest$n_clusters, 
    success = FALSE,
    all_results = results
  ))
}

# Main diagnostic function
main <- function() {
  cat("=== Cluster Consistency Check ===\n\n")
  
  for (dataset_name in names(DATASETS)) {
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n")
    cat("Dataset:", dataset_name, "\n")
    cat("━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━\n\n")
    
    dataset <- DATASETS[[dataset_name]]
    
    # Check Seurat object
    cat("1. Seurat Object Analysis:\n")
    seurat_info <- check_seurat_clusters(dataset$seurat_file)
    if (!is.null(seurat_info$error)) {
      cat("   ERROR:", seurat_info$error, "\n\n")
      next
    }
    
    cat("   - Total cells:", format(seurat_info$n_cells, big.mark = ","), "\n")
    cat("   - Default identities:", seurat_info$n_default_idents, "levels\n")
    cat("   - Seurat clusters:", seurat_info$n_seurat_clusters, "\n")
    cat("   - Available reductions:", paste(seurat_info$reductions, collapse = ", "), "\n\n")
    
    # Check DE results
    cat("2. DE Results Analysis:\n")
    de_clusters <- check_de_clusters(dataset$de_file)
    if (!is.null(de_clusters$error)) {
      cat("   ERROR:", de_clusters$error, "\n\n")
    } else {
      for (method in names(de_clusters)) {
        cat("   -", method, "clusters:", length(de_clusters[[method]]), 
            ":", paste(de_clusters[[method]][1:min(5, length(de_clusters[[method]]))], collapse = ", "),
            ifelse(length(de_clusters[[method]]) > 5, "...", ""), "\n")
      }
      cat("\n")
    }
    
    # Check enrichment results
    cat("3. Enrichment Results Analysis:\n")
    enrich_info <- check_enrichment_clusters(dataset$enrichment_file)
    if (!is.null(enrich_info$error)) {
      cat("   ERROR:", enrich_info$error, "\n\n")
    } else {
      cat("   - Total clusters:", enrich_info$n_clusters, "\n")
      cat("   - Clusters:", paste(enrich_info$all_clusters[1:min(5, length(enrich_info$all_clusters))], 
                                 collapse = ", "),
          ifelse(length(enrich_info$all_clusters) > 5, "...", ""), "\n")
      cat("   - Methods:", paste(enrich_info$methods, collapse = ", "), "\n\n")
    }
    
    # Check consistency
    cat("4. Consistency Check:\n")
    
    # Determine expected cluster count from DE/enrichment data
    expected_clusters <- NULL
    if (!is.null(de_clusters$MAST)) {
      expected_clusters <- length(de_clusters$MAST)
      cat("   - Expected clusters (from DE data):", expected_clusters, "\n")
    } else if (!is.null(enrich_info$n_clusters)) {
      expected_clusters <- enrich_info$n_clusters
      cat("   - Expected clusters (from enrichment):", expected_clusters, "\n")
    }
    
    if (!is.null(expected_clusters) && !is.na(seurat_info$n_seurat_clusters)) {
      if (seurat_info$n_seurat_clusters == expected_clusters) {
        cat("   ✓ Cluster counts match!\n")
      } else {
        cat("   ✗ MISMATCH: Seurat has", seurat_info$n_seurat_clusters, 
            "clusters but data expects", expected_clusters, "\n\n")
        
        # Ask user if they want to find matching resolution
        cat("   Would you like to find a resolution that gives", expected_clusters, "clusters? (y/n): ")
        response <- readline()
        
        if (tolower(response) == "y") {
          seurat_obj <- readRDS(dataset$seurat_file)
          result <- find_matching_resolution(seurat_obj, expected_clusters)
          
          if (result$success) {
            cat("\n   To fix this mismatch, re-run clustering with resolution =", result$resolution, "\n")
          } else {
            cat("\n   Suggested resolution:", result$resolution, 
                "(gives", result$n_clusters, "clusters)\n")
            cat("   You may need to try a custom resolution or check your preprocessing.\n")
          }
        }
      }
    }
    
    cat("\n")
  }
  
  cat("\n=== Diagnostic Complete ===\n")
}

# Run if called directly
if (!interactive()) {
  main()
} else {
  cat("Cluster consistency check script loaded.\n")
  cat("Run main() to check all datasets, or use individual functions:\n")
  cat("  - check_seurat_clusters(file): Check clusters in Seurat object\n")
  cat("  - check_de_clusters(file): Check clusters in DE results\n")
  cat("  - check_enrichment_clusters(file): Check clusters in enrichment results\n")
}