#!/usr/bin/env Rscript

# Test script to verify UMAP PC files contain different coordinates

library(SingleCellExperiment)

# Function to load and compare UMAP files
compare_umap_files <- function(file1, file2, dataset_name) {
  cat("\n=== Comparing UMAP files for", dataset_name, "===\n")
  cat("File 1:", file1, "\n")
  cat("File 2:", file2, "\n\n")
  
  # Check if files exist
  if (!file.exists(file1)) {
    cat("ERROR: File 1 does not exist\n")
    return(FALSE)
  }
  if (!file.exists(file2)) {
    cat("ERROR: File 2 does not exist\n")
    return(FALSE)
  }
  
  # Load files
  sce1 <- readRDS(file1)
  sce2 <- readRDS(file2)
  
  # Extract UMAP coordinates
  umap1 <- reducedDim(sce1, "UMAP")
  umap2 <- reducedDim(sce2, "UMAP")
  
  # Check dimensions
  cat("UMAP 1 dimensions:", dim(umap1), "\n")
  cat("UMAP 2 dimensions:", dim(umap2), "\n")
  
  # Check if dimensions match
  if (!identical(dim(umap1), dim(umap2))) {
    cat("WARNING: Dimensions don't match!\n")
    return(FALSE)
  }
  
  # Compare first 10 cells
  cat("\nFirst 5 cells from 30PC UMAP:\n")
  print(head(umap1, 5))
  
  cat("\nFirst 5 cells from 100PC UMAP:\n")
  print(head(umap2, 5))
  
  # Check if coordinates are identical
  coords_identical <- identical(umap1, umap2)
  cat("\nAre coordinates identical?", coords_identical, "\n")
  
  if (!coords_identical) {
    # Calculate correlation
    cor_umap1 <- cor(umap1[,1], umap2[,1])
    cor_umap2 <- cor(umap1[,2], umap2[,2])
    cat("Correlation UMAP_1:", round(cor_umap1, 4), "\n")
    cat("Correlation UMAP_2:", round(cor_umap2, 4), "\n")
    
    # Calculate mean absolute difference
    mean_diff_1 <- mean(abs(umap1[,1] - umap2[,1]))
    mean_diff_2 <- mean(abs(umap1[,2] - umap2[,2]))
    cat("Mean absolute difference UMAP_1:", round(mean_diff_1, 4), "\n")
    cat("Mean absolute difference UMAP_2:", round(mean_diff_2, 4), "\n")
  }
  
  # Check metadata
  if (!is.null(metadata(sce1)$pc_count) && !is.null(metadata(sce2)$pc_count)) {
    cat("\nMetadata PC counts:\n")
    cat("File 1 PC count:", metadata(sce1)$pc_count, "\n")
    cat("File 2 PC count:", metadata(sce2)$pc_count, "\n")
  }
  
  return(!coords_identical)
}

# Main execution
base_dir <- "inst/extdata/umap_data/"

# Test all dataset pairs
datasets <- c("iSCORE_PD", "iSCORE_PD_CRISPRi", "Full_Dataset")

for (dataset in datasets) {
  file_30pc <- file.path(base_dir, paste0(dataset, "_umap_data_30pc.rds"))
  file_100pc <- file.path(base_dir, paste0(dataset, "_umap_data_100pc.rds"))
  
  if (file.exists(file_30pc) && file.exists(file_100pc)) {
    different <- compare_umap_files(file_30pc, file_100pc, dataset)
    if (!different) {
      cat("\n*** WARNING: 30PC and 100PC files appear to have IDENTICAL coordinates! ***\n")
    }
  } else {
    cat("\nSkipping", dataset, "- files not found\n")
  }
}