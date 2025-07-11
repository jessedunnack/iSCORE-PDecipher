#!/usr/bin/env Rscript

# Script to check if UMAP files with different PC counts actually contain different coordinates

cat("=== UMAP Files Comparison Analysis ===\n\n")

# Find all UMAP files
umap_dir <- "inst/extdata/umap_data/"
all_files <- list.files(umap_dir, pattern = "umap_data.*\\.rds$", full.names = TRUE)

cat("Found UMAP files:\n")
for (f in all_files) {
  cat("  ", basename(f), "\n")
}
cat("\n")

# Function to safely extract UMAP coordinates
extract_umap_coords <- function(file_path) {
  tryCatch({
    obj <- readRDS(file_path)
    
    # Check if it's a SingleCellExperiment object
    if ("SingleCellExperiment" %in% class(obj)) {
      # Try to load SingleCellExperiment package
      if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        library(SingleCellExperiment)
        umap_coords <- reducedDim(obj, "UMAP")
        return(list(
          success = TRUE,
          coords = umap_coords,
          n_cells = ncol(obj),
          metadata = metadata(obj)
        ))
      } else {
        cat("  SingleCellExperiment package not available\n")
        return(list(success = FALSE, error = "SingleCellExperiment not available"))
      }
    } else {
      # Try to extract coordinates from other object types
      if (is.list(obj) && "UMAP" %in% names(obj)) {
        return(list(
          success = TRUE,
          coords = obj$UMAP,
          n_cells = nrow(obj$UMAP),
          metadata = obj$metadata %||% list()
        ))
      } else {
        return(list(success = FALSE, error = paste("Unknown object type:", class(obj))))
      }
    }
  }, error = function(e) {
    return(list(success = FALSE, error = e$message))
  })
}

# Compare specific dataset files
datasets_to_check <- c("iSCORE_PD", "iSCORE_PD_CRISPRi", "Full_Dataset")

for (dataset in datasets_to_check) {
  cat("=== Checking", dataset, "===\n")
  
  # Find files for this dataset
  dataset_files <- all_files[grepl(dataset, all_files)]
  pc_30_file <- dataset_files[grepl("30pc", dataset_files)]
  pc_100_file <- dataset_files[grepl("100pc", dataset_files)]
  pc_50_file <- dataset_files[grepl("50pc", dataset_files)]
  
  if (length(pc_30_file) == 0 || length(pc_100_file) == 0) {
    cat("  Missing 30PC or 100PC file for", dataset, "\n")
    next
  }
  
  cat("  30PC file:", basename(pc_30_file), "\n")
  cat("  100PC file:", basename(pc_100_file), "\n")
  if (length(pc_50_file) > 0) {
    cat("  50PC file:", basename(pc_50_file), "\n")
  }
  
  # Extract coordinates
  umap_30 <- extract_umap_coords(pc_30_file)
  umap_100 <- extract_umap_coords(pc_100_file)
  
  if (!umap_30$success) {
    cat("  ERROR loading 30PC:", umap_30$error, "\n")
    next
  }
  
  if (!umap_100$success) {
    cat("  ERROR loading 100PC:", umap_100$error, "\n") 
    next
  }
  
  cat("  30PC: ", nrow(umap_30$coords), " cells, ", ncol(umap_30$coords), " dimensions\n")
  cat("  100PC: ", nrow(umap_100$coords), " cells, ", ncol(umap_100$coords), " dimensions\n")
  
  # Check if coordinates are identical
  if (identical(dim(umap_30$coords), dim(umap_100$coords))) {
    coords_identical <- identical(umap_30$coords, umap_100$coords)
    
    if (coords_identical) {
      cat("  *** COORDINATES ARE IDENTICAL! ***\n")
    } else {
      cat("  Coordinates are DIFFERENT\n")
      
      # Calculate some differences
      diff_umap1 <- mean(abs(umap_30$coords[,1] - umap_100$coords[,1]))
      diff_umap2 <- mean(abs(umap_30$coords[,2] - umap_100$coords[,2]))
      
      cat("  Mean absolute difference UMAP_1:", round(diff_umap1, 4), "\n")
      cat("  Mean absolute difference UMAP_2:", round(diff_umap2, 4), "\n")
      
      # Show first few coordinates
      cat("  First 3 cells (30PC):\n")
      print(head(umap_30$coords, 3))
      cat("  First 3 cells (100PC):\n")
      print(head(umap_100$coords, 3))
    }
  } else {
    cat("  ERROR: Dimension mismatch!\n")
    cat("  30PC dims:", dim(umap_30$coords), "\n")
    cat("  100PC dims:", dim(umap_100$coords), "\n")
  }
  
  # Check metadata for PC count information
  if (!is.null(umap_30$metadata$pc_count)) {
    cat("  30PC metadata PC count:", umap_30$metadata$pc_count, "\n")
  }
  if (!is.null(umap_100$metadata$pc_count)) {
    cat("  100PC metadata PC count:", umap_100$metadata$pc_count, "\n")
  }
  
  cat("\n")
}

cat("=== Analysis Complete ===\n")