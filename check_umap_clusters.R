#!/usr/bin/env Rscript

# Script to check cluster consistency in UMAP data and marker files
# Compares with expected cluster counts from DE data

# Try to load packages, but continue if not available
suppressPackageStartupMessages({
  sce_available <- requireNamespace("SingleCellExperiment", quietly = TRUE)
  se_available <- requireNamespace("SummarizedExperiment", quietly = TRUE)
  
  if (sce_available) library(SingleCellExperiment)
  if (se_available) library(SummarizedExperiment)
})

# Expected clusters from DE data analysis
expected_clusters <- list(
  iSCORE_PD = 15,  # clusters 0-14
  iSCORE_PD_CRISPRi = 15,  # clusters 0-14
  Full_Dataset = 14  # clusters 0-13
)

# Path to UMAP data
umap_dir <- "inst/extdata/umap_data"

# Function to check clusters in SCE object
check_sce_clusters <- function(file_path, dataset_name) {
  cat("\n", strrep("=", 60), "\n")
  cat("Checking:", dataset_name, "\n")
  cat("File:", basename(file_path), "\n")
  cat(strrep("-", 60), "\n")
  
  if (!file.exists(file_path)) {
    cat("ERROR: File not found!\n")
    return(NULL)
  }
  
  tryCatch({
    # Load data
    data <- readRDS(file_path)
    
    # Check if it's a SingleCellExperiment
    if (sce_available && is(data, "SingleCellExperiment")) {
      # Get cluster information
      if ("cluster" %in% names(colData(data))) {
        clusters <- as.character(colData(data)$cluster)
        unique_clusters <- sort(unique(clusters))
        n_clusters <- length(unique_clusters)
        
        cat("Data type: SingleCellExperiment\n")
        cat("Total cells:", ncol(data), "\n")
        cat("Number of clusters:", n_clusters, "\n")
        cat("Cluster IDs:", paste(unique_clusters, collapse = ", "), "\n")
        
        # Check if numeric
        if (all(grepl("^[0-9]+$", unique_clusters))) {
          numeric_clusters <- as.numeric(unique_clusters)
          cat("Cluster range: ", min(numeric_clusters), "-", max(numeric_clusters), "\n", sep = "")
        }
        
        # Compare with expected
        expected <- expected_clusters[[dataset_name]]
        if (!is.null(expected)) {
          if (n_clusters == expected) {
            cat("✓ MATCHES expected cluster count (", expected, ")\n", sep = "")
          } else {
            cat("✗ MISMATCH: Expected", expected, "clusters, found", n_clusters, "\n")
          }
        }
        
        # Show cluster distribution
        cluster_table <- table(clusters)
        cat("\nCluster distribution:\n")
        for (i in seq_along(cluster_table)) {
          cat(sprintf("  Cluster %s: %d cells\n", names(cluster_table)[i], cluster_table[i]))
        }
        
      } else {
        cat("WARNING: No 'cluster' column found in colData\n")
        cat("Available columns:", paste(names(colData(data)), collapse = ", "), "\n")
      }
    } else {
      cat("Data type:", class(data)[1], "\n")
      cat("WARNING: Not a SingleCellExperiment object\n")
    }
    
  }, error = function(e) {
    cat("ERROR reading file:", conditionMessage(e), "\n")
  })
}

# Function to check marker files
check_marker_file <- function(file_path, dataset_name) {
  cat("\n", strrep("=", 60), "\n")
  cat("Checking markers:", dataset_name, "\n")
  cat("File:", basename(file_path), "\n")
  cat(strrep("-", 60), "\n")
  
  if (!file.exists(file_path)) {
    cat("ERROR: File not found!\n")
    return(NULL)
  }
  
  tryCatch({
    # Load data
    data <- readRDS(file_path)
    
    # Check structure
    if (is.list(data)) {
      cat("Data type: List\n")
      cat("Number of elements:", length(data), "\n")
      
      # If it's a named list, check names (likely cluster names)
      if (!is.null(names(data))) {
        cluster_names <- names(data)
        cat("Cluster names:", paste(cluster_names, collapse = ", "), "\n")
        
        # Check if numeric
        if (all(grepl("^cluster_[0-9]+$", cluster_names))) {
          # Extract numbers
          cluster_nums <- as.numeric(gsub("cluster_", "", cluster_names))
          cluster_nums <- sort(cluster_nums)
          cat("Cluster numbers:", paste(cluster_nums, collapse = ", "), "\n")
          cat("Cluster range: ", min(cluster_nums), "-", max(cluster_nums), "\n", sep = "")
          n_clusters <- length(unique(cluster_nums))
        } else {
          n_clusters <- length(cluster_names)
        }
        
        # Compare with expected
        expected <- expected_clusters[[dataset_name]]
        if (!is.null(expected)) {
          if (n_clusters == expected) {
            cat("✓ MATCHES expected cluster count (", expected, ")\n", sep = "")
          } else {
            cat("✗ MISMATCH: Expected", expected, "clusters, found", n_clusters, "\n")
          }
        }
        
        # Check content of first element
        if (length(data) > 0) {
          first_elem <- data[[1]]
          cat("\nFirst element structure:\n")
          cat("  Class:", class(first_elem)[1], "\n")
          if (is.data.frame(first_elem)) {
            cat("  Dimensions:", nrow(first_elem), "x", ncol(first_elem), "\n")
            cat("  Columns:", paste(names(first_elem), collapse = ", "), "\n")
          }
        }
      }
    } else if (is.data.frame(data)) {
      cat("Data type: data.frame\n")
      cat("Dimensions:", nrow(data), "x", ncol(data), "\n")
      cat("Columns:", paste(names(data), collapse = ", "), "\n")
      
      # Check for cluster column
      if ("cluster" %in% names(data)) {
        unique_clusters <- sort(unique(as.character(data$cluster)))
        cat("Unique clusters:", paste(unique_clusters, collapse = ", "), "\n")
      }
    }
    
  }, error = function(e) {
    cat("ERROR reading file:", conditionMessage(e), "\n")
  })
}

# Main analysis
cat("UMAP Data and Cluster Marker Analysis\n")
cat("=====================================\n")
cat("\nExpected cluster counts from DE analysis:\n")
cat("- iSCORE_PD: 15 clusters (0-14)\n")
cat("- iSCORE_PD_CRISPRi: 15 clusters (0-14)\n")
cat("- Full_Dataset: 14 clusters (0-13)\n")

# Check UMAP data files
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING UMAP DATA FILES\n")
cat(strrep("=", 70), "\n")

# iSCORE_PD
check_sce_clusters(file.path(umap_dir, "iSCORE_PD_umap_data.rds"), "iSCORE_PD")

# iSCORE_PD_CRISPRi
check_sce_clusters(file.path(umap_dir, "iSCORE_PD_CRISPRi_umap_data.rds"), "iSCORE_PD_CRISPRi")

# Full_Dataset
check_sce_clusters(file.path(umap_dir, "Full_Dataset_umap_data.rds"), "Full_Dataset")

# Check combined file
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING COMBINED UMAP FILE\n")
cat(strrep("=", 70), "\n")
combined_file <- file.path(umap_dir, "all_umap_data_combined.rds")
if (file.exists(combined_file)) {
  data <- readRDS(combined_file)
  if (is.list(data)) {
    cat("Combined file is a list with elements:\n")
    cat(paste("-", names(data)), sep = "\n")
    for (name in names(data)) {
      check_sce_clusters(data[[name]], name)
    }
  }
}

# Check marker files
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING CLUSTER MARKER FILES\n")
cat(strrep("=", 70), "\n")

# iSCORE_PD markers
check_marker_file(file.path(umap_dir, "iSCORE_PD_cluster_markers.rds"), "iSCORE_PD")

# iSCORE_PD_CRISPRi markers
check_marker_file(file.path(umap_dir, "iSCORE_PD_CRISPRi_cluster_markers.rds"), "iSCORE_PD_CRISPRi")

# Full_Dataset markers
check_marker_file(file.path(umap_dir, "Full_Dataset_cluster_markers.rds"), "Full_Dataset")

# Check summary file if exists
summary_file <- file.path(umap_dir, "umap_data_summary.csv")
if (file.exists(summary_file)) {
  cat("\n\n", strrep("=", 70), "\n", sep = "")
  cat("UMAP DATA SUMMARY FILE\n")
  cat(strrep("=", 70), "\n")
  summary_data <- read.csv(summary_file)
  print(summary_data)
}

cat("\n\nAnalysis complete.\n")