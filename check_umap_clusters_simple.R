#!/usr/bin/env Rscript

# Simple script to check cluster consistency without specialized packages
# Compares with expected cluster counts from DE data

# Expected clusters from DE data analysis
expected_clusters <- list(
  iSCORE_PD = 15,  # clusters 0-14
  iSCORE_PD_CRISPRi = 15,  # clusters 0-14
  Full_Dataset = 14  # clusters 0-13
)

# Path to UMAP data
umap_dir <- "inst/extdata/umap_data"

# Generic function to explore RDS file structure
explore_rds_file <- function(file_path, dataset_name) {
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
    
    # Get basic info
    cat("Object class:", class(data)[1], "\n")
    cat("Object type:", typeof(data), "\n")
    
    # If it's a list-like object
    if (is.list(data)) {
      cat("List length:", length(data), "\n")
      if (!is.null(names(data))) {
        cat("Named elements:", paste(names(data), collapse = ", "), "\n")
      }
      
      # For S4 objects, try to access slots
      if (isS4(data)) {
        cat("S4 object with slots:", paste(slotNames(data), collapse = ", "), "\n")
        
        # Try to find cluster information
        if ("colData" %in% slotNames(data)) {
          colData <- slot(data, "colData")
          if (is.data.frame(colData) || is.matrix(colData)) {
            cat("colData columns:", paste(colnames(colData), collapse = ", "), "\n")
            if ("cluster" %in% colnames(colData)) {
              clusters <- as.character(colData[, "cluster"])
              analyze_clusters(clusters, dataset_name)
            }
          }
        }
        
        # Check for cell names or barcodes
        if ("colnames" %in% slotNames(data)) {
          cell_names <- slot(data, "colnames")
          if (!is.null(cell_names)) {
            cat("Number of cells:", length(cell_names), "\n")
          }
        }
      } else if (is.data.frame(data)) {
        # Regular data frame
        cat("Data frame dimensions:", nrow(data), "x", ncol(data), "\n")
        cat("Columns:", paste(names(data), collapse = ", "), "\n")
        if ("cluster" %in% names(data)) {
          analyze_clusters(as.character(data$cluster), dataset_name)
        }
      } else if (is.list(data) && !is.null(names(data))) {
        # Named list - might be marker data
        if (all(grepl("cluster", names(data), ignore.case = TRUE))) {
          cat("Appears to be cluster marker data\n")
          analyze_cluster_names(names(data), dataset_name)
        }
      }
    }
    
    # Special handling for combined file
    if (basename(file_path) == "all_umap_data_combined.rds" && is.list(data)) {
      cat("\nThis appears to be a combined file. Checking each dataset:\n")
      for (name in names(data)) {
        cat("\n--- Subset:", name, "---\n")
        explore_rds_file(data[[name]], name)
      }
    }
    
  }, error = function(e) {
    cat("ERROR reading file:", conditionMessage(e), "\n")
  })
}

# Function to analyze clusters
analyze_clusters <- function(clusters, dataset_name) {
  unique_clusters <- sort(unique(clusters))
  n_clusters <- length(unique_clusters)
  
  cat("\nCluster Analysis:\n")
  cat("Total cells:", length(clusters), "\n")
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
  for (i in 1:min(5, length(cluster_table))) {
    cat(sprintf("  Cluster %s: %d cells\n", names(cluster_table)[i], cluster_table[i]))
  }
  if (length(cluster_table) > 5) {
    cat("  ... (showing first 5 clusters)\n")
  }
}

# Function to analyze cluster names
analyze_cluster_names <- function(cluster_names, dataset_name) {
  cat("Cluster names found:", length(cluster_names), "\n")
  
  # Extract numbers if in format "cluster_X"
  if (all(grepl("^cluster_[0-9]+$", cluster_names))) {
    cluster_nums <- as.numeric(gsub("cluster_", "", cluster_names))
    cluster_nums <- sort(cluster_nums)
    cat("Cluster numbers:", paste(cluster_nums, collapse = ", "), "\n")
    cat("Cluster range: ", min(cluster_nums), "-", max(cluster_nums), "\n", sep = "")
    n_clusters <- length(unique(cluster_nums))
  } else {
    n_clusters <- length(cluster_names)
    cat("Cluster names:", paste(head(cluster_names, 10), collapse = ", "))
    if (length(cluster_names) > 10) cat(", ...")
    cat("\n")
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
}

# Main analysis
cat("UMAP Data and Cluster Marker Analysis (Simple Version)\n")
cat("=====================================================\n")
cat("\nExpected cluster counts from DE analysis:\n")
cat("- iSCORE_PD: 15 clusters (0-14)\n")
cat("- iSCORE_PD_CRISPRi: 15 clusters (0-14)\n")
cat("- Full_Dataset: 14 clusters (0-13)\n")

# Check UMAP data files
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING UMAP DATA FILES\n")
cat(strrep("=", 70), "\n")

# iSCORE_PD
explore_rds_file(file.path(umap_dir, "iSCORE_PD_umap_data.rds"), "iSCORE_PD")

# iSCORE_PD_CRISPRi  
explore_rds_file(file.path(umap_dir, "iSCORE_PD_CRISPRi_umap_data.rds"), "iSCORE_PD_CRISPRi")

# Full_Dataset
explore_rds_file(file.path(umap_dir, "Full_Dataset_umap_data.rds"), "Full_Dataset")

# Check combined file
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING COMBINED UMAP FILE\n")
cat(strrep("=", 70), "\n")
explore_rds_file(file.path(umap_dir, "all_umap_data_combined.rds"), "all_combined")

# Check marker files
cat("\n\n", strrep("=", 70), "\n", sep = "")
cat("CHECKING CLUSTER MARKER FILES\n")
cat(strrep("=", 70), "\n")

# iSCORE_PD markers
explore_rds_file(file.path(umap_dir, "iSCORE_PD_cluster_markers.rds"), "iSCORE_PD")

# iSCORE_PD_CRISPRi markers
explore_rds_file(file.path(umap_dir, "iSCORE_PD_CRISPRi_cluster_markers.rds"), "iSCORE_PD_CRISPRi")

# Full_Dataset markers
explore_rds_file(file.path(umap_dir, "Full_Dataset_cluster_markers.rds"), "Full_Dataset")

# Check summary file if exists
summary_file <- file.path(umap_dir, "umap_data_summary.csv")
if (file.exists(summary_file)) {
  cat("\n\n", strrep("=", 70), "\n", sep = "")
  cat("UMAP DATA SUMMARY FILE\n")
  cat(strrep("=", 70), "\n")
  summary_data <- read.csv(summary_file, stringsAsFactors = FALSE)
  print(summary_data)
}

cat("\n\nAnalysis complete.\n")