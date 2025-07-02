#!/usr/bin/env Rscript

# Detailed analysis of cluster marker files

# Path to marker files
umap_dir <- "inst/extdata/umap_data"

# Function to analyze marker file
analyze_markers <- function(file_path, dataset_name, expected_clusters) {
  cat("\n", strrep("=", 70), "\n")
  cat("Analyzing markers for:", dataset_name, "\n")
  cat("File:", basename(file_path), "\n")
  cat("Expected clusters:", expected_clusters, "\n")
  cat(strrep("-", 70), "\n")
  
  if (!file.exists(file_path)) {
    cat("ERROR: File not found!\n")
    return(NULL)
  }
  
  # Load the data
  markers <- readRDS(file_path)
  
  # Check structure
  cat("Object class:", class(markers)[1], "\n")
  
  if (is.data.frame(markers)) {
    cat("Data frame with", nrow(markers), "rows and", ncol(markers), "columns\n")
    cat("Columns:", paste(names(markers), collapse=", "), "\n\n")
    
    # Analyze cluster column
    if ("cluster" %in% names(markers)) {
      clusters <- as.character(markers$cluster)
      unique_clusters <- sort(unique(clusters))
      n_unique <- length(unique_clusters)
      
      cat("Unique clusters found:", n_unique, "\n")
      cat("Cluster IDs:", paste(unique_clusters, collapse=", "), "\n")
      
      # Check if numeric
      if (all(grepl("^[0-9]+$", unique_clusters))) {
        numeric_clusters <- as.numeric(unique_clusters)
        cat("Numeric range:", min(numeric_clusters), "-", max(numeric_clusters), "\n")
      }
      
      # Compare with expected
      if (n_unique == expected_clusters) {
        cat("✓ MATCHES expected cluster count\n")
      } else {
        cat("✗ MISMATCH: Expected", expected_clusters, "clusters, found", n_unique, "\n")
      }
      
      # Show markers per cluster
      cat("\nMarkers per cluster:\n")
      cluster_counts <- table(clusters)
      for (cl in names(cluster_counts)) {
        cat("  Cluster", cl, ":", cluster_counts[cl], "markers\n")
      }
      
      # Show top markers for first few clusters
      cat("\nTop 3 markers per cluster (first 3 clusters):\n")
      for (cl in head(unique_clusters, 3)) {
        cluster_markers <- markers[markers$cluster == cl, ]
        # Sort by p_val_adj or avg_log2FC
        if ("p_val_adj" %in% names(cluster_markers)) {
          cluster_markers <- cluster_markers[order(cluster_markers$p_val_adj), ]
        } else if ("avg_log2FC" %in% names(cluster_markers)) {
          cluster_markers <- cluster_markers[order(-abs(cluster_markers$avg_log2FC)), ]
        }
        cat("\n  Cluster", cl, ":\n")
        top_genes <- head(cluster_markers$gene, 3)
        for (i in seq_along(top_genes)) {
          cat("    ", i, ". ", top_genes[i], sep="")
          if ("avg_log2FC" %in% names(cluster_markers)) {
            cat(" (log2FC=", round(cluster_markers$avg_log2FC[i], 2), ")", sep="")
          }
          if ("p_val_adj" %in% names(cluster_markers)) {
            cat(" (p.adj=", format(cluster_markers$p_val_adj[i], scientific=TRUE, digits=2), ")", sep="")
          }
          cat("\n")
        }
      }
    } else {
      cat("ERROR: No 'cluster' column found in data frame\n")
    }
  } else {
    cat("ERROR: Unexpected data structure (not a data frame)\n")
  }
}

# Main analysis
cat("Cluster Marker File Analysis\n")
cat("============================\n")
cat("\nExpected cluster counts from DE analysis:\n")
cat("- iSCORE_PD: 15 clusters (0-14)\n")
cat("- iSCORE_PD_CRISPRi: 15 clusters (0-14)\n")
cat("- Full_Dataset: 14 clusters (0-13)\n")

# Analyze each marker file
analyze_markers(file.path(umap_dir, "iSCORE_PD_cluster_markers.rds"), 
                "iSCORE_PD", 15)

analyze_markers(file.path(umap_dir, "iSCORE_PD_CRISPRi_cluster_markers.rds"), 
                "iSCORE_PD_CRISPRi", 15)

analyze_markers(file.path(umap_dir, "Full_Dataset_cluster_markers.rds"), 
                "Full_Dataset", 14)

# Check if there are other marker files
cat("\n\n", strrep("=", 70), "\n")
cat("Checking for other RDS files in directory:\n")
cat(strrep("-", 70), "\n")
all_files <- list.files(umap_dir, pattern = "\\.rds$", full.names = FALSE)
for (f in all_files) {
  cat("  -", f, "\n")
}

cat("\n\nAnalysis complete.\n")