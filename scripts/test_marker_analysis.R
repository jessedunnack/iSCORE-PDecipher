#!/usr/bin/env Rscript

# Test script to run marker-based cell type analysis
# This demonstrates the workflow with actual data

# Source required scripts
source("scripts/analyze_cluster_markers_celltype.R")
source("scripts/web_search_cell_types.R")

# Function to run quick test
test_marker_analysis <- function() {
  cat("\n=== Testing Marker-Based Cell Type Analysis ===\n")
  
  # Check if marker files exist
  marker_files <- list.files("inst/extdata/umap_data", 
                            pattern = "_cluster_markers\\.rds$", 
                            full.names = TRUE)
  
  if (length(marker_files) == 0) {
    cat("No marker files found. Please ensure marker data is available.\n")
    return(NULL)
  }
  
  cat("Found", length(marker_files), "marker files:\n")
  for (f in marker_files) {
    cat("  -", basename(f), "\n")
  }
  
  # Test with first available dataset
  test_file <- marker_files[1]
  dataset_name <- gsub("_cluster_markers\\.rds$", "", basename(test_file))
  
  cat("\nTesting with dataset:", dataset_name, "\n")
  
  # Load markers
  markers <- readRDS(test_file)
  cat("Loaded", nrow(markers), "markers\n")
  
  # Show structure
  cat("\nMarker data structure:\n")
  str(markers, max.level = 1)
  
  # Get clusters
  clusters <- sort(unique(as.character(markers$cluster)))
  cat("\nClusters found:", paste(clusters, collapse = ", "), "\n")
  
  # Analyze first cluster as example
  test_cluster <- clusters[1]
  cat("\n=== Analyzing Cluster", test_cluster, "===\n")
  
  # Get top markers
  cluster_markers <- markers %>% 
    filter(cluster == test_cluster) %>%
    filter(avg_log2FC > 0.25, p_val_adj < 0.05) %>%
    mutate(score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)) %>%
    arrange(desc(score)) %>%
    slice_head(n = 25)
  
  cat("\nTop 10 markers:\n")
  print(cluster_markers %>% 
        select(gene, avg_log2FC, p_val_adj, score) %>%
        slice_head(n = 10))
  
  # Generate search queries
  cat("\n=== Generating Web Search Queries ===\n")
  search_queries <- perform_marker_search(cluster_markers$gene, test_cluster)
  
  # Show example query
  cat("\nExample search query for cell type identification:\n")
  cat(search_queries$cell_type$query, "\n")
  
  # Interpret based on known markers
  interpretation <- interpret_search_results(NULL, cluster_markers$gene)
  
  cat("\n=== Preliminary Cell Type Inference ===\n")
  if (length(interpretation$cell_type_evidence) > 0) {
    cat("Evidence found:\n")
    for (type in names(interpretation$cell_type_evidence)) {
      cat("  ", type, ": ", 
          paste(interpretation$cell_type_evidence[[type]], collapse = ", "), "\n", sep = "")
    }
  } else {
    cat("No strong evidence from known marker patterns.\n")
    cat("Web searches needed for novel marker combinations.\n")
  }
  
  return(list(
    dataset = dataset_name,
    cluster = test_cluster,
    top_markers = cluster_markers,
    search_queries = search_queries,
    interpretation = interpretation
  ))
}

# Run test
if (!interactive()) {
  results <- test_marker_analysis()
  
  if (!is.null(results)) {
    cat("\n=== Next Steps ===\n")
    cat("1. Run the full analysis: source('scripts/analyze_cluster_markers_celltype.R'); main()\n")
    cat("2. Use WebSearch tool to execute the generated queries\n")
    cat("3. Review results in: results/cell_type_annotations/\n")
    cat("4. Create fine-to-coarse cluster mapping with Seurat object\n")
  }
}