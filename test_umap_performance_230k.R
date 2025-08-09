#!/usr/bin/env Rscript

# Test UMAP Performance with 230K+ Cell Dataset
# Tests caching implementation with real production data

cat("========================================\n")
cat("UMAP Performance Test for 230K+ Cells\n")
cat("========================================\n\n")

# Load required libraries
library(Seurat)
library(ggplot2)
library(SingleCellExperiment)

# Install tictoc if needed
if (!requireNamespace("tictoc", quietly = TRUE)) {
  install.packages("tictoc", repos = "https://cran.r-project.org")
}
library(tictoc)

# Source cache manager
source("inst/shiny/R/cache_manager.R")
source("inst/shiny/modules/umap_cache_integration.R")

# Initialize cache
cat("Initializing UMAP cache...\n")
GLOBAL_UMAP_CACHE <- CacheManager$new(
  max_size = 100,      # Larger cache for testing
  ttl_minutes = 240,   # 4 hour TTL
  verbose = TRUE
)

# File paths
ISCORE_PD_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/iSCORE-PD_final.rds"
ISCORE_PD_PLUS_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_final.rds"
METADATA_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_final_metadata.rds"

# Function to extract and test UMAP data
test_umap_extraction <- function(seurat_obj, dataset_name) {
  cat("\n----------------------------------------\n")
  cat("Testing:", dataset_name, "\n")
  cat("----------------------------------------\n")
  
  # Get basic info
  n_cells <- ncol(seurat_obj)
  n_genes <- nrow(seurat_obj)
  
  cat(sprintf("Dataset: %d cells, %d genes\n", n_cells, n_genes))
  
  # Extract UMAP coordinates
  cat("\nExtracting UMAP coordinates...\n")
  tic()
  
  if ("umap" %in% names(seurat_obj@reductions)) {
    umap_coords <- seurat_obj@reductions$umap@cell.embeddings
    clusters <- seurat_obj$seurat_clusters
  } else if ("UMAP" %in% names(seurat_obj@reductions)) {
    umap_coords <- seurat_obj@reductions$UMAP@cell.embeddings
    clusters <- seurat_obj$seurat_clusters
  } else {
    cat("ERROR: No UMAP reduction found!\n")
    return(NULL)
  }
  
  extraction_time <- toc(quiet = TRUE)
  cat(sprintf("UMAP extraction time: %.2f seconds\n", 
              extraction_time$toc - extraction_time$tic))
  
  # Create data frame for plotting
  umap_data <- data.frame(
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    cluster = paste0("cluster_", clusters)
  )
  
  # Test plot generation WITHOUT cache
  cat("\nGenerating UMAP plot (no cache)...\n")
  tic()
  
  p1 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.6) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = paste(dataset_name, "- First Render"))
  
  first_render_time <- toc(quiet = TRUE)
  cat(sprintf("First render time: %.2f seconds\n", 
              first_render_time$toc - first_render_time$tic))
  
  # Save plot to cache
  cache_key <- paste(dataset_name, "all_clusters", sep = "_")
  GLOBAL_UMAP_CACHE$set(cache_key, p1)
  
  # Test plot retrieval FROM cache
  cat("\nRetrieving UMAP plot (from cache)...\n")
  tic()
  
  cached_plot <- GLOBAL_UMAP_CACHE$get(cache_key)
  
  cache_retrieve_time <- toc(quiet = TRUE)
  cat(sprintf("Cache retrieve time: %.3f seconds\n", 
              cache_retrieve_time$toc - cache_retrieve_time$tic))
  
  # Calculate speedup
  speedup <- (first_render_time$toc - first_render_time$tic) / 
             (cache_retrieve_time$toc - cache_retrieve_time$tic)
  cat(sprintf("SPEEDUP: %.1fx faster with cache!\n", speedup))
  
  # Test with cluster subset (highlighting)
  cat("\nTesting cluster highlighting...\n")
  
  # Get unique clusters
  unique_clusters <- unique(umap_data$cluster)
  n_clusters <- length(unique_clusters)
  cat(sprintf("Found %d clusters\n", n_clusters))
  
  # Test highlighting cluster 0
  if ("cluster_0" %in% unique_clusters) {
    tic()
    
    umap_data$highlight <- ifelse(umap_data$cluster == "cluster_0", 
                                  "Cluster 0", "Other")
    
    p2 <- ggplot(umap_data, aes(x = UMAP1, y = UMAP2)) +
      geom_point(aes(color = highlight, alpha = highlight), size = 0.3) +
      scale_color_manual(values = c("Cluster 0" = "red", "Other" = "gray80")) +
      scale_alpha_manual(values = c("Cluster 0" = 0.8, "Other" = 0.2)) +
      theme_minimal() +
      labs(title = paste(dataset_name, "- Cluster 0 Highlighted"))
    
    highlight_time <- toc(quiet = TRUE)
    cat(sprintf("Highlight plot time: %.2f seconds\n", 
                highlight_time$toc - highlight_time$tic))
  }
  
  # Return results
  return(list(
    n_cells = n_cells,
    n_clusters = n_clusters,
    first_render = first_render_time$toc - first_render_time$tic,
    cache_retrieve = cache_retrieve_time$toc - cache_retrieve_time$tic,
    speedup = speedup,
    umap_data = umap_data
  ))
}

# Main testing
cat("\n========================================\n")
cat("Starting Performance Tests\n")
cat("========================================\n")

# Test smaller dataset first (iSCORE-PD only)
if (file.exists(ISCORE_PD_PATH)) {
  cat("\n1. Testing iSCORE-PD dataset (~210K cells)\n")
  cat("Loading dataset...\n")
  
  tic()
  seurat_pd <- readRDS(ISCORE_PD_PATH)
  load_time <- toc(quiet = TRUE)
  
  cat(sprintf("Load time: %.2f seconds\n", load_time$toc - load_time$tic))
  
  results_pd <- test_umap_extraction(seurat_pd, "iSCORE-PD")
  
  # Clean up memory
  rm(seurat_pd)
  gc()
}

# Test larger dataset (iSCORE-PD plus CRISPRi)
if (file.exists(ISCORE_PD_PLUS_PATH)) {
  cat("\n2. Testing iSCORE-PD_plus_CRISPRi dataset (230K+ cells)\n")
  cat("Loading dataset...\n")
  
  tic()
  seurat_plus <- readRDS(ISCORE_PD_PLUS_PATH)
  load_time <- toc(quiet = TRUE)
  
  cat(sprintf("Load time: %.2f seconds\n", load_time$toc - load_time$tic))
  
  results_plus <- test_umap_extraction(seurat_plus, "iSCORE-PD_plus_CRISPRi")
  
  # Test with metadata if available
  if (file.exists(METADATA_PATH)) {
    cat("\nLoading metadata...\n")
    metadata <- readRDS(METADATA_PATH)
    cat(sprintf("Metadata dimensions: %d rows, %d columns\n", 
                nrow(metadata), ncol(metadata)))
    
    # Check key metadata columns
    key_cols <- c("orig.ident", "nCount_RNA", "nFeature_RNA", 
                  "seurat_clusters", "cell_type", "condition")
    available_cols <- intersect(key_cols, colnames(metadata))
    cat("Available metadata columns:", paste(available_cols, collapse = ", "), "\n")
  }
  
  # Clean up memory
  rm(seurat_plus)
  gc()
}

# Final cache statistics
cat("\n========================================\n")
cat("Final Cache Statistics\n")
cat("========================================\n")

cache_stats <- GLOBAL_UMAP_CACHE$stats()
cat(sprintf("Cache usage: %d/%d plots\n", cache_stats$size, cache_stats$max_size))
cat(sprintf("Cache TTL: %d minutes\n", cache_stats$ttl_minutes))

if (cache_stats$size > 0) {
  cat("Cached keys:\n")
  for (key in cache_stats$keys) {
    cat("  -", key, "\n")
  }
}

# Summary
cat("\n========================================\n")
cat("Performance Summary\n")
cat("========================================\n")

if (exists("results_pd")) {
  cat(sprintf("\niSCORE-PD (~210K cells):\n"))
  cat(sprintf("  First render: %.2f seconds\n", results_pd$first_render))
  cat(sprintf("  Cache retrieve: %.3f seconds\n", results_pd$cache_retrieve))
  cat(sprintf("  Speedup: %.1fx\n", results_pd$speedup))
}

if (exists("results_plus")) {
  cat(sprintf("\niSCORE-PD_plus_CRISPRi (230K+ cells):\n"))
  cat(sprintf("  First render: %.2f seconds\n", results_plus$first_render))
  cat(sprintf("  Cache retrieve: %.3f seconds\n", results_plus$cache_retrieve))
  cat(sprintf("  Speedup: %.1fx\n", results_plus$speedup))
}

cat("\n✅ Testing complete!\n")

# Save results
results <- list(
  pd = if(exists("results_pd")) results_pd else NULL,
  plus = if(exists("results_plus")) results_plus else NULL,
  cache_stats = cache_stats,
  timestamp = Sys.time()
)

saveRDS(results, "test_results_230k_umap_performance.rds")
cat("\nResults saved to: test_results_230k_umap_performance.rds\n")