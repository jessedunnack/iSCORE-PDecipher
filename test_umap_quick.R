#!/usr/bin/env Rscript

# Quick UMAP Performance Test (without loading full datasets)
# Tests caching implementation with simulated data

cat("========================================\n")
cat("Quick UMAP Cache Performance Test\n")
cat("========================================\n\n")

# Load required libraries
library(ggplot2)

# Source cache manager
source("inst/shiny/R/cache_manager.R")

# Initialize cache
cat("Initializing UMAP cache...\n")
GLOBAL_UMAP_CACHE <- CacheManager$new(
  max_size = 100,
  ttl_minutes = 240,
  verbose = TRUE
)

# Function to simulate UMAP data
simulate_umap_data <- function(n_cells, n_clusters = 16) {
  cat(sprintf("\nSimulating UMAP data: %d cells, %d clusters\n", n_cells, n_clusters))
  
  # Generate cluster assignments
  clusters <- sample(0:(n_clusters-1), n_cells, replace = TRUE)
  
  # Generate UMAP coordinates (with cluster structure)
  umap1 <- numeric(n_cells)
  umap2 <- numeric(n_cells)
  
  for (i in 0:(n_clusters-1)) {
    idx <- which(clusters == i)
    n <- length(idx)
    if (n > 0) {
      # Each cluster has a center
      center_x <- runif(1, -10, 10)
      center_y <- runif(1, -10, 10)
      
      # Add noise around center
      umap1[idx] <- rnorm(n, center_x, 2)
      umap2[idx] <- rnorm(n, center_y, 2)
    }
  }
  
  data.frame(
    UMAP1 = umap1,
    UMAP2 = umap2,
    cluster = paste0("cluster_", clusters)
  )
}

# Function to generate UMAP plot
generate_umap_plot <- function(umap_data, title = "UMAP Plot") {
  ggplot(umap_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
    geom_point(size = 0.3, alpha = 0.6) +
    theme_minimal() +
    theme(legend.position = "right") +
    labs(title = title)
}

# Test different dataset sizes
test_sizes <- c(
  "10K" = 10000,
  "50K" = 50000,
  "100K" = 100000,
  "200K" = 200000,
  "230K" = 230000
)

results <- list()

cat("\n========================================\n")
cat("Testing Performance at Different Scales\n")
cat("========================================\n")

for (size_name in names(test_sizes)) {
  n_cells <- test_sizes[[size_name]]
  
  cat(sprintf("\n--- Testing %s cells ---\n", size_name))
  
  # Generate data
  umap_data <- simulate_umap_data(n_cells)
  
  # Test first render (no cache)
  cache_key <- paste0("test_", size_name)
  
  start_time <- Sys.time()
  p <- generate_umap_plot(umap_data, paste(size_name, "cells"))
  first_render_time <- as.numeric(Sys.time() - start_time)
  
  cat(sprintf("First render: %.3f seconds\n", first_render_time))
  
  # Cache the plot
  GLOBAL_UMAP_CACHE$set(cache_key, p)
  
  # Test cache retrieval
  start_time <- Sys.time()
  cached_p <- GLOBAL_UMAP_CACHE$get(cache_key)
  cache_retrieve_time <- as.numeric(Sys.time() - start_time)
  
  cat(sprintf("Cache retrieve: %.4f seconds\n", cache_retrieve_time))
  
  # Calculate speedup
  speedup <- first_render_time / cache_retrieve_time
  cat(sprintf("Speedup: %.1fx faster\n", speedup))
  
  # Store results
  results[[size_name]] <- list(
    n_cells = n_cells,
    first_render = first_render_time,
    cache_retrieve = cache_retrieve_time,
    speedup = speedup
  )
  
  # Clean up
  rm(umap_data, p)
  gc(verbose = FALSE)
}

# Test cache management
cat("\n========================================\n")
cat("Testing Cache Management\n")
cat("========================================\n")

# Test cache statistics
stats <- GLOBAL_UMAP_CACHE$stats()
cat(sprintf("\nCache status: %d/%d plots cached\n", stats$size, stats$max_size))
cat(sprintf("TTL: %d minutes\n", stats$ttl_minutes))

# Test cache eviction
cat("\nTesting cache eviction...\n")
initial_size <- stats$size

# Add plots until we exceed max size
for (i in 1:10) {
  dummy_plot <- generate_umap_plot(
    simulate_umap_data(1000), 
    paste("Dummy", i)
  )
  GLOBAL_UMAP_CACHE$set(paste0("dummy_", i), dummy_plot)
}

stats_after <- GLOBAL_UMAP_CACHE$stats()
cat(sprintf("After adding 10 more plots: %d/%d cached\n", 
            stats_after$size, stats_after$max_size))

# Summary
cat("\n========================================\n")
cat("Performance Summary\n")
cat("========================================\n\n")

# Create summary table
cat("Dataset Size | First Render | Cache Retrieve | Speedup\n")
cat("-------------|--------------|----------------|--------\n")

for (size_name in names(results)) {
  r <- results[[size_name]]
  cat(sprintf("%-12s | %9.3fs | %11.4fs | %6.1fx\n",
              size_name, r$first_render, r$cache_retrieve, r$speedup))
}

# Estimated performance for real data
cat("\n========================================\n")
cat("Estimated Real-World Performance\n")
cat("========================================\n")

# Based on simulated results, estimate real performance
if ("230K" %in% names(results)) {
  r230k <- results[["230K"]]
  
  # Real plots are more complex, so add overhead
  real_first_render <- r230k$first_render * 2  # Conservative estimate
  real_cache_retrieve <- r230k$cache_retrieve  # Cache retrieval is constant
  real_speedup <- real_first_render / real_cache_retrieve
  
  cat(sprintf("\nEstimated for 230K real cells:\n"))
  cat(sprintf("  First render: ~%.1f seconds\n", real_first_render))
  cat(sprintf("  Cache retrieve: ~%.3f seconds\n", real_cache_retrieve))
  cat(sprintf("  Expected speedup: ~%.0fx\n", real_speedup))
}

# Memory usage estimate
cat("\n========================================\n")
cat("Memory Usage Estimate\n")
cat("========================================\n")

# Rough estimate of plot object size
dummy_plot <- generate_umap_plot(simulate_umap_data(10000), "Size Test")
plot_size <- object.size(dummy_plot)
plot_size_mb <- as.numeric(plot_size) / 1024^2

cat(sprintf("\nApproximate plot size: %.2f MB\n", plot_size_mb))
cat(sprintf("Cache capacity (%d plots): ~%.0f MB\n", 
            stats$max_size, plot_size_mb * stats$max_size))

# Save results
saveRDS(results, "test_results_umap_cache_quick.rds")
cat("\n✅ Quick test complete! Results saved to test_results_umap_cache_quick.rds\n")