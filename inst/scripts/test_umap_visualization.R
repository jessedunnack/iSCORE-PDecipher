#!/usr/bin/env Rscript

# Test script for UMAP visualization with minimal SingleCellExperiment objects
# This demonstrates that dittoSeq works with our lightweight data structure

library(SingleCellExperiment)
library(dittoSeq)

# Test with synthetic data first
test_synthetic_umap <- function() {
  message("Testing UMAP visualization with synthetic data...")
  
  # Create synthetic UMAP coordinates
  n_cells <- 1000
  set.seed(42)
  
  # Generate clusters
  clusters <- factor(sample(0:5, n_cells, replace = TRUE))
  
  # Generate UMAP coordinates with cluster structure
  umap_coords <- matrix(0, nrow = n_cells, ncol = 2)
  colnames(umap_coords) <- c("UMAP_1", "UMAP_2")
  
  for (i in 1:6) {
    cluster_cells <- which(clusters == (i-1))
    n_cluster <- length(cluster_cells)
    
    # Create cluster-specific coordinates
    center_x <- cos(2 * pi * i / 6) * 10
    center_y <- sin(2 * pi * i / 6) * 10
    
    umap_coords[cluster_cells, 1] <- rnorm(n_cluster, center_x, 2)
    umap_coords[cluster_cells, 2] <- rnorm(n_cluster, center_y, 2)
  }
  
  # Create metadata
  cell_metadata <- DataFrame(
    seurat_clusters = clusters,
    cell_type = sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE),
    n_genes = rpois(n_cells, lambda = 2000)
  )
  
  # Create minimal SCE
  sce <- SingleCellExperiment(
    assays = list(counts = matrix(nrow = 0, ncol = n_cells)),
    colData = cell_metadata,
    reducedDims = list(UMAP = umap_coords)
  )
  
  # Test basic plot
  p1 <- dittoDimPlot(sce, 
                     var = "seurat_clusters",
                     reduction.use = "UMAP",
                     do.label = TRUE,
                     main = "Test: Synthetic UMAP by Cluster")
  
  # Test with different variable
  p2 <- dittoDimPlot(sce,
                     var = "cell_type",
                     reduction.use = "UMAP",
                     size = 0.5,
                     main = "Test: Synthetic UMAP by Cell Type")
  
  # Save test plots
  pdf("test_synthetic_umap.pdf", width = 12, height = 5)
  print(cowplot::plot_grid(p1, p2, ncol = 2))
  dev.off()
  
  message("✓ Synthetic data test completed successfully!")
  message("  Check test_synthetic_umap.pdf")
  
  return(sce)
}

# Test with actual extracted data
test_extracted_umap <- function(data_file) {
  if (!file.exists(data_file)) {
    message("Data file not found: ", data_file)
    return(NULL)
  }
  
  message("Testing UMAP visualization with extracted data...")
  message("Loading from: ", data_file)
  
  # Load data
  if (grepl("combined", data_file)) {
    sce_list <- readRDS(data_file)
    dataset_names <- names(sce_list)
    message("Found datasets: ", paste(dataset_names, collapse = ", "))
    
    # Test first dataset
    if (length(sce_list) > 0) {
      sce <- sce_list[[1]]
      dataset_name <- dataset_names[1]
    } else {
      message("No datasets in combined file")
      return(NULL)
    }
  } else {
    sce <- readRDS(data_file)
    dataset_name <- "Dataset"
  }
  
  # Check structure
  message("\nData structure:")
  message("- Cells: ", ncol(sce))
  message("- Reductions: ", paste(names(reducedDims(sce)), collapse = ", "))
  message("- Metadata columns: ", paste(names(colData(sce)), collapse = ", "))
  
  # Create plots
  plots <- list()
  
  # Plot 1: Coarse clusters
  if ("seurat_clusters" %in% names(colData(sce))) {
    plots$p1 <- dittoDimPlot(sce,
                            var = "seurat_clusters",
                            reduction.use = "UMAP",
                            do.label = TRUE,
                            main = paste(dataset_name, "- Coarse Clusters"))
  }
  
  # Plot 2: Fine clusters
  if ("seurat_clusters_fine" %in% names(colData(sce))) {
    plots$p2 <- dittoDimPlot(sce,
                            var = "seurat_clusters_fine",
                            reduction.use = "UMAP",
                            do.label = FALSE,  # Too many to label
                            legend.show = FALSE,  # Too many for legend
                            main = paste(dataset_name, "- Fine Clusters"))
  }
  
  # Plot 3: Perturbations (if available)
  if ("scMAGeCK_gene_assignment" %in% names(colData(sce))) {
    # Check number of unique perturbations
    n_perturbs <- length(unique(colData(sce)$scMAGeCK_gene_assignment))
    message("- Number of perturbations: ", n_perturbs)
    
    if (n_perturbs < 50) {
      plots$p3 <- dittoDimPlot(sce,
                              var = "scMAGeCK_gene_assignment",
                              reduction.use = "UMAP",
                              size = 0.3,
                              main = paste(dataset_name, "- Perturbations"))
    } else {
      # Too many perturbations - just show Non-Targeting vs others
      sce$is_control <- colData(sce)$scMAGeCK_gene_assignment == "Non-Targeting"
      plots$p3 <- dittoDimPlot(sce,
                              var = "is_control",
                              reduction.use = "UMAP",
                              size = 0.3,
                              main = paste(dataset_name, "- Control vs Perturbed"))
    }
  }
  
  # Save plots
  output_file <- paste0("test_extracted_", dataset_name, "_umap.pdf")
  pdf(output_file, width = 15, height = 5)
  if (length(plots) > 0) {
    print(cowplot::plot_grid(plotlist = plots, ncol = 3))
  }
  dev.off()
  
  message("✓ Extracted data test completed successfully!")
  message("  Check ", output_file)
  
  return(sce)
}

# Main test function
main <- function() {
  message("UMAP Visualization Test Script")
  message("==============================\n")
  
  # Test 1: Synthetic data
  sce_synthetic <- test_synthetic_umap()
  
  # Test 2: Extracted data (if available)
  data_dir <- "../../inst/extdata/umap_data"
  combined_file <- file.path(data_dir, "all_umap_data_combined.rds")
  
  if (file.exists(combined_file)) {
    sce_extracted <- test_extracted_umap(combined_file)
  } else {
    message("\nExtracted data not found at: ", combined_file)
    message("Run extract_umap_data.R first to generate the data.")
  }
  
  message("\nAll tests completed!")
  
  # Return test data for interactive exploration
  invisible(list(synthetic = sce_synthetic, extracted = sce_extracted))
}

# Run tests
if (!interactive()) {
  main()
} else {
  message("Test script loaded. Run main() to execute tests.")
}