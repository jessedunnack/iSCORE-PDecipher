#!/usr/bin/env Rscript

# Script to extract minimal UMAP visualization data from large Seurat objects
# Creates lightweight SingleCellExperiment objects for Shiny app visualization

library(Seurat)
library(SingleCellExperiment)
library(dplyr)
library(tidyr)

# Configuration
OUTPUT_DIR <- "E:/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data"

# Dataset configurations
DATASETS <- list(
  iSCORE_PD = list(
    file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/final_iSCORE-PD.rds",
    description = "iSCORE-PD dataset (mutations only)",
    expected_clusters = 15  # Based on DE analysis having cluster_0 to cluster_14
  ),
  iSCORE_PD_CRISPRi = list(
    file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
    description = "iSCORE-PD with CRISPRi perturbations",
    expected_clusters = NULL  # Will be determined from DE data
  ),
  Full_Dataset = list(
    file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_dataset.rds",
    description = "Complete dataset with CRISPRi and CRISPRa",
    expected_clusters = NULL  # Will be determined from DE data
  )
)

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("Created output directory: ", OUTPUT_DIR)
}

# Function to validate clusters and recalculate if needed
validate_and_update_clusters <- function(seurat_obj, expected_clusters = NULL, clustering_resolution = NULL) {
  current_clusters <- unique(seurat_obj@meta.data$seurat_clusters)
  n_current <- length(current_clusters)
  
  message(sprintf("  - Current clustering has %d clusters", n_current))
  
  # If expected clusters specified, check if they match
  if (!is.null(expected_clusters)) {
    if (n_current != expected_clusters) {
      message(sprintf("  - Expected %d clusters but found %d", expected_clusters, n_current))
      
      if (is.null(clustering_resolution)) {
        # Ask user for resolution
        message("  - Clustering resolution required for recalculation")
        cat("Current clusters do not match expected number.\n")
        cat(sprintf("Expected: %d, Current: %d\n", expected_clusters, n_current))
        clustering_resolution <- as.numeric(readline("Enter clustering resolution (e.g., 0.5): "))
      }
      
      message(sprintf("  - Recalculating clusters with resolution %.2f", clustering_resolution))
      
      # Recalculate clusters
      seurat_obj <- FindClusters(seurat_obj, resolution = clustering_resolution, verbose = FALSE)
      new_clusters <- unique(seurat_obj@meta.data$seurat_clusters)
      message(sprintf("  - New clustering has %d clusters", length(new_clusters)))
    }
  }
  
  return(seurat_obj)
}

# Function to calculate and save marker genes
calculate_cluster_markers <- function(seurat_obj, dataset_name, force_recalculate = FALSE) {
  markers_file <- file.path(OUTPUT_DIR, sprintf("%s_cluster_markers.rds", dataset_name))
  
  # Check if markers already exist
  if (file.exists(markers_file) && !force_recalculate) {
    message(sprintf("  - Loading existing markers from %s", markers_file))
    return(readRDS(markers_file))
  }
  
  message("  - Calculating cluster markers using MAST test...")
  
  # Set identity to seurat_clusters
  Idents(seurat_obj) <- "seurat_clusters"
  
  # Calculate markers using MAST (efficient and appropriate for droplet data)
  tryCatch({
    markers <- FindAllMarkers(
      seurat_obj,
      test.use = "MAST",
      logfc.threshold = 0.25,    # Reasonable LFC threshold
      min.pct = 0.1,             # Expressed in at least 10% of cells
      only.pos = TRUE,           # Only positive markers
      max.cells.per.ident = 500, # Downsample for efficiency
      verbose = FALSE
    )
    
    # Add cluster information and sort
    markers$cluster <- factor(markers$cluster)
    markers <- markers %>%
      arrange(cluster, desc(avg_log2FC)) %>%
      group_by(cluster) %>%
      slice_head(n = 50) %>%  # Keep top 50 markers per cluster
      ungroup()
    
    # Save markers
    saveRDS(markers, markers_file)
    message(sprintf("  - Saved %d markers to %s", nrow(markers), markers_file))
    
    # Print summary
    marker_summary <- markers %>%
      group_by(cluster) %>%
      summarise(n_markers = n(), 
                avg_logfc = round(mean(avg_log2FC), 2),
                .groups = 'drop')
    
    message("  - Marker summary by cluster:")
    print(marker_summary)
    
    return(markers)
    
  }, error = function(e) {
    warning(sprintf("Failed to calculate markers for %s: %s", dataset_name, e$message))
    return(NULL)
  })
}

# Function to extract minimal data from Seurat object with marker calculation
extract_minimal_umap_data <- function(seurat_obj, dataset_name, expected_clusters = NULL, 
                                      clustering_resolution = NULL, force_recalculate = FALSE) {
  message(sprintf("\nExtracting UMAP data for %s dataset...", dataset_name))
  
  # Get UMAP coordinates
  if (!"umap.cca" %in% names(seurat_obj@reductions)) {
    stop("umap.cca reduction not found in Seurat object")
  }
  
  umap_coords <- as.matrix(Embeddings(seurat_obj, reduction = "umap.cca"))
  message(sprintf("  - Extracted UMAP coordinates for %d cells", nrow(umap_coords)))
  
  # Validate and potentially recalculate clusters
  seurat_obj <- validate_and_update_clusters(seurat_obj, expected_clusters, clustering_resolution)
  
  # Calculate cluster markers
  markers <- calculate_cluster_markers(seurat_obj, dataset_name, force_recalculate)
  
  # Get metadata
  metadata_cols <- c(
    "seurat_clusters",
    "seurat_clusters_fine",
    "scMAGeCK_gene_assignment",
    "experiments",
    "orig.ident"
  )
  
  # Check which columns exist
  available_cols <- intersect(metadata_cols, colnames(seurat_obj@meta.data))
  
  if (!"seurat_clusters" %in% available_cols) {
    stop("seurat_clusters column not found in metadata")
  }
  
  # Extract metadata
  cell_metadata <- seurat_obj@meta.data[, available_cols, drop = FALSE]
  
  # Add dataset identifier
  cell_metadata$dataset <- dataset_name
  
  # Convert to DataFrame for SingleCellExperiment
  cell_metadata <- DataFrame(cell_metadata)
  
  message(sprintf("  - Extracted metadata with %d columns", ncol(cell_metadata)))
  
  # Create minimal SingleCellExperiment object
  # Use empty assay to minimize memory usage
  sce <- SingleCellExperiment(
    assays = list(counts = matrix(nrow = 0, ncol = nrow(umap_coords))),
    colData = cell_metadata,
    reducedDims = list(UMAP = umap_coords)
  )
  
  # Add dataset metadata including markers
  metadata(sce) <- list(
    dataset_name = dataset_name,
    n_cells = ncol(sce),
    n_clusters_coarse = length(unique(cell_metadata$seurat_clusters)),
    n_clusters_fine = if("seurat_clusters_fine" %in% available_cols) 
      length(unique(cell_metadata$seurat_clusters_fine)) else NA,
    extraction_date = Sys.Date(),
    umap_reduction = "umap.cca",
    has_markers = !is.null(markers),
    markers_file = if(!is.null(markers)) file.path(OUTPUT_DIR, sprintf("%s_cluster_markers.rds", dataset_name)) else NA
  )
  
  return(sce)
}

# Function to create summary statistics
create_summary_stats <- function(sce_list) {
  summary_df <- do.call(rbind, lapply(names(sce_list), function(dataset) {
    sce <- sce_list[[dataset]]
    meta <- metadata(sce)
    
    data.frame(
      dataset = dataset,
      n_cells = meta$n_cells,
      n_clusters_coarse = meta$n_clusters_coarse,
      n_clusters_fine = meta$n_clusters_fine,
      has_perturbation_info = "scMAGeCK_gene_assignment" %in% colnames(colData(sce)),
      has_markers = meta$has_markers %||% FALSE,
      extraction_date = as.character(meta$extraction_date)
    )
  }))
  
  return(summary_df)
}

# Utility function for Shiny app to load markers
load_markers_for_dataset <- function(dataset_name, output_dir = OUTPUT_DIR) {
  markers_file <- file.path(output_dir, sprintf("%s_cluster_markers.rds", dataset_name))
  
  if (file.exists(markers_file)) {
    return(readRDS(markers_file))
  } else {
    warning(sprintf("Markers file not found: %s", markers_file))
    return(NULL)
  }
}

# Main extraction process
main <- function() {
  message("Starting UMAP data extraction process...")
  message("=======================================")
  
  sce_list <- list()
  
  for (dataset_name in names(DATASETS)) {
    dataset_info <- DATASETS[[dataset_name]]
    seurat_file <- dataset_info$file
    
    if (!file.exists(seurat_file)) {
      warning(sprintf("Seurat file not found: %s", seurat_file))
      next
    }
    
    message(sprintf("\nProcessing %s dataset...", dataset_name))
    message(sprintf("Loading from: %s", seurat_file))
    
    # Load Seurat object
    seurat_obj <- tryCatch({
      readRDS(seurat_file)
    }, error = function(e) {
      warning(sprintf("Failed to load %s: %s", seurat_file, e$message))
      return(NULL)
    })
    
    if (is.null(seurat_obj)) next
    
    # Extract minimal data with expected clusters
    sce <- tryCatch({
      extract_minimal_umap_data(seurat_obj, dataset_name, 
                               expected_clusters = dataset_info$expected_clusters)
    }, error = function(e) {
      warning(sprintf("Failed to extract data from %s: %s", dataset_name, e$message))
      return(NULL)
    })
    
    if (!is.null(sce)) {
      # Save individual SCE object
      output_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data.rds", dataset_name))
      saveRDS(sce, output_file)
      
      # Check file size
      file_size_mb <- file.size(output_file) / 1024^2
      message(sprintf("  - Saved to: %s (%.1f MB)", output_file, file_size_mb))
      
      sce_list[[dataset_name]] <- sce
    }
    
    # Clear memory
    rm(seurat_obj)
    gc()
  }
  
  # Save combined data
  if (length(sce_list) > 0) {
    combined_file <- file.path(OUTPUT_DIR, "all_umap_data_combined.rds")
    saveRDS(sce_list, combined_file)
    message(sprintf("\nSaved combined data to: %s", combined_file))
    
    # Create summary statistics
    summary_stats <- create_summary_stats(sce_list)
    summary_file <- file.path(OUTPUT_DIR, "umap_data_summary.csv")
    write.csv(summary_stats, summary_file, row.names = FALSE)
    message(sprintf("Saved summary statistics to: %s", summary_file))
    
    # Print summary
    message("\nExtraction Summary:")
    message("===================")
    print(summary_stats)
  }
  
  message("\nExtraction process complete!")
}

# Run the script
if (!interactive()) {
  main()
} else {
  message("Script loaded. Run main() to extract UMAP data.")
}