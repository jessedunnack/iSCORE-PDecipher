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
    description = "iSCORE-PD dataset (mutations only)"
  ),
  iSCORE_PD_CRISPRi = list(
    file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
    description = "iSCORE-PD with CRISPRi perturbations"
  ),
  Full_Dataset = list(
    file = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa/full_dataset.rds",
    description = "Complete dataset with CRISPRi and CRISPRa"
  )
)

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  message("Created output directory: ", OUTPUT_DIR)
}

# Function to extract minimal data from Seurat object
extract_minimal_umap_data <- function(seurat_obj, dataset_name) {
  message(sprintf("\nExtracting UMAP data for %s dataset...", dataset_name))
  
  # Get UMAP coordinates
  if (!"umap.cca" %in% names(seurat_obj@reductions)) {
    stop("umap.cca reduction not found in Seurat object")
  }
  
  umap_coords <- as.matrix(Embeddings(seurat_obj, reduction = "umap.cca"))
  message(sprintf("  - Extracted UMAP coordinates for %d cells", nrow(umap_coords)))
  
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
  
  # Add dataset metadata
  metadata(sce) <- list(
    dataset_name = dataset_name,
    n_cells = ncol(sce),
    n_clusters_coarse = length(unique(cell_metadata$seurat_clusters)),
    n_clusters_fine = if("seurat_clusters_fine" %in% available_cols) 
      length(unique(cell_metadata$seurat_clusters_fine)) else NA,
    extraction_date = Sys.Date(),
    umap_reduction = "umap.cca"
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
      extraction_date = as.character(meta$extraction_date)
    )
  }))
  
  return(summary_df)
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
    
    # Extract minimal data
    sce <- tryCatch({
      extract_minimal_umap_data(seurat_obj, dataset_name)
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