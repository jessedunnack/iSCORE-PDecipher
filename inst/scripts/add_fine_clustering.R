#!/usr/bin/env Rscript
# Script to add fine clustering resolution to existing UMAP data files
# This will add a new clustering at default Seurat resolution to complement the coarse clustering

library(Seurat)
library(SingleCellExperiment)
library(SummarizedExperiment)
library(dplyr)
library(stringr)

# Function to add fine clustering to SCE object
add_fine_clustering_to_sce <- function(sce_path, seurat_path, output_path = NULL) {
  message("\n=== Processing: ", basename(sce_path), " ===")
  
  # Load the SingleCellExperiment object
  sce <- readRDS(sce_path)
  
  # Load the corresponding Seurat object
  message("Loading Seurat object from: ", seurat_path)
  seurat_obj <- readRDS(seurat_path)
  
  # Check if fine clustering already exists
  if ("seurat_clusters_fine" %in% colnames(colData(sce))) {
    message("Fine clustering already exists in this file. Skipping...")
    return(sce)
  }
  
  # Get current clustering info
  current_clusters <- unique(sce$seurat_clusters)
  message("Current coarse clustering has ", length(current_clusters), " clusters")
  
  # Find clusters at default resolution (0.8)
  message("Running FindClusters at resolution 0.8 (default)...")
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
  
  # Extract the fine clustering
  fine_clusters <- Idents(seurat_obj)
  n_fine <- length(unique(fine_clusters))
  message("Fine clustering generated ", n_fine, " clusters")
  
  # Add fine clustering to SCE object
  # Match cell barcodes between Seurat and SCE
  common_cells <- intersect(colnames(sce), names(fine_clusters))
  if (length(common_cells) == 0) {
    stop("No matching cells found between SCE and Seurat objects!")
  }
  
  message("Matching ", length(common_cells), " cells between objects")
  
  # Add fine clustering to SCE colData
  sce$seurat_clusters_fine <- NA
  sce$seurat_clusters_fine[match(names(fine_clusters), colnames(sce))] <- paste0("cluster_", fine_clusters)
  
  # Verify the addition
  fine_cluster_counts <- table(sce$seurat_clusters_fine)
  message("Fine cluster distribution:")
  print(fine_cluster_counts)
  
  # Save the updated SCE object
  if (is.null(output_path)) {
    output_path <- sce_path
  }
  
  message("Saving updated SCE object to: ", output_path)
  saveRDS(sce, output_path)
  
  return(sce)
}

# Function to calculate marker genes for fine clustering
calculate_fine_markers <- function(seurat_obj, sce_with_fine) {
  message("\n=== Calculating markers for fine clustering ===")
  
  # Set fine clusters as active ident
  Idents(seurat_obj) <- sce_with_fine$seurat_clusters_fine[match(colnames(seurat_obj), colnames(sce_with_fine))]
  
  # Run FindAllMarkers
  message("Running FindAllMarkers with MAST test...")
  markers <- FindAllMarkers(
    seurat_obj,
    test.use = "MAST",
    min.pct = 0.25,
    logfc.threshold = 0.5,
    min.diff.pct = 0.2,
    only.pos = FALSE,
    verbose = FALSE
  )
  
  # Add cluster column and filter
  markers <- markers %>%
    mutate(cluster = as.character(cluster)) %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    arrange(desc(avg_log2FC)) %>%
    slice_head(n = 50) %>%  # Keep top 50 per cluster for storage
    ungroup()
  
  message("Found ", nrow(markers), " significant markers across ", 
          length(unique(markers$cluster)), " clusters")
  
  return(markers)
}

# Main execution
main <- function() {
  # Set paths - adjust these based on your directory structure
  umap_data_dir <- "inst/extdata/umap_data"
  seurat_data_dir <- "../../iSCORE-PD_plus_CRISPRi"  # Adjust this path to where Seurat objects are stored
  
  # Define datasets to process
  datasets <- list(
    list(
      name = "iSCORE_PD_CRISPRi",
      sce_pattern = "iSCORE_PD_CRISPRi_umap_data.*\\.rds",
      seurat_file = "iSCORE-PD_plus_CRISPRi.rds"  # Adjust filename
    )
    # Add more datasets as needed
  )
  
  # Process each dataset
  for (dataset in datasets) {
    message("\n========================================")
    message("Processing dataset: ", dataset$name)
    message("========================================")
    
    # Find all SCE files for this dataset (different PC versions)
    sce_files <- list.files(umap_data_dir, pattern = dataset$sce_pattern, full.names = TRUE)
    
    if (length(sce_files) == 0) {
      warning("No SCE files found for dataset: ", dataset$name)
      next
    }
    
    # Path to Seurat object
    seurat_path <- file.path(seurat_data_dir, dataset$seurat_file)
    
    if (!file.exists(seurat_path)) {
      warning("Seurat object not found at: ", seurat_path)
      warning("Please update the path in the script")
      next
    }
    
    # Load Seurat object once
    seurat_obj <- readRDS(seurat_path)
    
    # Process first SCE file and calculate markers
    sce_updated <- NULL
    for (i in seq_along(sce_files)) {
      sce_file <- sce_files[i]
      
      if (i == 1) {
        # First file: add clustering and calculate markers
        sce_updated <- add_fine_clustering_to_sce(sce_file, seurat_path)
        
        # Calculate markers for fine clustering
        markers_fine <- calculate_fine_markers(seurat_obj, sce_updated)
        
        # Save markers
        markers_output <- file.path(umap_data_dir, 
                                   paste0(dataset$name, "_cluster_markers_fine.rds"))
        message("\nSaving fine clustering markers to: ", markers_output)
        saveRDS(markers_fine, markers_output)
        
      } else {
        # For other PC versions, just add the same fine clustering
        message("\n--- Adding fine clustering to: ", basename(sce_file), " ---")
        sce <- readRDS(sce_file)
        
        # Copy fine clustering from first file
        sce$seurat_clusters_fine <- sce_updated$seurat_clusters_fine[match(colnames(sce), colnames(sce_updated))]
        
        # Save
        saveRDS(sce, sce_file)
        message("Updated: ", basename(sce_file))
      }
    }
  }
  
  message("\n=== Fine clustering addition complete! ===")
  message("Don't forget to update the Shiny app to use the new clustering data")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
} else {
  message("Script loaded. Run main() to add fine clustering to UMAP data files.")
  message("Make sure to update the dataset paths in the script first!")
}