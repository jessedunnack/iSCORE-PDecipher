#!/usr/bin/env Rscript

# One-time script to generate UMAP embeddings with different PC counts
# This extracts ONLY the UMAP coordinates, not the full data processing

library(Seurat)
library(SingleCellExperiment)
library(dplyr)

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

# Function to extract UMAP for specific PC count
extract_umap_for_pcs <- function(seurat_obj, reduction_name, pc_count, dataset_name) {
  
  message(sprintf("  - Processing %d PCs...", pc_count))
  
  # Check if this UMAP already exists
  if (!reduction_name %in% names(seurat_obj@reductions)) {
    message(sprintf("    Creating %s reduction with %d PCs", reduction_name, pc_count))
    
    # Check if integrated.cca reduction exists
    if (!"integrated.cca" %in% names(seurat_obj@reductions)) {
      stop("integrated.cca reduction not found. Cannot compute UMAP.")
    }
    
    # Run UMAP with specified PCs
    seurat_obj <- RunUMAP(seurat_obj, 
                         reduction = "integrated.cca", 
                         dims = 1:pc_count,
                         reduction.name = reduction_name,
                         verbose = FALSE)
  } else {
    message(sprintf("    %s already exists, extracting coordinates", reduction_name))
  }
  
  # Extract UMAP coordinates
  umap_coords <- as.matrix(Embeddings(seurat_obj, reduction = reduction_name))
  
  # Get metadata (same for all PC versions)
  metadata_cols <- c(
    "seurat_clusters",
    "seurat_clusters_fine",
    "scMAGeCK_gene_assignment",
    "experiments",
    "orig.ident"
  )
  
  available_cols <- intersect(metadata_cols, colnames(seurat_obj@meta.data))
  cell_metadata <- seurat_obj@meta.data[, available_cols, drop = FALSE]
  cell_metadata$dataset <- dataset_name
  cell_metadata <- DataFrame(cell_metadata)
  
  # Create minimal SingleCellExperiment
  sce <- SingleCellExperiment(
    assays = list(counts = matrix(nrow = 0, ncol = nrow(umap_coords))),
    colData = cell_metadata,
    reducedDims = list(UMAP = umap_coords)
  )
  
  # Add metadata
  metadata(sce) <- list(
    dataset_name = dataset_name,
    n_cells = ncol(sce),
    n_clusters = length(unique(cell_metadata$seurat_clusters)),
    extraction_date = Sys.Date(),
    umap_reduction = reduction_name,
    pc_count = pc_count
  )
  
  # Save
  output_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data_%dpc.rds", dataset_name, pc_count))
  saveRDS(sce, output_file)
  
  file_size_mb <- file.size(output_file) / 1024^2
  message(sprintf("    Saved to: %s (%.1f MB)", basename(output_file), file_size_mb))
  
  return(seurat_obj)  # Return updated object with new reductions
}

# Main processing function
process_dataset_umaps <- function(dataset_name, dataset_info) {
  message(sprintf("\n=== Processing %s ===", dataset_name))
  message(sprintf("Loading from: %s", dataset_info$file))
  
  if (!file.exists(dataset_info$file)) {
    warning(sprintf("File not found: %s", dataset_info$file))
    return(NULL)
  }
  
  # Load Seurat object
  seurat_obj <- tryCatch({
    readRDS(dataset_info$file)
  }, error = function(e) {
    warning(sprintf("Failed to load %s: %s", dataset_info$file, e$message))
    return(NULL)
  })
  
  if (is.null(seurat_obj)) return(NULL)
  
  message(sprintf("Loaded object with %d cells", ncol(seurat_obj)))
  
  # Process each PC count
  pc_configs <- list(
    list(reduction = "umap.cca", pcs = 100),
    list(reduction = "umap.cca50", pcs = 50),
    list(reduction = "umap.cca30", pcs = 30)
  )
  
  for (config in pc_configs) {
    seurat_obj <- extract_umap_for_pcs(
      seurat_obj, 
      config$reduction, 
      config$pcs, 
      dataset_name
    )
  }
  
  # Clear memory
  rm(seurat_obj)
  gc()
  
  message(sprintf("✓ Completed %s\n", dataset_name))
}

# Interactive execution
if (interactive()) {
  cat("=== Generate Multi-PC UMAP Data ===\n\n")
  cat("This script will generate UMAP embeddings for 30, 50, and 100 PCs\n")
  cat("for each dataset. This is a one-time operation.\n\n")
  
  cat("Options:\n")
  cat("1. Process all datasets\n")
  cat("2. Process specific dataset\n")
  cat("3. Exit\n\n")
  
  choice <- readline("Enter choice (1-3): ")
  
  if (choice == "1") {
    # Process all datasets
    for (dataset_name in names(DATASETS)) {
      process_dataset_umaps(dataset_name, DATASETS[[dataset_name]])
    }
    cat("\n✓ All datasets processed!\n")
    
  } else if (choice == "2") {
    # Show dataset options
    cat("\nAvailable datasets:\n")
    for (i in seq_along(DATASETS)) {
      cat(sprintf("%d. %s\n", i, names(DATASETS)[i]))
    }
    
    dataset_idx <- as.numeric(readline("\nEnter dataset number: "))
    if (dataset_idx >= 1 && dataset_idx <= length(DATASETS)) {
      dataset_name <- names(DATASETS)[dataset_idx]
      process_dataset_umaps(dataset_name, DATASETS[[dataset_name]])
    } else {
      cat("Invalid selection\n")
    }
    
  } else {
    cat("Exiting...\n")
  }
  
} else {
  # If run as script, process all
  message("Processing all datasets...")
  for (dataset_name in names(DATASETS)) {
    process_dataset_umaps(dataset_name, DATASETS[[dataset_name]])
  }
  message("\n✓ All processing complete!")
}