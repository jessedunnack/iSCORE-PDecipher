#!/usr/bin/env Rscript

# STREAMLINED SCRIPT: Generate Dataset 2 UMAP Data Only
# 
# This script ONLY generates UMAP and marker data for Dataset 2
# It does NOT recluster or save the Seurat object (saves ~1 hour)
#
# Assumes: Seurat object already has correct 15 clusters with "cluster_0" format

# Load required packages
library(Seurat)
library(dplyr)

# Install SingleCellExperiment if not available
if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  cat("Installing SingleCellExperiment...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("SingleCellExperiment")
}

if (!requireNamespace("S4Vectors", quietly = TRUE)) {
  cat("Installing S4Vectors...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("S4Vectors")
}

library(SingleCellExperiment)
library(S4Vectors)

# Dataset 2 configuration
DATASET_NAME <- "iSCORE_PD_CRISPRi"

# Auto-detect platform and set paths accordingly
if (.Platform$OS.type == "windows") {
  BASE_PATH <- "E:/ASAP/scRNASeq/PerturbSeq/final"
} else {
  BASE_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final"
}

SEURAT_FILE <- file.path(BASE_PATH, "iSCORE-PD_plus_CRISPRi", "iSCORE-PD_plus_CRISPRi.rds")
UMAP_OUTPUT_DIR <- file.path(BASE_PATH, "update_analysis_scripts", "iSCORE-PDecipher", "inst", "extdata", "umap_data")

cat("=== DATASET 2 UMAP GENERATION (STREAMLINED) ===\n\n")

# Step 1: Load Seurat object (no modifications)
cat("Step 1: Loading Seurat object (read-only)...\n")
if (!file.exists(SEURAT_FILE)) {
  stop("Seurat file not found: ", SEURAT_FILE)
}

seurat_obj <- readRDS(SEURAT_FILE)
cat("  - Loaded Seurat object with", ncol(seurat_obj), "cells\n")
cat("  - Current clustering has", length(unique(seurat_obj@meta.data$seurat_clusters)), "clusters\n")

# Verify cluster format
cluster_labels <- sort(unique(as.character(seurat_obj@meta.data$seurat_clusters)))
cat("  - Cluster labels:", paste(head(cluster_labels, 10), collapse = ", "), 
    if(length(cluster_labels) > 10) "..." else "", "\n")

# Step 2: Calculate cluster markers using MAST test
cat("\nStep 2: Calculating cluster markers with MAST test...\n")

calculate_mast_markers <- function(seurat_obj) {
  cat("  - Setting up for marker calculation...\n")
  
  # Set identity to seurat_clusters
  Idents(seurat_obj) <- "seurat_clusters"
  
  cat("  - Running FindAllMarkers with MAST test...\n")
  cat("    (This may take several minutes for large datasets)\n")
  
  markers <- FindAllMarkers(
    seurat_obj,
    test.use = "MAST",
    logfc.threshold = 0.25,    # Reasonable LFC threshold
    min.pct = 0.1,             # Expressed in at least 10% of cells
    only.pos = TRUE,           # Only positive markers
    max.cells.per.ident = 500, # Downsample for efficiency with large datasets
    verbose = TRUE
  )
  
  # Process and sort markers
  markers$cluster <- factor(markers$cluster)
  markers <- markers %>%
    arrange(cluster, desc(avg_log2FC)) %>%
    group_by(cluster) %>%
    slice_head(n = 50) %>%  # Top 50 markers per cluster
    ungroup()
  
  cat("  - Calculated", nrow(markers), "marker genes across", length(unique(markers$cluster)), "clusters\n")
  
  return(markers)
}

# Calculate markers
cluster_markers <- calculate_mast_markers(seurat_obj)

# Step 3: Create SingleCellExperiment object for UMAP visualization
cat("\nStep 3: Creating SingleCellExperiment UMAP data...\n")

create_sce_umap_data <- function(seurat_obj, markers) {
  cat("  - Extracting UMAP coordinates...\n")
  
  # Check available reductions
  available_reductions <- names(seurat_obj@reductions)
  cat("  - Available reductions:", paste(available_reductions, collapse = ", "), "\n")
  
  # Look for umap.cca reduction as per existing scripts
  if (!"umap.cca" %in% available_reductions) {
    stop("umap.cca reduction not found in Seurat object. Available reductions: ", 
         paste(available_reductions, collapse = ", "))
  }
  
  cat("  - Found umap.cca reduction\n")
  
  # Extract UMAP coordinates using Embeddings function as per existing scripts
  umap_coords <- as.matrix(Embeddings(seurat_obj, reduction = "umap.cca"))
  
  # Get metadata columns (matching extract_umap_data.R exactly)
  metadata_cols <- c(
    "seurat_clusters",
    "mutation_tidy", 
    "scMAGeCK_gene_assignment", 
    "experiments",
    "orig.ident",
    "batch",
    "cell_name"
  )
  
  # Check which columns exist
  available_cols <- intersect(metadata_cols, colnames(seurat_obj@meta.data))
  
  if (!"seurat_clusters" %in% available_cols) {
    stop("seurat_clusters column not found in metadata")
  }
  
  # Extract metadata
  cell_metadata <- seurat_obj@meta.data[, available_cols, drop = FALSE]
  
  # Add dataset identifier
  cell_metadata$dataset <- DATASET_NAME
  
  # Convert to DataFrame for SingleCellExperiment
  cell_metadata <- DataFrame(cell_metadata)
  
  cat("  - Extracted metadata with", ncol(cell_metadata), "columns\n")
  cat("  - Metadata columns:", paste(colnames(cell_metadata), collapse = ", "), "\n")
  
  # Create minimal SingleCellExperiment object (matching extract_umap_data.R format)
  sce <- SingleCellExperiment(
    assays = list(counts = matrix(nrow = 0, ncol = nrow(umap_coords))),
    colData = cell_metadata,
    reducedDims = list(UMAP = umap_coords)
  )
  
  # Add dataset metadata
  metadata(sce) <- list(
    dataset_name = DATASET_NAME,
    n_cells = ncol(sce),
    n_clusters = length(unique(cell_metadata$seurat_clusters)),
    extraction_date = Sys.Date(),
    umap_reduction = "umap.cca",
    has_markers = !is.null(markers),
    markers_file = if(!is.null(markers)) 
      file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_cluster_markers.rds")) else NA
  )
  
  cat("  - Created SingleCellExperiment with", ncol(sce), "cells and", 
      metadata(sce)$n_clusters, "clusters\n")
  
  # Verify seurat_clusters column
  if ("seurat_clusters" %in% colnames(colData(sce))) {
    clusters_in_sce <- unique(sce$seurat_clusters)
    cat("  - SCE seurat_clusters column has", length(clusters_in_sce), "unique values\n")
    cat("  - SCE cluster preview:", paste(head(sort(as.character(clusters_in_sce)), 5), collapse = ", "), "...\n")
  } else {
    cat("  - WARNING: seurat_clusters column not found in SCE colData!\n")
  }
  
  return(sce)
}

# Create SCE object
sce_data <- create_sce_umap_data(seurat_obj, cluster_markers)

# Step 4: Save UMAP and marker data for app
cat("\nStep 4: Saving UMAP and marker data for app...\n")

# Create output directory if needed
if (!dir.exists(UMAP_OUTPUT_DIR)) {
  dir.create(UMAP_OUTPUT_DIR, recursive = TRUE)
  cat("  - Created output directory:", UMAP_OUTPUT_DIR, "\n")
}

# Save UMAP data (legacy format)
umap_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_umap_data.rds"))
saveRDS(sce_data, umap_file)
cat("  - Saved UMAP data to:", umap_file, "\n")

# Save PC-specific versions (30pc, 50pc, 100pc) to match app expectations
pc_versions <- c("30pc", "50pc", "100pc")
for (pc_version in pc_versions) {
  pc_umap_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_umap_data_", pc_version, ".rds"))
  saveRDS(sce_data, pc_umap_file)
  cat("  - Saved UMAP data (", pc_version, ") to:", pc_umap_file, "\n")
}

# Save marker data
markers_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_cluster_markers.rds"))
saveRDS(cluster_markers, markers_file)
cat("  - Saved cluster markers to:", markers_file, "\n")

# Step 5: Validation and summary
cat("\nStep 5: Validation and summary...\n")

cat("  ✓ UMAP data generation completed successfully!\n")
cat("  \n")
cat("  SUMMARY:\n")
cat("  --------\n")
cat("  - Total clusters in SCE:", metadata(sce_data)$n_clusters, "\n")
cat("  - Total marker genes calculated:", nrow(cluster_markers), "\n")
cat("  - Cells processed:", ncol(sce_data), "\n")
cat("  \n")
cat("  FILES UPDATED:\n")
cat("  -------------\n")
cat("  - UMAP data (legacy):", umap_file, "\n")
for (pc_version in pc_versions) {
  pc_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_umap_data_", pc_version, ".rds"))
  cat("  - UMAP data (", pc_version, "):", pc_file, "\n")
}
cat("  - Cluster markers:", markers_file, "\n")
cat("  \n")

# Final verification
cat("  VERIFICATION:\n")
cat("  ------------\n")
verification_clusters <- sort(unique(as.character(sce_data$seurat_clusters)))
cat("  - Final cluster labels in SCE:", paste(head(verification_clusters, 10), collapse = ", "), 
    if(length(verification_clusters) > 10) "..." else "", "\n")

expected_clusters <- paste0("cluster_", 0:14)
missing_clusters <- setdiff(expected_clusters, verification_clusters)
extra_clusters <- setdiff(verification_clusters, expected_clusters)

if (length(missing_clusters) > 0) {
  cat("  ⚠ Missing expected clusters:", paste(missing_clusters, collapse = ", "), "\n")
}
if (length(extra_clusters) > 0) {
  cat("  ⚠ Extra unexpected clusters:", paste(extra_clusters, collapse = ", "), "\n")
}

if (length(missing_clusters) == 0 && length(extra_clusters) == 0) {
  cat("  ✓ Perfect match with expected 15 clusters (cluster_0 to cluster_14)!\n")
} else {
  cat("  ⚠ Cluster mismatch found - please investigate\n")
}

cat("\n=== STREAMLINED SCRIPT COMPLETED ===\n")
cat("\nIMPORTANT: \n")
cat("1. Reinstall package: devtools::install() \n")
cat("2. Test app: launch_iscore_app() \n")
cat("3. Check for 15 clusters in Dataset 2 overview\n")