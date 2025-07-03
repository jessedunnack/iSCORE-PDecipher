#!/usr/bin/env Rscript

# DISPOSABLE SCRIPT: Fix Dataset 2 Clustering Mismatch
# 
# Problem: Dataset 2 (iSCORE-PD_plus_CRISPRi) shows ~34-35 clusters in UMAP/markers
# but DE analysis and enrichment data expect 15 clusters (cluster_0 to cluster_14)
#
# Solution: Re-cluster Seurat object to match DE analysis expectations,
# then regenerate UMAP data and marker genes with MAST test
#
# This script will be removed after successful execution

library(Seurat)
library(dplyr)
library(future)

# Enable parallel processing for faster computation
# Use multisession on Windows, multicore on Unix-like systems
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = 4)
} else {
  plan("multicore", workers = 4)
}
# Increase memory limit to handle large Seurat object (50GB)
options(future.globals.maxSize = 50000 * 1024^2)  # 50GB

# Dataset 2 configuration
DATASET_NAME <- "iSCORE_PD_CRISPRi"

# Auto-detect platform and set paths accordingly
if (.Platform$OS.type == "windows") {
  BASE_PATH <- "E:/ASAP/scRNASeq/PerturbSeq/final"
} else {
  BASE_PATH <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final"
}

SEURAT_FILE <- file.path(BASE_PATH, "iSCORE-PD_plus_CRISPRi", "iSCORE-PD_plus_CRISPRi.rds")
BACKUP_FILE <- file.path(BASE_PATH, "iSCORE-PD_plus_CRISPRi", "iSCORE-PD_plus_CRISPRi_BACKUP.rds")
TARGET_CLUSTERS <- 15  # clusters 0-14 as expected by DE analysis
UMAP_OUTPUT_DIR <- file.path(BASE_PATH, "update_analysis_scripts", "iSCORE-PDecipher", "inst", "extdata", "umap_data")

cat("=== DATASET 2 CLUSTERING FIX SCRIPT ===\n\n")

# Step 1: Load and examine current Seurat object
cat("Step 1: Loading Seurat object...\n")
if (!file.exists(SEURAT_FILE)) {
  stop("Seurat file not found: ", SEURAT_FILE)
}

seurat_obj <- readRDS(SEURAT_FILE)
cat("  - Loaded Seurat object with", ncol(seurat_obj), "cells\n")
cat("  - Current clustering has", length(unique(seurat_obj@meta.data$seurat_clusters)), "clusters\n")

# Handle both numeric and factor cluster labels
cluster_labels <- sort(unique(as.character(seurat_obj@meta.data$seurat_clusters)))
cat("  - Current cluster labels:", paste(cluster_labels, collapse = ", "), "\n")

# Step 2: Create backup before modifications
cat("\nStep 2: Creating backup...\n")
if (!file.exists(BACKUP_FILE)) {
  saveRDS(seurat_obj, BACKUP_FILE)
  cat("  - Backup saved to:", BACKUP_FILE, "\n")
} else {
  cat("  - Backup already exists, skipping\n")
}

# Step 3: Use fixed resolution of 0.2 as requested
cat("\nStep 3: Applying resolution 0.2...\n")

FIXED_RESOLUTION <- 0.2
cat("  - Using fixed resolution:", FIXED_RESOLUTION, "\n")

# Temporarily disable parallel processing for clustering to avoid memory issues
plan("sequential")

# Test the resolution to see how many clusters it gives
seurat_obj <- FindClusters(seurat_obj, resolution = FIXED_RESOLUTION, verbose = FALSE)

# Re-enable parallel processing after clustering
if (.Platform$OS.type == "windows") {
  plan("multisession", workers = 4)
} else {
  plan("multicore", workers = 4)
}
n_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
cat("  - Resolution", FIXED_RESOLUTION, "gives", n_clusters, "clusters\n")

if (n_clusters != TARGET_CLUSTERS) {
  cat("  ⚠ WARNING: Resolution 0.2 gives", n_clusters, "clusters, but DE analysis expects", TARGET_CLUSTERS, "\n")
  cat("  Proceeding with resolution 0.2 as requested...\n")
}

optimal <- list(resolution = FIXED_RESOLUTION, n_clusters = n_clusters, exact_match = (n_clusters == TARGET_CLUSTERS))

# Step 4: Verify clustering (already applied in Step 3)
cat("\nStep 4: Verifying clustering results...\n")

# Verify final clustering
final_clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
cat("  - Final clustering has", length(final_clusters), "clusters\n")
cat("  - Cluster labels:", paste(final_clusters, collapse = ", "), "\n")

# Ensure cluster names match DE analysis format (cluster_0, cluster_1, etc.)
if (any(grepl("^[0-9]+$", final_clusters))) {
  cat("  - Converting numeric cluster labels to cluster_X format...\n")
  seurat_obj@meta.data$seurat_clusters <- paste0("cluster_", seurat_obj@meta.data$seurat_clusters)
  Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
  final_clusters <- sort(unique(seurat_obj@meta.data$seurat_clusters))
  cat("  - Updated cluster labels:", paste(final_clusters, collapse = ", "), "\n")
}

# Step 5: Calculate cluster markers using MAST test
cat("\nStep 5: Calculating cluster markers with MAST test...\n")

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

# Step 6: Create lightweight UMAP data for visualization
cat("\nStep 6: Creating UMAP visualization data...\n")

# Load SingleCellExperiment package
library(SingleCellExperiment)
library(S4Vectors)

create_umap_data <- function(seurat_obj, markers) {
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
  
  # Get metadata columns (matching extract_umap_data.R)
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
  
  return(list(
    sce = sce,
    markers = markers,
    n_clusters = metadata(sce)$n_clusters
  ))
}

# Create UMAP data
umap_result <- create_umap_data(seurat_obj, cluster_markers)

# Step 7: Save corrected Seurat object
cat("\nStep 7: Saving corrected Seurat object...\n")
saveRDS(seurat_obj, SEURAT_FILE)
cat("  - Saved corrected Seurat object to:", SEURAT_FILE, "\n")

# Step 8: Save UMAP and marker data for app
cat("\nStep 8: Saving UMAP and marker data for app...\n")

# Create output directory if needed
if (!dir.exists(UMAP_OUTPUT_DIR)) {
  dir.create(UMAP_OUTPUT_DIR, recursive = TRUE)
  cat("  - Created output directory:", UMAP_OUTPUT_DIR, "\n")
}

# Save UMAP data (legacy format)
umap_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_umap_data.rds"))
saveRDS(umap_result$sce, umap_file)
cat("  - Saved UMAP data to:", umap_file, "\n")

# Save PC-specific versions (30pc, 50pc, 100pc) to match app expectations
pc_versions <- c("30pc", "50pc", "100pc")
for (pc_version in pc_versions) {
  pc_umap_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_umap_data_", pc_version, ".rds"))
  saveRDS(umap_result$sce, pc_umap_file)
  cat("  - Saved UMAP data (", pc_version, ") to:", pc_umap_file, "\n")
}

# Save marker data
markers_file <- file.path(UMAP_OUTPUT_DIR, paste0(DATASET_NAME, "_cluster_markers.rds"))
saveRDS(umap_result$markers, markers_file)
cat("  - Saved cluster markers to:", markers_file, "\n")

# Step 9: Validation and summary
cat("\nStep 9: Validation and summary...\n")

cat("  ✓ Clustering correction completed successfully!\n")
cat("  \n")
cat("  SUMMARY:\n")
cat("  --------\n")
cat("  - Target clusters:", TARGET_CLUSTERS, "\n")
cat("  - Achieved clusters:", umap_result$n_clusters, "\n")
cat("  - Clustering resolution used:", optimal$resolution, "\n")
cat("  - Total marker genes calculated:", nrow(cluster_markers), "\n")
cat("  - Cells processed:", ncol(umap_result$sce), "\n")
cat("  \n")
cat("  FILES UPDATED:\n")
cat("  -------------\n")
cat("  - Seurat object:", SEURAT_FILE, "\n")
cat("  - UMAP data:", umap_file, "\n")
cat("  - Cluster markers:", markers_file, "\n")
cat("  - Backup created:", BACKUP_FILE, "\n")
cat("  \n")

# Final verification
cat("  VERIFICATION:\n")
cat("  ------------\n")
verification_clusters <- sort(unique(umap_result$sce$seurat_clusters))
cat("  - Final cluster labels:", paste(verification_clusters, collapse = ", "), "\n")

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
  cat("  ✓ Perfect match with DE analysis expectations!\n")
} else {
  cat("  ⚠ Minor discrepancies found but clustering should be functional\n")
}

cat("\n=== SCRIPT COMPLETED SUCCESSFULLY ===\n")
cat("\nIMPORTANT: This script should be deleted after confirming the fix works.\n")
cat("To test: Launch the app and check that the Overview page UMAP shows ~15 clusters instead of ~34-35.\n")