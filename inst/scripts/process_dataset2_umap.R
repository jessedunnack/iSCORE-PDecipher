#!/usr/bin/env Rscript

# Interactive script to process UMAP data for dataset 2 (iSCORE_PD_plus_CRISPRi) only
# Run with: Rscript process_dataset2_umap.R [options]

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(dplyr)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-p", "--pc"), type="integer", default=30,
              help="Number of PCs for UMAP (default: 30)"),
  make_option(c("-m", "--markers"), type="logical", default=TRUE,
              help="Calculate cluster markers (default: TRUE)"),
  make_option(c("-f", "--force"), type="logical", default=FALSE,
              help="Force overwrite existing files without prompting (default: FALSE)"),
  make_option(c("-r", "--resolution"), type="numeric", default=NULL,
              help="Clustering resolution if reclustering needed"),
  make_option(c("-t", "--top-markers"), type="integer", default=10,
              help="Number of top markers to display per cluster (default: 10)"),
  make_option(c("-s", "--save-all-pc"), type="logical", default=FALSE,
              help="Save UMAP for all PC counts (30, 50, 100) (default: FALSE)"),
  make_option(c("-c", "--max-cells"), type="integer", default=500,
              help="Max cells per cluster for marker calculation (default: 500, use -1 for all cells)")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Process UMAP data for iSCORE_PD_plus_CRISPRi dataset")
opt <- parse_args(opt_parser)

# Configuration
DATASET_NAME <- "iSCORE_PD_CRISPRi"
SEURAT_FILE <- "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"
OUTPUT_DIR <- "E:/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data"
EXPECTED_CLUSTERS <- 15  # Based on DE analysis

cat("\n=== UMAP Processing for Dataset 2 (iSCORE_PD_plus_CRISPRi) ===\n")
cat(sprintf("PC count: %d\n", opt$pc))
cat(sprintf("Calculate markers: %s\n", opt$markers))
cat(sprintf("Force overwrite: %s\n", opt$force))
cat(sprintf("Save all PC versions: %s\n", opt$`save-all-pc`))
cat(sprintf("Max cells per cluster: %s\n", ifelse(opt$`max-cells` == -1, "all", opt$`max-cells`)))
cat("\n")

# Check if output directory exists
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", OUTPUT_DIR))
}

# Function to check file existence and prompt for overwrite
check_overwrite <- function(filepath, force = FALSE) {
  if (file.exists(filepath)) {
    if (!force) {
      cat(sprintf("\nFile already exists: %s\n", basename(filepath)))
      response <- readline("Overwrite? (y/n): ")
      return(tolower(response) == "y")
    }
    return(TRUE)
  }
  return(TRUE)
}

# Function to extract UMAP data for a specific PC count
extract_umap_for_pc <- function(seurat_obj, pc_count, dataset_name) {
  cat(sprintf("\nExtracting UMAP data for %d PCs...\n", pc_count))
  
  # Get the appropriate UMAP reduction
  umap_name <- ifelse(pc_count == 30, "umap.cca", 
                     ifelse(pc_count == 50, "umap.cca.50pc",
                           ifelse(pc_count == 100, "umap.cca.100pc", "umap.cca")))
  
  if (!umap_name %in% names(seurat_obj@reductions)) {
    warning(sprintf("UMAP reduction '%s' not found. Available: %s", 
                   umap_name, paste(names(seurat_obj@reductions), collapse=", ")))
    return(NULL)
  }
  
  # Extract UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, umap_name)
  cat(sprintf("  - Extracted UMAP coordinates: %d cells x 2 dimensions\n", nrow(umap_coords)))
  
  # Extract metadata
  metadata_cols <- c("seurat_clusters", "seurat_clusters_fine", 
                    "scMAGeCK_gene_assignment", "experiments", "orig.ident")
  available_cols <- intersect(metadata_cols, colnames(seurat_obj@meta.data))
  cell_metadata <- seurat_obj@meta.data[, available_cols, drop = FALSE]
  cell_metadata$dataset <- dataset_name
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment(
    assays = list(counts = matrix(nrow = 0, ncol = nrow(umap_coords))),
    colData = DataFrame(cell_metadata),
    reducedDims = list(UMAP = umap_coords)
  )
  
  # Add metadata
  metadata(sce) <- list(
    dataset_name = dataset_name,
    n_cells = ncol(sce),
    n_clusters = length(unique(cell_metadata$seurat_clusters)),
    extraction_date = Sys.Date(),
    pc_count = pc_count,
    umap_reduction = umap_name
  )
  
  return(sce)
}

# Function to calculate markers with progress
calculate_markers_with_progress <- function(seurat_obj, top_n = 10, max_cells = 500) {
  cat("\nCalculating cluster markers...\n")
  cat("This may take 5-10 minutes for ~200K cells\n")
  
  Idents(seurat_obj) <- "seurat_clusters"
  n_clusters <- length(unique(Idents(seurat_obj)))
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = n_clusters, style = 3)
  all_markers <- list()
  
  for (i in 0:(n_clusters-1)) {
    setTxtProgressBar(pb, i)
    cluster_name <- paste0("cluster_", i)
    
    markers <- tryCatch({
      FindMarkers(seurat_obj, 
                 ident.1 = as.character(i),
                 test.use = "MAST",
                 logfc.threshold = 0.5,
                 min.pct = 0.25,
                 min.diff.pct = 0.1,
                 only.pos = TRUE,
                 max.cells.per.ident = ifelse(max_cells == -1, Inf, max_cells),
                 verbose = FALSE)
    }, error = function(e) {
      warning(sprintf("Failed to find markers for cluster %d: %s", i, e$message))
      return(NULL)
    })
    
    if (!is.null(markers) && nrow(markers) > 0) {
      markers$gene <- rownames(markers)
      markers$cluster <- as.character(i)
      all_markers[[cluster_name]] <- markers
    }
  }
  
  close(pb)
  
  # Combine all markers
  combined_markers <- do.call(rbind, all_markers)
  combined_markers$cluster <- factor(combined_markers$cluster, 
                                   levels = as.character(0:(n_clusters-1)))
  
  # Print condensed summary
  cat("\n\n=== TOP MARKERS PER CLUSTER (Condensed View) ===\n")
  for (i in 0:(n_clusters-1)) {
    cluster_markers <- combined_markers %>%
      filter(cluster == as.character(i)) %>%
      arrange(desc(avg_log2FC)) %>%
      head(top_n)
    
    if (nrow(cluster_markers) > 0) {
      cat(sprintf("\nCluster %d: ", i))
      top_genes <- paste(cluster_markers$gene[1:min(5, nrow(cluster_markers))], collapse=", ")
      cat(sprintf("%s", top_genes))
      if (nrow(cluster_markers) > 5) cat("...")
    }
  }
  cat("\n\n")
  
  return(combined_markers)
}

# Main processing
tryCatch({
  # Check if Seurat file exists
  if (!file.exists(SEURAT_FILE)) {
    stop(sprintf("Seurat file not found: %s", SEURAT_FILE))
  }
  
  cat(sprintf("Loading Seurat object from: %s\n", SEURAT_FILE))
  cat("This may take a minute...\n")
  
  seurat_obj <- readRDS(SEURAT_FILE)
  
  cat(sprintf("\nLoaded dataset with %d cells and %d genes\n", 
              ncol(seurat_obj), nrow(seurat_obj)))
  
  # Check clusters
  current_clusters <- unique(seurat_obj@meta.data$seurat_clusters)
  n_clusters <- length(current_clusters)
  cat(sprintf("Current clustering: %d clusters (expected: %d)\n", n_clusters, EXPECTED_CLUSTERS))
  
  if (n_clusters != EXPECTED_CLUSTERS) {
    cat("\nWARNING: Cluster count mismatch!\n")
    if (!is.null(opt$resolution)) {
      cat(sprintf("Reclustering with resolution %.2f...\n", opt$resolution))
      seurat_obj <- FindClusters(seurat_obj, resolution = opt$resolution, verbose = FALSE)
      new_clusters <- length(unique(seurat_obj@meta.data$seurat_clusters))
      cat(sprintf("New clustering: %d clusters\n", new_clusters))
    } else {
      cat("Use --resolution flag to recluster if needed\n")
    }
  }
  
  # Process UMAP data
  if (opt$`save-all-pc`) {
    # Save all PC versions
    pc_counts <- c(30, 50, 100)
    cat("\nProcessing UMAP for all PC counts (30, 50, 100)...\n")
  } else {
    # Save only requested PC count
    pc_counts <- c(opt$pc)
  }
  
  saved_files <- c()
  
  for (pc in pc_counts) {
    # Extract UMAP data
    sce <- extract_umap_for_pc(seurat_obj, pc, DATASET_NAME)
    
    if (!is.null(sce)) {
      # Determine output filename
      if (pc == 30) {
        # Save both legacy and new format for 30 PC
        legacy_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data.rds", DATASET_NAME))
        if (check_overwrite(legacy_file, opt$force)) {
          saveRDS(sce, legacy_file)
          saved_files <- c(saved_files, legacy_file)
          cat(sprintf("  - Saved legacy format: %s\n", basename(legacy_file)))
        }
      }
      
      output_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data_%dpc.rds", DATASET_NAME, pc))
      if (check_overwrite(output_file, opt$force)) {
        saveRDS(sce, output_file)
        saved_files <- c(saved_files, output_file)
        file_size_mb <- file.size(output_file) / 1024^2
        cat(sprintf("  - Saved: %s (%.1f MB)\n", basename(output_file), file_size_mb))
      }
    }
  }
  
  # Calculate markers if requested
  if (opt$markers) {
    markers_file <- file.path(OUTPUT_DIR, sprintf("%s_cluster_markers.rds", DATASET_NAME))
    
    if (check_overwrite(markers_file, opt$force)) {
      markers <- calculate_markers_with_progress(seurat_obj, opt$`top-markers`, opt$`max-cells`)
      
      if (!is.null(markers)) {
        saveRDS(markers, markers_file)
        saved_files <- c(saved_files, markers_file)
        cat(sprintf("\nSaved markers to: %s\n", basename(markers_file)))
        
        # Update SCE objects with marker info
        for (file in saved_files) {
          if (grepl("umap_data", file)) {
            sce <- readRDS(file)
            metadata(sce)$has_markers <- TRUE
            metadata(sce)$markers_file <- markers_file
            saveRDS(sce, file)
          }
        }
      }
    }
  }
  
  # Summary
  cat("\n=== PROCESSING COMPLETE ===\n")
  cat(sprintf("Dataset: %s\n", DATASET_NAME))
  cat(sprintf("Cells: %d\n", ncol(seurat_obj)))
  cat(sprintf("Clusters: %d\n", length(unique(seurat_obj@meta.data$seurat_clusters))))
  cat(sprintf("Files saved: %d\n", length(saved_files)))
  
  if (length(saved_files) > 0) {
    cat("\nSaved files:\n")
    for (f in saved_files) {
      cat(sprintf("  - %s\n", basename(f)))
    }
  }
  
  cat("\nThe Shiny app should now be able to load the updated UMAP data.\n")
  
}, error = function(e) {
  cat(sprintf("\nERROR: %s\n", e$message))
  quit(status = 1)
})