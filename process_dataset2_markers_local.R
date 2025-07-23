#!/usr/bin/env Rscript

# Standalone script to process UMAP and markers for dataset 2
# Run this from: E:/scRNSeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher
# It will find the Seurat object in: ../iSCORE-PD_plus_CRISPRi

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
  make_option(c("-o", "--output-dir"), type="character", 
              default="E:/scRNSeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data",
              help="Output directory for results"),
  
  # FindAllMarkers parameters
  make_option("--cluster-markers", type="logical", default=TRUE,
              help="Calculate cluster markers using FindAllMarkers (default: TRUE)"),
  make_option("--cluster-test", type="character", default="MAST",
              help="Test to use for FindAllMarkers: MAST, wilcox, bimod, roc, t, negbinom, poisson, LR, DESeq2 (default: MAST)"),
  make_option("--cluster-logfc", type="numeric", default=0.5,
              help="Log fold change threshold for FindAllMarkers (default: 0.5)"),
  make_option("--cluster-min-pct", type="numeric", default=0.25,
              help="Min percent cells expressing for FindAllMarkers (default: 0.25)"),
  make_option("--cluster-min-diff-pct", type="numeric", default=0.2,
              help="Min difference in percent expressed for FindAllMarkers (default: 0.2)"),
  make_option("--cluster-only-pos", type="logical", default=FALSE,
              help="Only return positive markers in FindAllMarkers (default: FALSE)"),
  make_option("--cluster-max-cells", type="integer", default=500,
              help="Max cells per ident for FindAllMarkers, use -1 for Inf (default: 500)"),
  
  # FindMarkers parameters for group comparisons
  make_option("--group-markers", type="logical", default=TRUE,
              help="Calculate group comparison markers using FindMarkers (default: TRUE)"),
  make_option("--group-test", type="character", default="MAST",
              help="Test to use for FindMarkers: MAST, wilcox, bimod, roc, t, negbinom, poisson, LR, DESeq2 (default: MAST)"),
  make_option("--group-logfc", type="numeric", default=0.25,
              help="Log fold change threshold for FindMarkers (default: 0.25)"),
  make_option("--group-min-pct", type="numeric", default=0.1,
              help="Min percent cells expressing for FindMarkers (default: 0.1)"),
  make_option("--group-min-diff-pct", type="numeric", default=-Inf,
              help="Min difference in percent expressed for FindMarkers (default: -Inf)"),
  make_option("--group-only-pos", type="logical", default=FALSE,
              help="Only return positive markers in FindMarkers (default: FALSE)"),
  make_option("--group-max-cells", type="integer", default=500,
              help="Max cells per ident for FindMarkers, use -1 for Inf (default: 500)"),
  
  # Other options
  make_option("--force", type="logical", default=FALSE,
              help="Force overwrite existing files (default: FALSE)"),
  make_option("--verbose", type="logical", default=TRUE,
              help="Print detailed progress (default: TRUE)")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Process markers for iSCORE_PD_plus_CRISPRi with custom parameters")
opt <- parse_args(opt_parser)

# Configuration
DATASET_NAME <- "iSCORE_PD_CRISPRi"
SEURAT_FILE <- "../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"
OUTPUT_DIR <- opt$`output-dir`

cat("\n=== UMAP & Marker Processing for Dataset 2 ===\n")
cat("Working directory:", getwd(), "\n")
cat("Seurat file:", SEURAT_FILE, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

# Show parameter summary
cat("FindAllMarkers parameters:\n")
cat(sprintf("  Test: %s\n", opt$`cluster-test`))
cat(sprintf("  Log FC threshold: %.2f\n", opt$`cluster-logfc`))
cat(sprintf("  Min percent: %.2f\n", opt$`cluster-min-pct`))
cat(sprintf("  Min diff percent: %.2f\n", opt$`cluster-min-diff-pct`))
cat(sprintf("  Only positive: %s\n", opt$`cluster-only-pos`))
cat(sprintf("  Max cells: %s\n", ifelse(opt$`cluster-max-cells` == -1, "all", opt$`cluster-max-cells`)))

cat("\nFindMarkers parameters:\n")
cat(sprintf("  Test: %s\n", opt$`group-test`))
cat(sprintf("  Log FC threshold: %.2f\n", opt$`group-logfc`))
cat(sprintf("  Min percent: %.2f\n", opt$`group-min-pct`))
cat(sprintf("  Min diff percent: %s\n", 
            ifelse(is.infinite(opt$`group-min-diff-pct`), "-Inf", sprintf("%.2f", opt$`group-min-diff-pct`))))
cat(sprintf("  Only positive: %s\n", opt$`group-only-pos`))
cat(sprintf("  Max cells: %s\n", ifelse(opt$`group-max-cells` == -1, "all", opt$`group-max-cells`)))
cat("\n")

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
  cat("Created output directory:", OUTPUT_DIR, "\n")
}

# Function to check file overwrite
check_overwrite <- function(filepath, force = FALSE) {
  if (file.exists(filepath)) {
    if (!force) {
      cat(sprintf("\nFile exists: %s\n", basename(filepath)))
      response <- readline("Overwrite? (y/n): ")
      return(tolower(response) == "y")
    }
    return(TRUE)
  }
  return(TRUE)
}

# Function to extract UMAP data
extract_umap_data <- function(seurat_obj, dataset_name) {
  cat("\nExtracting UMAP data...\n")
  
  # Get UMAP coordinates
  umap_coords <- Embeddings(seurat_obj, "umap.cca")
  cat(sprintf("  - Extracted UMAP coordinates: %d cells x 2 dimensions\n", nrow(umap_coords)))
  
  # Extract metadata
  metadata_cols <- c("seurat_clusters", "mutation_tidy", "batch", 
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
  
  metadata(sce) <- list(
    dataset_name = dataset_name,
    n_cells = ncol(sce),
    n_clusters = length(unique(cell_metadata$seurat_clusters)),
    extraction_date = Sys.Date()
  )
  
  return(sce)
}

# Function to calculate cluster markers with custom parameters
calculate_cluster_markers_custom <- function(seurat_obj, opt) {
  cat("\n=== Calculating Cluster Markers with FindAllMarkers ===\n")
  cat(sprintf("Parameters: test=%s, logfc=%.2f, min.pct=%.2f, min.diff.pct=%.2f\n",
              opt$`cluster-test`, opt$`cluster-logfc`, opt$`cluster-min-pct`, opt$`cluster-min-diff-pct`))
  
  Idents(seurat_obj) <- "seurat_clusters"
  n_clusters <- length(unique(Idents(seurat_obj)))
  
  # Progress tracking
  pb <- txtProgressBar(min = 0, max = 1, style = 3)
  
  markers <- tryCatch({
    result <- FindAllMarkers(
      seurat_obj,
      test.use = opt$`cluster-test`,
      logfc.threshold = opt$`cluster-logfc`,
      min.pct = opt$`cluster-min-pct`,
      min.diff.pct = opt$`cluster-min-diff-pct`,
      only.pos = opt$`cluster-only-pos`,
      max.cells.per.ident = ifelse(opt$`cluster-max-cells` == -1, Inf, opt$`cluster-max-cells`),
      verbose = FALSE
    )
    setTxtProgressBar(pb, 1)
    result
  }, error = function(e) {
    setTxtProgressBar(pb, 1)
    warning(sprintf("Failed to find cluster markers: %s", e$message))
    return(NULL)
  })
  
  close(pb)
  
  if (!is.null(markers)) {
    markers$cluster <- factor(markers$cluster)
    cat(sprintf("\nFound %d marker genes across %d clusters\n", nrow(markers), n_clusters))
    
    # Show top markers per cluster
    cat("\nTop 5 markers per cluster:\n")
    for (i in 0:(n_clusters-1)) {
      cluster_markers <- markers %>%
        filter(cluster == as.character(i)) %>%
        arrange(desc(avg_log2FC)) %>%
        head(5)
      
      if (nrow(cluster_markers) > 0) {
        cat(sprintf("Cluster %d: %s\n", i, 
                   paste(cluster_markers$gene, collapse=", ")))
      }
    }
  }
  
  return(markers)
}

# Function to calculate group markers with custom parameters
calculate_group_markers_custom <- function(seurat_obj, opt) {
  cat("\n=== Calculating Group Comparison Markers with FindMarkers ===\n")
  cat(sprintf("Parameters: test=%s, logfc=%.2f, min.pct=%.2f\n",
              opt$`group-test`, opt$`group-logfc`, opt$`group-min-pct`))
  
  group_markers <- list()
  
  # 1. iSCORE-PD: eWT vs all mutants
  if ("mutation_tidy" %in% colnames(seurat_obj@meta.data)) {
    cat("\n1. iSCORE-PD comparison: eWT vs all mutants\n")
    
    seurat_obj@meta.data$iscore_group <- ifelse(
      seurat_obj@meta.data$mutation_tidy == "eWT", "eWT", "mutant"
    )
    
    ewt_cells <- sum(seurat_obj@meta.data$iscore_group == "eWT")
    mut_cells <- sum(seurat_obj@meta.data$iscore_group == "mutant")
    cat(sprintf("   - eWT cells: %d\n", ewt_cells))
    cat(sprintf("   - Mutant cells: %d\n", mut_cells))
    
    if (ewt_cells > 10 && mut_cells > 10) {
      Idents(seurat_obj) <- "iscore_group"
      
      markers <- tryCatch({
        FindMarkers(seurat_obj,
                   ident.1 = "mutant",
                   ident.2 = "eWT",
                   test.use = opt$`group-test`,
                   logfc.threshold = opt$`group-logfc`,
                   min.pct = opt$`group-min-pct`,
                   min.diff.pct = opt$`group-min-diff-pct`,
                   only.pos = opt$`group-only-pos`,
                   max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                   verbose = FALSE)
      }, error = function(e) {
        warning(sprintf("Failed: %s", e$message))
        return(NULL)
      })
      
      if (!is.null(markers)) {
        markers$gene <- rownames(markers)
        markers$comparison <- "iSCORE_mutants_vs_eWT"
        group_markers[["iSCORE_mutants_vs_eWT"]] <- markers
        cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
      }
    }
  }
  
  # 2. CRISPRi: Non-targeting vs perturbed
  if ("scMAGeCK_gene_assignment" %in% colnames(seurat_obj@meta.data)) {
    cat("\n2. CRISPRi comparison: Non-targeting vs perturbed cells\n")
    
    seurat_obj@meta.data$crispr_group <- ifelse(
      seurat_obj@meta.data$scMAGeCK_gene_assignment == "Non-Targeting", 
      "control", "perturbed"
    )
    
    ctrl_cells <- sum(seurat_obj@meta.data$crispr_group == "control")
    pert_cells <- sum(seurat_obj@meta.data$crispr_group == "perturbed")
    cat(sprintf("   - Control cells: %d\n", ctrl_cells))
    cat(sprintf("   - Perturbed cells: %d\n", pert_cells))
    
    if (ctrl_cells > 10 && pert_cells > 10) {
      Idents(seurat_obj) <- "crispr_group"
      
      markers <- tryCatch({
        FindMarkers(seurat_obj,
                   ident.1 = "perturbed",
                   ident.2 = "control",
                   test.use = opt$`group-test`,
                   logfc.threshold = opt$`group-logfc`,
                   min.pct = opt$`group-min-pct`,
                   min.diff.pct = opt$`group-min-diff-pct`,
                   only.pos = opt$`group-only-pos`,
                   max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                   verbose = FALSE)
      }, error = function(e) {
        warning(sprintf("Failed: %s", e$message))
        return(NULL)
      })
      
      if (!is.null(markers)) {
        markers$gene <- rownames(markers)
        markers$comparison <- "CRISPRi_perturbed_vs_control"
        group_markers[["CRISPRi_perturbed_vs_control"]] <- markers
        cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
      }
    }
  }
  
  # 3. Cross-method comparison
  if ("mutation_tidy" %in% colnames(seurat_obj@meta.data) && 
      "scMAGeCK_gene_assignment" %in% colnames(seurat_obj@meta.data)) {
    cat("\n3. Cross-method comparison: iSCORE-PD vs PerturbSeq cells\n")
    
    seurat_obj@meta.data$method_group <- ifelse(
      !is.na(seurat_obj@meta.data$mutation_tidy), "iSCORE",
      ifelse(!is.na(seurat_obj@meta.data$scMAGeCK_gene_assignment), "PerturbSeq", NA)
    )
    
    iscore_cells <- sum(seurat_obj@meta.data$method_group == "iSCORE", na.rm = TRUE)
    perturb_cells <- sum(seurat_obj@meta.data$method_group == "PerturbSeq", na.rm = TRUE)
    cat(sprintf("   - iSCORE cells: %d\n", iscore_cells))
    cat(sprintf("   - PerturbSeq cells: %d\n", perturb_cells))
    
    if (iscore_cells > 10 && perturb_cells > 10) {
      valid_cells <- !is.na(seurat_obj@meta.data$method_group)
      seurat_subset <- seurat_obj[, valid_cells]
      Idents(seurat_subset) <- "method_group"
      
      markers <- tryCatch({
        FindMarkers(seurat_subset,
                   ident.1 = "PerturbSeq",
                   ident.2 = "iSCORE",
                   test.use = opt$`group-test`,
                   logfc.threshold = opt$`group-logfc`,
                   min.pct = opt$`group-min-pct`,
                   min.diff.pct = opt$`group-min-diff-pct`,
                   only.pos = opt$`group-only-pos`,
                   max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                   verbose = FALSE)
      }, error = function(e) {
        warning(sprintf("Failed: %s", e$message))
        return(NULL)
      })
      
      if (!is.null(markers)) {
        markers$gene <- rownames(markers)
        markers$comparison <- "PerturbSeq_vs_iSCORE"
        group_markers[["PerturbSeq_vs_iSCORE"]] <- markers
        cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
      }
    }
  }
  
  if (length(group_markers) > 0) {
    return(do.call(rbind, group_markers))
  } else {
    return(NULL)
  }
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
  cat(sprintf("Current clustering: %d clusters\n", n_clusters))
  
  # Extract UMAP data
  sce <- extract_umap_data(seurat_obj, DATASET_NAME)
  
  # Save UMAP data
  umap_file <- file.path(OUTPUT_DIR, sprintf("%s_umap_data_%dpc.rds", DATASET_NAME, opt$pc))
  if (check_overwrite(umap_file, opt$force)) {
    saveRDS(sce, umap_file)
    cat(sprintf("\nSaved UMAP data to: %s\n", basename(umap_file)))
  }
  
  # Calculate cluster markers if requested
  if (opt$`cluster-markers`) {
    cluster_markers <- calculate_cluster_markers_custom(seurat_obj, opt)
    
    if (!is.null(cluster_markers)) {
      markers_file <- file.path(OUTPUT_DIR, sprintf("%s_cluster_markers.rds", DATASET_NAME))
      if (check_overwrite(markers_file, opt$force)) {
        saveRDS(cluster_markers, markers_file)
        cat(sprintf("\nSaved cluster markers to: %s\n", basename(markers_file)))
      }
    }
  }
  
  # Calculate group markers if requested
  if (opt$`group-markers`) {
    group_markers <- calculate_group_markers_custom(seurat_obj, opt)
    
    if (!is.null(group_markers)) {
      group_file <- file.path(OUTPUT_DIR, sprintf("%s_group_comparison_markers.rds", DATASET_NAME))
      if (check_overwrite(group_file, opt$force)) {
        saveRDS(group_markers, group_file)
        cat(sprintf("\nSaved group markers to: %s\n", basename(group_file)))
      }
    }
  }
  
  cat("\n=== PROCESSING COMPLETE ===\n")
  cat("Output files saved to:", OUTPUT_DIR, "\n")
  
}, error = function(e) {
  cat(sprintf("\nERROR: %s\n", e$message))
  quit(status = 1)
})