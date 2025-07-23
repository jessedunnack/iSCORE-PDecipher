#!/usr/bin/env Rscript

# Script to add gene expression data to UMAP SingleCellExperiment objects
# This enables gene expression visualization without loading the full Seurat object

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(Matrix)
  library(dplyr)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-s", "--seurat"), type="character", 
              help="Path to Seurat RDS file"),
  make_option(c("-u", "--umap"), type="character", 
              help="Path to UMAP SCE RDS file"),
  make_option(c("-o", "--output"), type="character", 
              help="Output path for updated SCE file"),
  make_option(c("-g", "--genes"), type="character", default=NULL,
              help="Comma-separated list of genes to include (optional)"),
  make_option(c("-m", "--markers"), type="character", default=NULL,
              help="Path to markers RDS file to include marker genes"),
  make_option(c("-t", "--top-markers"), type="integer", default=50,
              help="Number of top markers per cluster to include (default: 50)"),
  make_option(c("-p", "--pd-genes"), type="logical", default=TRUE,
              help="Include PD-relevant genes (default: TRUE)"),
  make_option(c("-n", "--normalize"), type="logical", default=TRUE,
              help="Use normalized data (default: TRUE)")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Add gene expression to UMAP SCE objects")
opt <- parse_args(opt_parser)

# Validate inputs
if (is.null(opt$seurat) || is.null(opt$umap) || is.null(opt$output)) {
  stop("Required arguments: --seurat, --umap, --output")
}

# PD-relevant genes to always include
PD_GENES <- c(
  # Core PD genes
  "LRRK2", "SNCA", "PINK1", "PRKN", "PARK2", "PARK7", "DJ1",
  "GBA", "ATP13A2", "VPS35", "FBXO7", "DNAJC6", "SYNJ1", "VPS13C",
  
  # Neuronal markers
  "TH", "SLC18A2", "DDC", "SLC6A3", "MAP2", "RBFOX3", "SYN1", 
  "TUBB3", "NCAM1", "GAP43", "STMN2",
  
  # Glial markers
  "GFAP", "AQP4", "S100B", "ALDH1L1", "CD68", "AIF1", "CX3CR1",
  "TMEM119", "OLIG2", "MBP", "PLP1",
  
  # Mitochondrial/stress
  "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ATP6",
  "HSP90AA1", "HSPA5", "HSPA8", "HSPB1",
  
  # Cell cycle
  "MKI67", "TOP2A", "PCNA", "MCM2", "CDK1"
)

cat("=== Adding Expression Data to UMAP SCE ===\n")
cat("Seurat file:", opt$seurat, "\n")
cat("UMAP file:", opt$umap, "\n")
cat("Output file:", opt$output, "\n\n")

# Load UMAP SCE
cat("Loading UMAP SCE object...\n")
sce <- readRDS(opt$umap)
n_cells <- ncol(sce)
cat(sprintf("  - Loaded SCE with %d cells\n", n_cells))

# Collect genes to include
genes_to_include <- character()

# 1. Add user-specified genes
if (!is.null(opt$genes)) {
  user_genes <- trimws(strsplit(opt$genes, ",")[[1]])
  genes_to_include <- c(genes_to_include, user_genes)
  cat(sprintf("  - User specified %d genes\n", length(user_genes)))
}

# 2. Add PD-relevant genes
if (opt$`pd-genes`) {
  genes_to_include <- c(genes_to_include, PD_GENES)
  cat(sprintf("  - Added %d PD-relevant genes\n", length(PD_GENES)))
}

# 3. Add top markers from each cluster
if (!is.null(opt$markers) && file.exists(opt$markers)) {
  cat("Loading marker genes...\n")
  markers <- readRDS(opt$markers)
  
  if (!is.null(markers) && nrow(markers) > 0) {
    top_markers <- markers %>%
      group_by(cluster) %>%
      arrange(desc(avg_log2FC)) %>%
      slice_head(n = opt$`top-markers`) %>%
      pull(gene) %>%
      unique()
    
    genes_to_include <- c(genes_to_include, top_markers)
    cat(sprintf("  - Added %d unique marker genes\n", length(top_markers)))
  }
}

# Remove duplicates
genes_to_include <- unique(genes_to_include)
cat(sprintf("\nTotal unique genes to include: %d\n", length(genes_to_include)))

# Load Seurat object (only the assay data)
cat("\nLoading expression data from Seurat object...\n")
cat("This may take a minute for large objects...\n")

seurat_obj <- readRDS(opt$seurat)

# Get available genes
available_genes <- rownames(seurat_obj)
genes_found <- intersect(genes_to_include, available_genes)
genes_missing <- setdiff(genes_to_include, available_genes)

cat(sprintf("  - Found %d/%d requested genes in Seurat object\n", 
            length(genes_found), length(genes_to_include)))

if (length(genes_missing) > 0 && length(genes_missing) < 20) {
  cat("  - Missing genes:", paste(head(genes_missing, 20), collapse=", "), "\n")
}

# Extract expression data
if (opt$normalize) {
  cat("Extracting normalized expression data...\n")
  expr_data <- GetAssayData(seurat_obj, slot = "data")[genes_found, , drop = FALSE]
} else {
  cat("Extracting raw count data...\n")
  expr_data <- GetAssayData(seurat_obj, slot = "counts")[genes_found, , drop = FALSE]
}

# Match cell order
cell_ids <- colnames(sce)
common_cells <- intersect(cell_ids, colnames(expr_data))

if (length(common_cells) != n_cells) {
  warning(sprintf("Cell mismatch: SCE has %d cells, found %d in Seurat", 
                  n_cells, length(common_cells)))
}

# Subset and reorder expression data to match SCE
expr_data <- expr_data[, cell_ids, drop = FALSE]

# Convert to sparse matrix if not already
if (!is(expr_data, "sparseMatrix")) {
  expr_data <- as(expr_data, "sparseMatrix")
}

cat(sprintf("\nExpression matrix: %d genes x %d cells\n", 
            nrow(expr_data), ncol(expr_data)))
cat(sprintf("Memory usage: %.1f MB\n", 
            object.size(expr_data) / 1024^2))

# Add expression data to SCE
assay(sce, "logcounts") <- expr_data

# Add metadata about included genes
metadata(sce)$expression_genes <- list(
  total = length(genes_found),
  pd_genes = intersect(genes_found, PD_GENES),
  marker_genes = if (!is.null(opt$markers)) intersect(genes_found, top_markers) else NULL,
  user_genes = if (!is.null(opt$genes)) intersect(genes_found, user_genes) else NULL
)

# Save updated SCE
cat("\nSaving updated SCE object...\n")
saveRDS(sce, opt$output)

# Summary
cat("\n=== Summary ===\n")
cat(sprintf("Original SCE size: %.1f MB\n", file.size(opt$umap) / 1024^2))
cat(sprintf("Updated SCE size: %.1f MB\n", file.size(opt$output) / 1024^2))
cat(sprintf("Added expression for %d genes\n", length(genes_found)))
cat("\nDone!\n")