#!/usr/bin/env Rscript

# Extract expression for ALL cluster markers plus key genes
# Run this in seuratv4 environment

library(Seurat)
library(Matrix)
library(dplyr)

# Key genes to always include
KEY_GENES <- c(
  # PD genes
  "LRRK2", "SNCA", "PINK1", "PRKN", "PARK2", "PARK7", "DJ1",
  "GBA", "ATP13A2", "VPS35", "FBXO7", "DNAJC6", "SYNJ1", "VPS13C",
  
  # Neuronal markers
  "TH", "SLC18A2", "DDC", "SLC6A3", "MAP2", "RBFOX3", "SYN1", 
  "TUBB3", "NCAM1", "GAP43", "STMN2", "ENO2", "SYP",
  
  # Glial markers
  "GFAP", "AQP4", "S100B", "ALDH1L1", "CD68", "AIF1", "CX3CR1",
  "TMEM119", "OLIG2", "MBP", "PLP1", "SOX10",
  
  # Mitochondrial/stress
  "MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ATP6",
  "HSP90AA1", "HSPA5", "HSPA8", "HSPB1",
  
  # Cell cycle
  "MKI67", "TOP2A", "PCNA", "MCM2", "CDK1",
  
  # Other important markers
  "APOE", "CLU", "C1QA", "C1QB", "TREM2", "TYROBP"
)

# Process each dataset
datasets <- list(
  "iSCORE_PD_CRISPRi" = list(
    seurat = "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
    markers = "inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds"
  )
)

for (dataset_name in names(datasets)) {
  cat("\n=== Processing", dataset_name, "===\n")
  
  paths <- datasets[[dataset_name]]
  
  if (!file.exists(paths$seurat)) {
    cat("Seurat file not found:", paths$seurat, "\n")
    next
  }
  
  # Load marker genes if available
  all_genes <- KEY_GENES
  
  if (file.exists(paths$markers)) {
    cat("Loading cluster markers...\n")
    markers <- readRDS(paths$markers)
    
    # Get all unique marker genes
    marker_genes <- unique(markers$gene)
    cat(sprintf("Found %d unique marker genes\n", length(marker_genes)))
    
    # Combine with key genes
    all_genes <- unique(c(all_genes, marker_genes))
    cat(sprintf("Total genes to extract: %d\n", length(all_genes)))
  }
  
  # Load Seurat object
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(paths$seurat)
  
  # Get available genes
  available_genes <- intersect(all_genes, rownames(seurat_obj))
  cat(sprintf("Found %d/%d requested genes in dataset\n", 
              length(available_genes), length(all_genes)))
  
  # Show some statistics
  cat("\nGene statistics:\n")
  cat(sprintf("  - Key PD/marker genes: %d\n", 
              length(intersect(KEY_GENES, available_genes))))
  cat(sprintf("  - Cluster marker genes: %d\n", 
              length(setdiff(available_genes, KEY_GENES))))
  
  # Extract normalized expression
  cat("\nExtracting expression data...\n")
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  
  # Convert to sparse matrix for efficiency
  if (!is(expr_matrix, "sparseMatrix")) {
    expr_matrix <- as(expr_matrix, "sparseMatrix")
  }
  
  # Check sparsity
  sparsity <- 1 - (nnzero(expr_matrix) / (nrow(expr_matrix) * ncol(expr_matrix)))
  cat(sprintf("Matrix sparsity: %.1f%%\n", sparsity * 100))
  
  # Save expression matrix
  output_file <- sprintf("inst/extdata/umap_data/%s_all_markers_expr.rds", dataset_name)
  cat(sprintf("\nSaving expression matrix to %s...\n", output_file))
  saveRDS(expr_matrix, output_file)
  
  cat(sprintf("Matrix size: %d genes x %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))
  cat(sprintf("File size: %.1f MB\n", file.size(output_file) / 1024^2))
  
  # Create gene info file
  gene_info <- data.frame(
    gene = rownames(expr_matrix),
    available = TRUE,
    is_key_gene = rownames(expr_matrix) %in% KEY_GENES,
    is_marker = rownames(expr_matrix) %in% marker_genes,
    stringsAsFactors = FALSE
  )
  
  # Save gene info
  info_file <- sprintf("inst/extdata/umap_data/%s_all_genes_info.rds", dataset_name)
  saveRDS(gene_info, info_file)
  cat(sprintf("Saved gene info to %s\n", info_file))
  
  # Clean up memory
  rm(seurat_obj)
  gc()
}

cat("\n=== All processing complete ===\n")