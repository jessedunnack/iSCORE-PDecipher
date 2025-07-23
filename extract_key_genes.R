#!/usr/bin/env Rscript

# Extract expression for key genes and save as lightweight matrix
# Run this in seuratv4 environment

library(Seurat)
library(Matrix)

# Key genes to extract
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

# Process iSCORE_PD_CRISPRi dataset
cat("\nProcessing iSCORE_PD_CRISPRi dataset...\n")
seurat_file <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"

if (file.exists(seurat_file)) {
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_file)
  
  # Get available genes
  available_genes <- intersect(KEY_GENES, rownames(seurat_obj))
  cat(sprintf("Found %d/%d key genes in dataset\n", length(available_genes), length(KEY_GENES)))
  
  # Extract normalized expression
  cat("Extracting expression data...\n")
  expr_matrix <- GetAssayData(seurat_obj, slot = "data")[available_genes, , drop = FALSE]
  
  # Convert to sparse matrix for efficiency
  if (!is(expr_matrix, "sparseMatrix")) {
    expr_matrix <- as(expr_matrix, "sparseMatrix")
  }
  
  # Save as lightweight RDS
  output_file <- "inst/extdata/umap_data/iSCORE_PD_CRISPRi_key_genes_expr.rds"
  cat(sprintf("Saving to %s...\n", output_file))
  saveRDS(expr_matrix, output_file)
  
  cat(sprintf("Matrix size: %d genes x %d cells\n", nrow(expr_matrix), ncol(expr_matrix)))
  cat(sprintf("File size: %.1f MB\n", file.size(output_file) / 1024^2))
  
  # Also save gene list
  gene_info <- data.frame(
    gene = rownames(expr_matrix),
    available = TRUE,
    stringsAsFactors = FALSE
  )
  
  # Add missing genes
  missing_genes <- setdiff(KEY_GENES, available_genes)
  if (length(missing_genes) > 0) {
    gene_info <- rbind(
      gene_info,
      data.frame(gene = missing_genes, available = FALSE, stringsAsFactors = FALSE)
    )
  }
  
  saveRDS(gene_info, "inst/extdata/umap_data/iSCORE_PD_CRISPRi_gene_info.rds")
  
  cat("\nMissing genes:\n")
  cat(paste(missing_genes, collapse = ", "), "\n")
  
} else {
  cat("Seurat file not found!\n")
}

cat("\nDone!\n")