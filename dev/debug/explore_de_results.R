# Script to explore the DE results RDS file structure
# For iSCORE-PDecipher data exploration

# Load required libraries
library(dplyr)

# Set path to your DE results file
# Update this path to match your local setup
de_path <- "path/to/your/dataset/full_DE_results.rds"

# Load the DE results
de_results <- readRDS(de_path)

# Explore the structure
cat("=== DE Results Structure ===\n")
cat("Top level structure:\n")
str(de_results, max.level = 1)

# The DE results contain three main components:
# 1. iSCORE_PD_MAST - Mutations from isogenic lines
# 2. CRISPRi_Mixscale - Knockdown experiments  
# 3. CRISPRa_Mixscale - Activation experiments (if present)

# Example: Explore MAST results
if ("iSCORE_PD_MAST" %in% names(de_results)) {
  cat("\n=== MAST Results ===\n")
  mast_genes <- names(de_results$iSCORE_PD_MAST)
  cat("Genes analyzed:", paste(head(mast_genes, 5), collapse = ", "), "...\n")
  cat("Total genes:", length(mast_genes), "\n")
  
  # Look at one gene's structure
  example_gene <- mast_genes[1]
  cat("\nExample gene:", example_gene, "\n")
  clusters <- names(de_results$iSCORE_PD_MAST[[example_gene]])
  cat("Clusters available:", paste(clusters, collapse = ", "), "\n")
  
  # Examine results for one gene/cluster
  if (length(clusters) > 0) {
    example_results <- de_results$iSCORE_PD_MAST[[example_gene]][[clusters[1]]]$results
    cat("\nDE results structure (first few rows):\n")
    print(head(example_results, 3))
    cat("\nColumns available:", paste(colnames(example_results), collapse = ", "), "\n")
  }
}

# Example: Explore MixScale results
if ("CRISPRi_Mixscale" %in% names(de_results)) {
  cat("\n=== CRISPRi Results ===\n")
  crispi_genes <- names(de_results$CRISPRi_Mixscale)
  cat("Genes perturbed:", paste(head(crispi_genes, 5), collapse = ", "), "...\n")
  cat("Total genes:", length(crispi_genes), "\n")
  
  # Look at experiment columns
  if (length(crispi_genes) > 0) {
    example_gene <- crispi_genes[1]
    clusters <- names(de_results$CRISPRi_Mixscale[[example_gene]])
    if (length(clusters) > 0) {
      example_results <- de_results$CRISPRi_Mixscale[[example_gene]][[clusters[1]]]$results
      
      # Find log2FC columns (experiment-specific)
      log2fc_cols <- grep("^log2FC_", colnames(example_results), value = TRUE)
      cat("\nExperiments available:", gsub("log2FC_", "", log2fc_cols), "\n")
    }
  }
}

# Quick function to extract significant DE genes
get_significant_de_genes <- function(de_results, gene, cluster, method = "MAST", pval_cutoff = 0.05) {
  if (method == "MAST" && "iSCORE_PD_MAST" %in% names(de_results)) {
    results <- de_results$iSCORE_PD_MAST[[gene]][[cluster]]$results
    sig_genes <- results[results$p_val_adj < pval_cutoff, ]
    return(sig_genes[order(sig_genes$p_val_adj), ])
  }
  # Add other methods as needed
}

# Example usage
cat("\n=== Example: Get significant DE genes ===\n")
# sig_genes <- get_significant_de_genes(de_results, "LRRK2", "cluster_0")
# print(head(sig_genes))