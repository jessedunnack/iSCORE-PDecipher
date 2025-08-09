#!/usr/bin/env Rscript

# DETAILED ANALYSIS OF "UNKNOWN" CLUSTERS
# Dig deeper into marker patterns

library(dplyr)

cat("=================================================================\n")
cat("DETAILED ANALYSIS OF UNKNOWN CLUSTERS\n")
cat("=================================================================\n\n")

# Load markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")

# Extended marker sets for investigation
SPECIAL_MARKERS <- list(
  # Neuroendocrine
  NEUROENDOCRINE = c("CGA", "CHGA", "CHGB", "CALCA", "CALCB", "SST", "VIP", "TPH1", "TPH2"),
  
  # Hypothalamic
  HYPOTHALAMIC = c("HCRT", "OXT", "AVP", "CRH", "TRH", "GHRH", "POMC", "NPY", "AGRP"),
  
  # Mesenchymal/Fibroblast
  MESENCHYMAL = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PRRX1", "PRRX2", "TWIST1", "TWIST2"),
  
  # Epithelial
  EPITHELIAL = c("EPCAM", "CDH1", "KRT8", "KRT18", "KRT19", "CLDN3", "CLDN4", "CLDN7"),
  
  # Immune/Microglia
  IMMUNE = c("PTPRC", "CD74", "HLA-DRA", "HLA-DRB1", "AIF1", "CX3CR1", "TMEM119", "P2RY12"),
  
  # Vascular
  VASCULAR = c("PECAM1", "CDH5", "VWF", "FLT1", "KDR", "TEK", "PDGFRB", "RGS5"),
  
  # Retinal pigment epithelium
  RPE = c("RPE65", "BEST1", "RLBP1", "RDH5", "RDH10", "TTR", "RBP1"),
  
  # Stress/Response
  STRESS_RESPONSE = c("FOS", "JUN", "ATF3", "ATF4", "DDIT3", "GADD45A", "CDKN1A", "GDF15"),
  
  # Extracellular matrix
  ECM = c("DCN", "LUM", "VCAN", "BGN", "FMOD", "COL1A1", "COL3A1", "FN1", "LAMB1"),
  
  # Progenitor subtypes
  RADIAL_GLIA = c("VIM", "FABP7", "BLBP", "SOX2", "PAX6", "HES1", "HES5", "HOPX"),
  INTERMEDIATE_PROG = c("EOMES", "NEUROG2", "NEUROD1", "NEUROD2", "NEUROD4", "NHLH1"),
  
  # Secretory/Transport
  SECRETORY = c("FOLR1", "CLIC6", "AQP1", "SLC4A10", "SLC12A2", "KCNJ13", "ATP1A1", "ATP1B1")
)

# Analyze each "Unknown" cluster
unknown_clusters <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14)

for (cl in unknown_clusters) {
  cl_name <- paste0("cluster_", cl)
  
  # Get markers
  cl_markers <- coarse_markers %>%
    filter(cluster == cl_name) %>%
    arrange(desc(avg_log2FC))
  
  cat(sprintf("\n=== CLUSTER %d ===\n", cl))
  cat(sprintf("Top marker: %s (%.2f FC, %.1f%% cells)\n", 
              cl_markers$gene[1], cl_markers$avg_log2FC[1], cl_markers$pct.1[1] * 100))
  
  # Check top 100 markers against special sets
  top100 <- head(cl_markers, 100)
  
  cat("\nSpecialized marker analysis (in top 100):\n")
  for (set_name in names(SPECIAL_MARKERS)) {
    found <- intersect(top100$gene, SPECIAL_MARKERS[[set_name]])
    if (length(found) > 0) {
      # Get their ranks
      ranks <- which(top100$gene %in% found)
      cat(sprintf("  %s: %d markers found\n", set_name, length(found)))
      for (i in 1:min(5, length(found))) {
        gene <- found[i]
        rank <- which(cl_markers$gene == gene)
        gene_data <- cl_markers[cl_markers$gene == gene, ]
        cat(sprintf("    - %s (rank %d): %.2f FC, %.1f%% cells\n", 
                    gene, rank, gene_data$avg_log2FC, gene_data$pct.1 * 100))
      }
    }
  }
  
  # Check for TH specifically (might not be in top 50)
  th_data <- cl_markers[cl_markers$gene == "TH", ]
  if (nrow(th_data) > 0) {
    cat(sprintf("\n  TH found at rank %d: %.2f FC, %.1f%% cells\n", 
                which(cl_markers$gene == "TH"), 
                th_data$avg_log2FC, th_data$pct.1 * 100))
  }
  
  # Suggest identity based on findings
  cat("\nPOSSIBLE IDENTITY: ")
  
  # Decision logic
  if (cl == 7 && "TTR" %in% head(cl_markers$gene, 5)) {
    cat("Choroid plexus epithelium\n")
  } else if (cl == 12 && "CGA" %in% head(cl_markers$gene, 5)) {
    cat("Neuroendocrine cells (CGA+, CALCA+, TPH1+)\n")
  } else if (cl == 14 && cl_markers$gene[1] == "HCRT") {
    cat("Hypothalamic neurons (Hypocretin/Orexin+)\n")
  } else if (cl == 3 && any(c("COL1A1", "PRRX1", "PRRX2") %in% head(cl_markers$gene, 10))) {
    cat("Mesenchymal/Fibroblast-like cells\n")
  } else if (cl == 8 && "DCN" %in% head(cl_markers$gene, 3)) {
    cat("Fibroblasts or ECM-producing cells\n")
  } else if (cl == 11 && cl_markers$gene[1] == "CRABP1") {
    cat("CRABP1+ progenitor or transitioning cells\n")
  } else if (cl == 6 && any(c("GDF15", "ATF3", "DDIT3") %in% head(cl_markers$gene, 10))) {
    cat("Stressed/damaged cells (stress response program)\n")
  } else {
    cat("Requires further investigation\n")
  }
}

cat("\n\n=== SUMMARY ===\n")
cat("Most clusters lack clear neuronal identity markers\n")
cat("Suggests many cells are:\n")
cat("- Progenitors at various stages\n")
cat("- Non-neuronal support cells\n")
cat("- Stressed or transitioning cells\n")
cat("- Specialized cell types (choroid plexus, neuroendocrine)\n")