#!/usr/bin/env Rscript

# RE-EVALUATE COARSE CLUSTER IDENTITIES USING ACTUAL MARKER FILES

library(dplyr)

cat("=================================================================\n")
cat("RE-CHECKING COARSE CLUSTER IDENTITIES WITH MARKER FILES\n")
cat("=================================================================\n\n")

# Load coarse markers
cat("1. Loading coarse cluster markers...\n")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")

# Load the identities from stress analysis
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")

cat("\nCurrent identity assignments:\n")
print(coarse_identities[, c("cluster", "identity", "n_cells")])

cat("\n\n2. Checking top markers for each cluster:\n")
cat("==========================================\n")

# Key markers to look for
DA_MARKERS <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", "ALDH1A1", "CALB1")
GLUT_MARKERS <- c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "TBR1", "SATB2")
GABA_MARKERS <- c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "SST", "PVALB", "VIP")
PROG_MARKERS <- c("SOX2", "NES", "PAX6", "HES1", "HES5", "TOP2A", "MKI67")
CHOROID_MARKERS <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13")

for (cl in 0:14) {
  cl_name <- paste0("cluster_", cl)
  
  if (cl_name %in% names(coarse_markers)) {
    markers <- coarse_markers[[cl_name]]
    
    cat(sprintf("\n--- CLUSTER %d (n=%d cells) ---\n", 
                cl, coarse_identities$n_cells[coarse_identities$cluster == cl]))
    cat(sprintf("Current assignment: %s\n", 
                coarse_identities$identity[coarse_identities$cluster == cl]))
    
    # Show top 20 markers
    top_markers <- head(markers, 20)
    cat("\nTop 20 markers:\n")
    print(top_markers[, c("gene", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")])
    
    # Check for key marker types
    cat("\nKey marker analysis:\n")
    
    # DA markers
    da_found <- intersect(top_markers$gene, DA_MARKERS)
    if (length(da_found) > 0) {
      cat(sprintf("  DA markers found: %s\n", paste(da_found, collapse=", ")))
    }
    
    # Glutamatergic
    glut_found <- intersect(top_markers$gene, GLUT_MARKERS)
    if (length(glut_found) > 0) {
      cat(sprintf("  Glutamatergic markers: %s\n", paste(glut_found, collapse=", ")))
    }
    
    # GABAergic
    gaba_found <- intersect(top_markers$gene, GABA_MARKERS)
    if (length(gaba_found) > 0) {
      cat(sprintf("  GABAergic markers: %s\n", paste(gaba_found, collapse=", ")))
    }
    
    # Progenitor
    prog_found <- intersect(top_markers$gene, PROG_MARKERS)
    if (length(prog_found) > 0) {
      cat(sprintf("  Progenitor markers: %s\n", paste(prog_found, collapse=", ")))
    }
    
    # Choroid
    choroid_found <- intersect(top_markers$gene, CHOROID_MARKERS)
    if (length(choroid_found) > 0) {
      cat(sprintf("  Choroid plexus markers: %s\n", paste(choroid_found, collapse=", ")))
    }
    
    # Check if TTR is in top 5 (potential choroid/stress marker)
    if ("TTR" %in% head(markers$gene, 5)) {
      cat("  *** TTR is in top 5 markers! ***\n")
    }
  }
}

cat("\n\n3. Summary of findings:\n")
cat("=======================\n")
cat("Based on marker analysis, clusters that may need re-assignment:\n")
cat("(Check console output above for detailed marker lists)\n")