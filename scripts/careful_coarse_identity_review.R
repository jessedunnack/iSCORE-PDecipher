#!/usr/bin/env Rscript

# CAREFUL RE-EVALUATION OF COARSE CLUSTER IDENTITIES
# Distinguishing neurons from progenitors

library(dplyr)

cat("=================================================================\n")
cat("CAREFUL COARSE CLUSTER IDENTITY REVIEW\n")
cat("=================================================================\n\n")

# Load coarse markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")

# Define marker categories more carefully
MATURE_NEURON_MARKERS <- c("MAP2", "RBFOX3", "SYN1", "SNAP25", "STMN2", "TUBB3", 
                          "NEFL", "NEFM", "NEFH", "ENO2", "SYT1")

PROGENITOR_MARKERS <- c("SOX2", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7", 
                       "BLBP", "SOX1", "NOTCH1", "ID1", "ID3", "ASCL1")

PROLIFERATION_MARKERS <- c("MKI67", "TOP2A", "PCNA", "CCND1", "CCND2", "CDC20", 
                          "UBE2C", "CENPF", "HIST1H3B")

# Specific neuronal subtypes
DA_SPECIFIC <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6")
GLUT_SPECIFIC <- c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B")
GABA_SPECIFIC <- c("GAD1", "GAD2", "SLC32A1")

# Non-neuronal
CHOROID_PLEXUS <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13", "ENPP2")
GLIA_MARKERS <- c("GFAP", "S100B", "ALDH1L1", "OLIG1", "OLIG2", "SOX10")

cat("Analyzing each cluster for neuronal maturity and identity...\n")
cat("==========================================================\n")

results <- data.frame(
  cluster = integer(),
  n_mature_neuron = integer(),
  n_progenitor = integer(),
  n_proliferation = integer(),
  has_DA = logical(),
  has_GLUT = logical(),
  has_GABA = logical(),
  has_choroid = logical(),
  suggested_identity = character(),
  stringsAsFactors = FALSE
)

for (cl in 0:14) {
  cl_name <- paste0("cluster_", cl)
  
  if (cl_name %in% names(coarse_markers)) {
    markers <- coarse_markers[[cl_name]]
    top50 <- head(markers, 50)  # Look at top 50 markers
    
    cat(sprintf("\n--- CLUSTER %d ---\n", cl))
    
    # Count marker types in top 50
    mature_count <- sum(top50$gene %in% MATURE_NEURON_MARKERS)
    prog_count <- sum(top50$gene %in% PROGENITOR_MARKERS)
    prolif_count <- sum(top50$gene %in% PROLIFERATION_MARKERS)
    
    # Check specific types
    has_da <- any(top50$gene %in% DA_SPECIFIC)
    has_glut <- any(top50$gene %in% GLUT_SPECIFIC)
    has_gaba <- any(top50$gene %in% GABA_SPECIFIC)
    has_choroid <- any(top50$gene %in% CHOROID_PLEXUS)
    
    # Print top 10 markers
    cat("Top 10 markers:\n")
    for (i in 1:min(10, nrow(markers))) {
      gene <- markers$gene[i]
      marker_type <- ""
      
      if (gene %in% MATURE_NEURON_MARKERS) marker_type <- "[MATURE]"
      if (gene %in% PROGENITOR_MARKERS) marker_type <- "[PROG]"
      if (gene %in% PROLIFERATION_MARKERS) marker_type <- "[PROLIF]"
      if (gene %in% DA_SPECIFIC) marker_type <- "[DA]"
      if (gene %in% GLUT_SPECIFIC) marker_type <- "[GLUT]"
      if (gene %in% GABA_SPECIFIC) marker_type <- "[GABA]"
      if (gene %in% CHOROID_PLEXUS) marker_type <- "[CHOROID]"
      
      cat(sprintf("  %2d. %-10s %.2f FC, %.1f%% cells %s\n", 
                  i, gene, markers$avg_log2FC[i], 
                  markers$pct.1[i] * 100, marker_type))
    }
    
    # Determine identity
    suggested_id <- "Unknown"
    
    # Strong choroid plexus signature?
    if (has_choroid && markers$gene[1] == "TTR") {
      suggested_id <- "Choroid plexus"
    }
    # Proliferating cells?
    else if (prolif_count >= 3 && markers$gene[1] %in% PROLIFERATION_MARKERS) {
      suggested_id <- "Proliferating cells"
    }
    # Progenitors?
    else if (prog_count > mature_count && prog_count >= 3) {
      suggested_id <- "Neural progenitors"
    }
    # Mature neurons - what type?
    else if (mature_count >= 2) {
      if (has_da && "TH" %in% head(markers$gene, 20)) {
        suggested_id <- "Dopaminergic neurons"
      } else if (has_glut) {
        suggested_id <- "Glutamatergic neurons"
      } else if (has_gaba) {
        suggested_id <- "GABAergic neurons"
      } else {
        suggested_id <- "Neurons (unspecified)"
      }
    }
    # Still progenitor-like?
    else if (prog_count >= 2) {
      suggested_id <- "Progenitor-like cells"
    }
    
    cat(sprintf("\nCounts: Mature=%d, Progenitor=%d, Proliferation=%d\n", 
                mature_count, prog_count, prolif_count))
    cat(sprintf("SUGGESTED: %s\n", suggested_id))
    
    results <- rbind(results, data.frame(
      cluster = cl,
      n_mature_neuron = mature_count,
      n_progenitor = prog_count,
      n_proliferation = prolif_count,
      has_DA = has_da,
      has_GLUT = has_glut,
      has_GABA = has_gaba,
      has_choroid = has_choroid,
      suggested_identity = suggested_id,
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n\n=== SUMMARY OF SUGGESTIONS ===\n")
cat("==============================\n")
print(results)

# Save suggestions
write.csv(results, "results/reclustered_analysis/coarse_identity_suggestions.csv", 
          row.names = FALSE)

cat("\n\nKey findings:\n")
cat("- Clusters with high progenitor markers should NOT be called neurons\n")
cat("- Only clusters with clear mature neuron markers should be called neurons\n")
cat("- TTR as top marker often indicates choroid plexus, not just stress\n")
cat("\nSaved to: results/reclustered_analysis/coarse_identity_suggestions.csv\n")