#!/usr/bin/env Rscript

# CAREFUL RE-EVALUATION OF COARSE CLUSTER IDENTITIES - FIXED
# Distinguishing neurons from progenitors

library(dplyr)

cat("=================================================================\n")
cat("CAREFUL COARSE CLUSTER IDENTITY REVIEW\n")
cat("=================================================================\n\n")

# Load coarse markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")

# Define marker categories more carefully
MATURE_NEURON_MARKERS <- c("MAP2", "RBFOX3", "SYN1", "SNAP25", "STMN2", "TUBB3", 
                          "NEFL", "NEFM", "NEFH", "ENO2", "SYT1", "UCHL1", "DCX")

PROGENITOR_MARKERS <- c("SOX2", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7", 
                       "BLBP", "SOX1", "NOTCH1", "ID1", "ID3", "ASCL1", "EOMES")

PROLIFERATION_MARKERS <- c("MKI67", "TOP2A", "PCNA", "CCND1", "CCND2", "CDC20", 
                          "UBE2C", "CENPF", "HIST1H3B", "HMGB2")

# Specific neuronal subtypes
DA_SPECIFIC <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", "ALDH1A1", "CALB1")
GLUT_SPECIFIC <- c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "GRM1", "GRM5")
GABA_SPECIFIC <- c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "DLX5", "DLX6")

# Non-neuronal
CHOROID_PLEXUS <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13", "ENPP2", "MSX1", "MSX2")
GLIA_MARKERS <- c("GFAP", "S100B", "ALDH1L1", "OLIG1", "OLIG2", "SOX10", "AQP4")

cat("Analyzing each cluster for neuronal maturity and identity...\n")
cat("==========================================================\n")

results <- data.frame(
  cluster = integer(),
  n_cells = integer(),
  n_mature_neuron = integer(),
  n_progenitor = integer(),
  n_proliferation = integer(),
  has_DA = logical(),
  has_GLUT = logical(),
  has_GABA = logical(),
  has_choroid = logical(),
  top_marker = character(),
  suggested_identity = character(),
  stringsAsFactors = FALSE
)

# Get cluster names
cluster_names <- unique(coarse_markers$cluster)

for (cl_name in cluster_names) {
  # Extract cluster number
  cl_num <- as.integer(gsub("cluster_", "", cl_name))
  
  # Get markers for this cluster
  cl_markers <- coarse_markers %>%
    filter(cluster == cl_name) %>%
    arrange(desc(avg_log2FC))
  
  # Get number of cells from previous analysis
  n_cells <- NA
  if (file.exists("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")) {
    prev_results <- read.csv("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")
    if (cl_num %in% prev_results$cluster) {
      n_cells <- prev_results$n_cells[prev_results$cluster == cl_num]
    }
  }
  
  cat(sprintf("\n--- CLUSTER %d (n=%s cells) ---\n", cl_num, 
              ifelse(is.na(n_cells), "?", as.character(n_cells))))
  
  # Look at top 50 markers
  top50 <- head(cl_markers, 50)
  
  # Count marker types in top 50
  mature_count <- sum(top50$gene %in% MATURE_NEURON_MARKERS)
  prog_count <- sum(top50$gene %in% PROGENITOR_MARKERS)
  prolif_count <- sum(top50$gene %in% PROLIFERATION_MARKERS)
  
  # Check specific types
  has_da <- any(top50$gene %in% DA_SPECIFIC)
  has_glut <- any(top50$gene %in% GLUT_SPECIFIC)
  has_gaba <- any(top50$gene %in% GABA_SPECIFIC)
  has_choroid <- any(top50$gene %in% CHOROID_PLEXUS)
  
  # Print top 15 markers
  cat("Top 15 markers:\n")
  for (i in 1:min(15, nrow(cl_markers))) {
    gene <- cl_markers$gene[i]
    marker_type <- ""
    
    if (gene %in% MATURE_NEURON_MARKERS) marker_type <- "[MATURE]"
    if (gene %in% PROGENITOR_MARKERS) marker_type <- "[PROG]"
    if (gene %in% PROLIFERATION_MARKERS) marker_type <- "[PROLIF]"
    if (gene %in% DA_SPECIFIC) marker_type <- "[DA]"
    if (gene %in% GLUT_SPECIFIC) marker_type <- "[GLUT]"
    if (gene %in% GABA_SPECIFIC) marker_type <- "[GABA]"
    if (gene %in% CHOROID_PLEXUS) marker_type <- "[CHOROID]"
    if (gene %in% GLIA_MARKERS) marker_type <- "[GLIA]"
    
    cat(sprintf("  %2d. %-12s %.2f FC, %.1f%% cells %s\n", 
                i, gene, cl_markers$avg_log2FC[i], 
                cl_markers$pct.1[i] * 100, marker_type))
  }
  
  # Determine identity
  suggested_id <- "Unknown"
  top_marker <- cl_markers$gene[1]
  
  # Strong choroid plexus signature?
  if (top_marker %in% CHOROID_PLEXUS || (has_choroid && sum(head(cl_markers$gene, 5) %in% CHOROID_PLEXUS) >= 2)) {
    suggested_id <- "Choroid plexus"
  }
  # Proliferating cells?
  else if (prolif_count >= 3 && top_marker %in% PROLIFERATION_MARKERS) {
    suggested_id <- "Proliferating cells"
  }
  # Progenitors?
  else if (prog_count > mature_count && prog_count >= 3) {
    if (prolif_count >= 2) {
      suggested_id <- "Proliferating progenitors"
    } else {
      suggested_id <- "Neural progenitors"
    }
  }
  # Mature neurons - what type?
  else if (mature_count >= 3 || (mature_count >= 2 && top50$gene[1] %in% MATURE_NEURON_MARKERS)) {
    if (has_da && any(head(cl_markers$gene, 20) %in% c("TH", "DDC"))) {
      suggested_id <- "Dopaminergic neurons"
    } else if (has_glut && any(head(cl_markers$gene, 20) %in% GLUT_SPECIFIC)) {
      suggested_id <- "Glutamatergic neurons"
    } else if (has_gaba && any(head(cl_markers$gene, 20) %in% GABA_SPECIFIC)) {
      suggested_id <- "GABAergic neurons"
    } else {
      suggested_id <- "Neurons (type unclear)"
    }
  }
  # Intermediate/transitioning cells?
  else if (mature_count >= 1 && prog_count >= 1) {
    suggested_id <- "Intermediate/transitioning cells"
  }
  # Still progenitor-like?
  else if (prog_count >= 2) {
    suggested_id <- "Progenitor-like cells"
  }
  # Check for glia
  else if (any(head(cl_markers$gene, 10) %in% GLIA_MARKERS)) {
    suggested_id <- "Glial cells"
  }
  
  cat(sprintf("\nCounts: Mature=%d, Progenitor=%d, Proliferation=%d\n", 
              mature_count, prog_count, prolif_count))
  cat(sprintf("Top marker: %s\n", top_marker))
  cat(sprintf("SUGGESTED: %s\n", suggested_id))
  
  results <- rbind(results, data.frame(
    cluster = cl_num,
    n_cells = ifelse(is.na(n_cells), 0, n_cells),
    n_mature_neuron = mature_count,
    n_progenitor = prog_count,
    n_proliferation = prolif_count,
    has_DA = has_da,
    has_GLUT = has_glut,
    has_GABA = has_gaba,
    has_choroid = has_choroid,
    top_marker = top_marker,
    suggested_identity = suggested_id,
    stringsAsFactors = FALSE
  ))
}

# Sort by cluster number
results <- results %>% arrange(cluster)

cat("\n\n=== SUMMARY OF SUGGESTIONS ===\n")
cat("==============================\n")
print(results[, c("cluster", "n_cells", "top_marker", "suggested_identity")])

# Save suggestions
write.csv(results, "results/reclustered_analysis/coarse_identity_suggestions.csv", 
          row.names = FALSE)

cat("\n\nKey findings:\n")
cat("- Clusters with high progenitor markers should NOT be called neurons\n")
cat("- Only clusters with clear mature neuron markers should be called neurons\n")
cat("- TTR/FOLR1 as top markers often indicates choroid plexus\n")
cat("- Check for proliferation markers (MKI67, TOP2A) for cycling cells\n")
cat("\nSaved to: results/reclustered_analysis/coarse_identity_suggestions.csv\n")