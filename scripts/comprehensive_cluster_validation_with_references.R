#!/usr/bin/env Rscript

# COMPREHENSIVE CLUSTER VALIDATION WITH REFERENCE DATASETS
# Using Kim 2021 protocol expectations, FOUNDIN-PD, macaque MPTP, and La Manno 2016

library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("COMPREHENSIVE CELL TYPE VALIDATION USING REFERENCE DATASETS\n")
cat("iPSC-derived midbrain neurons (Kim 2021 protocol)\n")
cat("=================================================================\n\n")

# Load marker data
cat("Loading marker data...\n")
markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")

# Define comprehensive marker sets based on reference datasets
# FOUNDIN-PD markers for dopaminergic neurons
DA_CORE_MARKERS <- c("TH", "DDC", "SLC6A3", "SLC18A2", "SLC18A1", "DRD2", "KCNJ6", "EN1", "EN2")
DA_TRANSCRIPTION_FACTORS <- c("LMX1A", "LMX1B", "FOXA2", "NR4A2", "PITX3", "OTX2")

# Vulnerability markers from macaque MPTP studies
DA_VULNERABLE_MARKERS <- c("SOX6", "ALDH1A1", "SNCG", "ATP13A2", "AGTR1", "RIT2")  # A9 ventral tier
DA_RESILIENT_MARKERS <- c("CALB1", "CALB2", "OTX2", "GRP", "CCK")  # A10 VTA

# MEIS1/MEIS2 network from MPTP studies
MEIS_NETWORK <- c("MEIS1", "MEIS2", "PBX1", "PBX3", "HOXA9", "HOXA10")

# Floor plate and midbrain markers (La Manno 2016)
FLOOR_PLATE_MARKERS <- c("CORIN", "LMX1A", "FOXA2", "ARX", "SHH", "WNT1", "MSX1", "LMO3")
MIDBRAIN_PROGENITORS <- c("OTX2", "GBX2", "EN1", "EN2", "PAX5", "PAX8", "WNT1", "FGF8")

# Neuronal maturation stages
NEURONAL_EARLY <- c("DCX", "TUBB3", "MAP2", "STMN2", "GAP43")
NEURONAL_MATURE <- c("SYN1", "SNAP25", "SYT1", "RBFOX3", "NEFL", "NEFM", "MAPT")

# Progenitor markers (should NOT be in neurons)
PROGENITOR_MARKERS <- c("SOX2", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7")

# Non-neuronal cell types (from La Manno 2016 human fetal midbrain)
RADIAL_GLIA <- c("VIM", "FABP7", "SLC1A3", "HES1", "HES5", "PAX6", "SOX2", "SOX9")
OLIGODENDROCYTE_LINEAGE <- c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "CSPG4", "NKX2-2")
VASCULAR <- c("CLDN5", "PECAM1", "VWF", "FLT1", "TAGLN", "ACTA2", "MYL9", "PDGFRB")
EPENDYMAL <- c("FOXJ1", "PIFO", "SPEF2", "CFAP126", "DNAH5")

# Choroid plexus (common contaminant)
CHOROID_PLEXUS <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "PRLR", "KCNJ13", "ENPP2", "MSX1", "MSX2")

# Other midbrain neurons (from Kim 2021 expectations)
GLUTAMATERGIC <- c("SLC17A6", "SLC17A7", "GRIN1", "GRIA1", "GRIK1")
GABAERGIC <- c("GAD1", "GAD2", "SLC32A1", "GABBR1", "GABRA1")
SEROTONERGIC <- c("TPH2", "SLC6A4", "HTR1A", "FEV", "PET1")

# Red nucleus and other midbrain
RED_NUCLEUS <- c("LHX9", "LMX1B", "POU4F1", "ISL1")

# Function to analyze each cluster
analyze_cluster <- function(cluster_id) {
  cluster_markers <- markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC))
  
  # Get top 50 markers
  top50 <- head(cluster_markers, 50)
  
  # Calculate scores for each cell type
  scores <- list()
  
  # Helper function
  calc_score <- function(marker_set, markers_df) {
    present <- marker_set[marker_set %in% markers_df$gene]
    if (length(present) == 0) return(list(score = 0, markers = character(0)))
    
    marker_scores <- markers_df %>%
      filter(gene %in% present) %>%
      summarise(
        score = sum(avg_log2FC),
        n_markers = n(),
        markers = paste(gene, collapse = ", ")
      )
    
    return(list(
      score = marker_scores$score,
      n_markers = marker_scores$n_markers,
      markers = marker_scores$markers
    ))
  }
  
  # Calculate all scores
  scores$da_core <- calc_score(DA_CORE_MARKERS, top50)
  scores$da_vulnerable <- calc_score(DA_VULNERABLE_MARKERS, top50)
  scores$da_resilient <- calc_score(DA_RESILIENT_MARKERS, top50)
  scores$floor_plate <- calc_score(FLOOR_PLATE_MARKERS, top50)
  scores$neuronal_early <- calc_score(NEURONAL_EARLY, top50)
  scores$neuronal_mature <- calc_score(NEURONAL_MATURE, top50)
  scores$progenitor <- calc_score(PROGENITOR_MARKERS, top50)
  scores$choroid_plexus <- calc_score(CHOROID_PLEXUS, top50)
  scores$glutamatergic <- calc_score(GLUTAMATERGIC, top50)
  scores$gabaergic <- calc_score(GABAERGIC, top50)
  scores$radial_glia <- calc_score(RADIAL_GLIA, top50)
  scores$oligodendrocyte <- calc_score(OLIGODENDROCYTE_LINEAGE, top50)
  scores$vascular <- calc_score(VASCULAR, top50)
  
  # MEIS network check
  meis_present <- MEIS_NETWORK[MEIS_NETWORK %in% top50$gene]
  
  return(list(
    cluster = cluster_id,
    top10_genes = paste(head(top50$gene, 10), collapse = ", "),
    scores = scores,
    meis_network = meis_present,
    n_cells = unique(cluster_markers$pct.1[1] * 100) # Approximate
  ))
}

# Analyze all 36 clusters
cat("\nANALYZING ALL 36 CLUSTERS:\n")
cat("==========================\n\n")

all_results <- list()

for (cluster in 0:35) {
  cat(sprintf("\n=== CLUSTER %d ===\n", cluster))
  
  result <- analyze_cluster(cluster)
  all_results[[as.character(cluster)]] <- result
  
  # Print top 10 genes
  cat("Top 10 markers:", result$top10_genes, "\n\n")
  
  # Decision logic based on Kim 2021 protocol expectations
  cell_type <- "Unknown"
  confidence <- "Low"
  subtype <- ""
  vulnerability <- "Unknown"
  
  # Extract scores
  da_core_score <- result$scores$da_core$score
  da_vuln_score <- result$scores$da_vulnerable$score
  da_resil_score <- result$scores$da_resilient$score
  floor_plate_score <- result$scores$floor_plate$score
  neuronal_early_score <- result$scores$neuronal_early$score
  neuronal_mature_score <- result$scores$neuronal_mature$score
  progenitor_score <- result$scores$progenitor$score
  cp_score <- result$scores$choroid_plexus$score
  
  # Decision tree
  if (cp_score > 20 && grepl("TTR", result$top10_genes)) {
    cell_type <- "Choroid Plexus"
    confidence <- "High"
    cat("CHOROID PLEXUS detected (TTR expression)\n")
  } else if (da_core_score > 5 && grepl("TH|DDC", result$top10_genes)) {
    cell_type <- "Dopaminergic Neurons"
    confidence <- ifelse(da_core_score > 10, "High", "Medium")
    
    # Determine subtype
    if (da_vuln_score > da_resil_score && length(result$meis_network) > 0) {
      subtype <- "A9-like (vulnerable)"
      vulnerability <- "High (SOX6+/ALDH1A1+/MEIS+)"
    } else if (da_resil_score > da_vuln_score) {
      subtype <- "A10-like (resilient)"
      vulnerability <- "Low (CALB1+)"
    } else {
      subtype <- "Unspecified"
      vulnerability <- "Medium"
    }
    
    cat(sprintf("DOPAMINERGIC NEURONS - %s\n", subtype))
    cat(sprintf("  Vulnerability: %s\n", vulnerability))
    cat("  DA markers:", result$scores$da_core$markers, "\n")
  } else if (progenitor_score > 5 && neuronal_early_score < 2) {
    cell_type <- "Neural Progenitors"
    confidence <- "High"
    
    if (floor_plate_score > 5) {
      subtype <- "Floor Plate Progenitors"
    } else {
      subtype <- "Midbrain Progenitors"
    }
    cat(sprintf("PROGENITORS - %s\n", subtype))
  } else if (neuronal_early_score > 10 || neuronal_mature_score > 5) {
    # Neurons but not DA
    if (floor_plate_score > 5) {
      cell_type <- "Floor Plate-Derived Neurons"
      cat("FLOOR PLATE-DERIVED NEURONS (non-DA)\n")
    } else if (result$scores$glutamatergic$score > 3) {
      cell_type <- "Glutamatergic Neurons"
      cat("GLUTAMATERGIC NEURONS\n")
    } else if (result$scores$gabaergic$score > 3) {
      cell_type <- "GABAergic Neurons"
      cat("GABAERGIC NEURONS\n")
    } else {
      cell_type <- "Neurons (unspecified)"
      cat("NEURONS (type unclear)\n")
    }
    confidence <- "Medium"
  } else if (result$scores$radial_glia$score > 10) {
    cell_type <- "Radial Glia"
    confidence <- "Medium"
    cat("RADIAL GLIA\n")
  } else if (result$scores$oligodendrocyte$score > 5) {
    cell_type <- "Oligodendrocyte Lineage"
    confidence <- "Medium"
    cat("OLIGODENDROCYTE LINEAGE\n")
  } else if (result$scores$vascular$score > 5) {
    cell_type <- "Vascular Cells"
    confidence <- "Medium"
    cat("VASCULAR CELLS\n")
  }
  
  # Print scores
  cat("\nMarker Scores:\n")
  if (da_core_score > 0) cat(sprintf("  DA core: %.2f\n", da_core_score))
  if (da_vuln_score > 0) cat(sprintf("  DA vulnerable (A9): %.2f\n", da_vuln_score))
  if (da_resil_score > 0) cat(sprintf("  DA resilient (A10): %.2f\n", da_resil_score))
  if (neuronal_early_score > 0) cat(sprintf("  Neuronal (early): %.2f\n", neuronal_early_score))
  if (neuronal_mature_score > 0) cat(sprintf("  Neuronal (mature): %.2f\n", neuronal_mature_score))
  if (progenitor_score > 0) cat(sprintf("  Progenitor: %.2f\n", progenitor_score))
  if (cp_score > 0) cat(sprintf("  Choroid plexus: %.2f\n", cp_score))
  
  if (length(result$meis_network) > 0) {
    cat("  MEIS network genes:", paste(result$meis_network, collapse = ", "), "\n")
  }
  
  # Store results
  all_results[[as.character(cluster)]]$cell_type <- cell_type
  all_results[[as.character(cluster)]]$confidence <- confidence
  all_results[[as.character(cluster)]]$subtype <- subtype
  all_results[[as.character(cluster)]]$vulnerability <- vulnerability
}

# Create summary table
cat("\n\n=== FINAL CLUSTER ASSIGNMENTS ===\n")
cat("=================================\n\n")

summary_df <- data.frame(
  cluster = 0:35,
  cell_type = sapply(0:35, function(x) all_results[[as.character(x)]]$cell_type),
  subtype = sapply(0:35, function(x) all_results[[as.character(x)]]$subtype),
  confidence = sapply(0:35, function(x) all_results[[as.character(x)]]$confidence),
  vulnerability = sapply(0:35, function(x) all_results[[as.character(x)]]$vulnerability),
  top_markers = sapply(0:35, function(x) {
    top5 <- strsplit(all_results[[as.character(x)]]$top10_genes, ", ")[[1]][1:5]
    paste(top5, collapse = ", ")
  })
)

# Print summary
print(summary_df)

# Save results
dir.create("results/comprehensive_validation", recursive = TRUE, showWarnings = FALSE)
write.csv(summary_df, "results/comprehensive_validation/cluster_assignments_with_references.csv", row.names = FALSE)
saveRDS(all_results, "results/comprehensive_validation/detailed_cluster_analysis.rds")

# Summary statistics
cat("\n\nSUMMARY STATISTICS:\n")
cat("===================\n")
cell_type_counts <- table(summary_df$cell_type)
print(cell_type_counts)

cat("\n\nDopaminergic neuron subtypes:\n")
da_clusters <- summary_df[summary_df$cell_type == "Dopaminergic Neurons", ]
if (nrow(da_clusters) > 0) {
  print(da_clusters[, c("cluster", "subtype", "vulnerability")])
}

cat("\n\nClusters requiring attention:\n")
low_conf <- summary_df[summary_df$confidence == "Low", ]
if (nrow(low_conf) > 0) {
  print(low_conf[, c("cluster", "cell_type", "top_markers")])
}

cat("\n\nAnalysis complete! Results saved to results/comprehensive_validation/\n")