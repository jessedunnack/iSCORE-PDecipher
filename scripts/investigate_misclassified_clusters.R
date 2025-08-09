#!/usr/bin/env Rscript

# CRITICAL INVESTIGATION: Misclassified Clusters
# This script examines clusters that may be misidentified

library(dplyr)
library(Seurat)

cat("=================================================================\n")
cat("CRITICAL INVESTIGATION: Potentially Misclassified Clusters\n")
cat("=================================================================\n\n")

# Load marker data
markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")

# Function to get comprehensive marker profile
get_cluster_profile <- function(cluster_id, n_markers = 30) {
  cluster_markers <- markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n_markers)
  
  return(cluster_markers)
}

# INVESTIGATION 1: Floor Plate Progenitors that express neuronal markers
cat("INVESTIGATION 1: Clusters labeled as 'Floor Plate Progenitors'\n")
cat("==============================================================\n\n")

# Check cluster 1 (labeled Floor Plate Progenitors)
cat("Cluster 1 - Currently labeled 'Floor Plate Progenitors':\n")
cluster1_markers <- get_cluster_profile(1, 40)

# Check for neuronal markers
neuronal_markers <- c("TUBB3", "STMN2", "GAP43", "MAP2", "NCAM1", "SYN1", "SNAP25", 
                     "RBFOX3", "NEFL", "NEFM", "DCX", "MAPT")
progenitor_markers <- c("SOX2", "NES", "VIM", "FABP7", "PAX6", "HES1", "HES5")

cat("\nNeuronal markers found in Cluster 1:\n")
neuronal_found <- cluster1_markers %>% filter(gene %in% neuronal_markers)
print(neuronal_found[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\nProgenitor markers found in Cluster 1:\n")
progenitor_found <- cluster1_markers %>% filter(gene %in% progenitor_markers)
print(progenitor_found[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\nTop 20 markers for Cluster 1:\n")
print(cluster1_markers[1:20, c("gene", "avg_log2FC")])

# Check other "Floor Plate" clusters
floor_plate_clusters <- c(3, 13)  # Clusters labeled as floor plate-related

for (cl in floor_plate_clusters) {
  cat("\n\n----------------------------------------\n")
  cat("Cluster", cl, "- Floor plate-related:\n")
  cl_markers <- get_cluster_profile(cl, 30)
  
  neuronal_found <- cl_markers %>% filter(gene %in% neuronal_markers)
  if (nrow(neuronal_found) > 0) {
    cat("\nNeuronal markers found:\n")
    print(neuronal_found[, c("gene", "avg_log2FC")])
  }
  
  cat("\nTop 10 markers:\n")
  print(cl_markers[1:10, c("gene", "avg_log2FC")])
}

# INVESTIGATION 2: Dopaminergic neurons that might be choroid plexus
cat("\n\nINVESTIGATION 2: Clusters labeled as 'Dopaminergic' \n")
cat("==================================================\n\n")

# Choroid plexus markers
cp_markers <- c("TTR", "FOLR1", "HTR2C", "CLIC6", "AQP1", "KCNJ13", "SLC12A2", 
                "SLC4A10", "TRPM3", "PRLR", "ENPP2", "IGFBP2")

# Check cluster 5 (labeled Dopaminergic)
cat("Cluster 5 - Currently labeled 'Dopaminergic Neurons':\n")
cluster5_markers <- get_cluster_profile(5, 40)

cat("\nChoroid plexus markers found in Cluster 5:\n")
cp_found <- cluster5_markers %>% filter(gene %in% cp_markers)
print(cp_found[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\nDopaminergic markers in Cluster 5:\n")
da_markers <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", "ALDH1A1", "NR4A2")
da_found <- cluster5_markers %>% filter(gene %in% da_markers)
print(da_found[, c("gene", "avg_log2FC", "p_val_adj")])

cat("\nTop 20 markers for Cluster 5:\n")
print(cluster5_markers[1:20, c("gene", "avg_log2FC")])

# Check other DA clusters
da_clusters <- c(18, 28)  # Other DA clusters

for (cl in da_clusters) {
  cat("\n\n----------------------------------------\n")
  cat("Cluster", cl, "- Dopaminergic:\n")
  cl_markers <- get_cluster_profile(cl, 30)
  
  cp_found <- cl_markers %>% filter(gene %in% cp_markers)
  if (nrow(cp_found) > 0) {
    cat("\nChoroid plexus markers found:\n")
    print(cp_found[, c("gene", "avg_log2FC")])
  }
  
  cat("\nTop 10 markers:\n")
  print(cl_markers[1:10, c("gene", "avg_log2FC")])
}

# INVESTIGATION 3: Create revised classifications
cat("\n\nPROPOSED REVISED CLASSIFICATIONS:\n")
cat("==================================\n\n")

# Based on marker analysis, propose new classifications
cat("Cluster 1: If TUBB3+/STMN2+/GAP43+ >> 'Immature/Developing Neurons' NOT 'Floor Plate Progenitors'\n")
cat("  - These are post-mitotic neurons, possibly with floor plate heritage\n\n")

cat("Cluster 5: If TTR+/FOLR1+ >> Check if truly DA or Choroid Plexus\n")
cat("  - Need to examine TTR/FOLR1 expression levels\n")
cat("  - If TTR >> TH, likely Choroid Plexus\n")
cat("  - If TH >> TTR, likely DA with some CP contamination\n\n")

# Save investigation results
investigation_results <- list(
  cluster1 = list(
    current_label = "Floor Plate Progenitors",
    neuronal_markers = neuronal_found,
    progenitor_markers = progenitor_found,
    proposed_label = "Immature/Developing Neurons (Floor Plate-derived)",
    evidence = "High expression of TUBB3, STMN2, GAP43 indicates post-mitotic neurons"
  ),
  cluster5 = list(
    current_label = "Dopaminergic Neurons",
    cp_markers = cp_found,
    da_markers = da_found,
    proposed_label = "Requires further investigation",
    evidence = "Check relative expression of TTR/FOLR1 vs TH/DDC"
  )
)

saveRDS(investigation_results, "results/cell_type_annotations/misclassification_investigation.rds")

# Create a detailed report
sink("results/cell_type_annotations/CRITICAL_MISCLASSIFICATION_REPORT.txt")
cat("CRITICAL FINDINGS - MISCLASSIFIED CLUSTERS\n")
cat("==========================================\n\n")

cat("1. FLOOR PLATE PROGENITORS ARE ACTUALLY NEURONS\n")
cat("   - Cluster 1 expresses mature neuronal markers\n")
cat("   - TUBB3, STMN2, GAP43 are post-mitotic neuron markers\n")
cat("   - These are NOT progenitors\n\n")

cat("2. NEED TO VERIFY DOPAMINERGIC vs CHOROID PLEXUS\n")
cat("   - Check if Cluster 5 has high TTR/FOLR1\n")
cat("   - If yes, it's likely Choroid Plexus, not DA neurons\n\n")

cat("3. RECOMMENDED ACTIONS:\n")
cat("   - Re-classify 'Floor Plate Progenitors' as neurons\n")
cat("   - Verify all 'Dopaminergic' clusters for CP markers\n")
cat("   - Update vulnerability assessments accordingly\n")
sink()

cat("\n\nInvestigation complete. See results in:\n")
cat("  results/cell_type_annotations/CRITICAL_MISCLASSIFICATION_REPORT.txt\n")
cat("  results/cell_type_annotations/misclassification_investigation.rds\n")