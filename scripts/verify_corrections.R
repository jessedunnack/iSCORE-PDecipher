#!/usr/bin/env Rscript

# Quick verification of corrections

library(Seurat)
library(dplyr)

cat("Loading corrected object...\n")
x <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_corrected.rds")

cat("\nChecking corrected cell types:\n")
cat("=============================\n")

# Check cluster 1 (was Floor Plate Progenitors)
cat("\nCluster 1:\n")
cat("  Original:", unique(x$validated_celltype[x$seurat_clusters_fine == 1])[1], "\n")
cat("  Corrected:", unique(x$corrected_celltype[x$seurat_clusters_fine == 1])[1], "\n")

# Check cluster 5 (was Dopaminergic Neurons)
cat("\nCluster 5:\n")
cat("  Original:", unique(x$validated_celltype[x$seurat_clusters_fine == 5])[1], "\n")
cat("  Corrected:", unique(x$corrected_celltype[x$seurat_clusters_fine == 5])[1], "\n")

# Summary of all changes
cat("\n\nSummary of all corrected clusters:\n")
changes <- x@meta.data %>%
  group_by(seurat_clusters_fine) %>%
  summarise(
    original = first(validated_celltype),
    corrected = first(corrected_celltype),
    changed = first(validated_celltype) != first(corrected_celltype)
  ) %>%
  filter(changed)

print(changes)