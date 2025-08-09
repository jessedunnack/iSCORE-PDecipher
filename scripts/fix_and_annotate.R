#\!/usr/bin/env Rscript
library(Seurat)
library(dplyr)

# Load and recluster
s <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds")

# Check neighbors
if (!"SCT_snn" %in% names(s@graphs)) {
  s <- FindNeighbors(s, dims = 1:30)
}

# Cluster at two resolutions
s <- FindClusters(s, resolution = 0.2)
s$seurat_clusters_coarse <- Idents(s)

s <- FindClusters(s, resolution = 0.8)
s$seurat_clusters_fine <- Idents(s)

# Save
saveRDS(s, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed.rds")

# Then annotate
source("scripts/add_celltype_annotations_to_seurat.R")
s_ann <- main("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed.rds")
