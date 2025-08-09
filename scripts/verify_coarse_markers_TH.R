#!/usr/bin/env Rscript

# VERIFY TH IN COARSE MARKERS - CHECK ALL CLUSTERS

library(dplyr)

cat("=== VERIFYING TH IN ALL COARSE CLUSTERS ===\n\n")

# Load coarse markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")

# Check the structure
cat("Structure of coarse markers data:\n")
cat("Number of rows:", nrow(coarse_markers), "\n")
cat("Columns:", paste(names(coarse_markers), collapse=", "), "\n")
cat("Unique clusters:", paste(sort(unique(coarse_markers$cluster)), collapse=", "), "\n\n")

# Search for TH in ALL data
th_all <- coarse_markers %>%
  filter(gene == "TH")

cat("ALL TH entries in coarse markers:\n")
print(th_all)

# Also check without modifying cluster names
coarse_markers_raw <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
th_raw <- coarse_markers_raw %>%
  filter(gene == "TH")

cat("\n\nTH in raw data (without cluster name modification):\n")
print(th_raw)

# Check if cluster_0 exists in different format
cluster0_variations <- c("0", "cluster_0", "Cluster_0", "cluster0")
cat("\n\nChecking cluster 0 variations:\n")
for (var in cluster0_variations) {
  count <- sum(coarse_markers_raw$cluster == var)
  cat(sprintf("'%s': %d rows\n", var, count))
}

# Get unique cluster names
cat("\n\nAll unique cluster names in raw data:\n")
print(unique(coarse_markers_raw$cluster))