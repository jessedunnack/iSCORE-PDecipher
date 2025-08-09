#!/usr/bin/env Rscript

# CHECK IF MARKER FILES ARE LIMITED TO TOP 50

library(dplyr)

cat("=== CHECKING MARKER FILE LIMITS ===\n\n")

# Load both marker files
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")

# Check coarse markers
cat("COARSE MARKERS FILE:\n")
cat("Total rows:", nrow(coarse_markers), "\n")

coarse_counts <- coarse_markers %>%
  group_by(cluster) %>%
  summarise(n_markers = n()) %>%
  arrange(cluster)

cat("\nMarkers per coarse cluster:\n")
print(coarse_counts)

# Check fine markers
cat("\n\nFINE MARKERS FILE:\n")
cat("Total rows:", nrow(fine_markers), "\n")

fine_counts <- fine_markers %>%
  group_by(cluster) %>%
  summarise(n_markers = n()) %>%
  arrange(as.numeric(as.character(cluster)))

cat("\nMarkers per fine cluster:\n")
print(fine_counts)

# Check if all are exactly 50
cat("\n\nANALYSIS:\n")
cat("Coarse clusters - all have exactly 50?", all(coarse_counts$n_markers == 50), "\n")
cat("Fine clusters - all have exactly 50?", all(fine_counts$n_markers == 50), "\n")

# Look for TH position in clusters where it might be just outside top 50
cat("\n\nIMPLICATION: If TH is ranked 51st or lower in cluster 0, it won't appear in this file!\n")
cat("This could explain why you see TH in cluster 0 but I don't - it might be just outside the top 50.\n")