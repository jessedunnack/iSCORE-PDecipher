#!/usr/bin/env Rscript

# CHECK COARSE CLUSTER 0 FOR TH

library(dplyr)

cat("=== CHECKING COARSE CLUSTER 0 FOR TH ===\n\n")

# Load coarse markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# Get ALL markers for cluster 0
cluster0_all <- coarse_markers %>%
  filter(cluster == "0") %>%
  arrange(desc(avg_log2FC))

cat("Total markers in cluster 0:", nrow(cluster0_all), "\n\n")

# Search for TH
if ("TH" %in% cluster0_all$gene) {
  th_position <- which(cluster0_all$gene == "TH")
  th_data <- cluster0_all[th_position, ]
  cat("TH FOUND in cluster 0!\n")
  cat("Position in ranked list:", th_position, "\n")
  print(th_data)
} else {
  cat("TH NOT FOUND in cluster 0 markers\n")
}

# Show top 50 markers
cat("\n\nTop 50 markers in cluster 0:\n")
print(head(cluster0_all[, c("gene", "avg_log2FC", "pct.1", "p_val_adj")], 50))

# Search for other DA markers in cluster 0
cat("\n\nSearching for other DA markers in cluster 0:\n")
da_markers <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", 
                "LMX1A", "FOXA2", "NR4A2", "PITX3", "EN1", "EN2")

da_found <- cluster0_all %>%
  filter(gene %in% da_markers)

if (nrow(da_found) > 0) {
  cat("DA markers found:\n")
  print(da_found[, c("gene", "avg_log2FC", "pct.1", "p_val_adj")])
} else {
  cat("No DA markers found in cluster 0\n")
}