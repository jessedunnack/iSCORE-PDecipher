#!/usr/bin/env Rscript

# Quick check of TH expression in clusters 1 and 5

library(dplyr)

cat("=== CHECKING TH EXPRESSION IN CLUSTERS 1 AND 5 ===\n\n")

# Load fine cluster markers
fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")

# Check cluster 1
cat("CLUSTER 1 markers:\n")
cluster1_markers <- fine_markers %>%
  filter(cluster == "1") %>%
  arrange(desc(avg_log2FC)) %>%
  head(30)

if ("TH" %in% cluster1_markers$gene) {
  th_data <- cluster1_markers %>% filter(gene == "TH")
  cat("TH FOUND in cluster 1:\n")
  print(th_data)
} else {
  cat("TH NOT FOUND in top 30 markers of cluster 1\n")
}

cat("\nTop 10 markers for cluster 1:\n")
print(head(cluster1_markers$gene, 10))

# Check cluster 5
cat("\n\nCLUSTER 5 markers:\n")
cluster5_markers <- fine_markers %>%
  filter(cluster == "5") %>%
  arrange(desc(avg_log2FC)) %>%
  head(30)

if ("TH" %in% cluster5_markers$gene) {
  th_data <- cluster5_markers %>% filter(gene == "TH")
  cat("TH FOUND in cluster 5:\n")
  print(th_data)
} else {
  cat("TH NOT FOUND in top 30 markers of cluster 5\n")
}

cat("\nTop 10 markers for cluster 5:\n")
print(head(cluster5_markers$gene, 10))

# Search for TH in all clusters
cat("\n\nSearching for TH in ALL clusters:\n")
th_all <- fine_markers %>% 
  filter(gene == "TH") %>%
  arrange(desc(avg_log2FC))

print(th_all)