#!/usr/bin/env Rscript

# SYSTEMATIC COARSE CLUSTER REVIEW
# For user approval of each cluster assignment

library(dplyr)

cat("=================================================================\n")
cat("SYSTEMATIC COARSE CLUSTER REVIEW\n")
cat("15 coarse clusters (0-14) for review\n")
cat("=================================================================\n\n")

# Load coarse cluster markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# Function to display cluster info
display_cluster <- function(cluster_id) {
  cat(sprintf("\n========== COARSE CLUSTER %d ==========\n", cluster_id))
  
  # Get top markers
  top_markers <- coarse_markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(15)
  
  cat("\nTop 15 marker genes:\n")
  for (i in 1:nrow(top_markers)) {
    cat(sprintf("%2d. %-12s (FC=%.2f, pct.1=%.3f)\n", 
                i, 
                top_markers$gene[i], 
                top_markers$avg_log2FC[i],
                top_markers$pct.1[i]))
  }
  
  cat("\n")
}

# Review each coarse cluster
for (i in 0:14) {
  display_cluster(i)
  cat("PROPOSED LABEL: [Awaiting systematic review]\n")
  cat("USER FEEDBACK: [Pending]\n")
  cat("\n-------------------------------------------\n")
}