#!/usr/bin/env Rscript

# INVESTIGATE TH EXPRESSION ACROSS ALL CLUSTERS
# Understanding why we detect so few dopaminergic neurons

library(dplyr)
library(ggplot2)
library(tidyr)

cat("=================================================================\n")
cat("INVESTIGATING TH AND DOPAMINERGIC MARKER EXPRESSION\n")
cat("=================================================================\n\n")

# Load marker data
fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# Key dopaminergic markers to investigate
DA_MARKERS <- c("TH", "DDC", "SLC6A3", "SLC18A2", "LMX1A", "FOXA2", 
                "NR4A2", "PITX3", "EN1", "EN2", "DRD2", "KCNJ6",
                "ALDH1A1", "SOX6", "CALB1", "SNCG", "ATP13A2")

cat("1. Searching for TH expression across all clusters...\n")
cat("====================================================\n\n")

# Find TH in fine clusters
th_fine <- fine_markers %>%
  filter(gene == "TH") %>%
  arrange(desc(avg_log2FC))

cat("TH expression in fine clusters:\n")
if (nrow(th_fine) > 0) {
  print(th_fine)
  cat(sprintf("\nTH is expressed in %d out of 36 fine clusters\n", nrow(th_fine)))
  cat("\nTop 5 clusters by TH expression (avg_log2FC):\n")
  print(head(th_fine[, c("cluster", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")], 5))
} else {
  cat("WARNING: TH not found in any fine cluster markers!\n")
}

# Find TH in coarse clusters
th_coarse <- coarse_markers %>%
  filter(gene == "TH") %>%
  arrange(desc(avg_log2FC))

cat("\n\nTH expression in coarse clusters:\n")
if (nrow(th_coarse) > 0) {
  print(th_coarse)
} else {
  cat("WARNING: TH not found in any coarse cluster markers!\n")
}

# 2. Check all DA markers across clusters
cat("\n\n2. Checking all dopaminergic markers...\n")
cat("======================================\n\n")

# Create matrix of DA marker expression
da_expression_matrix <- matrix(0, nrow = 36, ncol = length(DA_MARKERS))
rownames(da_expression_matrix) <- paste0("Cluster_", 0:35)
colnames(da_expression_matrix) <- DA_MARKERS

for (i in 0:35) {
  cluster_markers <- fine_markers %>%
    filter(cluster == as.character(i))
  
  for (j in 1:length(DA_MARKERS)) {
    marker <- DA_MARKERS[j]
    marker_data <- cluster_markers %>% filter(gene == marker)
    if (nrow(marker_data) > 0) {
      da_expression_matrix[i+1, j] <- marker_data$avg_log2FC[1]
    }
  }
}

# Find clusters with any DA marker expression
clusters_with_da <- which(rowSums(da_expression_matrix > 0) > 0) - 1
cat(sprintf("Clusters with ANY dopaminergic marker expression: %d out of 36\n", 
            length(clusters_with_da)))
cat("Cluster IDs:", paste(clusters_with_da, collapse = ", "), "\n")

# 3. Look for clusters with multiple DA markers
cat("\n\n3. Clusters with multiple dopaminergic markers...\n")
cat("================================================\n\n")

da_marker_counts <- data.frame(
  cluster = 0:35,
  n_da_markers = rowSums(da_expression_matrix > 0),
  total_da_score = rowSums(da_expression_matrix),
  has_TH = da_expression_matrix[, "TH"] > 0,
  has_DDC = da_expression_matrix[, "DDC"] > 0,
  TH_score = da_expression_matrix[, "TH"],
  DDC_score = da_expression_matrix[, "DDC"]
)

# Sort by number of DA markers
da_marker_counts <- da_marker_counts %>%
  arrange(desc(n_da_markers), desc(total_da_score))

cat("Top 10 clusters by number of DA markers expressed:\n")
print(head(da_marker_counts, 10))

# 4. Investigate clusters with TH or DDC
cat("\n\n4. Detailed analysis of TH/DDC positive clusters...\n")
cat("==================================================\n\n")

th_or_ddc_clusters <- da_marker_counts %>%
  filter(has_TH | has_DDC)

for (idx in 1:nrow(th_or_ddc_clusters)) {
  cluster_id <- th_or_ddc_clusters$cluster[idx]
  cat(sprintf("\n--- CLUSTER %d ---\n", cluster_id))
  
  # Get all DA markers for this cluster
  cluster_da_markers <- DA_MARKERS[da_expression_matrix[cluster_id + 1, ] > 0]
  cat("DA markers present:", paste(cluster_da_markers, collapse = ", "), "\n")
  
  # Get top 20 markers
  top_markers <- fine_markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(20)
  
  cat("Top 10 overall markers:", paste(head(top_markers$gene, 10), collapse = ", "), "\n")
  
  # Check for neuronal markers
  neuronal_markers <- c("TUBB3", "MAP2", "STMN2", "GAP43", "SYN1", "SNAP25")
  neuronal_present <- top_markers$gene[top_markers$gene %in% neuronal_markers]
  if (length(neuronal_present) > 0) {
    cat("Neuronal markers present:", paste(neuronal_present, collapse = ", "), "\n")
  }
}

# 5. Check if TH might be below detection threshold
cat("\n\n5. Searching for TH in extended marker lists...\n")
cat("==============================================\n\n")

# Function to search deeper in marker lists
search_gene_deeply <- function(gene_name, markers_df, top_n = 200) {
  found_in <- c()
  
  for (cluster_id in 0:35) {
    extended_markers <- markers_df %>%
      filter(cluster == as.character(cluster_id)) %>%
      head(top_n)
    
    if (gene_name %in% extended_markers$gene) {
      rank <- which(extended_markers$gene == gene_name)
      fc <- extended_markers$avg_log2FC[rank]
      found_in <- c(found_in, sprintf("Cluster_%d(rank=%d,FC=%.2f)", cluster_id, rank, fc))
    }
  }
  
  return(found_in)
}

th_deep_search <- search_gene_deeply("TH", fine_markers, top_n = 200)
cat("TH found in extended search (top 200 genes per cluster):\n")
if (length(th_deep_search) > 0) {
  cat(paste(th_deep_search, collapse = "\n"), "\n")
} else {
  cat("TH not found even in top 200 genes of any cluster\n")
}

# 6. Create visualization
cat("\n\n6. Creating visualizations...\n")

# Plot DA marker counts
p1 <- ggplot(da_marker_counts, aes(x = factor(cluster), y = n_da_markers)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_point(aes(y = ifelse(has_TH, n_da_markers + 0.5, NA)), 
             color = "red", size = 3, shape = 16) +
  labs(x = "Cluster", y = "Number of DA markers", 
       title = "Dopaminergic Marker Count by Cluster",
       subtitle = "Red dots indicate TH expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Heatmap of DA marker expression
# Commented out pheatmap code - requires pheatmap package
# library(pheatmap)
# pdf("results/refined_analysis/da_marker_heatmap.pdf", width = 10, height = 8)
# pheatmap(
#   t(da_expression_matrix),
#   cluster_cols = TRUE,
#   cluster_rows = FALSE,
#   color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
#   main = "Dopaminergic Marker Expression Heatmap",
#   fontsize = 8,
#   cellwidth = 15,
#   cellheight = 12,
#   show_colnames = TRUE,
#   angle_col = 45
# )
# dev.off()

ggsave("results/refined_analysis/da_marker_counts.pdf", p1, width = 12, height = 6)

# 7. Summary report
cat("\n\n7. SUMMARY REPORT\n")
cat("=================\n\n")

cat(sprintf("Total clusters analyzed: 36\n"))
cat(sprintf("Clusters with TH expression: %d (%.1f%%)\n", 
            sum(da_marker_counts$has_TH), 
            100 * sum(da_marker_counts$has_TH) / 36))
cat(sprintf("Clusters with DDC expression: %d (%.1f%%)\n", 
            sum(da_marker_counts$has_DDC), 
            100 * sum(da_marker_counts$has_DDC) / 36))
cat(sprintf("Clusters with any DA marker: %d (%.1f%%)\n", 
            length(clusters_with_da), 
            100 * length(clusters_with_da) / 36))

# Expected vs observed
cat("\n\nEXPECTATION vs REALITY:\n")
cat("- Expected: Multiple dopaminergic neuron clusters from Kim 2021 protocol\n")
cat("- Observed: Very limited TH expression across clusters\n")
cat("\nPossible explanations:\n")
cat("1. Differentiation efficiency: Not all cells became dopaminergic\n")
cat("2. Heterogeneity: DA neurons may be a small fraction mixed with other types\n")
cat("3. Maturation state: Cells may be immature and not yet expressing TH strongly\n")
cat("4. Technical issues: TH expression might be below detection threshold\n")
cat("5. Cell stress: Differentiation or experimental conditions may affect TH expression\n")

# Save detailed results
write.csv(da_marker_counts, 
          "results/refined_analysis/da_marker_expression_summary.csv", 
          row.names = FALSE)

cat("\n\nAnalysis complete. Results saved to results/refined_analysis/\n")