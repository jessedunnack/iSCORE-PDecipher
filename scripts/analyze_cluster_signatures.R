# Analyze which clusters show significant signatures for each gene
# Focus on cluster identity and cell type associations

# Load previous results
flat_results <- read.csv("fisher_test_results_flat.csv", stringsAsFactors = FALSE)

# Sort by gene pair and p-value
flat_results <- flat_results[order(flat_results$gene_pair, flat_results$fisher_p), ]

cat("=== CLUSTER-SPECIFIC SIGNATURES BY GENE ===\n\n")

# Analyze each gene pair
gene_pairs <- unique(flat_results$gene_pair)

for (gene_pair in gene_pairs) {
  cat("\n", gene_pair, "\n", sep="")
  cat(rep("=", nchar(gene_pair)), "\n", sep="")
  
  # Get results for this gene pair
  gene_data <- flat_results[flat_results$gene_pair == gene_pair, ]
  
  # Only show significant results
  sig_data <- gene_data[gene_data$fisher_p < 0.05, ]
  
  if (nrow(sig_data) > 0) {
    # Group by cluster
    clusters <- unique(sig_data$cluster)
    
    cat("Significant in", length(clusters), "clusters:\n\n")
    
    for (cluster in clusters) {
      cluster_data <- sig_data[sig_data$cluster == cluster, ]
      
      # Find best p-value for this cluster
      best_idx <- which.min(cluster_data$fisher_p)
      best_row <- cluster_data[best_idx, ]
      
      cat("  ", cluster, ":\n", sep="")
      cat("    Best p-value: ", format(best_row$fisher_p, scientific=TRUE), 
          " (", best_row$experiment, ")\n", sep="")
      cat("    Overlap: ", best_row$total_overlap, " genes (", 
          best_row$same_direction, " same, ", 
          best_row$opposite_direction, " opposite)\n", sep="")
      cat("    Odds ratio: ", round(best_row$odds_ratio, 2), "\n", sep="")
      
      # Show all significant experiments for this cluster
      if (nrow(cluster_data) > 1) {
        cat("    Also significant in: ")
        other_exps <- setdiff(cluster_data$experiment, best_row$experiment)
        cat(paste(other_exps, collapse=", "), "\n")
      }
      cat("\n")
    }
    
    # Summary statistics
    cat("  Summary:\n")
    cat("    Total significant comparisons: ", nrow(sig_data), "\n", sep="")
    cat("    Clusters involved: ", paste(sort(clusters), collapse=", "), "\n", sep="")
    cat("    Strongest cluster: ", sig_data$cluster[1], 
        " (p=", format(sig_data$fisher_p[1], scientific=TRUE), ")\n", sep="")
  } else {
    cat("  No significant overlaps found\n")
  }
}

# Create cluster summary table
cat("\n\n=== CLUSTER SUMMARY TABLE ===\n\n")

library(dplyr)
cluster_summary <- flat_results %>%
  filter(fisher_p < 0.05) %>%
  group_by(cluster) %>%
  summarise(
    n_significant = n(),
    n_genes = n_distinct(gene_pair),
    best_p_value = min(fisher_p),
    best_gene = gene_pair[which.min(fisher_p)],
    avg_overlap = round(mean(total_overlap)),
    max_overlap = max(total_overlap),
    .groups = 'drop'
  ) %>%
  arrange(best_p_value)

print(cluster_summary)

# Identify cluster patterns
cat("\n\n=== CLUSTER PATTERNS ===\n\n")

# Which clusters show the most PD gene signatures?
top_clusters <- cluster_summary %>%
  filter(n_genes >= 5) %>%
  arrange(desc(n_genes))

cat("Clusters with signatures from 5+ PD genes:\n")
print(top_clusters)

# Create heatmap data
cat("\n\nCreating cluster-gene heatmap data...\n")

# Create matrix of -log10(p-values)
gene_cluster_matrix <- flat_results %>%
  filter(fisher_p < 0.05) %>%
  group_by(gene_pair, cluster) %>%
  summarise(
    best_p = min(fisher_p),
    .groups = 'drop'
  ) %>%
  mutate(neg_log_p = -log10(best_p)) %>%
  select(-best_p) %>%
  pivot_wider(names_from = cluster, values_from = neg_log_p, values_fill = 0)

# Save for visualization
write.csv(gene_cluster_matrix, "gene_cluster_significance_matrix.csv", row.names = FALSE)

# Top cluster-specific signatures
cat("\n\n=== TOP CLUSTER-SPECIFIC SIGNATURES ===\n\n")

for (cluster in unique(flat_results$cluster)) {
  cluster_data <- flat_results[flat_results$cluster == cluster & flat_results$fisher_p < 0.05, ]
  
  if (nrow(cluster_data) > 0) {
    cat("\n", cluster, ":\n", sep="")
    
    # Top 3 signatures in this cluster
    top3 <- head(cluster_data[order(cluster_data$fisher_p), ], 3)
    
    for (i in 1:nrow(top3)) {
      cat("  ", i, ". ", top3$gene_pair[i], 
          " (p=", format(top3$fisher_p[i], scientific=TRUE), 
          ", ", top3$total_overlap[i], " genes)\n", sep="")
    }
  }
}

cat("\n\nAnalysis complete. Files created:\n")
cat("- gene_cluster_significance_matrix.csv\n")