# PD Signature Cluster Analysis
# Purpose: Analyze pathway signatures across all clusters
# Date: 2025-07-19

library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)

# Configuration
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_cluster"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
all_data <- readRDS(data_file)

# Add gene column if needed
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# PD keyword search
pd_terms <- c("mitochondr", "lysosom", "autophagy", "synap", "dopamin", 
              "proteasom", "synuclein", "oxidative", "microglia", "calcium",
              "apoptosis", "ER stress", "unfolded protein", "ubiquitin")
pd_pattern <- paste(pd_terms, collapse = "|")

# Filter for PD-relevant and significant
pd_data <- all_data %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(pd_pattern, Description, ignore.case = TRUE))

# Get all clusters
all_clusters <- sort(unique(pd_data$cluster))
cat("Found", length(all_clusters), "clusters:", paste(all_clusters, collapse = ", "), "\n\n")

# Function to analyze pathways in a cluster
analyze_cluster_pathways <- function(data, cluster_id) {
  cat("Analyzing", cluster_id, "...\n")
  
  # Filter data for this cluster
  cluster_data <- data %>% filter(cluster == cluster_id)
  
  if (nrow(cluster_data) == 0) {
    return(list(
      cluster = cluster_id,
      top_pathways = data.frame(),
      method_breakdown = data.frame(
        cluster = cluster_id,
        n_mast = 0,
        n_mixscale = 0,
        n_genes = 0
      )
    ))
  }
  
  # Get top pathways
  top_pathways <- cluster_data %>%
    group_by(Description, enrichment_type) %>%
    summarise(
      n_genes = n_distinct(gene),
      genes = paste(unique(gene)[1:min(3, length(unique(gene)))], collapse = ", "),
      mean_neg_log_p = mean(-log10(p.adjust)),
      methods = paste(unique(method), collapse = "+"),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_neg_log_p)) %>%
    head(20)
  
  # Method breakdown
  method_stats <- data.frame(
    cluster = cluster_id,
    n_mast = sum(cluster_data$method == "MAST"),
    n_mixscale = sum(cluster_data$method == "MixScale"),
    n_genes = n_distinct(cluster_data$gene),
    n_pathways = n_distinct(cluster_data$Description)
  )
  
  return(list(
    cluster = cluster_id,
    top_pathways = top_pathways,
    method_breakdown = method_stats
  ))
}

# Analyze all clusters
cluster_results <- list()
for (cluster in all_clusters) {
  cluster_results[[cluster]] <- analyze_cluster_pathways(pd_data, cluster)
}

# Combine method breakdowns
all_cluster_stats <- bind_rows(lapply(cluster_results, function(x) x$method_breakdown))
write.csv(all_cluster_stats, file.path(output_dir, "cluster_method_breakdown.csv"), row.names = FALSE)

# Create cluster × pathway matrix for top pathways
# Get top 30 pathways across all clusters
all_top_pathways <- pd_data %>%
  group_by(Description) %>%
  summarise(
    total_occurrences = n(),
    n_clusters = n_distinct(cluster),
    mean_neg_log_p = mean(-log10(p.adjust)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_clusters), desc(mean_neg_log_p)) %>%
  head(30) %>%
  pull(Description)

# Create matrix
pathway_cluster_matrix <- pd_data %>%
  filter(Description %in% all_top_pathways) %>%
  group_by(Description, cluster) %>%
  summarise(
    score = mean(-log10(p.adjust)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = cluster, values_from = score, values_fill = 0)

# Convert to matrix for heatmap
mat <- as.matrix(pathway_cluster_matrix[,-1])
rownames(mat) <- substr(pathway_cluster_matrix$Description, 1, 60)

# Create heatmap
pdf(file.path(output_dir, "cluster_pathway_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  border_color = NA,
  cellwidth = 30,
  cellheight = 12,
  fontsize = 10,
  main = "PD Pathway Enrichment Across Clusters",
  angle_col = 0
)
dev.off()

# Analyze cluster-specific vs ubiquitous pathways
pathway_cluster_counts <- pd_data %>%
  group_by(Description) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    clusters = paste(unique(cluster), collapse = ", "),
    .groups = "drop"
  )

# Cluster-specific pathways (only in 1-2 clusters)
cluster_specific <- pathway_cluster_counts %>%
  filter(n_clusters <= 2) %>%
  arrange(desc(n_clusters))

# Ubiquitous pathways (in 8+ clusters)
ubiquitous <- pathway_cluster_counts %>%
  filter(n_clusters >= 8) %>%
  arrange(desc(n_clusters))

write.csv(cluster_specific, file.path(output_dir, "cluster_specific_pathways.csv"), row.names = FALSE)
write.csv(ubiquitous, file.path(output_dir, "ubiquitous_pathways.csv"), row.names = FALSE)

# Create visualization of cluster patterns
cluster_plot_data <- all_cluster_stats %>%
  pivot_longer(cols = c(n_mast, n_mixscale), 
               names_to = "Method", 
               values_to = "Count") %>%
  mutate(Method = case_when(
    Method == "n_mast" ~ "Mutation\niSCORE-PD",
    Method == "n_mixscale" ~ "CRISPRi\nPerturbation"
  ))

p1 <- ggplot(cluster_plot_data, aes(x = cluster, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4", 
                              "CRISPRi\nPerturbation" = "#ff7f0e")) +
  labs(
    title = "PD Pathway Distribution Across Clusters by Method",
    x = "Cluster",
    y = "Number of Enriched Pathways"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "cluster_method_distribution.pdf"), p1, width = 12, height = 8)

# Gene diversity per cluster
p2 <- ggplot(all_cluster_stats, aes(x = cluster, y = n_genes)) +
  geom_bar(stat = "identity", fill = "#2ca02c") +
  geom_text(aes(label = n_genes), vjust = -0.5) +
  labs(
    title = "Gene Diversity Across Clusters",
    subtitle = "Number of unique genes with PD pathway enrichment",
    x = "Cluster",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(file.path(output_dir, "cluster_gene_diversity.pdf"), p2, width = 10, height = 6)

# Identify cluster signatures
cat("\n=== CLUSTER SIGNATURE ANALYSIS ===\n")

# Find clusters with strongest enrichment
top_enriched_clusters <- all_cluster_stats %>%
  arrange(desc(n_pathways)) %>%
  head(5)

cat("\nClusters with most PD pathway enrichment:\n")
print(top_enriched_clusters)

# Method preference by cluster
method_preference <- all_cluster_stats %>%
  mutate(
    mast_ratio = n_mast / (n_mast + n_mixscale),
    preference = case_when(
      mast_ratio > 0.7 ~ "Mutation-dominant",
      mast_ratio < 0.3 ~ "CRISPRi-dominant",
      TRUE ~ "Balanced"
    )
  ) %>%
  select(cluster, n_mast, n_mixscale, mast_ratio, preference)

cat("\nMethod preference by cluster:\n")
print(method_preference)

# Save detailed cluster results
for (cluster in names(cluster_results)) {
  result <- cluster_results[[cluster]]
  if (nrow(result$top_pathways) > 0) {
    write.csv(result$top_pathways, 
              file.path(output_dir, paste0(cluster, "_top_pathways.csv")), 
              row.names = FALSE)
  }
}

cat("\nCluster-specific pathways found:", nrow(cluster_specific), "\n")
cat("Ubiquitous pathways (8+ clusters):", nrow(ubiquitous), "\n")

# Example cluster-specific pathways
cat("\nExample cluster-specific pathways:\n")
print(head(cluster_specific, 10))

cat("\nExample ubiquitous pathways:\n")
print(head(ubiquitous, 10))

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")