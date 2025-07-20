# Comprehensive Visualization Suite with Fixes
# Purpose: Create visualizations with gene harmonization and proper clustering
# Date: 2025-07-20

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
library(stringr)

# Helper functions
harmonize_gene_names <- function(gene_name) {
  gene_map <- c(
    "PRKN" = "PARK2",
    "VPS13C_A444P" = "VPS13C",
    "VPS13C_W395C" = "VPS13C", 
    "SNCA_A30P" = "SNCA",
    "SNCA_A53T" = "SNCA",
    "SNCA_Triplication" = "SNCA"
  )
  if (gene_name %in% names(gene_map)) {
    return(gene_map[gene_name])
  }
  return(gene_name)
}

natural_sort_clusters <- function(cluster_vec) {
  numeric_part <- as.numeric(gsub("cluster_", "", cluster_vec))
  cluster_vec[order(numeric_part)]
}

# Configuration
results_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
output_dir <- file.path(results_dir, "visualizations/comprehensive_fixed")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load harmonized data
harmonized_summary <- read.csv(file.path(results_dir, "by_gene_fixed/harmonized_gene_summary.csv"), stringsAsFactors = FALSE)
cluster_summary <- read.csv(file.path(results_dir, "by_cluster_fixed/cluster_summary_sorted.csv"), stringsAsFactors = FALSE)

# Load original data files with fixes
mast_top <- read.csv(file.path(results_dir, "mast_top_fast.csv"), stringsAsFactors = FALSE)
mixscale_top <- read.csv(file.path(results_dir, "mixscale_top_fast.csv"), stringsAsFactors = FALSE)
convergent_top <- read.csv(file.path(results_dir, "convergent_top_fast.csv"), stringsAsFactors = FALSE)

# Set theme
theme_set(theme_minimal(base_size = 14))

# ==============================================================================
# 1. HARMONIZED GENE PROFILES
# ==============================================================================
cat("Creating harmonized gene profiles...\n")

gene_profile_data <- harmonized_summary %>%
  select(Gene = gene, 
         `Mutation-Only` = n_mast_only,
         `CRISPRi-Only` = n_mixscale_only,
         `Convergent` = n_convergent) %>%
  pivot_longer(cols = -Gene, names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("Mutation-Only", "CRISPRi-Only", "Convergent")))

# Order genes by total pathways
gene_order <- harmonized_summary %>%
  arrange(desc(n_total)) %>%
  pull(gene)

gene_profile_data$Gene <- factor(gene_profile_data$Gene, levels = gene_order)

p_harmonized <- ggplot(gene_profile_data, aes(x = Gene, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "Harmonized Gene Signature Profiles",
       subtitle = "Gene variants merged (e.g., PRKN→PARK2, VPS13C variants→VPS13C, SNCA variants→SNCA)",
       y = "Number of Enriched Pathways") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 18),
        legend.position = "top") +
  coord_flip()

ggsave(file.path(output_dir, "01_harmonized_gene_profiles.pdf"), p_harmonized, width = 10, height = 10)

# ==============================================================================
# 2. CLUSTERED GENE HEATMAP
# ==============================================================================
cat("Creating clustered gene heatmap...\n")

# Prepare matrix
gene_matrix <- harmonized_summary %>%
  select(gene, `Mutation-Only` = n_mast_only, 
         `CRISPRi-Only` = n_mixscale_only, 
         `Convergent` = n_convergent)

gene_mat <- as.matrix(gene_matrix[,-1])
rownames(gene_mat) <- gene_matrix$gene

# Create clustered heatmap with annotations
pdf(file.path(output_dir, "02_gene_heatmap_clustered.pdf"), width = 10, height = 10)

# Add row annotations showing total pathways
row_annotation <- data.frame(
  Total = harmonized_summary$n_total,
  row.names = harmonized_summary$gene
)

pheatmap(
  gene_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "complete",
  display_numbers = TRUE,
  number_format = "%.0f",
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  border_color = "gray80",
  main = "PD Pathway Distribution - Hierarchically Clustered",
  annotation_row = row_annotation,
  annotation_colors = list(Total = colorRampPalette(c("lightblue", "darkblue"))(100)),
  fontsize = 10,
  cellwidth = 50,
  cellheight = 20
)

dev.off()

# ==============================================================================
# 3. CONVERGENCE ANALYSIS WITH HARMONIZED GENES
# ==============================================================================
cat("Creating convergence analysis...\n")

# Show which genes have strongest convergence
convergence_data <- harmonized_summary %>%
  mutate(
    convergence_ratio = n_convergent / n_total,
    convergence_strength = n_convergent * convergence_ratio
  ) %>%
  arrange(desc(convergence_strength))

p_convergence <- ggplot(convergence_data, aes(x = reorder(gene, convergence_strength), 
                                              y = convergence_strength)) +
  geom_segment(aes(xend = gene, yend = 0), size = 2, color = "gray50") +
  geom_point(aes(size = n_total, color = convergence_ratio), alpha = 0.8) +
  scale_color_viridis(name = "Convergence\nRatio", limits = c(0, 1)) +
  scale_size_continuous(name = "Total\nPathways", range = c(4, 12)) +
  coord_flip() +
  labs(title = "Gene Convergence Strength Analysis",
       subtitle = "Combining pathway count and convergence ratio",
       x = "",
       y = "Convergence Strength (count × ratio)") +
  theme(plot.title = element_text(face = "bold", size = 18))

ggsave(file.path(output_dir, "03_convergence_strength_harmonized.pdf"), p_convergence, width = 10, height = 8)

# ==============================================================================
# 4. CLUSTER ANALYSIS WITH NATURAL SORTING
# ==============================================================================
cat("Creating cluster analysis with natural sorting...\n")

# Ensure clusters are properly sorted
cluster_summary$cluster <- factor(cluster_summary$cluster, 
                                 levels = natural_sort_clusters(unique(cluster_summary$cluster)))

# Create cluster visualization
cluster_method_data <- cluster_summary %>%
  select(cluster, `Mutation\niSCORE-PD` = n_mast, `CRISPRi\nPerturbation` = n_mixscale) %>%
  pivot_longer(cols = -cluster, names_to = "Method", values_to = "Count")

p_clusters <- ggplot(cluster_method_data, aes(x = cluster, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4", 
                              "CRISPRi\nPerturbation" = "#ff7f0e")) +
  labs(title = "PD Pathway Distribution Across Clusters",
       subtitle = "Properly sorted: cluster_0 through cluster_14",
       x = "Cluster",
       y = "Number of Enriched Pathways") +
  theme(plot.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom")

ggsave(file.path(output_dir, "04_cluster_distribution_sorted.pdf"), p_clusters, width = 14, height = 8)

# ==============================================================================
# 5. SUMMARY STATISTICS TABLE
# ==============================================================================
cat("Creating summary statistics...\n")

# Calculate key statistics
total_genes_original <- length(unique(c(
  "ATP13A2", "DNAJC6", "FBXO7", "GBA", "LRRK2", "PARK7", "PINK1", 
  "PRKN", "PARK2", "SNCA", "SNCA_A30P", "SNCA_A53T", "SYNJ1", 
  "VPS13C", "VPS13C_A444P", "VPS13C_W395C"
)))

total_genes_harmonized <- nrow(harmonized_summary)

stats_summary <- data.frame(
  Metric = c(
    "Original gene variants",
    "Harmonized genes",
    "Total clusters",
    "Total pathways (all methods)",
    "Gene with most convergence",
    "Gene with highest ratio",
    "Most enriched cluster"
  ),
  Value = c(
    total_genes_original,
    total_genes_harmonized,
    length(unique(cluster_summary$cluster)),
    sum(harmonized_summary$n_total),
    convergence_data$gene[which.max(convergence_data$n_convergent)],
    convergence_data$gene[which.max(convergence_data$convergence_ratio)],
    cluster_summary$cluster[which.max(cluster_summary$n_pathways)]
  )
)

# Create visual table
p_stats <- ggplot(data = NULL) +
  theme_void() +
  annotate("text", x = 0.5, y = 0.95, label = "Analysis Summary", 
           size = 8, fontface = "bold") +
  annotate("text", x = 0.5, y = 0.85, 
           label = paste(capture.output(print(stats_summary, row.names = FALSE)), collapse = "\n"),
           size = 5, family = "mono") +
  xlim(0, 1) + ylim(0, 1)

ggsave(file.path(output_dir, "05_summary_statistics.pdf"), p_stats, width = 10, height = 8)

# ==============================================================================
# 6. VARIANT MAPPING VISUALIZATION
# ==============================================================================
cat("Creating variant mapping visualization...\n")

# Show how variants map to harmonized genes
variant_mapping <- data.frame(
  Original = c("PRKN", "VPS13C_A444P", "VPS13C_W395C", "SNCA_A30P", "SNCA_A53T"),
  Harmonized = c("PARK2", "VPS13C", "VPS13C", "SNCA", "SNCA"),
  Type = c("Name change", "Variant merge", "Variant merge", "Variant merge", "Variant merge")
)

# Add pathway counts
for (i in 1:nrow(variant_mapping)) {
  harm_gene <- variant_mapping$Harmonized[i]
  gene_data <- harmonized_summary[harmonized_summary$gene == harm_gene, ]
  if (nrow(gene_data) > 0) {
    variant_mapping$Pathways[i] <- gene_data$n_total[1]
  } else {
    variant_mapping$Pathways[i] <- 0
  }
}

p_mapping <- ggplot(variant_mapping, aes(x = Original, y = Harmonized)) +
  geom_point(aes(size = Pathways, color = Type), alpha = 0.7) +
  geom_text(aes(label = Pathways), vjust = -1.5) +
  scale_size_continuous(range = c(5, 15)) +
  scale_color_brewer(palette = "Set1") +
  labs(title = "Gene Variant Harmonization Mapping",
       subtitle = "How MAST variants map to CRISPRi gene names",
       x = "Original Gene Name (MAST)",
       y = "Harmonized Gene Name (CRISPRi)") +
  theme(plot.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(output_dir, "06_variant_mapping.pdf"), p_mapping, width = 10, height = 8)

cat("\nComprehensive fixed visualizations complete!\n")
cat("Created 6 improved figures in:", output_dir, "\n\n")

# Print summary
cat("Key improvements implemented:\n")
cat("1. Gene harmonization: PRKN→PARK2, VPS13C variants→VPS13C, SNCA variants→SNCA\n")
cat("2. Hierarchical clustering in heatmaps\n")
cat("3. Natural cluster sorting (0,1,2...9,10,11,12,13,14)\n")
cat("4. Clear variant mapping visualization\n")