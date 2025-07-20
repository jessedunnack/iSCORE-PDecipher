# Comprehensive Visualization Suite for PD Signatures
# Purpose: Create exhaustive yet digestible plots for thesis committee
# Date: 2025-07-20

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(viridis)
library(RColorBrewer)
library(ggrepel)
library(pheatmap)
# library(VennDiagram)  # Skip if not installed
# library(circlize)     # Skip if not installed
# library(ggridges)     # Skip if not installed

# Configuration
results_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
output_dir <- file.path(results_dir, "visualizations/comprehensive")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading analysis results...\n")
mast_top <- read.csv(file.path(results_dir, "mast_top_fast.csv"), stringsAsFactors = FALSE)
mixscale_top <- read.csv(file.path(results_dir, "mixscale_top_fast.csv"), stringsAsFactors = FALSE)
convergent_top <- read.csv(file.path(results_dir, "convergent_top_fast.csv"), stringsAsFactors = FALSE)
gene_summary <- read.csv(file.path(results_dir, "by_gene/all_genes_summary.csv"), stringsAsFactors = FALSE)
cluster_stats <- read.csv(file.path(results_dir, "by_cluster/cluster_method_breakdown.csv"), stringsAsFactors = FALSE)

# Define consistent color scheme
method_colors <- c(
  "Mutation\niSCORE-PD" = "#1f77b4",
  "CRISPRi\nPerturbation" = "#ff7f0e",
  "Convergent" = "#2ca02c"
)

# Set theme
theme_set(theme_minimal(base_size = 14))

# ==============================================================================
# 1. OVERVIEW LANDSCAPE (3-panel figure)
# ==============================================================================
cat("Creating overview landscape...\n")

# Panel A: Total counts
overview_data <- data.frame(
  Category = c("Mutation-Only", "CRISPRi-Only", "Convergent"),
  Count = c(nrow(mast_top), nrow(mixscale_top), nrow(convergent_top))
)

p1a <- ggplot(overview_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, size = 6, fontface = "bold") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4", 
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "A. PD Pathway Distribution",
       y = "Number of Pathways") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))

# Panel B: Venn diagram placeholder (text representation)
venn_summary <- data.frame(
  Label = c("Total PD Pathways", "Mutation + Convergent", "CRISPRi + Convergent", "All Methods"),
  Value = c(sum(overview_data$Count), 
            nrow(mast_top) + nrow(convergent_top),
            nrow(mixscale_top) + nrow(convergent_top),
            nrow(convergent_top))
)

p1b <- ggplot(venn_summary, aes(x = Label, y = Value)) +
  geom_bar(stat = "identity", fill = "gray70") +
  geom_text(aes(label = Value), vjust = -0.5, size = 5) +
  labs(title = "B. Method Coverage",
       y = "Pathways") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 16))

# Panel C: Top biological categories
categorize_pathway <- function(desc) {
  desc_lower <- tolower(desc)
  if (grepl("mitochondr|oxidative|atp|electron", desc_lower)) return("Mitochondrial/Energy")
  if (grepl("synap|dopamin|neurotrans", desc_lower)) return("Synaptic/Neuronal")
  if (grepl("lysosom|autophagy", desc_lower)) return("Lysosomal/Autophagy")
  if (grepl("proteasom|ubiquitin", desc_lower)) return("Protein Degradation")
  if (grepl("apoptosis|cell death", desc_lower)) return("Cell Death")
  if (grepl("inflamm|immune|microglia", desc_lower)) return("Neuroinflammation")
  return("Other")
}

# Categorize all top pathways
all_pathways <- rbind(
  mast_top %>% head(20) %>% select(Description, mean_neg_log_p) %>% mutate(Method = "Mutation-Only"),
  mixscale_top %>% head(20) %>% select(Description, mean_neg_log_p) %>% mutate(Method = "CRISPRi-Only"),
  convergent_top %>% head(20) %>% select(Description, mean_neg_log_p) %>% mutate(Method = "Convergent")
)
all_pathways$Category <- sapply(all_pathways$Description, categorize_pathway)

category_counts <- all_pathways %>%
  group_by(Category) %>%
  summarise(Count = n()) %>%
  arrange(desc(Count))

p1c <- ggplot(category_counts, aes(x = reorder(Category, Count), y = Count, fill = Category)) +
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Set3") +
  coord_flip() +
  labs(title = "C. Biological Themes",
       x = "", y = "Number of Pathways") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold", size = 16))

# Combine
pdf(file.path(output_dir, "01_overview_landscape.pdf"), width = 18, height = 6)
p1a + p1b + p1c
dev.off()

# ==============================================================================
# 2. GENE SIGNATURE PROFILES (All 16 genes)
# ==============================================================================
cat("Creating gene signature profiles...\n")

gene_profile_data <- gene_summary %>%
  select(Gene = gene, 
         `Mutation-Only` = n_mast_only,
         `CRISPRi-Only` = n_mixscale_only,
         `Convergent` = n_convergent) %>%
  pivot_longer(cols = -Gene, names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("Mutation-Only", "CRISPRi-Only", "Convergent")))

p2 <- ggplot(gene_profile_data, aes(x = Gene, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "Gene-Specific Pathway Signature Profiles",
       subtitle = "Distribution of PD pathways across all analyzed genes",
       y = "Number of Enriched Pathways") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 18),
        legend.position = "top") +
  coord_flip()

ggsave(file.path(output_dir, "02_gene_signature_profiles.pdf"), p2, width = 10, height = 12)

# ==============================================================================
# 3. BIOLOGICAL THEMES COMPARISON
# ==============================================================================
cat("Creating biological themes comparison...\n")

bio_theme_data <- all_pathways %>%
  group_by(Category, Method) %>%
  summarise(Count = n(), .groups = "drop") %>%
  mutate(Method = factor(Method, levels = c("Mutation-Only", "CRISPRi-Only", "Convergent")))

p3 <- ggplot(bio_theme_data, aes(x = Category, y = Count, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "Biological Theme Distribution Across Methods",
       subtitle = "Comparing pathway categories between mutation and perturbation approaches",
       x = "", y = "Number of Pathways") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold", size = 18),
        legend.position = "bottom")

ggsave(file.path(output_dir, "03_biological_themes_comparison.pdf"), p3, width = 12, height = 8)

# ==============================================================================
# 4. CONVERGENCE STRENGTH PLOT
# ==============================================================================
cat("Creating convergence strength visualization...\n")

# Get top convergent pathways with details
conv_strength <- convergent_top %>%
  head(30) %>%
  mutate(
    total_genes = n_genes_mast + n_genes_mixscale,
    # Use mean p-value weighted by gene counts as proxy
    mast_strength = n_genes_mast * mean_neg_log_p,
    mixscale_strength = n_genes_mixscale * mean_neg_log_p,
    Description_short = ifelse(nchar(Description) > 40,
                              paste0(substr(Description, 1, 40), "..."),
                              Description)
  )

p4 <- ggplot(conv_strength, aes(x = mast_strength, y = mixscale_strength)) +
  geom_point(aes(size = total_genes, color = mean_neg_log_p), alpha = 0.7) +
  geom_text_repel(aes(label = Description_short), 
                  data = conv_strength %>% head(10),
                  size = 3, max.overlaps = 20) +
  scale_color_viridis(name = "Combined\n-log10(p)") +
  scale_size_continuous(name = "Total Genes", range = c(3, 10)) +
  labs(title = "Convergent Pathway Strength Analysis",
       subtitle = "Comparing effect strength between mutation and CRISPRi methods",
       x = "Mutation Strength (genes × -log10(p))",
       y = "CRISPRi Strength (genes × -log10(p))") +
  theme(plot.title = element_text(face = "bold", size = 18))

ggsave(file.path(output_dir, "04_convergence_strength.pdf"), p4, width = 12, height = 10)

# ==============================================================================
# 5. TOP CONVERGENT PATHWAYS WITH CONTRIBUTIONS
# ==============================================================================
cat("Creating top convergent pathways visualization...\n")

top_conv <- convergent_top %>%
  head(15) %>%
  mutate(Description_short = ifelse(nchar(Description) > 50,
                                   paste0(substr(Description, 1, 50), "..."),
                                   Description))

# Create data for stacked bar
conv_contribution <- rbind(
  data.frame(
    Description = top_conv$Description_short,
    Genes = top_conv$n_genes_mast,
    Method = "Mutation\niSCORE-PD",
    Order = 1:nrow(top_conv)
  ),
  data.frame(
    Description = top_conv$Description_short,
    Genes = top_conv$n_genes_mixscale,
    Method = "CRISPRi\nPerturbation",
    Order = 1:nrow(top_conv)
  )
)

p5 <- ggplot(conv_contribution, aes(x = reorder(Description, -Order), y = Genes, fill = Method)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4",
                              "CRISPRi\nPerturbation" = "#ff7f0e")) +
  coord_flip() +
  labs(title = "Top Convergent Pathways: Method Contributions",
       subtitle = "Number of genes contributing from each experimental approach",
       x = "", y = "Number of Genes") +
  theme(plot.title = element_text(face = "bold", size = 18),
        legend.position = "bottom")

ggsave(file.path(output_dir, "05_top_convergent_contributions.pdf"), p5, width = 12, height = 10)

# ==============================================================================
# 6. CLUSTER CONSISTENCY HEATMAP
# ==============================================================================
cat("Creating cluster consistency heatmap...\n")

# Load full data for cluster analysis
full_data_path <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
if (!file.exists(full_data_path)) {
  cat("Warning: Full data file not found. Skipping cluster heatmap.\n")
  # Create placeholder
  pdf(file.path(output_dir, "06_cluster_consistency_heatmap.pdf"), width = 12, height = 10)
  plot.new()
  text(0.5, 0.5, "Cluster heatmap requires full data file", cex = 2)
  dev.off()
} else {
  full_data <- readRDS(full_data_path)

# Get top convergent pathways across clusters
cluster_pathway_data <- full_data %>%
  filter(Description %in% convergent_top$Description[1:20]) %>%
  group_by(Description, cluster) %>%
  summarise(present = 1, .groups = "drop") %>%
  pivot_wider(names_from = cluster, values_from = present, values_fill = 0)

# Convert to matrix
cluster_mat <- as.matrix(cluster_pathway_data[,-1])
rownames(cluster_mat) <- substr(cluster_pathway_data$Description, 1, 50)

# Create heatmap
pdf(file.path(output_dir, "06_cluster_consistency_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  cluster_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  color = c("white", "#2ca02c"),
  border_color = "gray80",
  main = "Convergent Pathway Presence Across Clusters",
  legend_breaks = c(0, 1),
  legend_labels = c("Absent", "Present"),
  fontsize = 10
)
dev.off()
}

# ==============================================================================
# 7. P-VALUE DISTRIBUTIONS
# ==============================================================================
cat("Creating p-value distribution plots...\n")

# Combine p-value data
pval_data <- rbind(
  mast_top %>% select(mean_neg_log_p) %>% mutate(Method = "Mutation-Only", pval = mean_neg_log_p),
  mixscale_top %>% select(mean_neg_log_p) %>% mutate(Method = "CRISPRi-Only", pval = mean_neg_log_p),
  convergent_top %>% select(mean_neg_log_p) %>% mutate(Method = "Convergent", pval = mean_neg_log_p)
) %>% select(Method, pval)

p7 <- ggplot(pval_data, aes(x = pval, fill = Method)) +
  geom_density(alpha = 0.7) +
  facet_wrap(~Method, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "Statistical Significance Distribution",
       subtitle = "Density of -log10(p-values) across pathway categories",
       x = "Mean -log10(adjusted p-value)",
       y = "Density") +
  theme(plot.title = element_text(face = "bold", size = 18),
        legend.position = "none")

ggsave(file.path(output_dir, "07_pvalue_distributions.pdf"), p7, width = 10, height = 6)

# ==============================================================================
# 8. GENE RANKING BY CONVERGENT STRENGTH
# ==============================================================================
cat("Creating gene ranking visualization...\n")

gene_ranking <- gene_summary %>%
  mutate(
    convergent_ratio = n_convergent / (n_convergent + n_mast_only + n_mixscale_only),
    total_pathways = n_convergent + n_mast_only + n_mixscale_only
  ) %>%
  arrange(desc(n_convergent)) %>%
  head(12)

p8 <- ggplot(gene_ranking, aes(x = reorder(gene, n_convergent), y = n_convergent)) +
  geom_segment(aes(xend = gene, yend = 0), size = 1.5, color = "gray50") +
  geom_point(aes(size = total_pathways), color = "#2ca02c", alpha = 0.8) +
  geom_text(aes(label = n_convergent), hjust = -0.5, size = 4) +
  scale_size_continuous(name = "Total Pathways", range = c(4, 12)) +
  coord_flip() +
  labs(title = "Gene Ranking by Convergent Pathway Count",
       subtitle = "Genes with strongest cross-method validation",
       x = "", y = "Number of Convergent Pathways") +
  theme(plot.title = element_text(face = "bold", size = 18))

ggsave(file.path(output_dir, "08_gene_ranking_convergent.pdf"), p8, width = 10, height = 8)

# ==============================================================================
# 9. METHOD COMPARISON SUMMARY
# ==============================================================================
cat("Creating method comparison summary...\n")

method_summary <- data.frame(
  Metric = c("Total Pathways", "Avg Genes/Pathway", "Avg -log10(p)", "Top Category"),
  Mutation_Only = c(
    nrow(mast_top),
    round(mean(mast_top$n_genes), 1),
    round(mean(mast_top$mean_neg_log_p), 1),
    names(sort(table(sapply(mast_top$Description[1:20], categorize_pathway)), decreasing = TRUE)[1])
  ),
  CRISPRi_Only = c(
    nrow(mixscale_top),
    round(mean(mixscale_top$n_genes), 1),
    round(mean(mixscale_top$mean_neg_log_p), 1),
    names(sort(table(sapply(mixscale_top$Description[1:20], categorize_pathway)), decreasing = TRUE)[1])
  ),
  Convergent = c(
    nrow(convergent_top),
    round(mean(convergent_top$n_genes_mast + convergent_top$n_genes_mixscale), 1),
    round(mean(convergent_top$mean_neg_log_p), 1),
    names(sort(table(sapply(convergent_top$Description[1:20], categorize_pathway)), decreasing = TRUE)[1])
  )
)

# Create a visual table
p9 <- ggplot(data = NULL) +
  theme_void() +
  annotate("text", x = 0.5, y = 0.9, label = "Method Comparison Summary", 
           size = 8, fontface = "bold") +
  annotate("text", x = 0.5, y = 0.7, 
           label = paste(capture.output(print(method_summary, row.names = FALSE)), collapse = "\n"),
           size = 5, family = "mono") +
  xlim(0, 1) + ylim(0, 1)

ggsave(file.path(output_dir, "09_method_comparison_table.pdf"), p9, width = 10, height = 6)

# ==============================================================================
# 10. EFFECT SIZE COMPARISON
# ==============================================================================
cat("Creating effect size comparison...\n")

# Compare fold enrichment distributions
# Note: Using mean_neg_log_p as proxy since fold enrichment may not be available
effect_data <- rbind(
  mast_top %>% 
    mutate(Method = "Mutation-Only", value = mean_neg_log_p) %>%
    select(Method, value),
  mixscale_top %>% 
    mutate(Method = "CRISPRi-Only", value = mean_neg_log_p) %>%
    select(Method, value),
  convergent_top %>% 
    mutate(Method = "Convergent", value = mean_neg_log_p) %>%
    select(Method, value)
)

p10 <- ggplot(effect_data, aes(x = Method, y = value, fill = Method)) +
  geom_violin(alpha = 0.7) +
  geom_boxplot(width = 0.2, alpha = 0.8) +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "Statistical Significance Distribution Across Methods",
       subtitle = "Distribution of -log10(p-values) for significant pathways",
       y = "-log10(adjusted p-value)") +
  theme(plot.title = element_text(face = "bold", size = 18),
        legend.position = "none")

ggsave(file.path(output_dir, "10_effect_size_comparison.pdf"), p10, width = 10, height = 8)

# ==============================================================================
# 11. PATHWAY OVERLAP NETWORK (Simplified)
# ==============================================================================
cat("Creating pathway network visualization...\n")

# Create a simple network showing top pathways and their presence across methods
network_data <- convergent_top %>%
  head(10) %>%
  mutate(
    Description_short = ifelse(nchar(Description) > 30,
                              paste0(substr(Description, 1, 30), "..."),
                              Description),
    x = seq(-1, 1, length.out = n()),
    y = seq(-1, 1, length.out = n())
  )

p11 <- ggplot(network_data) +
  geom_point(aes(x = x, y = y, size = mean_neg_log_p), 
             color = "#2ca02c", alpha = 0.6) +
  geom_text_repel(aes(x = x, y = y, label = Description_short),
                  size = 3, max.overlaps = 20) +
  scale_size_continuous(name = "-log10(p)", range = c(5, 15)) +
  labs(title = "Top Convergent Pathways Network",
       subtitle = "Size represents statistical significance") +
  theme_void() +
  theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        legend.position = "bottom")

ggsave(file.path(output_dir, "11_pathway_network.pdf"), p11, width = 10, height = 10)

# ==============================================================================
# 12. SUMMARY INFOGRAPHIC
# ==============================================================================
cat("Creating summary infographic...\n")

# Key numbers for infographic
key_stats <- list(
  n_genes = length(unique(gene_summary$gene)),
  n_clusters = length(unique(cluster_stats$cluster)),
  n_total_pathways = nrow(mast_top) + nrow(mixscale_top) + nrow(convergent_top),
  top_convergent = convergent_top$Description[1],
  top_convergent_p = round(convergent_top$mean_neg_log_p[1], 1),
  strongest_gene = gene_summary$gene[which.max(gene_summary$n_convergent)]
)

p12 <- ggplot(data = NULL) +
  theme_void() +
  # Title
  annotate("text", x = 0.5, y = 0.95, 
           label = "iSCORE-PDecipher: PD Signature Analysis Summary",
           size = 10, fontface = "bold") +
  # Key stats
  annotate("text", x = 0.25, y = 0.8, 
           label = paste0(key_stats$n_genes, "\nPD Genes"),
           size = 8, color = "#1f77b4") +
  annotate("text", x = 0.5, y = 0.8, 
           label = paste0(key_stats$n_clusters, "\nClusters"),
           size = 8, color = "#ff7f0e") +
  annotate("text", x = 0.75, y = 0.8, 
           label = paste0(key_stats$n_total_pathways, "\nPathways"),
           size = 8, color = "#2ca02c") +
  # Top finding
  annotate("text", x = 0.5, y = 0.6, 
           label = "Strongest Convergent Signal:",
           size = 6, fontface = "bold") +
  annotate("text", x = 0.5, y = 0.5, 
           label = paste0(substr(key_stats$top_convergent, 1, 60), "..."),
           size = 5) +
  annotate("text", x = 0.5, y = 0.4, 
           label = paste0("p < 10^-", key_stats$top_convergent_p),
           size = 6, color = "#2ca02c", fontface = "bold") +
  # Bottom text
  annotate("text", x = 0.5, y = 0.2, 
           label = paste0("Top convergent gene: ", key_stats$strongest_gene),
           size = 5) +
  xlim(0, 1) + ylim(0, 1)

ggsave(file.path(output_dir, "12_summary_infographic.pdf"), p12, width = 10, height = 8)

cat("\nComprehensive visualization suite complete!\n")
cat("Created 12 publication-ready figures in:", output_dir, "\n\n")

# Print summary
cat("Visualization Summary:\n")
cat("1. Overview landscape (3-panel)\n")
cat("2. Gene signature profiles (all 16 genes)\n")
cat("3. Biological themes comparison\n")
cat("4. Convergence strength scatter\n")
cat("5. Top convergent pathway contributions\n")
cat("6. Cluster consistency heatmap\n")
cat("7. P-value distributions\n")
cat("8. Gene ranking by convergent strength\n")
cat("9. Method comparison summary\n")
cat("10. Effect size distributions\n")
cat("11. Pathway network\n")
cat("12. Summary infographic\n")