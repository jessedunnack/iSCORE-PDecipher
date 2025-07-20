# PD Signature Comprehensive Visualization
# Purpose: Create comprehensive multi-panel visualizations combining all analyses
# Date: 2025-07-19

library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(gridExtra)
library(RColorBrewer)

# Helper function for word wrapping
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) {
    if (is.na(x)) return('')
    stringr::str_wrap(x, width = width)
  })
}


# Configuration
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
results_base <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
output_dir <- file.path(results_base, "comprehensive")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load pre-computed results
cat("Loading analysis results...\n")

# Gene summaries
gene_summary <- read.csv(file.path(results_base, "by_gene/all_genes_summary.csv"), stringsAsFactors = FALSE)

# Cluster summaries
cluster_stats <- read.csv(file.path(results_base, "by_cluster/cluster_method_breakdown.csv"), stringsAsFactors = FALSE)

# Top signatures
mast_top <- read.csv(file.path(results_base, "mast_top_fast.csv"), stringsAsFactors = FALSE)
mixscale_top <- read.csv(file.path(results_base, "mixscale_top_fast.csv"), stringsAsFactors = FALSE)
convergent_top <- read.csv(file.path(results_base, "convergent_top_fast.csv"), stringsAsFactors = FALSE)

# Load full data for detailed analysis
all_data <- readRDS(data_file)
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Filter for PD-relevant
pd_terms <- c("mitochondr", "lysosom", "autophagy", "synap", "dopamin", 
              "proteasom", "synuclein", "oxidative", "microglia", "calcium",
              "apoptosis", "ER stress", "unfolded protein", "ubiquitin")
pd_pattern <- paste(pd_terms, collapse = "|")

pd_data <- all_data %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(pd_pattern, Description, ignore.case = TRUE))

# ==============================================================================
# 1. Gene × Pathway Matrices for Each Category
# ==============================================================================

create_gene_pathway_matrix <- function(data, pathways, method_filter = NULL, title = "") {
  # Filter by method if specified
  if (!is.null(method_filter)) {
    data <- data %>% filter(method == method_filter)
  }
  
  # Create matrix
  matrix_data <- data %>%
    filter(Description %in% pathways) %>%
    group_by(gene, Description) %>%
    summarise(
      score = mean(-log10(p.adjust)),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = gene, values_from = score, values_fill = 0)
  
  if (nrow(matrix_data) == 0) return(NULL)
  
  # Convert to matrix
  mat <- as.matrix(matrix_data[,-1])
  rownames(mat) <- wrap_text(matrix_data$Description, width = 50)
  
  # Create heatmap
  if (ncol(mat) > 1 && nrow(mat) > 1) {
    p <- pheatmap(
      mat,
      cluster_rows = TRUE,
      cluster_cols = TRUE,
      color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
      border_color = NA,
      main = title,
      silent = TRUE
    )
    return(p$gtable)
  }
  return(NULL)
}

# Get top pathways for each category
top_mast_pathways <- head(mast_top$Description, 20)
top_mixscale_pathways <- head(mixscale_top$Description, 20)
top_convergent_pathways <- head(convergent_top$Description, 20)

# Create matrices
cat("Creating gene × pathway matrices...\n")

pdf(file.path(output_dir, "gene_pathway_matrices.pdf"), width = 14, height = 10)

# Mutation-only matrix
p1 <- create_gene_pathway_matrix(pd_data, top_mast_pathways, "MAST", 
                                "Top Mutation - iSCORE-PD Only Pathways Across Genes")
if (!is.null(p1)) grid.arrange(p1)

# CRISPRi-only matrix
p2 <- create_gene_pathway_matrix(pd_data, top_mixscale_pathways, "MixScale",
                                "Top CRISPRi Perturbation Only Pathways Across Genes")
if (!is.null(p2)) grid.arrange(p2)

# Convergent matrix
p3 <- create_gene_pathway_matrix(pd_data, top_convergent_pathways, NULL,
                                "Top Convergent Pathways Across Genes")
if (!is.null(p3)) grid.arrange(p3)

dev.off()

# ==============================================================================
# 2. Comprehensive Summary Figure
# ==============================================================================

cat("Creating comprehensive summary figure...\n")

# Panel A: Gene signature distribution
panel_a_data <- gene_summary %>%
  select(gene, `Mutation Only` = n_mast_only, 
         `CRISPRi Only` = n_mixscale_only, 
         `Convergent` = n_convergent) %>%
  pivot_longer(cols = -gene, names_to = "Category", values_to = "Count") %>%
  mutate(Category = factor(Category, levels = c("Mutation Only", "CRISPRi Only", "Convergent")))

panel_a <- ggplot(panel_a_data, aes(x = gene, y = Count, fill = Category)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("Mutation Only" = "#1f77b4", 
                              "CRISPRi Only" = "#ff7f0e", 
                              "Convergent" = "#2ca02c")) +
  labs(title = "A. PD Pathway Distribution by Gene", x = "", y = "Number of Pathways") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    plot.title = element_text(face = "bold")
  )

# Panel B: Cluster enrichment patterns
panel_b <- ggplot(cluster_stats, aes(x = cluster, y = n_pathways)) +
  geom_bar(stat = "identity", fill = "#d62728") +
  geom_text(aes(label = n_pathways), vjust = -0.5, size = 3) +
  labs(title = "B. Total PD Pathways by Cluster", x = "Cluster", y = "Number of Pathways") +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

# Panel C: Top convergent pathways
panel_c_data <- convergent_top %>%
  head(10) %>%
  mutate(
    total_genes = n_genes_mast + n_genes_mixscale,
    Description_short = ifelse(nchar(Description) > 40,
                              paste0(wrap_text(Description, width = 40), "..."),
                              Description)
  )

panel_c <- ggplot(panel_c_data, aes(x = reorder(Description_short, total_genes), y = total_genes)) +
  geom_bar(stat = "identity", fill = "#2ca02c") +
  coord_flip() +
  labs(title = "C. Top Convergent Pathways", x = "", y = "Total Genes") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold"))

# Panel D: Method comparison
method_comp_data <- data.frame(
  Method = c("Mutation\niSCORE-PD", "CRISPRi\nPerturbation"),
  Unique_Pathways = c(nrow(mast_top), nrow(mixscale_top)),
  Avg_Genes = c(mean(mast_top$n_genes), mean(mixscale_top$n_genes))
)

panel_d1 <- ggplot(method_comp_data, aes(x = Method, y = Unique_Pathways, fill = Method)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4", 
                              "CRISPRi\nPerturbation" = "#ff7f0e")) +
  labs(title = "D. Method Comparison", y = "Unique Pathways") +
  theme_minimal(base_size = 10) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

# Combine all panels
pdf(file.path(output_dir, "comprehensive_summary_figure.pdf"), width = 16, height = 12)
grid.arrange(panel_a, panel_b, panel_c, panel_d1, ncol = 2)
dev.off()

# ==============================================================================
# 3. Gene Trend Analysis
# ==============================================================================

cat("Analyzing trends across genes...\n")

# Identify genes with different patterns
gene_patterns <- gene_summary %>%
  mutate(
    pattern = case_when(
      n_mast_only > n_mixscale_only & n_mast_only > n_convergent ~ "Mutation-dominant",
      n_mixscale_only > n_mast_only & n_mixscale_only > n_convergent ~ "CRISPRi-dominant",
      n_convergent > n_mast_only & n_convergent > n_mixscale_only ~ "Convergent-dominant",
      TRUE ~ "Balanced"
    )
  )

# Create pattern visualization
pattern_summary <- gene_patterns %>%
  group_by(pattern) %>%
  summarise(
    n_genes = n(),
    genes = paste(gene, collapse = ", "),
    .groups = "drop"
  )

p_patterns <- ggplot(pattern_summary, aes(x = pattern, y = n_genes, fill = pattern)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = n_genes), vjust = -0.5) +
  scale_fill_brewer(palette = "Set1") +
  labs(
    title = "Gene Classification by Dominant Signature Pattern",
    x = "Pattern",
    y = "Number of Genes"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "gene_pattern_classification.pdf"), p_patterns, width = 10, height = 6)

# ==============================================================================
# 4. Pathway Category Analysis
# ==============================================================================

# Categorize pathways
categorize_pathway <- function(desc) {
  desc_lower <- tolower(desc)
  if (grepl("mitochondr|oxidative|atp|electron", desc_lower)) return("Mitochondrial")
  if (grepl("synap|dopamin|neurotrans", desc_lower)) return("Synaptic")
  if (grepl("lysosom|autophagy|endocyt", desc_lower)) return("Lysosomal")
  if (grepl("proteasom|ubiquitin|protein degradation", desc_lower)) return("Protein Degradation")
  if (grepl("apoptosis|cell death", desc_lower)) return("Cell Death")
  return("Other")
}

# Apply categorization
pd_data$category <- sapply(pd_data$Description, categorize_pathway)

# Category distribution by method
category_method <- pd_data %>%
  group_by(category, method) %>%
  summarise(
    n_pathways = n_distinct(Description),
    .groups = "drop"
  ) %>%
  mutate(method = ifelse(method == "MAST", "Mutation\niSCORE-PD", "CRISPRi\nPerturbation"))

p_categories <- ggplot(category_method, aes(x = category, y = n_pathways, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4", 
                              "CRISPRi\nPerturbation" = "#ff7f0e")) +
  labs(
    title = "PD Pathway Categories by Method",
    x = "Pathway Category",
    y = "Number of Unique Pathways"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(output_dir, "pathway_categories_by_method.pdf"), p_categories, width = 12, height = 8)

# ==============================================================================
# 5. Summary Statistics
# ==============================================================================

cat("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")

# Overall statistics
total_genes <- length(unique(pd_data$gene))
total_clusters <- length(unique(pd_data$cluster))
total_pathways <- length(unique(pd_data$Description))

cat("\nDataset Overview:\n")
cat("- Total PD genes analyzed:", total_genes, "\n")
cat("- Total clusters:", total_clusters, "\n")
cat("- Total PD-relevant pathways:", total_pathways, "\n")

# Method breakdown
method_summary <- pd_data %>%
  group_by(method) %>%
  summarise(
    n_enrichments = n(),
    n_pathways = n_distinct(Description),
    n_genes = n_distinct(gene),
    n_clusters = n_distinct(cluster),
    .groups = "drop"
  )

cat("\nMethod Summary:\n")
print(method_summary)

# Save comprehensive summary
summary_report <- list(
  overview = data.frame(
    total_genes = total_genes,
    total_clusters = total_clusters,
    total_pathways = total_pathways,
    mutation_only_pathways = nrow(mast_top),
    crispri_only_pathways = nrow(mixscale_top),
    convergent_pathways = nrow(convergent_top)
  ),
  gene_patterns = pattern_summary,
  method_stats = method_summary,
  top_convergent = head(convergent_top, 5)
)

saveRDS(summary_report, file.path(output_dir, "comprehensive_summary.rds"))

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")

# Create a final summary table for the manuscript
manuscript_table <- data.frame(
  Category = c("Total PD Genes", "Total Clusters", "Mutation-only Pathways", 
               "CRISPRi-only Pathways", "Convergent Pathways", 
               "Strongest Convergent Signal"),
  Value = c(
    total_genes,
    total_clusters,
    nrow(mast_top),
    nrow(mixscale_top),
    nrow(convergent_top),
    paste0(convergent_top$Description[1], " (p<1e-", 
           round(convergent_top$mean_neg_log_p[1]), ")")
  )
)

write.csv(manuscript_table, file.path(output_dir, "manuscript_summary_table.csv"), row.names = FALSE)

cat("\nManuscript summary table created!\n")
