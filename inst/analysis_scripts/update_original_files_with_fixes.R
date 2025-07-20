# Update Original Files with Fixes
# Purpose: Replace original files with fixed versions (gene harmonization, clustering, sorting)
# Date: 2025-07-20

library(dplyr)
library(ggplot2)
library(tidyr)
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

# ==============================================================================
# 1. UPDATE BY-GENE ANALYSIS
# ==============================================================================
cat("Updating by-gene analysis with harmonization...\n")

# Configuration - use ORIGINAL directories
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_gene"

# Load data
all_data <- readRDS(data_file)

# Add gene column if needed
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Add harmonized gene column
all_data$gene_harmonized <- sapply(all_data$gene, harmonize_gene_names)

# PD keyword filtering
pd_terms <- c("mitochondr", "lysosom", "autophagy", "synap", "dopamin", 
              "proteasom", "synuclein", "oxidative", "microglia", "calcium",
              "apoptosis", "ER stress", "unfolded protein", "ubiquitin")
pd_pattern <- paste(pd_terms, collapse = "|")

pd_data <- all_data %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(pd_pattern, Description, ignore.case = TRUE))

# Get unique harmonized genes
unique_harmonized_genes <- unique(pd_data$gene_harmonized)
cat("Found", length(unique_harmonized_genes), "unique genes after harmonization\n")

# Analyze each harmonized gene group
gene_results_harmonized <- list()

for (harm_gene in unique_harmonized_genes) {
  cat("  Processing:", harm_gene, "\n")
  
  # Get all variants for this harmonized gene
  gene_variants <- unique(pd_data$gene[pd_data$gene_harmonized == harm_gene])
  
  # Filter data for this harmonized gene group
  gene_data <- pd_data %>% filter(gene_harmonized == harm_gene)
  
  # Separate by method
  mast_data <- gene_data %>% filter(method == "MAST")
  mixscale_data <- gene_data %>% filter(method == "MixScale")
  
  # Find unique and convergent terms
  mast_terms <- unique(mast_data$Description)
  mixscale_terms <- unique(mixscale_data$Description)
  convergent_terms <- intersect(mast_terms, mixscale_terms)
  
  # Get top pathways for each category
  get_top_pathways <- function(data, terms, n_top = 10) {
    if (length(terms) == 0) return(data.frame())
    
    data %>%
      filter(Description %in% terms) %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        n_clusters = n_distinct(cluster),
        mean_neg_log_p = mean(-log10(p.adjust)),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_neg_log_p)) %>%
      head(n_top)
  }
  
  # Store results
  mast_only_terms <- setdiff(mast_terms, mixscale_terms)
  mixscale_only_terms <- setdiff(mixscale_terms, mast_terms)
  
  # Save CSV files for this gene
  if (length(mast_only_terms) > 0) {
    mast_only_top <- get_top_pathways(mast_data, mast_only_terms, 50)
    write.csv(mast_only_top, 
              file.path(output_dir, paste0(harm_gene, "_mast_only_pathways.csv")), 
              row.names = FALSE)
  }
  
  if (length(mixscale_only_terms) > 0) {
    mixscale_only_top <- get_top_pathways(mixscale_data, mixscale_only_terms, 50)
    write.csv(mixscale_only_top, 
              file.path(output_dir, paste0(harm_gene, "_mixscale_only_pathways.csv")), 
              row.names = FALSE)
  }
  
  if (length(convergent_terms) > 0) {
    convergent_data <- pd_data %>%
      filter(gene_harmonized == harm_gene, Description %in% convergent_terms) %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        n_clusters_mast = n_distinct(cluster[method == "MAST"]),
        n_clusters_mixscale = n_distinct(cluster[method == "MixScale"]),
        mean_neg_log_p = mean(-log10(p.adjust)),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_neg_log_p)) %>%
      head(50)
    
    write.csv(convergent_data, 
              file.path(output_dir, paste0(harm_gene, "_convergent_pathways.csv")), 
              row.names = FALSE)
  }
  
  # Store summary
  gene_results_harmonized[[harm_gene]] <- list(
    gene = harm_gene,
    variants = gene_variants,
    n_mast_only = length(mast_only_terms),
    n_mixscale_only = length(mixscale_only_terms),
    n_convergent = length(convergent_terms),
    n_total = length(unique(c(mast_terms, mixscale_terms)))
  )
}

# Create updated summary table
summary_data <- bind_rows(lapply(gene_results_harmonized, function(x) {
  data.frame(
    gene = x$gene,
    n_mast_only = x$n_mast_only,
    n_mixscale_only = x$n_mixscale_only,
    n_convergent = x$n_convergent,
    n_total = x$n_total,
    top_mast_only = "See CSV files",
    top_mixscale_only = "See CSV files",
    top_convergent = "See CSV files",
    stringsAsFactors = FALSE
  )
}))

# Replace the original summary file
write.csv(summary_data, file.path(output_dir, "all_genes_summary.csv"), row.names = FALSE)

# Update individual gene plots
plots_dir <- file.path(output_dir, "plots")

for (harm_gene in names(gene_results_harmonized)) {
  result <- gene_results_harmonized[[harm_gene]]
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Category = c("Mutation\niSCORE-PD\nOnly", "CRISPRi\nPerturbation\nOnly", "Convergent"),
    Count = c(result$n_mast_only, result$n_mixscale_only, result$n_convergent),
    stringsAsFactors = FALSE
  )
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c(
      "Mutation\niSCORE-PD\nOnly" = "#1f77b4",
      "CRISPRi\nPerturbation\nOnly" = "#ff7f0e",
      "Convergent" = "#2ca02c"
    )) +
    labs(
      title = paste("PD Pathway Signatures for", harm_gene),
      subtitle = paste("Includes:", paste(result$variants, collapse = ", "), 
                      "\nTotal pathways:", result$n_total),
      y = "Number of Pathways"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(size = 12)
    ) +
    ylim(0, max(plot_data$Count) * 1.2)
  
  # Save with harmonized gene name
  ggsave(file.path(plots_dir, paste0(harm_gene, "_signature_summary.pdf")), 
         p, width = 8, height = 6)
}

# Update the all genes heatmap with clustering
cat("Updating all genes heatmap with clustering...\n")

# Create heatmap data
heatmap_data <- summary_data %>%
  select(gene, `Mutation Only` = n_mast_only, 
         `CRISPRi Only` = n_mixscale_only, 
         `Convergent` = n_convergent)

# Convert to matrix
heatmap_mat <- as.matrix(heatmap_data[,-1])
rownames(heatmap_mat) <- heatmap_data$gene

# Create CLUSTERED heatmap
pdf(file.path(output_dir, "all_genes_pathway_heatmap.pdf"), width = 12, height = 8)

# Create heatmap with clustering
pheatmap(
  heatmap_mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.0f",
  number_color = "white",
  color = colorRampPalette(c("lightblue", "yellow", "orange", "red"))(50),
  border_color = NA,
  main = "PD Pathway Distribution Across All Genes (Hierarchically Clustered)",
  fontsize = 10,
  cellwidth = 60,
  cellheight = 20
)

dev.off()

# ==============================================================================
# 2. UPDATE CLUSTER ANALYSIS WITH NATURAL SORTING
# ==============================================================================
cat("\nUpdating cluster analysis with natural sorting...\n")

cluster_output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_cluster"

# Get all clusters and sort naturally
all_clusters <- unique(pd_data$cluster)
all_clusters_sorted <- natural_sort_clusters(all_clusters)

cat("Clusters in natural order:", paste(all_clusters_sorted, collapse = ", "), "\n")

# Analyze clusters
cluster_results <- list()
for (cluster in all_clusters_sorted) {
  cluster_data <- pd_data %>% filter(cluster == !!cluster)
  
  if (nrow(cluster_data) > 0) {
    # Get top pathways for this cluster
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
    
    # Save top pathways
    write.csv(top_pathways, 
              file.path(cluster_output_dir, paste0(cluster, "_top_pathways.csv")), 
              row.names = FALSE)
    
    # Store summary
    cluster_results[[cluster]] <- data.frame(
      cluster = cluster,
      n_pathways = n_distinct(cluster_data$Description),
      n_mast = sum(cluster_data$method == "MAST"),
      n_mixscale = sum(cluster_data$method == "MixScale"),
      n_genes = n_distinct(cluster_data$gene),
      stringsAsFactors = FALSE
    )
  }
}

# Combine and sort results
cluster_stats <- bind_rows(cluster_results)
cluster_stats$cluster <- factor(cluster_stats$cluster, levels = all_clusters_sorted)
cluster_stats <- cluster_stats %>% arrange(cluster)

# Replace original file
write.csv(cluster_stats, file.path(cluster_output_dir, "cluster_method_breakdown.csv"), row.names = FALSE)

# Update cluster visualizations
cluster_plot_data <- cluster_stats %>%
  pivot_longer(cols = c(n_mast, n_mixscale), 
               names_to = "Method", 
               values_to = "Count") %>%
  mutate(Method = case_when(
    Method == "n_mast" ~ "Mutation\niSCORE-PD",
    Method == "n_mixscale" ~ "CRISPRi\nPerturbation"
  ))

p_cluster <- ggplot(cluster_plot_data, aes(x = cluster, y = Count, fill = Method)) +
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

ggsave(file.path(cluster_output_dir, "cluster_method_distribution.pdf"), p_cluster, width = 12, height = 8)

# Update cluster heatmap with natural sorting
top_pathways_all <- pd_data %>%
  group_by(Description) %>%
  summarise(
    n_clusters = n_distinct(cluster),
    mean_neg_log_p = mean(-log10(p.adjust)),
    .groups = "drop"
  ) %>%
  arrange(desc(n_clusters), desc(mean_neg_log_p)) %>%
  head(30) %>%
  pull(Description)

# Create cluster x pathway matrix
cluster_pathway_data <- pd_data %>%
  filter(Description %in% top_pathways_all) %>%
  group_by(Description, cluster) %>%
  summarise(
    score = mean(-log10(p.adjust)),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = cluster, values_from = score, values_fill = 0)

# Ensure clusters are in natural order
cluster_order <- all_clusters_sorted[all_clusters_sorted %in% names(cluster_pathway_data)]
cluster_pathway_data <- cluster_pathway_data %>%
  select(Description, all_of(cluster_order))

# Convert to matrix
mat <- as.matrix(cluster_pathway_data[,-1])
rownames(mat) <- substr(cluster_pathway_data$Description, 1, 60)

# Create clustered heatmap
pdf(file.path(cluster_output_dir, "cluster_pathway_heatmap.pdf"), width = 12, height = 10)
pheatmap(
  mat,
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # Keep natural cluster order
  color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
  border_color = NA,
  cellwidth = 30,
  cellheight = 12,
  fontsize = 10,
  main = "PD Pathway Enrichment Across Clusters (Pathways Clustered)",
  angle_col = 0
)
dev.off()

cat("\n=== ALL FILES UPDATED IN PLACE ===\n")
cat("Updated files in:\n")
cat("- By gene:", output_dir, "\n")
cat("- By cluster:", cluster_output_dir, "\n")
cat("\nKey improvements:\n")
cat("1. Gene harmonization applied (PRKN→PARK2, etc.)\n")
cat("2. All heatmaps now clustered\n")
cat("3. Clusters sorted naturally\n")