# Fix Gene Harmonization, Clustering, and Sorting Issues
# Purpose: Comprehensive fixes for visualization improvements
# Date: 2025-07-20

library(dplyr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(stringr)

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

# Gene harmonization function
harmonize_gene_names <- function(gene_name) {
  # Map MAST names to CRISPRi names for convergent analysis
  gene_map <- c(
    "PRKN" = "PARK2",
    "VPS13C_A444P" = "VPS13C",
    "VPS13C_W395C" = "VPS13C", 
    "SNCA_A30P" = "SNCA",
    "SNCA_A53T" = "SNCA",
    "SNCA_Triplication" = "SNCA"
  )
  
  # Return mapped name if exists, otherwise return original
  if (gene_name %in% names(gene_map)) {
    return(gene_map[gene_name])
  }
  return(gene_name)
}

# Natural sorting function for clusters
natural_sort_clusters <- function(cluster_vec) {
  # Extract numeric part from cluster names
  numeric_part <- as.numeric(gsub("cluster_", "", cluster_vec))
  # Return sorted vector
  cluster_vec[order(numeric_part)]
}

# Function to create properly clustered heatmap
create_clustered_heatmap <- function(data_matrix, title, filename, 
                                   cluster_rows = TRUE, cluster_cols = TRUE,
                                   width = 14, height = 8) {
  pdf(filename, width = width, height = height)
  
  # Ensure matrix has proper row and column names
  if (is.null(rownames(data_matrix))) {
    rownames(data_matrix) <- paste0("Row", 1:nrow(data_matrix))
  }
  if (is.null(colnames(data_matrix))) {
    colnames(data_matrix) <- paste0("Col", 1:ncol(data_matrix))
  }
  
  # Create heatmap with clustering
  pheatmap(
    data_matrix,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    clustering_distance_rows = "euclidean",
    clustering_distance_cols = "euclidean",
    clustering_method = "complete",
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
    border_color = NA,
    main = title,
    fontsize = 10,
    cellwidth = 20,
    cellheight = 12
  )
  
  dev.off()
}

# ==============================================================================
# 1. FIX BY-GENE ANALYSIS WITH HARMONIZATION
# ==============================================================================

cat("Fixing gene-by-gene analysis with harmonization...\n")

# Configuration
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_gene_fixed"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "plots"), showWarnings = FALSE)

# Load data
all_data <- readRDS(data_file)

# Add gene column if needed
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Add harmonized gene column for cross-method comparison
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
  cat("\nAnalyzing harmonized gene group:", harm_gene, "\n")
  
  # Get all variants for this harmonized gene
  gene_variants <- unique(pd_data$gene[pd_data$gene_harmonized == harm_gene])
  cat("  Includes variants:", paste(gene_variants, collapse = ", "), "\n")
  
  # Filter data for this harmonized gene group
  gene_data <- pd_data %>% filter(gene_harmonized == harm_gene)
  
  # Separate by method
  mast_data <- gene_data %>% filter(method == "MAST")
  mixscale_data <- gene_data %>% filter(method == "MixScale")
  
  # Find unique and convergent terms
  mast_terms <- unique(mast_data$Description)
  mixscale_terms <- unique(mixscale_data$Description)
  convergent_terms <- intersect(mast_terms, mixscale_terms)
  
  cat("  MAST pathways:", length(mast_terms), "\n")
  cat("  MixScale pathways:", length(mixscale_terms), "\n")
  cat("  Convergent pathways:", length(convergent_terms), "\n")
  
  # Store results
  gene_results_harmonized[[harm_gene]] <- list(
    gene = harm_gene,
    variants = gene_variants,
    n_mast_only = length(setdiff(mast_terms, mixscale_terms)),
    n_mixscale_only = length(setdiff(mixscale_terms, mast_terms)),
    n_convergent = length(convergent_terms),
    n_total = length(unique(c(mast_terms, mixscale_terms)))
  )
}

# Create summary table
harmonized_summary <- bind_rows(lapply(gene_results_harmonized, function(x) {
  data.frame(
    gene = x$gene,
    variants = paste(x$variants, collapse = "; "),
    n_mast_only = x$n_mast_only,
    n_mixscale_only = x$n_mixscale_only,
    n_convergent = x$n_convergent,
    n_total = x$n_total,
    stringsAsFactors = FALSE
  )
}))

write.csv(harmonized_summary, file.path(output_dir, "harmonized_gene_summary.csv"), row.names = FALSE)

# Create visualization for each harmonized gene
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
  
  ggsave(file.path(output_dir, "plots", paste0(harm_gene, "_harmonized_summary.pdf")), 
         p, width = 8, height = 6)
}

# ==============================================================================
# 2. CREATE CLUSTERED HEATMAP WITH HARMONIZED GENES
# ==============================================================================

cat("\nCreating clustered heatmap with harmonized genes...\n")

# Create matrix for heatmap
heatmap_data <- harmonized_summary %>%
  select(gene, `Mutation Only` = n_mast_only, 
         `CRISPRi Only` = n_mixscale_only, 
         `Convergent` = n_convergent)

# Convert to data frame with row names
heatmap_df <- as.data.frame(heatmap_data[,-1])
rownames(heatmap_df) <- heatmap_data$gene

# Convert to matrix
heatmap_matrix <- as.matrix(heatmap_df)

# Create clustered heatmap
create_clustered_heatmap(
  heatmap_matrix,
  title = "PD Pathway Distribution - Harmonized Genes (Clustered)",
  filename = file.path(output_dir, "harmonized_genes_clustered_heatmap.pdf"),
  cluster_rows = TRUE,
  cluster_cols = FALSE,  # Don't cluster the three categories
  width = 10,
  height = 12
)

# Also create a version with text annotations
pdf(file.path(output_dir, "harmonized_genes_heatmap_annotated.pdf"), width = 10, height = 12)

# Scale data for better visualization
scaled_matrix <- t(scale(t(heatmap_matrix)))

pheatmap(
  scaled_matrix,
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  display_numbers = heatmap_matrix,  # Show actual counts
  number_format = "%.0f",
  number_color = "black",
  color = colorRampPalette(c("blue", "white", "red"))(50),
  border_color = NA,
  main = "PD Pathway Distribution - Harmonized Genes (Z-score scaled)",
  fontsize = 10,
  cellwidth = 60,
  cellheight = 20
)

dev.off()

# ==============================================================================
# 3. FIX CLUSTER ANALYSIS WITH NATURAL SORTING
# ==============================================================================

cat("\nFixing cluster analysis with natural sorting...\n")

cluster_output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_cluster_fixed"
dir.create(cluster_output_dir, showWarnings = FALSE, recursive = TRUE)

# Get all clusters and sort naturally
all_clusters <- unique(pd_data$cluster)
all_clusters_sorted <- natural_sort_clusters(all_clusters)

cat("Clusters in natural order:", paste(all_clusters_sorted, collapse = ", "), "\n")

# Analyze clusters with proper sorting
cluster_results <- list()
for (cluster in all_clusters_sorted) {
  cluster_data <- pd_data %>% filter(cluster == !!cluster)
  
  if (nrow(cluster_data) > 0) {
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

# Combine results
cluster_summary <- bind_rows(cluster_results)

# Ensure natural sorting in output
cluster_summary$cluster <- factor(cluster_summary$cluster, levels = all_clusters_sorted)
cluster_summary <- cluster_summary %>% arrange(cluster)

write.csv(cluster_summary, file.path(cluster_output_dir, "cluster_summary_sorted.csv"), row.names = FALSE)

# Create visualization with proper sorting
p_cluster <- ggplot(cluster_summary, aes(x = cluster, y = n_pathways)) +
  geom_bar(stat = "identity", fill = "#d62728") +
  geom_text(aes(label = n_pathways), vjust = -0.5, size = 3) +
  labs(
    title = "PD Pathways by Cluster (Natural Sort)",
    x = "Cluster",
    y = "Number of Pathways"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  )

ggsave(file.path(cluster_output_dir, "cluster_pathways_sorted.pdf"), p_cluster, width = 12, height = 6)

# Create cluster × gene matrix with proper sorting
top_genes <- names(head(sort(table(pd_data$gene_harmonized), decreasing = TRUE), 15))

cluster_gene_data <- pd_data %>%
  filter(gene_harmonized %in% top_genes) %>%
  group_by(cluster, gene_harmonized) %>%
  summarise(n_pathways = n_distinct(Description), .groups = "drop") %>%
  pivot_wider(names_from = gene_harmonized, values_from = n_pathways, values_fill = 0)

# Ensure clusters are sorted
cluster_gene_data$cluster <- factor(cluster_gene_data$cluster, levels = all_clusters_sorted)
cluster_gene_data <- cluster_gene_data %>% arrange(cluster)

# Convert to matrix
cluster_gene_matrix <- as.matrix(cluster_gene_data[,-1])
rownames(cluster_gene_matrix) <- cluster_gene_data$cluster

# Create clustered heatmap
create_clustered_heatmap(
  cluster_gene_matrix,
  title = "Cluster × Gene Pathway Counts (Genes Clustered)",
  filename = file.path(cluster_output_dir, "cluster_gene_heatmap_clustered.pdf"),
  cluster_rows = FALSE,  # Keep cluster order
  cluster_cols = TRUE,   # Cluster genes
  width = 14,
  height = 10
)

cat("\n=== FIXES COMPLETE ===\n")
cat("Results saved to:\n")
cat("- Gene harmonization:", output_dir, "\n")
cat("- Cluster analysis:", cluster_output_dir, "\n")
cat("\nKey improvements:\n")
cat("1. Gene variants properly harmonized (PRKN→PARK2, VPS13C variants→VPS13C, etc.)\n")
cat("2. All heatmaps now include hierarchical clustering\n")
cat("3. Clusters sorted naturally (0,1,2...9,10,11...14)\n")