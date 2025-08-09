#!/usr/bin/env Rscript

# GENERATE COARSE CLUSTER DOTPLOT WITH IMPROVED BIOLOGICAL ORDERING
# Implements proper developmental trajectory from early progenitors to DA neurons
# with row clustering within developmental groups

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

cat("=================================================================\n")
cat("COARSE DOTPLOT WITH IMPROVED BIOLOGICAL ORDERING\n")
cat("=================================================================\n\n")

# 1. Load data
cat("1. Loading data...\n")

# Load Seurat object
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")
cat("  - Loaded Seurat object\n")

# Original markers from user's dotplot code
original_markers <- c("CALB1","SOX6","FOXA2","SHH","NTN1","CORIN","HES1","FOXP2","TFF3","MAOB",
                     "SPARCL1","NMU","ERBB4","RGCC","APCDD1","CRABP1","DCC","COL1A1","MKI67",
                     "HTR2C","SYN1","SNCA","NRXN1","STMN2","GAP43","SNAP25","CARTPT","CNTNAP5",
                     "LMX1A","LMX1B","NR4A2","KCNJ6","TH","MYT11")
cat(sprintf("  - Loaded %d original markers\n", length(original_markers)))

# Load selected markers
coarse_genes <- readLines("results/dotplot_markers/coarse_genes.txt")
cat(sprintf("  - Loaded %d selected markers\n", length(coarse_genes)))

# Combine gene lists (union to avoid duplicates)
all_genes <- union(original_markers, coarse_genes)
cat(sprintf("  - Total unique markers: %d\n", length(all_genes)))

coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

# 2. Define IMPROVED developmental trajectory ordering
cat("\n2. Setting up improved biological trajectory ordering...\n")

# Define developmental stage groups with proper biological progression
# Order will be REVERSED for display (bottom to top: early prog → DA neurons)
dev_stages <- list(
  "DA Neurons" = c(0, 14),                    # Dopaminergic & Hypothalamic neurons (TOP)
  "Neuroblasts/Immature" = c(6, 5, 9, 12, 13), # Stressed, Unidentified, PTGDS+, Neuroendocrine, RBP4+
  "Dividing Cells" = c(10),                   # Proliferating cells
  "Choroid Plexus" = c(7),                    # Choroid plexus epithelium
  "ECM/Fibroblasts" = c(3, 8),                # Mesenchymal fibroblasts & ECM fibroblasts
  "Mature Progenitors" = c(2, 11),            # PTPRZ1+ & CRABP1+ progenitors
  "Early Progenitors" = c(4, 1)               # Uncommitted & Intermediate progenitors (BOTTOM)
)

# Reverse the order for bottom-to-top display
dev_stages <- rev(dev_stages)
cluster_order <- unlist(dev_stages)

cat("  Developmental trajectory (bottom to top):\n")
for (stage_name in names(dev_stages)) {
  cat(sprintf("    %s: clusters %s\n", stage_name, 
              paste(dev_stages[[stage_name]], collapse = ", ")))
}

# 3. Get expression data for row clustering
cat("\n3. Extracting expression data for clustering...\n")

# Set up Seurat object
Idents(seurat_obj) <- "seurat_clusters_coarse"
DefaultAssay(seurat_obj) <- "SCT"

# Check gene availability
available_genes <- intersect(all_genes, rownames(seurat_obj))
cat(sprintf("  Using %d available genes (of %d requested)\n", length(available_genes), length(all_genes)))

# Get expression matrix using DotPlot
p_data <- DotPlot(
  seurat_obj, 
  features = available_genes,
  assay = 'SCT',
  cluster.idents = FALSE
)

plot_data <- p_data$data
expr_matrix <- plot_data %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  column_to_rownames("features.plot") %>%
  as.matrix()

# 4. Cluster rows (cell types) within each developmental group
cat("\n4. Clustering cell types within developmental groups...\n")

clustered_order <- numeric()

for (stage_name in names(dev_stages)) {
  stage_clusters <- dev_stages[[stage_name]]
  cat(sprintf("\n  %s (clusters: %s)\n", stage_name, paste(stage_clusters, collapse = ", ")))
  
  if (length(stage_clusters) > 1) {
    # Get expression for these clusters
    stage_expr <- t(expr_matrix[, as.character(stage_clusters), drop = FALSE])
    
    # Perform hierarchical clustering on clusters
    cluster_dist <- dist(stage_expr, method = "euclidean")
    cluster_hclust <- hclust(cluster_dist, method = "ward.D2")
    
    # Get reordered clusters
    clusters_reordered <- stage_clusters[cluster_hclust$order]
    clustered_order <- c(clustered_order, clusters_reordered)
    
    cat(sprintf("    Reordered by similarity: %s\n", paste(clusters_reordered, collapse = ", ")))
  } else {
    clustered_order <- c(clustered_order, stage_clusters)
    cat("    Single cluster - no reordering needed\n")
  }
}

# 5. Cluster genes within each developmental stage
cat("\n5. Clustering genes within developmental stages...\n")

gene_order <- character()
stage_boundaries <- numeric()
genes_already_added <- character()

for (stage_name in names(dev_stages)) {
  stage_clusters <- as.character(dev_stages[[stage_name]])
  cat(sprintf("\n  %s:\n", stage_name))
  
  # Get genes that are expressed in at least one cluster of this stage
  stage_expr <- expr_matrix[, stage_clusters, drop = FALSE]
  
  # Filter for genes with meaningful expression in this stage
  max_expr <- apply(stage_expr, 1, max, na.rm = TRUE)
  stage_genes <- names(max_expr[max_expr > -0.5])
  
  # Remove genes that have already been added
  stage_genes <- setdiff(stage_genes, genes_already_added)
  
  if (length(stage_genes) > 1) {
    # Perform hierarchical clustering
    stage_expr_subset <- stage_expr[stage_genes, , drop = FALSE]
    gene_dist <- dist(stage_expr_subset, method = "euclidean")
    gene_clust <- hclust(gene_dist, method = "ward.D2")
    
    # Get ordered genes
    genes_ordered <- stage_genes[gene_clust$order]
    gene_order <- c(gene_order, genes_ordered)
    genes_already_added <- c(genes_already_added, genes_ordered)
    
    cat(sprintf("    Clustered %d genes\n", length(genes_ordered)))
  } else if (length(stage_genes) == 1) {
    gene_order <- c(gene_order, stage_genes)
    genes_already_added <- c(genes_already_added, stage_genes)
  }
  
  # Store boundary position
  if (length(gene_order) > 0) {
    stage_boundaries <- c(stage_boundaries, length(gene_order))
  }
}

# Add any remaining genes
remaining_genes <- setdiff(available_genes, genes_already_added)
if (length(remaining_genes) > 0) {
  gene_order <- c(gene_order, remaining_genes)
  cat(sprintf("\n  Added %d remaining genes\n", length(remaining_genes)))
}

# 6. Create enhanced dotplot with improved ordering
cat("\n6. Creating enhanced dotplot...\n")

# Update Seurat object with clustered order
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(clustered_order)
)

# Get new plot data with reordered clusters
p_data <- DotPlot(
  seurat_obj, 
  features = available_genes,
  assay = 'SCT',
  cluster.idents = FALSE
)
plot_data <- p_data$data

# Reorder plot data
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)

# Create identity labels
identity_labels <- paste0("C", coarse_identities$cluster, ": ", coarse_identities$identity)
names(identity_labels) <- as.character(coarse_identities$cluster)

plot_data$id <- factor(
  identity_labels[as.character(plot_data$id)],
  levels = identity_labels[as.character(clustered_order)]
)

# Create the plot
p_enhanced <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,  # Reverse to get blue-white-red
    name = "Scaled\nExpression",
    limits = c(-2.5, 2.5),
    oob = scales::squish
  ) +
  scale_size_continuous(
    name = "Percent\nExpressed",
    range = c(0, 8),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    plot.margin = margin(10, 10, 40, 10)  # Extra bottom margin for labels
  ) +
  labs(
    x = "Marker Genes (Clustered within Developmental Stages)",
    y = "Cell Type",
    title = "Coarse Cluster Markers - Biological Trajectory",
    subtitle = "Early Progenitors → Mature Progenitors → Support Cells → Dividing → Neuroblasts → Neurons"
  )

# Add vertical lines to separate gene stages
if (length(stage_boundaries) > 1) {
  for (i in 1:(length(stage_boundaries)-1)) {
    p_enhanced <- p_enhanced +
      geom_vline(xintercept = stage_boundaries[i] + 0.5, 
                 linetype = "dashed", color = "grey40", linewidth = 0.5)
  }
}

# Add horizontal lines to separate developmental stages
stage_positions <- cumsum(table(factor(clustered_order, levels = clustered_order)))
stage_group_positions <- numeric()
current_pos <- 0

for (stage_name in names(dev_stages)) {
  n_clusters <- length(dev_stages[[stage_name]])
  if (n_clusters > 0) {
    current_pos <- current_pos + n_clusters
    if (current_pos < length(clustered_order)) {
      stage_group_positions <- c(stage_group_positions, current_pos)
    }
  }
}

for (pos in stage_group_positions) {
  p_enhanced <- p_enhanced +
    geom_hline(yintercept = pos + 0.5, 
               linetype = "solid", color = "black", linewidth = 0.5)
}

# 7. Save plots
cat("\n7. Saving plots...\n")

dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

pdf("results/dotplots/dotplot_coarse_clusters_improved_order.pdf", width = 14, height = 10)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_coarse_clusters_improved_order.pdf\n")

png("results/dotplots/dotplot_coarse_clusters_improved_order.png", width = 1400, height = 1000, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_coarse_clusters_improved_order.png\n")

# 8. Create summary
cat("\n8. Creating summary statistics...\n")

# Save ordering information
ordering_summary <- data.frame(
  plot_order = seq_along(clustered_order),
  cluster = clustered_order,
  developmental_stage = rep(names(dev_stages), sapply(dev_stages, length)),
  stringsAsFactors = FALSE
)

# Add identity information
ordering_summary <- ordering_summary %>%
  left_join(coarse_identities, by = "cluster")

write.csv(ordering_summary, 
          "results/dotplots/coarse_improved_order_summary.csv", 
          row.names = FALSE)

cat("\nFinal cluster ordering (bottom to top):\n")
print(ordering_summary[, c("plot_order", "cluster", "identity", "developmental_stage")])

cat("\n=== IMPROVED COARSE DOTPLOT COMPLETE ===\n")
cat("Biological trajectory implemented:\n")
cat("- Bottom: Early progenitors (uncommitted → intermediate)\n")
cat("- Middle: Mature progenitors → ECM/Support cells → Dividing cells\n")
cat("- Top: Neuroblasts → Mature neurons (including DA neurons)\n")
cat("- Clusters within each stage ordered by expression similarity\n")