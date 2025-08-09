#!/usr/bin/env Rscript

# GENERATE COARSE CLUSTER DOTPLOT WITH CUSTOM USER-SPECIFIED ORDERING
# Uses exact order requested by user for Y-axis display

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

cat("=================================================================\n")
cat("COARSE DOTPLOT WITH CUSTOM USER-SPECIFIED ORDERING\n")
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

# 2. Define CUSTOM ordering based on user preference
cat("\n2. Setting up custom user-specified ordering...\n")

# User's preferred order (bottom to top)
# Based on clarified order without duplicates
custom_order <- c(
  1,   # Progenitors_Intermediate
  8,   # Fibroblasts_ECM
  4,   # Progenitors_Uncommitted
  3,   # Mesenchymal_Fibroblasts
  5,   # Cells_Unidentified
  9,   # Cells_PTGDS+
  7,   # Choroid_Plexus
  10,  # Cells_Proliferating
  6,   # Cells_Stressed
  13,  # Cells_RBP4+
  11,  # Progenitors_CRABP1+
  2,   # Progenitors_PTPRZ1+
  12,  # Cells_Neuroendocrine
  14,  # Neurons_Hypothalamic_HCRT
  0    # Neurons_Dopaminergic
)

# Define groups for visual separation (optional)
custom_groups <- list(
  "Early Progenitors" = c(1, 4),
  "Support/Structural" = c(8, 3),
  "Unidentified/Transitional" = c(5, 9),
  "Specialized Non-neuronal" = c(12, 7),
  "Stressed/Proliferating" = c(6, 10, 13),
  "Late Progenitors" = c(11, 2),
  "Mature Neurons" = c(14, 0)
)

cat("  Custom order (bottom to top):\n")
for (i in seq_along(custom_order)) {
  cluster <- custom_order[i]
  info <- coarse_identities[coarse_identities$cluster == cluster, ]
  cat(sprintf("    %2d. C%d: %s\n", i, cluster, info$identity))
}

# 3. Get expression data
cat("\n3. Extracting expression data...\n")

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

# 4. Cluster genes within custom groups
cat("\n4. Clustering genes within custom groups...\n")

gene_order <- character()
group_boundaries <- numeric()
genes_already_added <- character()

for (group_name in names(custom_groups)) {
  group_clusters <- as.character(custom_groups[[group_name]])
  cat(sprintf("\n  %s:\n", group_name))
  
  # Get genes that are expressed in at least one cluster of this group
  group_expr <- expr_matrix[, group_clusters, drop = FALSE]
  
  # Filter for genes with meaningful expression in this group
  max_expr <- apply(group_expr, 1, max, na.rm = TRUE)
  group_genes <- names(max_expr[max_expr > -0.5])
  
  # Remove genes that have already been added
  group_genes <- setdiff(group_genes, genes_already_added)
  
  if (length(group_genes) > 1) {
    # Perform hierarchical clustering
    group_expr_subset <- group_expr[group_genes, , drop = FALSE]
    gene_dist <- dist(group_expr_subset, method = "euclidean")
    gene_clust <- hclust(gene_dist, method = "ward.D2")
    
    # Get ordered genes
    genes_ordered <- group_genes[gene_clust$order]
    gene_order <- c(gene_order, genes_ordered)
    genes_already_added <- c(genes_already_added, genes_ordered)
    
    cat(sprintf("    Clustered %d genes\n", length(genes_ordered)))
  } else if (length(group_genes) == 1) {
    gene_order <- c(gene_order, group_genes)
    genes_already_added <- c(genes_already_added, group_genes)
  }
  
  # Store boundary position
  if (length(gene_order) > 0) {
    group_boundaries <- c(group_boundaries, length(gene_order))
  }
}

# Add any remaining genes
remaining_genes <- setdiff(available_genes, genes_already_added)
if (length(remaining_genes) > 0) {
  gene_order <- c(gene_order, remaining_genes)
  cat(sprintf("\n  Added %d remaining genes\n", length(remaining_genes)))
}

# 5. Create enhanced dotplot with custom ordering
cat("\n5. Creating enhanced dotplot...\n")

# Update Seurat object with custom order
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(custom_order)
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
  levels = identity_labels[as.character(custom_order)]
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
    plot.margin = margin(10, 10, 40, 10)
  ) +
  labs(
    x = "Marker Genes (Clustered within Groups)",
    y = "Cell Type",
    title = "Coarse Cluster Markers - Custom Order",
    subtitle = "User-specified ordering with grouped gene clustering"
  )

# Add vertical lines to separate gene groups
if (length(group_boundaries) > 1) {
  for (i in 1:(length(group_boundaries)-1)) {
    p_enhanced <- p_enhanced +
      geom_vline(xintercept = group_boundaries[i] + 0.5, 
                 linetype = "dashed", color = "grey40", linewidth = 0.5)
  }
}

# Add horizontal lines to separate custom groups
group_positions <- numeric()
current_pos <- 0

for (group_name in names(custom_groups)) {
  n_clusters <- length(custom_groups[[group_name]])
  if (n_clusters > 0) {
    current_pos <- current_pos + n_clusters
    if (current_pos < length(custom_order)) {
      group_positions <- c(group_positions, current_pos)
    }
  }
}

for (pos in group_positions) {
  p_enhanced <- p_enhanced +
    geom_hline(yintercept = pos + 0.5, 
               linetype = "solid", color = "black", linewidth = 0.5)
}

# 6. Save plots
cat("\n6. Saving plots...\n")

dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

pdf("results/dotplots/dotplot_coarse_clusters_custom_order.pdf", width = 14, height = 10)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_coarse_clusters_custom_order.pdf\n")

png("results/dotplots/dotplot_coarse_clusters_custom_order.png", width = 1400, height = 1000, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_coarse_clusters_custom_order.png\n")

# 7. Create summary
cat("\n7. Creating summary statistics...\n")

# Save ordering information
ordering_summary <- data.frame(
  plot_order = seq_along(custom_order),
  cluster = custom_order,
  stringsAsFactors = FALSE
)

# Add group information
ordering_summary$group <- NA
for (group_name in names(custom_groups)) {
  clusters_in_group <- custom_groups[[group_name]]
  ordering_summary$group[ordering_summary$cluster %in% clusters_in_group] <- group_name
}

# Add identity information
ordering_summary <- ordering_summary %>%
  left_join(coarse_identities, by = "cluster")

write.csv(ordering_summary, 
          "results/dotplots/coarse_custom_order_summary.csv", 
          row.names = FALSE)

cat("\nFinal custom cluster ordering (bottom to top):\n")
print(ordering_summary[, c("plot_order", "cluster", "identity", "group")])

cat("\n=== CUSTOM ORDER COARSE DOTPLOT COMPLETE ===\n")
cat("Features:\n")
cat("- User-specified cluster ordering implemented\n")
cat("- Genes clustered within custom groups\n")
cat("- Visual separators between groups\n")
cat("- RdBu color palette and all original markers included\n")