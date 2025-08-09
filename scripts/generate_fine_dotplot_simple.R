#!/usr/bin/env Rscript

# GENERATE FINE CLUSTER DOTPLOT (SIMPLIFIED VERSION)
# Creates a dotplot showing fine cluster distinctions

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

cat("=================================================================\n")
cat("GENERATING FINE CLUSTER DOTPLOT (SIMPLIFIED)\n")
cat("=================================================================\n\n")

# 1. Load data
cat("1. Loading data...\n")

# Load Seurat object
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")
cat("  - Loaded Seurat object\n")

# Load marker genes and metadata
fine_genes <- readLines("results/dotplot_markers/fine_genes.txt")
fine_order <- read.csv("results/dotplot_markers/fine_cluster_order.csv")$cluster
fine_identities <- read.csv("results/reclustered_analysis/fine_cluster_identities_FINAL.csv")
fine_to_coarse <- read.csv("results/reclustered_analysis/fine_to_coarse_mapping.csv")
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

cat(sprintf("  - Loaded %d marker genes\n", length(fine_genes)))
cat(sprintf("  - %d fine clusters to plot\n", length(fine_order)))

# 2. Set up cluster identities and order
cat("\n2. Setting up fine cluster identities...\n")

# Create detailed labels with coarse parent info
coarse_info <- coarse_identities %>%
  select(cluster, coarse_identity = identity, cell_type_broad)

fine_labels <- fine_identities %>%
  left_join(coarse_info, 
            by = c("coarse_cluster" = "cluster")) %>%
  mutate(
    # Label for plot with parent info
    label = paste0("F", fine_cluster, " (C", coarse_cluster, "): ", 
                   gsub("Progenitors_|Neurons_|Cells_|Fibroblasts_|Mesenchymal_", "", fine_identity))
  )

# Create label mapping
identity_labels <- fine_labels$label
names(identity_labels) <- as.character(fine_labels$fine_cluster)

# Set active ident to fine clusters
Idents(seurat_obj) <- "seurat_clusters_fine"

# Reorder clusters according to our specified order
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(fine_order)
)

# 3. Check gene availability
cat("\n3. Checking gene availability...\n")
DefaultAssay(seurat_obj) <- "SCT"
available_genes <- rownames(seurat_obj)
genes_to_plot <- intersect(fine_genes, available_genes)
missing_genes <- setdiff(fine_genes, available_genes)

if (length(missing_genes) > 0) {
  cat(sprintf("  WARNING: %d genes not found\n", length(missing_genes)))
}
cat(sprintf("  Plotting %d available genes\n", length(genes_to_plot)))

# 4. Create dotplot
cat("\n4. Creating fine cluster dotplot...\n")

# Create base dotplot
p <- DotPlot(
  seurat_obj, 
  features = genes_to_plot,
  assay = 'SCT',
  dot.scale = 5,  # Smaller dots for more clusters
  scale.by = "size",
  cluster.idents = FALSE
) + 
  # Diverging color scale
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = "Avg Expression"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2)
  ) +
  labs(
    x = "Marker Genes",
    y = "Fine Cluster Identity",
    title = "Fine Cluster Marker Expression"
  )

# 5. Add cluster identity labels
cat("\n5. Adding cluster labels...\n")

# Get plot data and map cluster numbers to identities
plot_data <- p$data
plot_data$id <- factor(
  identity_labels[as.character(plot_data$id)],
  levels = identity_labels[as.character(fine_order)]
)

# Recreate plot with labeled identities
p_labeled <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = "Scaled\nExpression"
  ) +
  scale_size_continuous(
    name = "Percent\nExpressed",
    range = c(0, 5),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = "Marker Genes",
    y = "Fine Cluster Identity",
    title = "Fine Cluster Marker Expression"
  )

# Add separators between coarse parent groups
cat("  - Adding parent group separators...\n")

# Find positions where coarse parent changes
separator_positions <- fine_labels %>%
  arrange(match(fine_cluster, fine_order)) %>%
  mutate(position = seq_len(n())) %>%
  group_by(coarse_cluster) %>%
  summarise(max_pos = max(position), .groups = 'drop') %>%
  filter(max_pos < length(fine_order)) %>%
  pull(max_pos)

# Add horizontal lines at separator positions
for (pos in separator_positions) {
  p_labeled <- p_labeled +
    geom_hline(yintercept = pos + 0.5, linetype = "dashed", 
               color = "grey50", linewidth = 0.5)
}

# 6. Save plots
cat("\n6. Saving plots...\n")

# Create output directory if needed
dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

# Save as PDF
pdf("results/dotplots/dotplot_fine_clusters_simple.pdf", width = 16, height = 12)
print(p_labeled)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_fine_clusters_simple.pdf\n")

# Save as PNG
png("results/dotplots/dotplot_fine_clusters_simple.png", width = 1600, height = 1200, res = 150)
print(p_labeled)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_fine_clusters_simple.png\n")

# 7. Create summary statistics
cat("\n7. Creating summary statistics...\n")

# Summarize by coarse parent
parent_summary <- fine_labels %>%
  group_by(coarse_cluster, coarse_identity, cell_type_broad) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    fine_clusters = paste(fine_cluster, collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(match(coarse_cluster, unique(fine_to_coarse$coarse_cluster[match(fine_order, fine_to_coarse$fine_cluster)])))

cat("\nFine clusters by coarse parent:\n")
print(parent_summary)

# Save summary
write.csv(parent_summary, "results/dotplots/fine_cluster_parent_summary.csv", row.names = FALSE)

cat("\n=== FINE DOTPLOT COMPLETE ===\n")
cat("Outputs saved to results/dotplots/\n")
cat("\nVisualization features:\n")
cat("- Fine clusters grouped by coarse parent order\n")
cat("- Dashed lines separate coarse parent groups\n")
cat("- Diverging color scale (blue-white-red)\n")
cat("- Mix of specific and shared markers\n")