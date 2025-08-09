#!/usr/bin/env Rscript

# GENERATE FINE CLUSTER DOTPLOT
# Creates a dotplot showing fine cluster distinctions with shared parent markers

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

cat("=================================================================\n")
cat("GENERATING FINE CLUSTER DOTPLOT\n")
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
# First, rename the identity column in coarse_identities to avoid confusion
coarse_info <- coarse_identities %>%
  select(cluster, coarse_identity = identity)

fine_labels <- fine_identities %>%
  left_join(coarse_info, 
            by = c("coarse_cluster" = "cluster")) %>%
  mutate(
    # Short label for plot
    label_short = paste0("F", fine_cluster, ": ", 
                        gsub("Progenitors_|Neurons_|Cells_|Fibroblasts_|Mesenchymal_", "", fine_identity)),
    # Full label with parent
    label_full = paste0("F", fine_cluster, " (C", coarse_cluster, "): ", fine_identity)
  )

# Create label mapping
identity_labels <- fine_labels$label_short
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
  cat(sprintf("  WARNING: %d genes not found: %s\n", 
              length(missing_genes), paste(head(missing_genes, 10), collapse = ", ")))
}
cat(sprintf("  Plotting %d available genes\n", length(genes_to_plot)))

# 4. Create dotplot
cat("\n4. Creating fine cluster dotplot...\n")

# Create base dotplot
p <- DotPlot(
  seurat_obj, 
  features = genes_to_plot,
  assay = 'SCT',
  dot.scale = 6,  # Slightly smaller dots for more clusters
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
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey90", size = 0.2)
  ) +
  labs(
    x = "Marker Genes",
    y = "Fine Cluster Identity",
    title = "Fine Cluster Marker Expression"
  )

# 5. Add cluster identity labels and parent grouping
cat("\n5. Adding cluster labels and parent grouping...\n")

# Get plot data
plot_data <- p$data

# Map cluster numbers to identities
plot_data$id <- factor(
  identity_labels[as.character(plot_data$id)],
  levels = identity_labels[as.character(fine_order)]
)

# Add coarse parent info for coloring
fine_coarse_map <- fine_to_coarse %>%
  select(fine_cluster, coarse_cluster) %>%
  distinct()

plot_data <- plot_data %>%
  left_join(
    data.frame(
      id = identity_labels,
      fine_cluster = as.numeric(names(identity_labels))
    ),
    by = "id"
  ) %>%
  left_join(fine_coarse_map, by = "fine_cluster") %>%
  left_join(coarse_identities[, c("cluster", "cell_type_broad")], 
            by = c("coarse_cluster" = "cluster"))

# 6. Create enhanced plot with parent grouping
cat("\n6. Creating enhanced plot with visual grouping...\n")

# Define colors for broad cell types
cell_type_colors <- c(
  "Progenitors" = "#E8F4E8",
  "Non-neuronal" = "#FFF4E8", 
  "Unknown" = "#F0F0F0",
  "Neurons" = "#E8E8FF"
)

# Get y-axis positions for each fine cluster
y_positions <- seq_along(levels(plot_data$id))
names(y_positions) <- levels(plot_data$id)

# Find boundaries between coarse parent groups
coarse_boundaries <- fine_labels %>%
  arrange(match(fine_cluster, fine_order)) %>%
  mutate(position = seq_len(n())) %>%
  group_by(coarse_cluster) %>%
  summarise(
    start = as.numeric(min(position) - 0.5),
    end = as.numeric(max(position) + 0.5),
    .groups = 'drop'
  ) %>%
  left_join(coarse_identities[, c("cluster", "identity", "cell_type_broad")], 
            by = c("coarse_cluster" = "cluster"))

# Ensure numeric values
coarse_boundaries$start <- as.numeric(coarse_boundaries$start)
coarse_boundaries$end <- as.numeric(coarse_boundaries$end)

# Create the enhanced plot
p_enhanced <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  # Add background rectangles for coarse groups
  geom_rect(
    data = coarse_boundaries,
    aes(xmin = -Inf, xmax = Inf, ymin = start, ymax = end, fill = cell_type_broad),
    alpha = 0.1,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  scale_fill_manual(
    values = cell_type_colors,
    name = "Cell Type\nCategory",
    guide = guide_legend(order = 1),
    na.value = "grey90"
  ) +
  # Add the dots
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
    range = c(0, 6),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 8),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = "Marker Genes",
    y = "Fine Cluster Identity",
    title = "Fine Cluster Marker Expression Grouped by Coarse Parent"
  )

# Add separator lines between coarse groups
for (i in 1:(nrow(coarse_boundaries)-1)) {
  p_enhanced <- p_enhanced +
    geom_hline(yintercept = coarse_boundaries$end[i], 
               linetype = "dashed", color = "grey40", linewidth = 0.5)
}

# 7. Save plots
cat("\n7. Saving plots...\n")

# Create output directory if needed
dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

# Save as PDF (larger size for many clusters)
pdf("results/dotplots/dotplot_fine_clusters.pdf", width = 16, height = 12)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_fine_clusters.pdf\n")

# Save as PNG
png("results/dotplots/dotplot_fine_clusters.png", width = 1600, height = 1200, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_fine_clusters.png\n")

# 8. Create summary statistics
cat("\n8. Creating summary statistics...\n")

# Summarize markers by coarse parent
marker_summary <- fine_labels %>%
  group_by(coarse_cluster, coarse_identity) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    fine_clusters = paste(fine_cluster, collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(match(coarse_cluster, unique(fine_to_coarse$coarse_cluster[match(fine_order, fine_to_coarse$fine_cluster)])))

cat("\nFine clusters by coarse parent:\n")
print(marker_summary)

# Save summaries
write.csv(marker_summary, "results/dotplots/fine_cluster_summary.csv", row.names = FALSE)

# Create gene usage summary
selected_markers <- read.csv("results/dotplot_markers/selected_markers_fine.csv")
gene_usage <- selected_markers %>%
  group_by(gene) %>%
  summarise(
    n_clusters = n_distinct(fine_cluster),
    clusters = paste(unique(fine_cluster), collapse = ", "),
    .groups = 'drop'
  ) %>%
  arrange(desc(n_clusters))

cat("\nMost frequently used marker genes:\n")
print(head(gene_usage, 20))

write.csv(gene_usage, "results/dotplots/fine_marker_gene_usage.csv", row.names = FALSE)

cat("\n=== FINE DOTPLOT COMPLETE ===\n")
cat("Outputs saved to results/dotplots/\n")
cat("\nVisualization features:\n")
cat("- Fine clusters grouped by coarse parent\n")
cat("- Background shading shows cell type categories\n")
cat("- Mix of specific and shared markers\n")
cat("- Diverging color scale (blue-white-red)\n")
cat("- Visual blocks for related clusters\n")