#!/usr/bin/env Rscript

# GENERATE COARSE CLUSTER DOTPLOT
# Creates a dotplot with optimized marker genes and diverging color scale

library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

cat("=================================================================\n")
cat("GENERATING COARSE CLUSTER DOTPLOT\n")
cat("=================================================================\n\n")

# 1. Load data
cat("1. Loading data...\n")

# Load Seurat object
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")
cat("  - Loaded Seurat object\n")

# Load marker genes
coarse_genes <- readLines("results/dotplot_markers/coarse_genes.txt")
coarse_order <- read.csv("results/dotplot_markers/coarse_cluster_order.csv")$cluster
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

cat(sprintf("  - Loaded %d marker genes\n", length(coarse_genes)))
cat(sprintf("  - Cluster order: %s\n", paste(coarse_order, collapse = ", ")))

# 2. Set up cluster identities and order
cat("\n2. Setting up cluster identities...\n")

# Create identity labels
identity_labels <- coarse_identities$identity
names(identity_labels) <- as.character(coarse_identities$cluster)

# Add cell type broad categories for visual grouping
identity_labels_formatted <- paste0(
  "C", coarse_identities$cluster, ": ",
  coarse_identities$identity
)
names(identity_labels_formatted) <- as.character(coarse_identities$cluster)

# Set active ident to coarse clusters
Idents(seurat_obj) <- "seurat_clusters_coarse"

# Reorder clusters according to our specified order
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(coarse_order)
)

# 3. Check gene availability
cat("\n3. Checking gene availability...\n")
DefaultAssay(seurat_obj) <- "SCT"
available_genes <- rownames(seurat_obj)
genes_to_plot <- intersect(coarse_genes, available_genes)
missing_genes <- setdiff(coarse_genes, available_genes)

if (length(missing_genes) > 0) {
  cat(sprintf("  WARNING: %d genes not found in data: %s\n", 
              length(missing_genes), paste(missing_genes, collapse = ", ")))
}
cat(sprintf("  Plotting %d available genes\n", length(genes_to_plot)))

# 4. Create dotplot
cat("\n4. Creating dotplot...\n")

# Create base dotplot
p <- DotPlot(
  seurat_obj, 
  features = genes_to_plot,
  assay = 'SCT',
  dot.scale = 8,
  scale.by = "size",
  cluster.idents = FALSE  # Keep our custom order
) + 
  # Use diverging color scale
  scale_color_gradient2(
    low = "blue", 
    mid = "white", 
    high = "red",
    midpoint = 0,
    name = "Avg Expression"
  ) +
  # Rotate x-axis labels
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.background = element_rect(fill = "white", colour = "black"),
    panel.grid.major = element_line(colour = "grey90", size = 0.2)
  ) +
  labs(
    x = "Marker Genes",
    y = "Cluster Identity",
    title = "Coarse Cluster Marker Expression"
  )

# 5. Add cluster identity labels
cat("\n5. Adding cluster identity labels...\n")

# Get the plot data to modify y-axis labels
plot_data <- p$data

# Map cluster numbers to identities
plot_data$id <- factor(
  identity_labels_formatted[as.character(plot_data$id)],
  levels = identity_labels_formatted[as.character(coarse_order)]
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
    range = c(0, 8),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 9),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", size = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold")
  ) +
  labs(
    x = "Marker Genes",
    y = "Cell Type",
    title = "Coarse Cluster Marker Expression"
  )

# Add horizontal lines to separate cell type categories
cat("  - Adding category separators...\n")

# Find positions for separators
progenitor_end <- which(coarse_order == 11)
non_neuronal_end <- which(coarse_order == 13)

p_final <- p_labeled +
  geom_hline(yintercept = progenitor_end + 0.5, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = non_neuronal_end + 0.5, linetype = "dashed", color = "grey50")

# 6. Save plots
cat("\n6. Saving plots...\n")

# Create output directory
dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

# Save as PDF (vector format, best for publications)
pdf("results/dotplots/dotplot_coarse_clusters.pdf", width = 12, height = 8)
print(p_final)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_coarse_clusters.pdf\n")

# Save as PNG (raster format, for presentations)
png("results/dotplots/dotplot_coarse_clusters.png", width = 1200, height = 800, res = 150)
print(p_final)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_coarse_clusters.png\n")

# 7. Create summary statistics
cat("\n7. Creating summary statistics...\n")

# Calculate average expression and percent expressed per cluster
expr_summary <- plot_data %>%
  group_by(id) %>%
  summarise(
    n_genes = n(),
    avg_pct_expressed = mean(pct.exp),
    avg_scaled_expr = mean(avg.exp.scaled),
    .groups = 'drop'
  ) %>%
  arrange(match(id, levels(id)))

cat("\nExpression summary by cluster:\n")
print(expr_summary)

# Save summary
write.csv(expr_summary, "results/dotplots/coarse_dotplot_summary.csv", row.names = FALSE)

cat("\n=== COARSE DOTPLOT COMPLETE ===\n")
cat("Outputs saved to results/dotplots/\n")
cat("\nVisualization shows:\n")
cat("- Progenitor clusters at bottom\n")
cat("- Non-neuronal/Unknown in middle\n")
cat("- Neuronal clusters at top\n")
cat("- Diverging color scale (blue-white-red)\n")
cat("- Most definitive markers per cluster\n")