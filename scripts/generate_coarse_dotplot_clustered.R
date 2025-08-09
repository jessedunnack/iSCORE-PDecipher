#!/usr/bin/env Rscript

# GENERATE COARSE CLUSTER DOTPLOT WITH DEVELOPMENTAL ORDERING AND GENE CLUSTERING
# Orders clusters by developmental trajectory and clusters genes within stages

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

cat("=================================================================\n")
cat("COARSE DOTPLOT WITH DEVELOPMENTAL ORDERING AND GENE CLUSTERING\n")
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

# 2. Define developmental trajectory ordering
cat("\n2. Setting up developmental trajectory ordering...\n")

# Define developmental stage groups with biologically meaningful order
dev_stages <- list(
  "Early Progenitors" = c(4, 1, 2),      # Uncommitted → Intermediate → PTPRZ1+
  "Late Progenitors" = c(11),            # CRABP1+ (weak TH)
  "Dividing Cells" = c(10),              # Proliferating
  "Immature Neurons" = c(0),             # Dopaminergic
  "Mature Neurons" = c(14),              # Hypothalamic HCRT
  "Other Cell Types" = c(3, 8, 6, 7, 12, 5, 9, 13)  # Non-neuronal
)

# Create ordered cluster list
cluster_order <- unlist(dev_stages)
cat("  Cluster order:", paste(cluster_order, collapse = ", "), "\n")

# 3. Set up Seurat object for plotting
cat("\n3. Preparing data for dotplot...\n")

Idents(seurat_obj) <- "seurat_clusters_coarse"
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(cluster_order)
)

# Check gene availability
DefaultAssay(seurat_obj) <- "SCT"
available_genes <- intersect(all_genes, rownames(seurat_obj))
cat(sprintf("  Using %d available genes (of %d requested)\n", length(available_genes), length(all_genes)))

# 4. Get expression data from DotPlot
cat("\n4. Extracting expression data...\n")

# Create initial dotplot to get data
p_data <- DotPlot(
  seurat_obj, 
  features = available_genes,
  assay = 'SCT',
  cluster.idents = FALSE
)

# Extract the data
plot_data <- p_data$data

# Create expression matrix for clustering
expr_matrix <- plot_data %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  column_to_rownames("features.plot") %>%
  as.matrix()

cat(sprintf("  Expression matrix: %d genes x %d clusters\n", nrow(expr_matrix), ncol(expr_matrix)))

# 5. Cluster genes within each developmental stage
cat("\n5. Clustering genes within developmental stages...\n")

gene_order <- character()
stage_boundaries <- numeric()
genes_already_added <- character()  # Track which genes have been added

for (stage_name in names(dev_stages)) {
  stage_clusters <- as.character(dev_stages[[stage_name]])
  cat(sprintf("\n  %s (clusters: %s)\n", stage_name, paste(stage_clusters, collapse = ", ")))
  
  # Get genes that are expressed in at least one cluster of this stage
  stage_expr <- expr_matrix[, stage_clusters, drop = FALSE]
  
  # Filter for genes with meaningful expression in this stage
  max_expr <- apply(stage_expr, 1, max, na.rm = TRUE)
  stage_genes <- names(max_expr[max_expr > -0.5])  # Adjust threshold as needed
  
  # Remove genes that have already been added
  stage_genes <- setdiff(stage_genes, genes_already_added)
  
  if (length(stage_genes) > 1) {
    # Perform hierarchical clustering
    stage_expr_subset <- stage_expr[stage_genes, , drop = FALSE]
    
    # Calculate distance and cluster
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
    cat(sprintf("    Single gene: %s\n", stage_genes))
  }
  
  # Store boundary position
  if (length(gene_order) > 0) {
    stage_boundaries <- c(stage_boundaries, length(gene_order))
  }
}

# Add any remaining genes not assigned to stages
remaining_genes <- setdiff(available_genes, genes_already_added)
if (length(remaining_genes) > 0) {
  gene_order <- c(gene_order, remaining_genes)
  cat(sprintf("\n  Added %d remaining genes\n", length(remaining_genes)))
}

# 6. Create enhanced dotplot with ordered genes
cat("\n6. Creating enhanced dotplot...\n")

# Reorder plot data
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)

# Create identity labels
identity_labels <- paste0("C", coarse_identities$cluster, ": ", coarse_identities$identity)
names(identity_labels) <- as.character(coarse_identities$cluster)

plot_data$id <- factor(
  identity_labels[as.character(plot_data$id)],
  levels = identity_labels[as.character(cluster_order)]
)

# Create the plot
p_enhanced <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_color_distiller(
    palette = "RdBu",
    direction = -1,  # Reverse to get blue-white-red
    name = "Scaled\nExpression",
    limits = c(-2.5, 2.5),  # Set reasonable limits
    oob = scales::squish  # Squish values outside limits
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
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  ) +
  labs(
    x = "Marker Genes (Clustered within Developmental Stages)",
    y = "Cell Type",
    title = "Coarse Cluster Markers - Developmental Trajectory",
    subtitle = "Early Progenitors → Late Progenitors → Dividing → Neurons → Other"
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
stage_positions <- cumsum(sapply(dev_stages, length))
if (length(stage_positions) > 1) {
  for (i in 1:(length(stage_positions)-1)) {
    p_enhanced <- p_enhanced +
      geom_hline(yintercept = stage_positions[i] + 0.5, 
                 linetype = "solid", color = "black", linewidth = 0.5)
  }
}

# 7. Add gene group labels
cat("\n7. Adding gene group labels...\n")

# Calculate positions for stage labels
stage_label_positions <- numeric()
stage_label_names <- character()
current_pos <- 0

for (i in 1:length(dev_stages)) {
  stage_name <- names(dev_stages)[i]
  if (i <= length(stage_boundaries) && !is.na(stage_boundaries[i])) {
    stage_start <- ifelse(i == 1, 1, stage_boundaries[i-1] + 1)
    stage_end <- stage_boundaries[i]
    
    if (stage_end > stage_start) {
      # Calculate center position for this stage
      stage_center <- (stage_start + stage_end) / 2
      stage_label_positions <- c(stage_label_positions, stage_center)
      stage_label_names <- c(stage_label_names, stage_name)
      cat(sprintf("  %s: position %.1f\n", stage_name, stage_center))
    }
  }
}

# Add stage labels as a second x-axis
if (length(stage_label_positions) > 0) {
  p_enhanced <- p_enhanced +
    annotate("text", 
             x = stage_label_positions, 
             y = -0.5,  # Below the x-axis
             label = stage_label_names,
             angle = 45, 
             hjust = 1,
             vjust = 1,
             size = 3,
             fontface = "bold") +
    coord_cartesian(clip = "off")  # Allow drawing outside plot area
}

# 8. Save plots
cat("\n8. Saving plots...\n")

dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

# Save main plot
pdf("results/dotplots/dotplot_coarse_clusters_dev_trajectory.pdf", width = 14, height = 10)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_coarse_clusters_dev_trajectory.pdf\n")

png("results/dotplots/dotplot_coarse_clusters_dev_trajectory.png", width = 1400, height = 1000, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_coarse_clusters_dev_trajectory.png\n")

# 9. Create summary
cat("\n9. Creating summary statistics...\n")

# Gene clustering summary
gene_stage_assignment <- data.frame(
  gene = gene_order,
  position = seq_along(gene_order),
  stage = NA
)

# Assign stages to genes
current_pos <- 0
for (i in 1:length(dev_stages)) {
  stage_name <- names(dev_stages)[i]
  stage_end <- stage_boundaries[i]
  if (!is.na(stage_end) && stage_end > current_pos) {
    gene_stage_assignment$stage[(current_pos + 1):stage_end] <- stage_name
    current_pos <- stage_end
  }
}

# Save gene ordering
write.csv(gene_stage_assignment, 
          "results/dotplots/coarse_gene_clustering_summary.csv", 
          row.names = FALSE)

# Print summary
stage_summary <- gene_stage_assignment %>%
  filter(!is.na(stage)) %>%
  group_by(stage) %>%
  summarise(
    n_genes = n(),
    genes = paste(head(gene, 5), collapse = ", ")
  )

cat("\nGene distribution by developmental stage:\n")
print(stage_summary)

cat("\n=== ENHANCED COARSE DOTPLOT COMPLETE ===\n")
cat("Features:\n")
cat("- Clusters ordered by developmental trajectory\n")
cat("- Genes clustered within each developmental stage\n")
cat("- Visual separators between stages\n")
cat("- Clear progression from progenitors to neurons\n")