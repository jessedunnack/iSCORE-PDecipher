#!/usr/bin/env Rscript

# GENERATE FINE CLUSTER DOTPLOT WITH IMPROVED BIOLOGICAL ORDERING
# Orders fine clusters based on coarse parent developmental trajectory
# with clustering within each coarse group

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

cat("=================================================================\n")
cat("FINE DOTPLOT WITH IMPROVED BIOLOGICAL ORDERING\n")
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

# Load selected fine cluster markers
fine_genes <- readLines("results/dotplot_markers/fine_genes.txt")
cat(sprintf("  - Loaded %d selected markers\n", length(fine_genes)))

# Combine gene lists (union to avoid duplicates)
all_genes <- union(original_markers, fine_genes)
cat(sprintf("  - Total unique markers: %d\n", length(all_genes)))

# Load metadata
fine_identities <- read.csv("results/reclustered_analysis/fine_cluster_identities_FINAL.csv")
fine_to_coarse <- read.csv("results/reclustered_analysis/fine_to_coarse_mapping.csv")
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

# 2. Define improved developmental trajectory for coarse clusters
cat("\n2. Setting up improved biological trajectory...\n")

# Same ordering as improved coarse dotplot
dev_stages <- list(
  "DA Neurons" = c(0, 14),
  "Neuroblasts/Immature" = c(6, 5, 9, 12, 13),
  "Dividing Cells" = c(10),
  "Choroid Plexus" = c(7),
  "ECM/Fibroblasts" = c(3, 8),
  "Mature Progenitors" = c(2, 11),
  "Early Progenitors" = c(4, 1)
)

# Reverse for bottom-to-top display
dev_stages <- rev(dev_stages)

# Get coarse cluster order from improved analysis
# Based on expression similarity within groups
coarse_order_improved <- c(4, 1, 2, 11, 3, 8, 7, 10, 12, 5, 13, 6, 9, 0, 14)

cat("  Coarse cluster order (bottom to top):\n")
for (i in seq_along(coarse_order_improved)) {
  coarse_cl <- coarse_order_improved[i]
  coarse_info <- coarse_identities[coarse_identities$cluster == coarse_cl, ]
  cat(sprintf("    %2d. C%d: %s\n", i, coarse_cl, coarse_info$identity))
}

# 3. Order fine clusters within each coarse parent
cat("\n3. Ordering fine clusters within coarse parents...\n")

# Function to order fine clusters by expression similarity
order_fine_clusters_by_similarity <- function(seurat_obj, fine_clusters) {
  if (length(fine_clusters) <= 1) {
    return(fine_clusters)
  }
  
  # Get expression data for these fine clusters
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Calculate mean expression for each cluster
  cluster_means <- list()
  for (fc in fine_clusters) {
    cells <- which(seurat_obj$seurat_clusters_fine == as.character(fc))
    if (length(cells) > 0) {
      # Use a subset of highly variable genes for clustering
      var_genes <- intersect(VariableFeatures(seurat_obj)[1:500], rownames(expr_data))
      cluster_means[[as.character(fc)]] <- rowMeans(expr_data[var_genes, cells, drop = FALSE])
    }
  }
  
  if (length(cluster_means) > 1) {
    # Convert to matrix
    mean_matrix <- do.call(cbind, cluster_means)
    
    # Cluster
    cluster_dist <- dist(t(mean_matrix), method = "euclidean")
    cluster_hclust <- hclust(cluster_dist, method = "ward.D2")
    
    # Return ordered clusters
    ordered_names <- colnames(mean_matrix)[cluster_hclust$order]
    return(as.numeric(ordered_names))
  } else {
    return(fine_clusters)
  }
}

# Create ordered fine cluster list
fine_cluster_order <- numeric()

for (coarse_cl in coarse_order_improved) {
  # Get fine clusters for this coarse cluster
  fine_cls <- fine_to_coarse$fine_cluster[fine_to_coarse$coarse_cluster == coarse_cl]
  
  if (length(fine_cls) > 0) {
    # Order them by similarity
    fine_cls_ordered <- order_fine_clusters_by_similarity(seurat_obj, fine_cls)
    fine_cluster_order <- c(fine_cluster_order, fine_cls_ordered)
    
    coarse_info <- coarse_identities[coarse_identities$cluster == coarse_cl, ]
    cat(sprintf("  C%d (%s): fine clusters %s\n", 
                coarse_cl, coarse_info$identity,
                paste(fine_cls_ordered, collapse = ", ")))
  }
}

# 4. Set up Seurat object for plotting
cat("\n4. Preparing data for dotplot...\n")

Idents(seurat_obj) <- "seurat_clusters_fine"
seurat_obj@active.ident <- factor(
  seurat_obj@active.ident,
  levels = as.character(fine_cluster_order)
)

# Check gene availability
DefaultAssay(seurat_obj) <- "SCT"
available_genes <- intersect(all_genes, rownames(seurat_obj))
cat(sprintf("  Using %d available genes (of %d requested)\n", length(available_genes), length(all_genes)))

# 5. Get expression data
cat("\n5. Extracting expression data...\n")

p_data <- DotPlot(
  seurat_obj, 
  features = available_genes,
  assay = 'SCT',
  cluster.idents = FALSE
)

plot_data <- p_data$data

# Create expression matrix
expr_matrix <- plot_data %>%
  select(features.plot, id, avg.exp.scaled) %>%
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>%
  column_to_rownames("features.plot") %>%
  as.matrix()

# 6. Cluster genes within each coarse parent group
cat("\n6. Clustering genes within coarse parent groups...\n")

gene_order <- character()
coarse_boundaries <- numeric()
genes_already_added <- character()

for (coarse_cl in coarse_order_improved) {
  # Get fine clusters for this coarse cluster
  fine_cls <- fine_to_coarse$fine_cluster[fine_to_coarse$coarse_cluster == coarse_cl]
  fine_cls <- intersect(fine_cluster_order, fine_cls)  # Maintain order
  
  if (length(fine_cls) > 0) {
    coarse_info <- coarse_identities[coarse_identities$cluster == coarse_cl, ]
    cat(sprintf("\n  Coarse %d (%s):\n", coarse_cl, coarse_info$identity))
    
    # Get expression for these fine clusters
    fine_cls_chr <- as.character(fine_cls)
    available_cols <- intersect(fine_cls_chr, colnames(expr_matrix))
    
    if (length(available_cols) > 0) {
      coarse_expr <- expr_matrix[, available_cols, drop = FALSE]
      
      # Filter for genes with meaningful expression
      max_expr <- apply(coarse_expr, 1, max, na.rm = TRUE)
      coarse_genes <- names(max_expr[max_expr > -0.5])
      
      # Remove genes that have already been added
      coarse_genes <- setdiff(coarse_genes, genes_already_added)
      
      if (length(coarse_genes) > 1) {
        # Perform hierarchical clustering
        coarse_expr_subset <- coarse_expr[coarse_genes, , drop = FALSE]
        gene_dist <- dist(coarse_expr_subset, method = "euclidean")
        gene_clust <- hclust(gene_dist, method = "ward.D2")
        
        # Get ordered genes
        genes_ordered <- coarse_genes[gene_clust$order]
        gene_order <- c(gene_order, genes_ordered)
        genes_already_added <- c(genes_already_added, genes_ordered)
        
        cat(sprintf("    Clustered %d genes for %d fine clusters\n", 
                    length(genes_ordered), length(available_cols)))
      } else if (length(coarse_genes) == 1) {
        gene_order <- c(gene_order, coarse_genes)
        genes_already_added <- c(genes_already_added, coarse_genes)
      }
    }
    
    # Store boundary
    if (length(gene_order) > 0) {
      coarse_boundaries <- c(coarse_boundaries, length(gene_order))
    }
  }
}

# Add remaining genes
remaining_genes <- setdiff(available_genes, genes_already_added)
if (length(remaining_genes) > 0) {
  gene_order <- c(gene_order, remaining_genes)
  cat(sprintf("\n  Added %d remaining genes\n", length(remaining_genes)))
}

# 7. Create enhanced dotplot
cat("\n7. Creating enhanced dotplot...\n")

# Create labels with hierarchical information
fine_labels <- fine_identities %>%
  left_join(coarse_identities, by = c("coarse_cluster" = "cluster")) %>%
  mutate(
    # Shorten labels for readability
    short_fine = gsub("Progenitors_|Neurons_|Cells_|Fibroblasts_|Mesenchymal_", "", fine_identity),
    short_coarse = gsub("Progenitors_|Neurons_|Cells_|Fibroblasts_|Mesenchymal_", "", identity),
    label = paste0("F", fine_cluster, ": ", short_fine, " (C", coarse_cluster, ")")
  )

identity_labels <- fine_labels$label
names(identity_labels) <- as.character(fine_labels$fine_cluster)

# Reorder plot data
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)
plot_data$id <- factor(
  identity_labels[as.character(plot_data$id)],
  levels = identity_labels[as.character(fine_cluster_order)]
)

# Create plot
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
    range = c(0, 5),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 12),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey90", linewidth = 0.2),
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  ) +
  labs(
    x = "Marker Genes (Clustered within Coarse Parents)",
    y = "Fine Cluster Identity",
    title = "Fine Cluster Markers - Biological Trajectory",
    subtitle = "Ordered by developmental progression within coarse cell types"
  )

# Add vertical lines between coarse parent gene groups
if (length(coarse_boundaries) > 1) {
  for (i in 1:(length(coarse_boundaries)-1)) {
    p_enhanced <- p_enhanced +
      geom_vline(xintercept = coarse_boundaries[i] + 0.5, 
                 linetype = "dashed", color = "grey40", linewidth = 0.5)
  }
}

# Add horizontal lines between coarse parents
fine_positions <- fine_labels %>%
  arrange(match(fine_cluster, fine_cluster_order)) %>%
  mutate(position = seq_len(n())) %>%
  group_by(coarse_cluster) %>%
  summarise(max_pos = max(position), .groups = 'drop')

# Order by coarse cluster appearance
coarse_boundaries_y <- numeric()
current_pos <- 0
for (coarse_cl in coarse_order_improved) {
  n_fine <- sum(fine_to_coarse$coarse_cluster == coarse_cl)
  if (n_fine > 0) {
    current_pos <- current_pos + n_fine
    if (current_pos < length(fine_cluster_order)) {
      coarse_boundaries_y <- c(coarse_boundaries_y, current_pos)
    }
  }
}

for (boundary in coarse_boundaries_y) {
  p_enhanced <- p_enhanced +
    geom_hline(yintercept = boundary + 0.5, 
               linetype = "solid", color = "black", linewidth = 0.5)
}

# 8. Save plots
cat("\n8. Saving plots...\n")

dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

pdf("results/dotplots/dotplot_fine_clusters_improved_order.pdf", width = 18, height = 14)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_fine_clusters_improved_order.pdf\n")

png("results/dotplots/dotplot_fine_clusters_improved_order.png", width = 1800, height = 1400, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_fine_clusters_improved_order.png\n")

# 9. Create summary
cat("\n9. Creating summary statistics...\n")

# Fine cluster ordering summary
fine_order_summary <- fine_labels %>%
  filter(fine_cluster %in% fine_cluster_order) %>%
  arrange(match(fine_cluster, fine_cluster_order)) %>%
  mutate(plot_order = seq_len(n())) %>%
  select(plot_order, fine_cluster, coarse_cluster, identity, fine_identity, n_cells) %>%
  rename(coarse_identity = identity)

# Add developmental stage
stage_mapping <- data.frame(
  coarse_cluster = unlist(dev_stages),
  developmental_stage = rep(names(dev_stages), sapply(dev_stages, length))
)

fine_order_summary <- fine_order_summary %>%
  left_join(stage_mapping, by = "coarse_cluster")

write.csv(fine_order_summary, 
          "results/dotplots/fine_improved_order_summary.csv", 
          row.names = FALSE)

# Print summary
dev_stage_summary <- fine_order_summary %>%
  group_by(developmental_stage) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    .groups = 'drop'
  ) %>%
  mutate(
    developmental_stage = factor(developmental_stage, levels = names(dev_stages))
  ) %>%
  arrange(developmental_stage)

cat("\nFine clusters by developmental stage:\n")
print(dev_stage_summary)

cat("\n=== IMPROVED FINE DOTPLOT COMPLETE ===\n")
cat("Features:\n")
cat("- Fine clusters ordered by biological trajectory\n")
cat("- Grouped by coarse parent cell types\n")
cat("- Genes clustered within each coarse parent group\n")
cat("- Visual separators between cell types\n")
cat("- Bottom to top: Early prog → Mature prog → Support → Dividing → Neuroblasts → Neurons\n")