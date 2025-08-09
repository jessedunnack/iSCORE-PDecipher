#!/usr/bin/env Rscript

# GENERATE FINE CLUSTER DOTPLOT WITH DEVELOPMENTAL ORDERING AND GENE CLUSTERING
# Orders fine clusters by maturity within coarse parents and clusters genes within groups

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(RColorBrewer)

cat("=================================================================\n")
cat("FINE DOTPLOT WITH DEVELOPMENTAL ORDERING AND GENE CLUSTERING\n")
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

# 2. Define developmental trajectory for coarse clusters
cat("\n2. Setting up developmental trajectory...\n")

dev_stages <- list(
  "Early Progenitors" = c(4, 1, 2),
  "Late Progenitors" = c(11),
  "Dividing Cells" = c(10),
  "Immature Neurons" = c(0),
  "Mature Neurons" = c(14),
  "Other Cell Types" = c(3, 8, 6, 7, 12, 5, 9, 13)
)

coarse_order <- unlist(dev_stages)

# 3. Order fine clusters within each coarse parent
cat("\n3. Ordering fine clusters by maturity within coarse parents...\n")

# For DA neurons (coarse cluster 0), we'll order by neuronal maturity markers
order_fine_clusters_by_maturity <- function(seurat_obj, fine_clusters, coarse_id) {
  
  if (coarse_id == 0) {  # DA neurons
    # Calculate maturity score based on key markers
    maturity_markers <- c("TH", "SLC18A2", "SLC6A3", "CALB1", "KCNJ6")
    maturity_scores <- numeric(length(fine_clusters))
    
    for (i in seq_along(fine_clusters)) {
      fc <- fine_clusters[i]
      cells <- which(seurat_obj$seurat_clusters_fine == as.character(fc))
      
      if (length(cells) > 0) {
        expr_data <- GetAssayData(seurat_obj, slot = "data")
        available_markers <- intersect(maturity_markers, rownames(expr_data))
        
        if (length(available_markers) > 0) {
          marker_expr <- expr_data[available_markers, cells, drop = FALSE]
          maturity_scores[i] <- mean(colMeans(as.matrix(marker_expr)))
        }
      }
    }
    
    # Order by maturity score (low to high)
    return(fine_clusters[order(maturity_scores)])
    
  } else if (coarse_id %in% c(1, 2, 4, 11)) {  # Progenitors
    # Order by differentiation markers
    diff_markers <- c("DCX", "NEUROD1", "NEUROD2", "STMN2")
    diff_scores <- numeric(length(fine_clusters))
    
    for (i in seq_along(fine_clusters)) {
      fc <- fine_clusters[i]
      cells <- which(seurat_obj$seurat_clusters_fine == as.character(fc))
      
      if (length(cells) > 0) {
        expr_data <- GetAssayData(seurat_obj, slot = "data")
        available_markers <- intersect(diff_markers, rownames(expr_data))
        
        if (length(available_markers) > 0) {
          marker_expr <- expr_data[available_markers, cells, drop = FALSE]
          diff_scores[i] <- mean(colMeans(as.matrix(marker_expr)))
        }
      }
    }
    
    return(fine_clusters[order(diff_scores)])
    
  } else {
    # For other clusters, maintain original order
    return(sort(fine_clusters))
  }
}

# Create ordered fine cluster list
fine_cluster_order <- numeric()

for (coarse_cl in coarse_order) {
  # Get fine clusters for this coarse cluster
  fine_cls <- fine_to_coarse$fine_cluster[fine_to_coarse$coarse_cluster == coarse_cl]
  
  if (length(fine_cls) > 0) {
    # Order them by maturity/differentiation
    fine_cls_ordered <- order_fine_clusters_by_maturity(seurat_obj, fine_cls, coarse_cl)
    fine_cluster_order <- c(fine_cluster_order, fine_cls_ordered)
    
    cat(sprintf("  Coarse %d: fine clusters %s\n", 
                coarse_cl, paste(fine_cls_ordered, collapse = ", ")))
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
genes_already_added <- character()  # Track which genes have been added

for (coarse_cl in coarse_order) {
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

# Create labels
fine_labels <- fine_identities %>%
  left_join(coarse_identities %>% select(cluster, coarse_identity = identity), 
            by = c("coarse_cluster" = "cluster")) %>%
  mutate(
    label = paste0("F", fine_cluster, " (C", coarse_cluster, "): ", 
                   gsub("Progenitors_|Neurons_|Cells_|Fibroblasts_|Mesenchymal_", "", 
                        fine_identity))
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
    limits = c(-2.5, 2.5),  # Set reasonable limits
    oob = scales::squish  # Squish values outside limits
  ) +
  scale_size_continuous(
    name = "Percent\nExpressed",
    range = c(0, 5),
    limits = c(0, 100)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5),  # Smaller for more genes
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
    title = "Fine Cluster Markers - Developmental Trajectory",
    subtitle = "Ordered by maturity within each cell type"
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

for (i in 1:(nrow(fine_positions)-1)) {
  p_enhanced <- p_enhanced +
    geom_hline(yintercept = fine_positions$max_pos[i] + 0.5, 
               linetype = "solid", color = "black", linewidth = 0.5)
}

# 8. Save plots
cat("\n8. Saving plots...\n")

dir.create("results/dotplots", showWarnings = FALSE, recursive = TRUE)

pdf("results/dotplots/dotplot_fine_clusters_dev_trajectory.pdf", width = 18, height = 14)
print(p_enhanced)
dev.off()
cat("  - Saved PDF: results/dotplots/dotplot_fine_clusters_dev_trajectory.pdf\n")

png("results/dotplots/dotplot_fine_clusters_dev_trajectory.png", width = 1800, height = 1400, res = 150)
print(p_enhanced)
dev.off()
cat("  - Saved PNG: results/dotplots/dotplot_fine_clusters_dev_trajectory.png\n")

# 9. Create summary
cat("\n9. Creating summary statistics...\n")

# Fine cluster ordering summary
fine_order_summary <- fine_labels %>%
  filter(fine_cluster %in% fine_cluster_order) %>%
  arrange(match(fine_cluster, fine_cluster_order)) %>%
  mutate(plot_order = seq_len(n())) %>%
  select(plot_order, fine_cluster, coarse_cluster, coarse_identity, fine_identity, n_cells)

write.csv(fine_order_summary, 
          "results/dotplots/fine_cluster_dev_order_summary.csv", 
          row.names = FALSE)

# Print summary
dev_stage_summary <- fine_order_summary %>%
  mutate(dev_stage = case_when(
    coarse_cluster %in% c(4, 1, 2) ~ "Early Progenitors",
    coarse_cluster == 11 ~ "Late Progenitors",
    coarse_cluster == 10 ~ "Dividing Cells",
    coarse_cluster == 0 ~ "Immature Neurons",
    coarse_cluster == 14 ~ "Mature Neurons",
    TRUE ~ "Other Cell Types"
  )) %>%
  group_by(dev_stage) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    .groups = 'drop'
  )

cat("\nFine clusters by developmental stage:\n")
print(dev_stage_summary)

cat("\n=== ENHANCED FINE DOTPLOT COMPLETE ===\n")
cat("Features:\n")
cat("- Fine clusters ordered by maturity within coarse parents\n")
cat("- Genes clustered separately within each coarse parent group\n")
cat("- Visual separators between cell types\n")
cat("- Developmental trajectory preserved\n")