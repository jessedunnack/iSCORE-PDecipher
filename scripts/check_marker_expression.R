#!/usr/bin/env Rscript

# Check marker expression in key clusters

library(Seurat)
library(dplyr)

cat("Loading Seurat object...\n")
x <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")

# Key clusters to check based on critical misclassifications document
clusters_to_check <- list(
  "1" = "Floor Plate Progenitors (but has neuronal markers!)",
  "4" = "Choroid Plexus",
  "5" = "Dopaminergic Neurons"
)

# Define marker sets
neuronal_markers <- c("TUBB3", "STMN2", "GAP43", "MAP2", "MAPT")
progenitor_markers <- c("SOX2", "NES", "VIM", "FABP7")
da_markers <- c("TH", "DDC", "SLC6A3", "SLC18A2")
cp_markers <- c("TTR", "FOLR1", "HTR2C", "CLIC6")

# Check which assay to use
if ("SCT" %in% names(x@assays)) {
  assay_to_use <- "SCT"
} else if ("RNA" %in% names(x@assays)) {
  assay_to_use <- "RNA"
} else {
  stop("No RNA or SCT assay found")
}

cat("\nUsing assay:", assay_to_use, "\n")

# Calculate average expression
cat("\nCalculating average expression by cluster...\n")
avg_expr <- AverageExpression(x, 
                              features = c(neuronal_markers, progenitor_markers, da_markers, cp_markers),
                              group.by = "seurat_clusters_fine", 
                              assays = assay_to_use)[[assay_to_use]]

# Check each cluster
for (cluster in names(clusters_to_check)) {
  cat("\n==================================\n")
  cat("Cluster", cluster, "-", clusters_to_check[[cluster]], "\n")
  cat("==================================\n")
  
  # Get cell type label
  cluster_cells <- which(x$seurat_clusters_fine == cluster)
  celltype <- unique(x$celltypes_fine[cluster_cells])[1]
  cat("Current label:", celltype, "\n")
  cat("Number of cells:", length(cluster_cells), "\n\n")
  
  # Check marker expression
  if (cluster %in% colnames(avg_expr)) {
    cluster_expr <- avg_expr[, cluster]
    
    # Neuronal markers
    cat("Neuronal markers:\n")
    for (marker in neuronal_markers) {
      if (marker %in% names(cluster_expr)) {
        cat(sprintf("  %s: %.2f\n", marker, cluster_expr[marker]))
      }
    }
    
    # Progenitor markers
    cat("\nProgenitor markers:\n")
    for (marker in progenitor_markers) {
      if (marker %in% names(cluster_expr)) {
        cat(sprintf("  %s: %.2f\n", marker, cluster_expr[marker]))
      }
    }
    
    # DA markers (for cluster 5)
    if (cluster == "5") {
      cat("\nDopaminergic markers:\n")
      for (marker in da_markers) {
        if (marker %in% names(cluster_expr)) {
          cat(sprintf("  %s: %.2f\n", marker, cluster_expr[marker]))
        }
      }
    }
    
    # CP markers (for clusters 4 and 5)
    if (cluster %in% c("4", "5")) {
      cat("\nChoroid plexus markers:\n")
      for (marker in cp_markers) {
        if (marker %in% names(cluster_expr)) {
          cat(sprintf("  %s: %.2f\n", marker, cluster_expr[marker]))
        }
      }
    }
    
    # Top 10 genes
    top_genes <- sort(cluster_expr, decreasing = TRUE)[1:10]
    cat("\nTop 10 expressed genes:\n")
    for (i in 1:10) {
      cat(sprintf("  %d. %s: %.2f\n", i, names(top_genes)[i], top_genes[i]))
    }
  }
}

# Create a quick visualization
cat("\n\nCreating marker heatmap...\n")
library(pheatmap)

# Select clusters to visualize
clusters_to_plot <- c("0", "1", "3", "4", "5", "13", "18", "28")
markers_to_plot <- c("TUBB3", "STMN2", "SOX2", "NES", "TH", "DDC", "TTR", "FOLR1")

# Get available markers
available_markers <- markers_to_plot[markers_to_plot %in% rownames(avg_expr)]
available_clusters <- clusters_to_plot[clusters_to_plot %in% colnames(avg_expr)]

if (length(available_markers) > 0 && length(available_clusters) > 0) {
  expr_subset <- avg_expr[available_markers, available_clusters]
  
  # Add cell type annotations
  cluster_labels <- sapply(available_clusters, function(cl) {
    cells <- which(x$seurat_clusters_fine == cl)
    if (length(cells) > 0) {
      label <- unique(x$celltypes_fine[cells])[1]
      paste0("C", cl, ": ", substr(label, 1, 20))
    } else {
      paste0("C", cl)
    }
  })
  
  colnames(expr_subset) <- cluster_labels
  
  # Save heatmap
  pdf("results/metadata_validation/marker_expression_heatmap.pdf", width = 10, height = 6)
  pheatmap(log2(expr_subset + 1),
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           main = "Key Marker Expression by Cluster",
           color = colorRampPalette(c("white", "yellow", "red"))(100))
  dev.off()
  
  cat("Heatmap saved to: results/metadata_validation/marker_expression_heatmap.pdf\n")
}

cat("\nDone!\n")