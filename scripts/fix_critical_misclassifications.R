#!/usr/bin/env Rscript

# FIX CRITICAL MISCLASSIFICATIONS
# Based on marker expression analysis

library(Seurat)
library(dplyr)

cat("=================================================================\n")
cat("FIXING CRITICAL CELL TYPE MISCLASSIFICATIONS\n")
cat("=================================================================\n\n")

# Load object
x <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")

# Create new corrected celltype columns
x$corrected_celltype <- x$validated_celltype
x$corrected_celltype_major <- x$validated_celltype_major

# Fix Cluster 1: Floor Plate Progenitors -> Floor Plate-Derived Neurons
cat("FIXING Cluster 1: Floor Plate Progenitors -> Floor Plate-Derived Neurons\n")
cluster1_cells <- WhichCells(x, expression = seurat_clusters_fine == 1)
x$corrected_celltype[cluster1_cells] <- "Floor Plate-Derived Neurons"
x$corrected_celltype_major[cluster1_cells] <- "Neurons"
cat("  Updated", length(cluster1_cells), "cells\n")

# Fix Cluster 5: Dopaminergic Neurons -> Choroid Plexus
cat("\nFIXING Cluster 5: Dopaminergic Neurons -> Choroid Plexus\n")
cluster5_cells <- WhichCells(x, expression = seurat_clusters_fine == 5)
x$corrected_celltype[cluster5_cells] <- "Choroid Plexus (TTR-high)"
x$corrected_celltype_major[cluster5_cells] <- "Choroid Plexus"
cat("  Updated", length(cluster5_cells), "cells\n")

# Check other Floor Plate Progenitor clusters for neuronal markers
cat("\nChecking other clusters labeled as progenitors...\n")

# Get average expression
if ("SCT" %in% names(x@assays)) {
  avg_expr <- AverageExpression(x, 
                                features = c("TUBB3", "STMN2", "GAP43", "SOX2", "NES"),
                                group.by = "seurat_clusters_fine", 
                                assays = "SCT")$SCT
} else {
  avg_expr <- AverageExpression(x, 
                                features = c("TUBB3", "STMN2", "GAP43", "SOX2", "NES"),
                                group.by = "seurat_clusters_fine", 
                                assays = "RNA")$RNA
}

# Check each cluster
for (cluster in 0:35) {
  cluster_cells <- WhichCells(x, expression = seurat_clusters_fine == cluster)
  if (length(cluster_cells) > 0) {
    current_label <- unique(x$validated_celltype[cluster_cells])[1]
    
    # If labeled as progenitor but has neuronal markers
    if (!is.na(current_label) && grepl("Progenitor", current_label) && 
        cluster %in% colnames(avg_expr)) {
      
      neuronal_score <- sum(avg_expr[c("TUBB3", "STMN2", "GAP43"), as.character(cluster)])
      progenitor_score <- sum(avg_expr[c("SOX2", "NES"), as.character(cluster)])
      
      if (neuronal_score > 10 && neuronal_score > progenitor_score * 10) {
        cat("\n  Cluster", cluster, "(", current_label, "):")
        cat("\n    Neuronal score:", round(neuronal_score, 2))
        cat("\n    Progenitor score:", round(progenitor_score, 2))
        cat("\n    -> Reclassifying as neurons\n")
        
        x$corrected_celltype[cluster_cells] <- gsub("Progenitors", "Neurons", current_label)
        x$corrected_celltype_major[cluster_cells] <- "Neurons"
      }
    }
  }
}

# Check for other potential DA clusters
cat("\nLooking for true dopaminergic neurons...\n")
if ("TH" %in% rownames(x)) {
  assay_name <- ifelse("SCT" %in% names(x@assays), "SCT", "RNA")
  th_expr <- AverageExpression(x, features = "TH", group.by = "seurat_clusters_fine", 
                               assays = assay_name)[[assay_name]]
  
  # Find clusters with high TH - check if TH row exists
  if ("TH" %in% rownames(th_expr)) {
    th_values <- th_expr["TH",]
    high_th_clusters <- names(which(th_values > 5))
    
    for (cluster in high_th_clusters) {
      cluster_cells <- WhichCells(x, expression = seurat_clusters_fine == cluster)
      if (length(cluster_cells) > 0) {
        current_label <- unique(x$corrected_celltype[cluster_cells])[1]
        if (!grepl("Dopaminergic", current_label)) {
          cat("  Cluster", cluster, "has high TH expression (", 
              round(th_values[cluster], 2), ") - currently labeled as:", current_label, "\n")
        }
      }
    }
  } else {
    cat("  TH not found in expression matrix\n")
  }
}

# Summary of changes
cat("\n\nSUMMARY OF CORRECTIONS:\n")
cat("=======================\n")

# Compare original vs corrected
comparison <- x@meta.data %>%
  group_by(seurat_clusters_fine) %>%
  summarise(
    n_cells = n(),
    original = first(validated_celltype),
    corrected = first(corrected_celltype),
    changed = first(validated_celltype) != first(corrected_celltype)
  ) %>%
  filter(changed)

print(comparison)

# Save corrected object
cat("\nSaving corrected object...\n")
saveRDS(x, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_corrected.rds")

# Also save a summary table
write.csv(comparison, "results/metadata_validation/celltype_corrections_summary.csv", row.names = FALSE)

cat("\nCorrected object saved as: iSCORE-PD_plus_CRISPRi_corrected.rds\n")
cat("Summary saved to: results/metadata_validation/celltype_corrections_summary.csv\n")

# Create before/after visualization
cat("\nCreating visualization of corrections...\n")

tryCatch({
  # Sample 50,000 cells for visualization
  set.seed(42)
  cells_to_plot <- sample(colnames(x), min(50000, ncol(x)))
  
  # Check if we have UMAP coordinates
  if ("umap" %in% names(x@reductions) || "umap.cca" %in% names(x@reductions)) {
    reduction_to_use <- ifelse("umap.cca" %in% names(x@reductions), "umap.cca", "umap")
    
    p1 <- DimPlot(x[, cells_to_plot], reduction = reduction_to_use,
                  group.by = "validated_celltype", label = TRUE, repel = TRUE) +
      ggtitle("Original Classifications") +
      theme(legend.position = "none")
    
    p2 <- DimPlot(x[, cells_to_plot], reduction = reduction_to_use,
                  group.by = "corrected_celltype", label = TRUE, repel = TRUE) +
      ggtitle("Corrected Classifications") +
      theme(legend.position = "none")
    
    library(patchwork)
    p_combined <- p1 + p2
    ggsave("results/metadata_validation/celltype_corrections_comparison.pdf", p_combined, width = 20, height = 10)
  } else {
    cat("  No UMAP reduction found, skipping visualization\n")
  }
}, error = function(e) {
  cat("  Error creating visualization:", e$message, "\n")
})

cat("\nVisualization saved to: results/metadata_validation/celltype_corrections_comparison.pdf\n")
cat("\nDONE! Please use the corrected object for all downstream analyses.\n")