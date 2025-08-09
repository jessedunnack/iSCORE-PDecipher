#!/usr/bin/env Rscript

# WORKING Complete Cell Type Refinement Analysis
# Fixed version that properly handles the metadata assignment

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)

cat("================================================================\n")
cat("Complete Cell Type Refinement Analysis with Vulnerability Status\n")
cat("================================================================\n\n")

# Step 1: Load data
cat("Step 1: Loading data...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds")
current_annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv")
marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
all_markers_df <- readRDS(marker_file)

cat("- Loaded", ncol(seurat_obj), "cells\n")
cat("- Found", length(unique(seurat_obj$seurat_clusters_fine)), "fine clusters\n")

# Step 2: Process markers and assess vulnerability
cat("\nStep 2: Processing markers and assessing vulnerability...\n")

# Simplified vulnerability assessment based on key markers
vulnerability_status <- character(36)
da_subtypes <- character(36)
refined_types <- current_annotations$cell_type
evidence <- character(36)

for (i in 0:35) {
  cluster_markers <- all_markers_df %>%
    filter(cluster == as.character(i),
           avg_log2FC > 0.5,
           p_val_adj < 0.05) %>%
    pull(gene)
  
  # Default to original annotation
  idx <- which(current_annotations$fine_cluster == i)
  
  # Check for dopaminergic markers
  has_TH <- "TH" %in% cluster_markers
  has_DDC <- "DDC" %in% cluster_markers
  has_SOX6 <- "SOX6" %in% cluster_markers
  has_ALDH1A1 <- "ALDH1A1" %in% cluster_markers
  has_CALB1 <- "CALB1" %in% cluster_markers
  has_KCNJ6 <- "KCNJ6" %in% cluster_markers
  has_MEIS1 <- "MEIS1" %in% cluster_markers
  has_MEIS2 <- "MEIS2" %in% cluster_markers
  
  # Determine vulnerability
  if (has_SOX6 || has_ALDH1A1 || has_KCNJ6 || has_MEIS1 || has_MEIS2) {
    vulnerability_status[i+1] <- "Vulnerable"
  } else if (has_CALB1 || ("FOXP2" %in% cluster_markers)) {
    vulnerability_status[i+1] <- "Resilient"
  } else {
    vulnerability_status[i+1] <- "Unknown"
  }
  
  # Refine dopaminergic classifications
  if (has_TH && has_DDC) {
    if (has_SOX6 && has_ALDH1A1 && !has_CALB1) {
      refined_types[idx] <- "Dopaminergic Neurons (A9-like)"
      da_subtypes[i+1] <- "A9 SNc (vulnerable)"
      evidence[i+1] <- "TH+/DDC+/SOX6+/ALDH1A1+"
    } else if (has_CALB1 && !has_SOX6) {
      refined_types[idx] <- "Dopaminergic Neurons (A10-like)"
      da_subtypes[i+1] <- "A10 VTA (resilient)"
      evidence[i+1] <- "TH+/DDC+/CALB1+/SOX6-"
    } else if (has_TH) {
      da_subtypes[i+1] <- "DA (unspecified)"
      evidence[i+1] <- "TH+/DDC+"
    }
  }
}

# Update annotations
current_annotations$refined_type <- refined_types
current_annotations$vulnerability <- vulnerability_status
current_annotations$da_subtype <- da_subtypes
current_annotations$evidence <- evidence

# Step 3: Apply to Seurat object - FIXED VERSION
cat("\nStep 3: Updating Seurat object...\n")

# Initialize new metadata columns
seurat_obj$refined_celltype <- as.character(seurat_obj$celltypes_fine_major)  # Start with original
seurat_obj$vulnerability_status <- "Unknown"
seurat_obj$da_subtype <- ""

# Apply annotations cluster by cluster
for (i in 0:35) {
  cells_in_cluster <- WhichCells(seurat_obj, expression = seurat_clusters_fine == i)
  
  if (length(cells_in_cluster) > 0) {
    idx <- which(current_annotations$fine_cluster == i)
    
    seurat_obj@meta.data[cells_in_cluster, "refined_celltype"] <- current_annotations$refined_type[idx]
    seurat_obj@meta.data[cells_in_cluster, "vulnerability_status"] <- vulnerability_status[i+1]
    seurat_obj@meta.data[cells_in_cluster, "da_subtype"] <- da_subtypes[i+1]
  }
}

# Print summary
cat("\nSummary of assignments:\n")
cat("- Vulnerability:", table(seurat_obj$vulnerability_status), "\n")
cat("- DA subtypes:", sum(seurat_obj$da_subtype != ""), "cells\n")

# Step 4: Create visualizations
cat("\nStep 4: Creating visualizations...\n")
dir.create("results/cell_type_annotations/final_refined_plots", recursive = TRUE, showWarnings = FALSE)

# 1. Vulnerability UMAP
p1 <- DimPlot(seurat_obj,
              group.by = "vulnerability_status",
              cols = c("Vulnerable" = "red", 
                      "Resilient" = "blue",
                      "Unknown" = "gray90"),
              pt.size = 0.5) +
  ggtitle("Vulnerability Status (Based on MPTP Markers)") +
  theme(legend.position = "right")

ggsave("results/cell_type_annotations/final_refined_plots/umap_vulnerability.pdf", p1, width = 10, height = 8)

# 2. Refined cell types
p2 <- DimPlot(seurat_obj,
              group.by = "refined_celltype",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.3) +
  ggtitle("Refined Cell Type Classifications") +
  NoLegend()

ggsave("results/cell_type_annotations/final_refined_plots/umap_refined_celltypes.pdf", p2, width = 12, height = 10)

# 3. DA subtypes if any
da_cells <- WhichCells(seurat_obj, expression = da_subtype != "")
if (length(da_cells) > 0) {
  p3 <- DimPlot(seurat_obj,
                cells = da_cells,
                group.by = "da_subtype",
                pt.size = 1.5) +
    ggtitle(paste("Dopaminergic Subtypes (n =", length(da_cells), ")"))
  
  ggsave("results/cell_type_annotations/final_refined_plots/umap_da_subtypes.pdf", p3, width = 10, height = 8)
}

# 4. Feature plots
markers_to_plot <- c("SOX6", "ALDH1A1", "KCNJ6", "CALB1", "MEIS1", "MEIS2")
available <- markers_to_plot[markers_to_plot %in% rownames(seurat_obj)]

if (length(available) > 0) {
  p4 <- FeaturePlot(seurat_obj, 
                   features = available,
                   ncol = 3,
                   pt.size = 0.1)
  
  ggsave("results/cell_type_annotations/final_refined_plots/vulnerability_markers.pdf",
         p4, width = 15, height = 10)
}

# 5. Summary by cell type
summary_plot <- seurat_obj@meta.data %>%
  group_by(refined_celltype, vulnerability_status) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ggplot(aes(x = refined_celltype, y = count, fill = vulnerability_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Vulnerable" = "red", "Resilient" = "blue", "Unknown" = "gray90")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Vulnerability by Cell Type", y = "Proportion") +
  coord_flip()

ggsave("results/cell_type_annotations/final_refined_plots/vulnerability_by_celltype.pdf", summary_plot, width = 10, height = 8)

# Step 5: Save results
cat("\nStep 5: Saving results...\n")

# Save annotations
write.csv(current_annotations, 
          "results/cell_type_annotations/final_refined_annotations.csv",
          row.names = FALSE)

# Save Seurat object
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_final_refined.rds")

# Create summary report
sink("results/cell_type_annotations/final_summary_report.txt")
cat("FINAL CELL TYPE REFINEMENT SUMMARY\n")
cat("==================================\n\n")
cat("Total cells:", ncol(seurat_obj), "\n")
cat("Total clusters:", length(unique(seurat_obj$seurat_clusters_fine)), "\n\n")

cat("Vulnerability Assessment:\n")
print(table(seurat_obj$vulnerability_status))
cat("\n")

cat("Dopaminergic Subtypes:\n")
da_summary <- table(seurat_obj$da_subtype[seurat_obj$da_subtype != ""])
if (length(da_summary) > 0) print(da_summary)
cat("\n")

cat("Refined Classifications:\n")
for (i in 1:nrow(current_annotations)) {
  if (current_annotations$refined_type[i] != current_annotations$cell_type[i]) {
    cat("  Cluster", current_annotations$fine_cluster[i], ":",
        current_annotations$cell_type[i], "->", 
        current_annotations$refined_type[i], "\n")
  }
}
sink()

cat("\n=== COMPLETE! ===\n")
cat("Results saved to: results/cell_type_annotations/final_refined_plots/\n")
cat("Seurat object: iSCORE-PD_plus_CRISPRi_final_refined.rds\n")
cat("Report: results/cell_type_annotations/final_summary_report.txt\n")

# Return success
invisible(TRUE)