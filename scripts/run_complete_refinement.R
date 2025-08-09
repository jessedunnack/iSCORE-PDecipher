#!/usr/bin/env Rscript

# Complete refinement pipeline that works

library(Seurat)
library(dplyr)
library(ggplot2)

cat("=== Complete Refinement Pipeline ===\n\n")

# Step 1: Load the annotated Seurat object
cat("Loading Seurat object...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds")
cat("Loaded", ncol(seurat_obj), "cells with", length(unique(seurat_obj$seurat_clusters_fine)), "fine clusters\n")

# Step 2: Load vulnerability assessment
cat("\nLoading vulnerability assessment...\n")
vulnerability <- readRDS("results/pos_neg_markers/vulnerability_assessment.rds")

# Step 3: Create simple vulnerability annotations
cat("Adding vulnerability status to Seurat object...\n")

# Initialize all cells as Unknown
seurat_obj$vulnerability_status <- "Unknown"

# Map vulnerability to cells based on cluster
for (cluster_id in names(vulnerability)) {
  cluster_num <- as.numeric(cluster_id)
  if (cluster_num <= 35) {  # Only process clusters 0-35
    cells_in_cluster <- WhichCells(seurat_obj, expression = seurat_clusters_fine == cluster_num)
    vuln_status <- vulnerability[[cluster_id]]$vulnerability
    seurat_obj$vulnerability_status[cells_in_cluster] <- vuln_status
  }
}

# Print summary
cat("\nVulnerability Summary:\n")
print(table(seurat_obj$vulnerability_status))

# Step 4: Create visualizations
cat("\nCreating visualizations...\n")
dir.create("results/cell_type_annotations/refined_plots", recursive = TRUE, showWarnings = FALSE)

# UMAP with vulnerability status
p1 <- DimPlot(seurat_obj,
              group.by = "vulnerability_status",
              cols = c("Vulnerable" = "red", 
                      "Resilient" = "blue",
                      "Unknown" = "gray"),
              pt.size = 0.3) +
  ggtitle("Vulnerability Status Based on MPTP Markers") +
  theme(legend.position = "right")

ggsave("results/cell_type_annotations/refined_plots/umap_vulnerability.pdf", p1, width = 10, height = 8)

# UMAP with original cell types colored by vulnerability
p2 <- DimPlot(seurat_obj,
              group.by = "celltypes_fine_major",
              split.by = "vulnerability_status",
              pt.size = 0.3) +
  ggtitle("Cell Types Split by Vulnerability Status")

ggsave("results/cell_type_annotations/refined_plots/umap_celltypes_by_vulnerability.pdf", p2, width = 18, height = 6)

# Feature plots for key vulnerability markers
vulnerability_markers <- c("SOX6", "ALDH1A1", "KCNJ6", "MEIS1", "MEIS2", "CALB1")
vulnerability_markers <- vulnerability_markers[vulnerability_markers %in% rownames(seurat_obj)]

if (length(vulnerability_markers) > 0) {
  p3 <- FeaturePlot(seurat_obj, 
                   features = vulnerability_markers,
                   ncol = 3,
                   pt.size = 0.1)
  
  ggsave("results/cell_type_annotations/refined_plots/vulnerability_markers_features.pdf",
         p3, width = 15, height = 10)
  cat("Created feature plots for:", paste(vulnerability_markers, collapse = ", "), "\n")
}

# Summary by cell type
vuln_by_type <- seurat_obj@meta.data %>%
  group_by(celltypes_fine_major, vulnerability_status) %>%
  summarise(count = n(), .groups = 'drop')

p4 <- ggplot(vuln_by_type, 
            aes(x = celltypes_fine_major, 
                y = count, 
                fill = vulnerability_status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Vulnerable" = "red", 
                             "Resilient" = "blue",
                             "Unknown" = "gray")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Vulnerability Proportion by Cell Type",
       x = "Cell Type",
       y = "Proportion") +
  coord_flip()

ggsave("results/cell_type_annotations/refined_plots/vulnerability_by_celltype.pdf", p4, width = 10, height = 8)

# Save the updated Seurat object
cat("\nSaving Seurat object with vulnerability annotations...\n")
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_vulnerability_annotated.rds")

cat("\n=== Complete! ===\n")
cat("Visualizations saved to: results/cell_type_annotations/refined_plots/\n")
cat("Annotated object saved as: iSCORE-PD_plus_CRISPRi_vulnerability_annotated.rds\n")

# Return summary
summary_data <- list(
  vulnerability_counts = table(seurat_obj$vulnerability_status),
  vulnerable_clusters = names(vulnerability)[sapply(vulnerability, function(x) x$vulnerability == "Vulnerable")],
  resilient_clusters = names(vulnerability)[sapply(vulnerability, function(x) x$vulnerability == "Resilient")]
)

print(summary_data)