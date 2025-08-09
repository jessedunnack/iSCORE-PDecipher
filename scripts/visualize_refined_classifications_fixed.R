#!/usr/bin/env Rscript

# Fixed visualization script

library(Seurat)
library(ggplot2)
library(dplyr)

# Load data
seurat_path <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds"
seurat_obj <- readRDS(seurat_path)
refined <- read.csv("results/cell_type_annotations/refined_classifications.csv")

# Fix the indexing issue - subtract 1 from fine_cluster to match 0-based Seurat clustering
refined$fine_cluster_zero <- refined$fine_cluster - 1

# Create lookups using 0-based indexing
refined_lookup <- setNames(refined$refined_type, as.character(refined$fine_cluster_zero))
vulnerability_lookup <- setNames(refined$vulnerability_status, as.character(refined$fine_cluster_zero))
da_subtype_lookup <- setNames(refined$da_subtype, as.character(refined$fine_cluster_zero))

# Add to Seurat object
seurat_obj$refined_celltype <- refined_lookup[as.character(seurat_obj$seurat_clusters_fine)]
seurat_obj$vulnerability <- vulnerability_lookup[as.character(seurat_obj$seurat_clusters_fine)]
seurat_obj$da_subtype <- da_subtype_lookup[as.character(seurat_obj$seurat_clusters_fine)]

# Handle NAs
seurat_obj$refined_celltype[is.na(seurat_obj$refined_celltype)] <- "Unknown"
seurat_obj$vulnerability[is.na(seurat_obj$vulnerability)] <- "Unknown"
seurat_obj$da_subtype[is.na(seurat_obj$da_subtype)] <- ""

# Create output directory
dir.create("results/cell_type_annotations/refined_plots", recursive = TRUE, showWarnings = FALSE)

# Create plots
cat("Creating UMAP plots...\n")

# 1. Refined cell types
p1 <- DimPlot(seurat_obj, 
              group.by = "refined_celltype",
              label = TRUE,
              repel = TRUE,
              pt.size = 0.3) +
  ggtitle("Refined Cell Type Classifications") +
  theme(legend.position = "right")

ggsave("results/cell_type_annotations/refined_plots/umap_refined_types.pdf", p1, width = 12, height = 10)

# 2. Vulnerability status
p2 <- DimPlot(seurat_obj,
              group.by = "vulnerability",
              cols = c("Vulnerable" = "red", 
                      "Resilient" = "blue",
                      "Unknown" = "gray"),
              pt.size = 0.3) +
  ggtitle("Vulnerability Status (MPTP Model)")

ggsave("results/cell_type_annotations/refined_plots/umap_vulnerability.pdf", p2, width = 10, height = 8)

# 3. DA subtypes only (if any exist)
da_cells <- WhichCells(seurat_obj, expression = da_subtype != "")
if (length(da_cells) > 0) {
  p3 <- DimPlot(seurat_obj,
                cells = da_cells,
                group.by = "da_subtype",
                pt.size = 1) +
    ggtitle("Dopaminergic Neuron Subtypes") +
    theme(legend.position = "right")
  
  ggsave("results/cell_type_annotations/refined_plots/umap_da_subtypes.pdf", p3, width = 10, height = 8)
}

# 4. Summary plots
cat("Creating summary plots...\n")

# Confidence comparison
conf_data <- refined %>%
  select(fine_cluster, original_confidence, refined_confidence) %>%
  pivot_longer(cols = c(original_confidence, refined_confidence),
              names_to = "type",
              values_to = "confidence") %>%
  mutate(type = gsub("_confidence", "", type))

p4 <- ggplot(conf_data, aes(x = type, fill = confidence)) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = c("High" = "darkgreen", 
                             "Medium" = "orange", 
                             "Low" = "red")) +
  labs(title = "Confidence Level Comparison",
       x = "Classification Type",
       y = "Proportion") +
  theme_minimal()

ggsave("results/cell_type_annotations/refined_plots/summary_confidence.pdf", p4, width = 8, height = 6)

# Vulnerability by cell type
vuln_summary <- refined %>%
  group_by(refined_type, vulnerability_status) %>%
  summarise(count = n(), .groups = 'drop') %>%
  filter(!is.na(vulnerability_status))

p5 <- ggplot(vuln_summary, 
            aes(x = refined_type, 
                y = count, 
                fill = vulnerability_status)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("Vulnerable" = "red", 
                             "Resilient" = "blue",
                             "Unknown" = "gray")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Vulnerability Status by Cell Type",
       x = "Cell Type",
       y = "Number of Clusters")

ggsave("results/cell_type_annotations/refined_plots/summary_vulnerability_by_type.pdf", p5, width = 10, height = 8)

# 5. Feature plots for key markers
cat("Creating feature plots...\n")
vulnerability_markers <- c("SOX6", "ALDH1A1", "KCNJ6", "MEIS1", "MEIS2", "CALB1")
vulnerability_markers <- vulnerability_markers[vulnerability_markers %in% rownames(seurat_obj)]

if (length(vulnerability_markers) > 0) {
  p6 <- FeaturePlot(seurat_obj, 
                   features = vulnerability_markers,
                   ncol = 3,
                   pt.size = 0.1)
  
  ggsave("results/cell_type_annotations/refined_plots/vulnerability_markers_features.pdf",
         p6, width = 15, height = 10)
}

cat("\nVisualization complete! Plots saved to:\n")
cat("  results/cell_type_annotations/refined_plots/\n")

# Print summary
cat("\nSummary:\n")
cat("- Total clusters with vulnerability assessment:", sum(!is.na(seurat_obj$vulnerability)), "\n")
cat("- Vulnerable clusters:", sum(seurat_obj$vulnerability == "Vulnerable", na.rm = TRUE), "\n")
cat("- Resilient clusters:", sum(seurat_obj$vulnerability == "Resilient", na.rm = TRUE), "\n")
cat("- DA neurons identified:", sum(seurat_obj$da_subtype != "", na.rm = TRUE), "\n")