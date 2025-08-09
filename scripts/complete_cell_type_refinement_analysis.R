#!/usr/bin/env Rscript

# Complete Cell Type Refinement Analysis Pipeline
# This script combines all three steps into one comprehensive analysis

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)

cat("================================================================\n")
cat("Complete Cell Type Refinement Analysis with Vulnerability Status\n")
cat("================================================================\n\n")

# Step 1: Load and prepare data
cat("Step 1: Loading data...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds")
current_annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv")

cat("- Loaded", ncol(seurat_obj), "cells\n")
cat("- Found", length(unique(seurat_obj$seurat_clusters_fine)), "fine clusters\n")
cat("- Annotations available for clusters 0-", max(current_annotations$fine_cluster), "\n")

# Step 2: Extract positive and negative markers
cat("\nStep 2: Extracting positive and negative markers...\n")

# Load pre-calculated markers
marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
all_markers_df <- readRDS(marker_file)

# Process markers for each cluster
all_pos_neg_markers <- list()
clusters_to_process <- 0:35  # Only process clusters with annotations

for (cluster_id in clusters_to_process) {
  cat("\r  Processing cluster", cluster_id, "...    ")
  
  # Get positive markers from pre-calculated results
  pos_markers <- all_markers_df %>%
    filter(cluster == as.character(cluster_id),
           avg_log2FC > 0.5,
           p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC)) %>%
    head(50)
  
  # For negative markers, find genes with low expression in this cluster
  # but high expression in other clusters
  other_cluster_markers <- all_markers_df %>%
    filter(cluster != as.character(cluster_id),
           avg_log2FC > 1.5,
           p_val_adj < 0.01) %>%
    group_by(gene) %>%
    summarise(n_clusters = n(), .groups = 'drop') %>%
    filter(n_clusters >= 3) %>%  # Gene must be marker in at least 3 other clusters
    pull(gene)
  
  # Exclude this cluster's top markers from negative list
  neg_genes <- setdiff(other_cluster_markers, pos_markers$gene[1:10])
  
  neg_markers <- data.frame(
    gene = head(neg_genes, 50),
    cluster = cluster_id,
    direction = "negative"
  )
  
  all_pos_neg_markers[[as.character(cluster_id)]] <- list(
    positive = pos_markers,
    negative = neg_markers
  )
}

cat("\n  Marker extraction complete!\n")

# Step 3: Apply vulnerability assessment
cat("\nStep 3: Assessing vulnerability patterns...\n")

# Define vulnerability markers
vulnerability_markers <- list(
  vulnerable_da = list(
    positive = c("SOX6", "ALDH1A1", "KCNJ6", "GIRK2", "GRM5"),
    negative = c("CALB1", "CALB2", "OTX2")
  ),
  resilient_da = list(
    positive = c("CALB1", "CALB2", "OTX2", "FOXP2", "SORCS3"),
    negative = c("SOX6", "KCNJ6")
  ),
  vulnerable_glut = list(
    positive = c("MEIS1", "MEIS2", "SLC17A6"),
    negative = c("LMX1A", "FOXP2")
  )
)

vulnerability_assessment <- list()

for (cluster_id in names(all_pos_neg_markers)) {
  pos_genes <- all_pos_neg_markers[[cluster_id]]$positive$gene
  neg_genes <- all_pos_neg_markers[[cluster_id]]$negative$gene
  
  scores <- list()
  for (pattern_name in names(vulnerability_markers)) {
    pattern <- vulnerability_markers[[pattern_name]]
    pos_score <- sum(pattern$positive %in% pos_genes)
    neg_score <- sum(pattern$negative %in% neg_genes)
    scores[[pattern_name]] <- pos_score + neg_score
  }
  
  best_match <- names(which.max(unlist(scores)))
  vulnerability_assessment[[cluster_id]] <- ifelse(grepl("vulnerable", best_match), "Vulnerable", "Resilient")
}

# Step 4: Apply classification refinement rules
cat("\nStep 4: Refining cell type classifications...\n")

# Check which key markers exist in the dataset
available_genes <- rownames(seurat_obj)
key_markers <- c("TH", "DDC", "SOX6", "ALDH1A1", "KCNJ6", "CALB1", "OTX2", 
                "FOXA2", "LMX1A", "NR4A2", "SLC6A3", "SLC18A2", 
                "MEIS1", "MEIS2", "SLC17A6", "POU3F2", "NKX2.1")
available_key_markers <- key_markers[key_markers %in% available_genes]
cat("  Available key markers:", length(available_key_markers), "of", length(key_markers), "\n")

# Define simplified rules based on available markers
refined_annotations <- current_annotations
refined_annotations$refined_type <- refined_annotations$cell_type
refined_annotations$refined_confidence <- refined_annotations$confidence
refined_annotations$da_subtype <- ""
refined_annotations$vulnerability <- ""
refined_annotations$evidence <- ""

for (i in 1:nrow(refined_annotations)) {
  cluster_id <- as.character(refined_annotations$fine_cluster[i])
  
  if (cluster_id %in% names(all_pos_neg_markers)) {
    pos_genes <- all_pos_neg_markers[[cluster_id]]$positive$gene
    neg_genes <- all_pos_neg_markers[[cluster_id]]$negative$gene
    
    # Check for dopaminergic neurons
    if ("TH" %in% pos_genes && "TH" %in% available_genes) {
      if ("SOX6" %in% pos_genes && "ALDH1A1" %in% pos_genes && !"CALB1" %in% pos_genes) {
        refined_annotations$refined_type[i] <- "Dopaminergic Neurons (A9-like)"
        refined_annotations$da_subtype[i] <- "A9 SNc (vulnerable)"
        refined_annotations$refined_confidence[i] <- "High"
        refined_annotations$evidence[i] <- "TH+/SOX6+/ALDH1A1+/CALB1-"
      } else if ("CALB1" %in% pos_genes && !"SOX6" %in% pos_genes) {
        refined_annotations$refined_type[i] <- "Dopaminergic Neurons (A10-like)"
        refined_annotations$da_subtype[i] <- "A10 VTA (resilient)"
        refined_annotations$refined_confidence[i] <- "High"
        refined_annotations$evidence[i] <- "TH+/CALB1+/SOX6-"
      } else {
        refined_annotations$da_subtype[i] <- "Unspecified DA"
      }
    }
    
    # Check for floor plate refinements
    if ("FOXA2" %in% pos_genes && "FOXA2" %in% available_genes) {
      if ("SHH" %in% pos_genes && !"TH" %in% pos_genes) {
        refined_annotations$refined_type[i] <- "Early Floor Plate Progenitors"
        refined_annotations$refined_confidence[i] <- "High"
        refined_annotations$evidence[i] <- paste(refined_annotations$evidence[i], "FOXA2+/SHH+/TH-")
      } else if ("NR4A2" %in% pos_genes && !"TH" %in% pos_genes) {
        refined_annotations$refined_type[i] <- "Late Floor Plate Progenitors"
        refined_annotations$refined_confidence[i] <- "High"
        refined_annotations$evidence[i] <- paste(refined_annotations$evidence[i], "FOXA2+/NR4A2+/TH-")
      }
    }
    
    # Add vulnerability
    refined_annotations$vulnerability[i] <- vulnerability_assessment[[cluster_id]]
  }
}

# Count changes
n_changed <- sum(refined_annotations$refined_type != refined_annotations$cell_type)
n_improved_conf <- sum(refined_annotations$refined_confidence == "High" & refined_annotations$confidence != "High")
cat("  Classifications changed:", n_changed, "\n")
cat("  Confidence improved:", n_improved_conf, "\n")

# Step 5: Update Seurat object with refined annotations
cat("\nStep 5: Updating Seurat object...\n")

# Create lookup tables
refined_lookup <- setNames(refined_annotations$refined_type, as.character(refined_annotations$fine_cluster))
vulnerability_lookup <- setNames(refined_annotations$vulnerability, as.character(refined_annotations$fine_cluster))
da_subtype_lookup <- setNames(refined_annotations$da_subtype, as.character(refined_annotations$fine_cluster))
confidence_lookup <- setNames(refined_annotations$refined_confidence, as.character(refined_annotations$fine_cluster))

# Apply to Seurat object
seurat_obj$refined_celltype <- refined_lookup[as.character(seurat_obj$seurat_clusters_fine)]
seurat_obj$vulnerability_status <- vulnerability_lookup[as.character(seurat_obj$seurat_clusters_fine)]
seurat_obj$da_subtype <- da_subtype_lookup[as.character(seurat_obj$seurat_clusters_fine)]
seurat_obj$refined_confidence <- confidence_lookup[as.character(seurat_obj$seurat_clusters_fine)]

# Handle NAs
seurat_obj$refined_celltype[is.na(seurat_obj$refined_celltype)] <- "Unknown"
seurat_obj$vulnerability_status[is.na(seurat_obj$vulnerability_status)] <- "Unknown"
seurat_obj$da_subtype[is.na(seurat_obj$da_subtype)] <- ""

# Step 6: Create comprehensive visualizations
cat("\nStep 6: Creating visualizations...\n")
dir.create("results/cell_type_annotations/complete_refined_plots", recursive = TRUE, showWarnings = FALSE)

# 1. Original vs refined cell types
p1a <- DimPlot(seurat_obj, group.by = "celltypes_fine_major", label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Original Cell Type Classifications") + NoLegend()

p1b <- DimPlot(seurat_obj, group.by = "refined_celltype", label = TRUE, repel = TRUE, pt.size = 0.3) +
  ggtitle("Refined Cell Type Classifications") + NoLegend()

p1 <- p1a + p1b
ggsave("results/cell_type_annotations/complete_refined_plots/comparison_original_vs_refined.pdf", p1, width = 20, height = 10)

# 2. Vulnerability status
p2 <- DimPlot(seurat_obj,
              group.by = "vulnerability_status",
              cols = c("Vulnerable" = "red", "Resilient" = "blue", "Unknown" = "gray90"),
              pt.size = 0.5) +
  ggtitle("Vulnerability Status (MPTP Model)") +
  theme(legend.position = "right")

ggsave("results/cell_type_annotations/complete_refined_plots/umap_vulnerability.pdf", p2, width = 10, height = 8)

# 3. Dopaminergic subtypes
da_cells <- WhichCells(seurat_obj, expression = da_subtype != "")
if (length(da_cells) > 0) {
  p3 <- DimPlot(seurat_obj,
                cells = da_cells,
                group.by = "da_subtype",
                pt.size = 1.5) +
    ggtitle(paste("Dopaminergic Neuron Subtypes (n =", length(da_cells), "cells)")) +
    theme(legend.position = "right")
  
  ggsave("results/cell_type_annotations/complete_refined_plots/umap_da_subtypes.pdf", p3, width = 10, height = 8)
}

# 4. Feature plots for vulnerability markers
vulnerability_markers <- c("SOX6", "ALDH1A1", "KCNJ6", "CALB1", "MEIS1", "MEIS2", "TH", "DDC", "OTX2")
available_markers <- vulnerability_markers[vulnerability_markers %in% rownames(seurat_obj)]

if (length(available_markers) >= 4) {
  p4 <- FeaturePlot(seurat_obj, 
                   features = available_markers[1:min(6, length(available_markers))],
                   ncol = 3,
                   pt.size = 0.1)
  
  ggsave("results/cell_type_annotations/complete_refined_plots/vulnerability_markers_features.pdf",
         p4, width = 15, height = 10)
}

# 5. Summary statistics visualization
summary_data <- data.frame(
  Category = c("Cell Types", "Vulnerability", "DA Subtypes", "High Confidence"),
  Original = c(
    length(unique(seurat_obj$celltypes_fine_major)),
    sum(seurat_obj$vulnerability_status == "Vulnerable", na.rm = TRUE),
    sum(refined_annotations$cell_type == "Dopaminergic Neurons"),
    sum(refined_annotations$confidence == "High")
  ),
  Refined = c(
    length(unique(seurat_obj$refined_celltype)),
    sum(seurat_obj$vulnerability_status == "Vulnerable", na.rm = TRUE),
    sum(seurat_obj$da_subtype != "", na.rm = TRUE),
    sum(refined_annotations$refined_confidence == "High")
  )
)

# 6. Create summary table
cat("\n=== ANALYSIS SUMMARY ===\n")
print(summary_data)

# 7. Save all results
cat("\nStep 7: Saving results...\n")

# Save refined annotations
write.csv(refined_annotations, 
          "results/cell_type_annotations/complete_refined_annotations.csv",
          row.names = FALSE)

# Save Seurat object with all annotations
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fully_refined.rds")

# Save marker analysis
saveRDS(all_pos_neg_markers, "results/cell_type_annotations/complete_pos_neg_markers.rds")

# Create final summary report
sink("results/cell_type_annotations/refinement_summary_report.txt")
cat("Cell Type Refinement Analysis Summary\n")
cat("=====================================\n\n")
cat("Total cells analyzed:", ncol(seurat_obj), "\n")
cat("Total clusters:", length(unique(seurat_obj$seurat_clusters_fine)), "\n\n")

cat("Vulnerability Assessment:\n")
print(table(seurat_obj$vulnerability_status))
cat("\n")

cat("Dopaminergic Subtypes Identified:\n")
da_summary <- table(seurat_obj$da_subtype[seurat_obj$da_subtype != ""])
if (length(da_summary) > 0) print(da_summary)
cat("\n")

cat("Cell Type Changes:\n")
for (i in 1:nrow(refined_annotations)) {
  if (refined_annotations$refined_type[i] != refined_annotations$cell_type[i]) {
    cat("  Cluster", refined_annotations$fine_cluster[i], ":",
        refined_annotations$cell_type[i], "->", refined_annotations$refined_type[i],
        "(", refined_annotations$evidence[i], ")\n")
  }
}
sink()

cat("\n=== COMPLETE! ===\n")
cat("All results saved to: results/cell_type_annotations/complete_refined_plots/\n")
cat("Refined Seurat object: iSCORE-PD_plus_CRISPRi_fully_refined.rds\n")
cat("Summary report: results/cell_type_annotations/refinement_summary_report.txt\n")

# Return results
invisible(list(
  annotations = refined_annotations,
  markers = all_pos_neg_markers,
  vulnerability = vulnerability_assessment,
  seurat_obj = seurat_obj
))