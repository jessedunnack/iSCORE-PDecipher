#!/usr/bin/env Rscript

# CRITICAL VALIDATION: Check all metadata labels are correctly assigned to clusters
# Multiple independent checks to ensure data integrity

library(Seurat)
library(dplyr)
library(ggplot2)
library(pheatmap)

cat("=====================================================================\n")
cat("CRITICAL METADATA VALIDATION: Ensuring Labels Match Correct Clusters\n")
cat("=====================================================================\n\n")

# Load the Seurat object
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")

# Create output directory for validation plots
dir.create("results/metadata_validation", recursive = TRUE, showWarnings = FALSE)

# ====================
# VALIDATION CHECK 1: Cluster consistency across metadata columns
# ====================
cat("CHECK 1: Verifying cluster assignments are consistent\n")
cat("---------------------------------------------------\n")

# Get all cluster-related metadata columns
cluster_cols <- grep("cluster", colnames(seurat_obj@meta.data), value = TRUE)
celltype_cols <- grep("celltype", colnames(seurat_obj@meta.data), value = TRUE)

# Check unique values in seurat_clusters_fine
unique_clusters <- sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_fine))))
cat("Unique clusters in seurat_clusters_fine:", paste(unique_clusters, collapse=", "), "\n")
cat("Total clusters:", length(unique_clusters), "\n\n")

# Create a cross-tabulation to check alignment
cat("Cross-tabulation of cluster assignments:\n")
if ("seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  cluster_comparison <- table(
    seurat_clusters = seurat_obj$seurat_clusters,
    seurat_clusters_fine = seurat_obj$seurat_clusters_fine
  )
  print(head(cluster_comparison, 10))
}

# ====================
# VALIDATION CHECK 2: Marker gene expression matches cell type labels
# ====================
cat("\n\nCHECK 2: Validating marker expression matches cell type labels\n")
cat("--------------------------------------------------------------\n")

# Define key marker genes for validation
validation_markers <- list(
  Dopaminergic = c("TH", "DDC", "SLC6A3", "SLC18A2"),
  Choroid_Plexus = c("TTR", "FOLR1", "HTR2C", "CLIC6"),
  Neurons = c("TUBB3", "MAP2", "STMN2", "SYN1"),
  Progenitors = c("SOX2", "NES", "VIM", "FABP7"),
  Floor_Plate = c("FOXA2", "CORIN", "ARX", "LMX1A"),
  Oligodendrocytes = c("OLIG2", "SOX10", "MBP", "PLP1"),
  Vascular = c("TAGLN", "ACTA2", "MYL9", "PECAM1"),
  Proliferating = c("MKI67", "TOP2A", "PCNA")
)

# Calculate average expression by cluster
cat("Calculating average marker expression per cluster...\n")
avg_expression_by_cluster <- AverageExpression(
  seurat_obj,
  features = unlist(validation_markers),
  group.by = "seurat_clusters_fine",
  assays = "RNA"
)$RNA

# Create heatmap of marker expression
pdf("results/metadata_validation/marker_expression_by_cluster.pdf", width = 12, height = 10)
pheatmap(
  log2(avg_expression_by_cluster + 1),
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Marker Expression by Cluster",
  fontsize_row = 8,
  fontsize_col = 10
)
dev.off()

# ====================
# VALIDATION CHECK 3: Cell type assignment validation
# ====================
cat("\n\nCHECK 3: Checking cell type assignments match marker expression\n")
cat("---------------------------------------------------------------\n")

# For each major cell type, check which clusters are labeled as such
celltype_summary <- seurat_obj@meta.data %>%
  group_by(seurat_clusters_fine, validated_celltype_major) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  arrange(seurat_clusters_fine)

cat("\nCell type assignments per cluster:\n")
print(celltype_summary)

# Check specific marker expression in assigned cell types
validation_results <- list()

# Check Dopaminergic assignments
if (any(grepl("Dopaminergic", seurat_obj$validated_celltype_major))) {
  da_clusters <- unique(seurat_obj$seurat_clusters_fine[grepl("Dopaminergic", seurat_obj$validated_celltype_major)])
  
  cat("\n\nDopaminergic neuron validation:\n")
  cat("Clusters labeled as Dopaminergic:", paste(da_clusters, collapse=", "), "\n")
  
  # Check TH expression in these clusters
  th_expr <- avg_expression_by_cluster["TH", as.character(da_clusters)]
  cat("TH expression in DA clusters:", paste(round(th_expr, 2), collapse=", "), "\n")
  
  # Flag if TH is low in "DA" clusters
  if (any(th_expr < 1)) {
    cat("WARNING: Low TH expression in some DA clusters!\n")
  }
}

# Check Choroid Plexus assignments
if (any(grepl("Choroid", seurat_obj$validated_celltype_major))) {
  cp_clusters <- unique(seurat_obj$seurat_clusters_fine[grepl("Choroid", seurat_obj$validated_celltype_major)])
  
  cat("\n\nChoroid Plexus validation:\n")
  cat("Clusters labeled as Choroid Plexus:", paste(cp_clusters, collapse=", "), "\n")
  
  # Check TTR expression
  ttr_expr <- avg_expression_by_cluster["TTR", as.character(cp_clusters)]
  cat("TTR expression in CP clusters:", paste(round(ttr_expr, 2), collapse=", "), "\n")
}

# ====================
# VALIDATION CHECK 4: Visual validation with feature plots
# ====================
cat("\n\nCHECK 4: Creating feature plots for visual validation\n")
cat("-----------------------------------------------------\n")

# Create feature plots for key markers
key_markers_to_plot <- c("TH", "TTR", "TUBB3", "SOX2", "MKI67", "FOXA2")
available_markers <- key_markers_to_plot[key_markers_to_plot %in% rownames(seurat_obj)]

p1 <- FeaturePlot(seurat_obj, features = available_markers, ncol = 3, pt.size = 0.1)
ggsave("results/metadata_validation/key_markers_feature_plot.pdf", p1, width = 15, height = 10)

# Create side-by-side DimPlots
p2 <- DimPlot(seurat_obj, group.by = "seurat_clusters_fine", label = TRUE, repel = TRUE) +
  ggtitle("Cluster Numbers")

p3 <- DimPlot(seurat_obj, group.by = "validated_celltype_major", label = TRUE, repel = TRUE) +
  ggtitle("Cell Type Labels") +
  theme(legend.position = "bottom")

p_combined <- p2 + p3
ggsave("results/metadata_validation/clusters_vs_celltypes.pdf", p_combined, width = 20, height = 10)

# ====================
# VALIDATION CHECK 5: Specific cluster-marker validation
# ====================
cat("\n\nCHECK 5: Detailed cluster-by-cluster validation\n")
cat("-----------------------------------------------\n")

# For each cluster, check if its label matches its top markers
problematic_clusters <- c()

for (cluster in unique_clusters) {
  cluster_label <- unique(seurat_obj$validated_celltype[seurat_obj$seurat_clusters_fine == cluster])[1]
  
  # Get top expressed genes in this cluster
  cluster_expr <- avg_expression_by_cluster[, as.character(cluster)]
  top_genes <- names(sort(cluster_expr, decreasing = TRUE)[1:10])
  
  # Check for mismatches
  issues <- c()
  
  # If labeled as DA but no TH
  if (grepl("Dopaminergic", cluster_label) && !("TH" %in% top_genes[1:20])) {
    issues <- c(issues, "Labeled as Dopaminergic but TH not in top 20 genes")
  }
  
  # If labeled as Choroid Plexus but no TTR
  if (grepl("Choroid", cluster_label) && !("TTR" %in% top_genes[1:10])) {
    issues <- c(issues, "Labeled as Choroid Plexus but TTR not in top 10 genes")
  }
  
  # If labeled as Neurons but has progenitor markers
  if (grepl("Neuron", cluster_label) && any(c("SOX2", "NES") %in% top_genes[1:10])) {
    issues <- c(issues, "Labeled as Neurons but has progenitor markers in top 10")
  }
  
  # If labeled as Progenitors but has neuronal markers
  if (grepl("Progenitor", cluster_label) && any(c("TUBB3", "MAP2", "STMN2") %in% top_genes[1:10])) {
    issues <- c(issues, "Labeled as Progenitors but has neuronal markers in top 10")
  }
  
  if (length(issues) > 0) {
    problematic_clusters <- c(problematic_clusters, cluster)
    cat("\nCluster", cluster, "(", cluster_label, "):\n")
    cat("  Issues:", paste(issues, collapse = "; "), "\n")
    cat("  Top 5 markers:", paste(top_genes[1:5], collapse = ", "), "\n")
  }
}

# ====================
# VALIDATION CHECK 6: Metadata alignment check
# ====================
cat("\n\nCHECK 6: Checking metadata alignment\n")
cat("-------------------------------------\n")

# Sample 100 random cells and verify their metadata
set.seed(42)
sample_cells <- sample(colnames(seurat_obj), 100)

sample_check <- seurat_obj@meta.data[sample_cells, c("seurat_clusters_fine", 
                                                     "celltypes_fine_major",
                                                     "validated_celltype_major",
                                                     "vulnerability_status")]

cat("Sample of cell metadata (first 10):\n")
print(head(sample_check, 10))

# Check for any cells with mismatched metadata
mismatches <- seurat_obj@meta.data %>%
  filter(seurat_clusters_fine == 5 & !grepl("Dopaminergic", validated_celltype)) %>%
  nrow()

if (mismatches > 0) {
  cat("\nWARNING: Found", mismatches, "cells in cluster 5 not labeled as Dopaminergic!\n")
}

# ====================
# FINAL SUMMARY
# ====================
cat("\n\n=====================\n")
cat("VALIDATION SUMMARY\n")
cat("=====================\n")

if (length(problematic_clusters) > 0) {
  cat("\nPROBLEMATIC CLUSTERS FOUND:\n")
  cat("Clusters with potential mislabeling:", paste(problematic_clusters, collapse = ", "), "\n")
  cat("\nACTION REQUIRED: Review these clusters manually!\n")
} else {
  cat("\nAll clusters passed basic validation checks.\n")
}

# Save detailed validation report
sink("results/metadata_validation/VALIDATION_REPORT.txt")
cat("METADATA VALIDATION REPORT\n")
cat("=========================\n\n")
cat("Date:", Sys.Date(), "\n")
cat("Object:", "iSCORE-PD_plus_CRISPRi_validated.rds\n\n")

cat("Clusters checked:", length(unique_clusters), "\n")
cat("Total cells:", ncol(seurat_obj), "\n\n")

if (length(problematic_clusters) > 0) {
  cat("PROBLEMATIC CLUSTERS:\n")
  for (cluster in problematic_clusters) {
    cat("\nCluster", cluster, "\n")
    cluster_label <- unique(seurat_obj$validated_celltype[seurat_obj$seurat_clusters_fine == cluster])[1]
    cat("  Current label:", cluster_label, "\n")
    cluster_expr <- avg_expression_by_cluster[, as.character(cluster)]
    top_genes <- names(sort(cluster_expr, decreasing = TRUE)[1:10])
    cat("  Top 10 markers:", paste(top_genes, collapse = ", "), "\n")
  }
}
sink()

cat("\nValidation complete! Check results in: results/metadata_validation/\n")