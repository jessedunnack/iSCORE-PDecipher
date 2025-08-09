#!/usr/bin/env Rscript

# ESTABLISH CLUSTER HIERARCHY BY RUNNING FindClusters AT TWO RESOLUTIONS
# This will create proper fine-to-coarse cluster mapping

library(Seurat)
library(dplyr)

cat("=================================================================\n")
cat("ESTABLISHING CLUSTER HIERARCHY\n")
cat("Running FindClusters at two resolutions\n")
cat("=================================================================\n\n")

# Load the Seurat object
cat("Loading Seurat object (this may take a few minutes)...\n")
start_time <- Sys.time()
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")
load_time <- difftime(Sys.time(), start_time, units = "mins")
cat("Object loaded in", round(load_time, 2), "minutes\n")
cat("Total cells:", ncol(seurat_obj), "\n\n")

# Check if we already have cluster columns
existing_clusters <- grep("cluster", colnames(seurat_obj@meta.data), value = TRUE)
cat("Existing cluster columns:", paste(existing_clusters, collapse = ", "), "\n\n")

# Run FindNeighbors if not already done
if (!"SCT_snn" %in% names(seurat_obj@graphs)) {
  cat("Running FindNeighbors...\n")
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:30)
} else {
  cat("Using existing neighbor graph\n")
}

# 1. Find COARSE clusters (lower resolution)
cat("\n1. Finding COARSE clusters (resolution = 0.2)...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, verbose = FALSE)
seurat_obj$seurat_clusters_coarse <- Idents(seurat_obj)
n_coarse <- length(unique(seurat_obj$seurat_clusters_coarse))
cat("Found", n_coarse, "coarse clusters\n")

# 2. Find FINE clusters (higher resolution)
cat("\n2. Finding FINE clusters (resolution = 0.8)...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.8, verbose = FALSE)
seurat_obj$seurat_clusters_fine <- Idents(seurat_obj)
n_fine <- length(unique(seurat_obj$seurat_clusters_fine))
cat("Found", n_fine, "fine clusters\n")

# 3. Create fine-to-coarse mapping
cat("\n3. Creating fine-to-coarse cluster mapping...\n")
cluster_mapping <- seurat_obj@meta.data %>%
  group_by(seurat_clusters_fine, seurat_clusters_coarse) %>%
  summarise(n_cells = n(), .groups = 'drop') %>%
  arrange(seurat_clusters_fine, desc(n_cells))

# Get the dominant coarse cluster for each fine cluster
fine_to_coarse_map <- cluster_mapping %>%
  group_by(seurat_clusters_fine) %>%
  slice_max(n_cells, n = 1) %>%
  select(fine_cluster = seurat_clusters_fine, 
         coarse_cluster = seurat_clusters_coarse,
         n_cells)

cat("\nFine-to-Coarse Cluster Mapping:\n")
print(as.data.frame(fine_to_coarse_map))

# 4. Save the mapping
cat("\n4. Saving cluster hierarchy information...\n")
dir.create("results/cluster_hierarchy", recursive = TRUE, showWarnings = FALSE)

# Save mapping table
write.csv(fine_to_coarse_map, 
          "results/cluster_hierarchy/fine_to_coarse_mapping.csv", 
          row.names = FALSE)

# Save full mapping details
write.csv(cluster_mapping, 
          "results/cluster_hierarchy/detailed_cluster_mapping.csv", 
          row.names = FALSE)

# 5. Create summary by coarse cluster
cat("\n5. Creating summary by coarse cluster...\n")
coarse_summary <- seurat_obj@meta.data %>%
  group_by(seurat_clusters_coarse) %>%
  summarise(
    n_cells = n(),
    fine_clusters = paste(sort(unique(seurat_clusters_fine)), collapse = ", "),
    n_fine_clusters = n_distinct(seurat_clusters_fine),
    .groups = 'drop'
  )

cat("\nCoarse Cluster Summary:\n")
print(coarse_summary)
write.csv(coarse_summary, 
          "results/cluster_hierarchy/coarse_cluster_summary.csv", 
          row.names = FALSE)

# 6. Save metadata with cluster assignments
cat("\n6. Saving metadata with cluster hierarchy...\n")
metadata_subset <- seurat_obj@meta.data %>%
  select(seurat_clusters_coarse, seurat_clusters_fine, 
         all_of(intersect(c("validated_celltype", "corrected_celltype", 
                           "mutation_tidy", "scMAGeCK_gene_assignment"), 
                         colnames(seurat_obj@meta.data))))

saveRDS(metadata_subset, "results/cluster_hierarchy/cluster_metadata.rds")

# Also save as CSV for easy viewing (first 1000 rows as example)
write.csv(head(metadata_subset, 1000), 
          "results/cluster_hierarchy/cluster_metadata_sample.csv", 
          row.names = TRUE)

# 7. Save the updated Seurat object
cat("\n7. Saving Seurat object with cluster hierarchy...\n")
save_start <- Sys.time()
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_with_hierarchy.rds")
save_time <- difftime(Sys.time(), save_start, units = "mins")
cat("Object saved in", round(save_time, 2), "minutes\n")

# 8. Create visualization plots
cat("\n8. Creating visualization plots...\n")
library(ggplot2)
library(patchwork)

# Plot both clustering resolutions
if ("umap" %in% names(seurat_obj@reductions) || "umap.cca" %in% names(seurat_obj@reductions)) {
  reduction_to_use <- ifelse("umap.cca" %in% names(seurat_obj@reductions), "umap.cca", "umap")
  
  # Sample cells for faster plotting
  set.seed(42)
  cells_to_plot <- sample(colnames(seurat_obj), min(50000, ncol(seurat_obj)))
  
  p1 <- DimPlot(seurat_obj[, cells_to_plot], 
                reduction = reduction_to_use,
                group.by = "seurat_clusters_coarse", 
                label = TRUE, 
                repel = TRUE) +
    ggtitle(paste("Coarse Clusters (n =", n_coarse, ")"))
  
  p2 <- DimPlot(seurat_obj[, cells_to_plot], 
                reduction = reduction_to_use,
                group.by = "seurat_clusters_fine", 
                label = TRUE, 
                repel = TRUE) +
    ggtitle(paste("Fine Clusters (n =", n_fine, ")"))
  
  p_combined <- p1 + p2
  ggsave("results/cluster_hierarchy/cluster_hierarchy_umap.pdf", 
         p_combined, width = 16, height = 8)
  
  cat("UMAP plots saved\n")
}

# Final summary
cat("\n\n=== CLUSTER HIERARCHY ESTABLISHED ===\n")
cat("=====================================\n")
cat("Coarse clusters (res=0.2):", n_coarse, "\n")
cat("Fine clusters (res=0.8):", n_fine, "\n")
cat("\nFiles saved:\n")
cat("- results/cluster_hierarchy/fine_to_coarse_mapping.csv\n")
cat("- results/cluster_hierarchy/coarse_cluster_summary.csv\n")
cat("- results/cluster_hierarchy/cluster_metadata.rds\n")
cat("- results/cluster_hierarchy/cluster_hierarchy_umap.pdf\n")
cat("- ../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_with_hierarchy.rds\n")

cat("\nTotal runtime:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n")
cat("\nDONE! You can now proceed with hierarchical cluster analysis.\n")