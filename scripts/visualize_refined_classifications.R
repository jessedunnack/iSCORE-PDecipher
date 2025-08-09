#!/usr/bin/env Rscript

# Visualize refined cell type classifications with vulnerability status

library(Seurat)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

# Function to create enhanced UMAP plots
create_refined_umaps <- function(seurat_obj, refined_annotations) {
  
  # Add refined annotations to Seurat object
  refined_lookup <- setNames(refined_annotations$refined_type, 
                            as.character(refined_annotations$fine_cluster))
  vulnerability_lookup <- setNames(refined_annotations$vulnerability_status,
                                  as.character(refined_annotations$fine_cluster))
  da_subtype_lookup <- setNames(refined_annotations$da_subtype,
                               as.character(refined_annotations$fine_cluster))
  
  seurat_obj$refined_celltype <- refined_lookup[as.character(seurat_obj$seurat_clusters_fine)]
  seurat_obj$vulnerability <- vulnerability_lookup[as.character(seurat_obj$seurat_clusters_fine)]
  seurat_obj$da_subtype <- da_subtype_lookup[as.character(seurat_obj$seurat_clusters_fine)]
  
  # Create plots
  plots <- list()
  
  # 1. Refined cell types
  plots$refined_types <- DimPlot(seurat_obj, 
                                group.by = "refined_celltype",
                                label = TRUE,
                                repel = TRUE,
                                pt.size = 0.3) +
    ggtitle("Refined Cell Type Classifications") +
    theme(legend.position = "right")
  
  # 2. Vulnerability status
  plots$vulnerability <- DimPlot(seurat_obj,
                               group.by = "vulnerability",
                               cols = c("Vulnerable" = "red", 
                                       "Resilient" = "blue",
                                       "Unknown" = "gray"),
                               pt.size = 0.3) +
    ggtitle("Vulnerability Status (MPTP Model)")
  
  # 3. DA subtypes only
  da_cells <- WhichCells(seurat_obj, expression = da_subtype != "")
  if (length(da_cells) > 0) {
    plots$da_subtypes <- DimPlot(seurat_obj,
                                cells = da_cells,
                                group.by = "da_subtype",
                                pt.size = 1) +
      ggtitle("Dopaminergic Neuron Subtypes") +
      theme(legend.position = "right")
  }
  
  return(plots)
}

# Function to create marker heatmap
create_marker_heatmap <- function(seurat_obj, key_markers) {
  
  # Average expression by cluster
  avg_exp <- AverageExpression(seurat_obj, 
                              features = key_markers,
                              group.by = "seurat_clusters_fine",
                              assays = "RNA")$RNA
  
  # Scale for visualization
  avg_exp_scaled <- t(scale(t(avg_exp)))
  
  # Get cluster annotations
  cluster_order <- as.character(0:35)
  
  # Create annotation for columns (clusters)
  refined <- read.csv("results/cell_type_annotations/refined_classifications.csv")
  
  col_ann <- HeatmapAnnotation(
    CellType = refined$refined_type[match(cluster_order, refined$fine_cluster)],
    Vulnerability = refined$vulnerability_status[match(cluster_order, refined$fine_cluster)],
    Confidence = refined$refined_confidence[match(cluster_order, refined$fine_cluster)],
    col = list(
      Vulnerability = c("Vulnerable" = "red", "Resilient" = "blue", "Unknown" = "gray"),
      Confidence = c("High" = "darkgreen", "Medium" = "orange", "Low" = "red")
    )
  )
  
  # Create row annotation for genes
  gene_categories <- list(
    DA_markers = c("TH", "DDC", "SLC6A3", "SLC18A2", "NR4A2"),
    Vulnerability = c("SOX6", "ALDH1A1", "KCNJ6", "MEIS1", "MEIS2"),
    Resilience = c("CALB1", "FOXP2", "SORCS3", "RELN", "TMEFF2"),
    Regional = c("OTX2", "GBX2", "EN1", "EN2", "HOXD10", "HOXD11"),
    FloorPlate = c("FOXA2", "LMX1A", "LMX1B", "SHH", "CORIN")
  )
  
  gene_cat_vector <- character(length(key_markers))
  for (cat in names(gene_categories)) {
    gene_cat_vector[key_markers %in% gene_categories[[cat]]] <- cat
  }
  
  row_ann <- HeatmapAnnotation(
    Category = gene_cat_vector,
    which = "row",
    col = list(
      Category = c(
        "DA_markers" = "purple",
        "Vulnerability" = "red",
        "Resilience" = "blue",
        "Regional" = "green",
        "FloorPlate" = "orange"
      )
    )
  )
  
  # Create heatmap
  ht <- Heatmap(
    avg_exp_scaled,
    name = "Scaled\nExpression",
    column_title = "Refined Cluster Classifications",
    row_title = "Key Markers",
    top_annotation = col_ann,
    left_annotation = row_ann,
    column_order = cluster_order,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
  )
  
  return(ht)
}

# Function to create summary plots
create_summary_plots <- function(refined_annotations) {
  
  plots <- list()
  
  # 1. Confidence comparison
  conf_data <- refined_annotations %>%
    select(fine_cluster, original_confidence, refined_confidence) %>%
    pivot_longer(cols = c(original_confidence, refined_confidence),
                names_to = "type",
                values_to = "confidence") %>%
    mutate(type = gsub("_confidence", "", type))
  
  plots$confidence <- ggplot(conf_data, aes(x = type, fill = confidence)) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = c("High" = "darkgreen", 
                               "Medium" = "orange", 
                               "Low" = "red")) +
    labs(title = "Confidence Level Comparison",
         x = "Classification Type",
         y = "Proportion") +
    theme_minimal()
  
  # 2. Vulnerability by cell type
  vuln_summary <- refined_annotations %>%
    group_by(refined_type, vulnerability_status) %>%
    summarise(count = n(), .groups = 'drop') %>%
    filter(!is.na(vulnerability_status))
  
  plots$vulnerability_by_type <- ggplot(vuln_summary, 
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
  
  return(plots)
}

# Main visualization function
main <- function() {
  cat("Creating Visualizations for Refined Classifications\n")
  cat("=================================================\n\n")
  
  # Load data
  seurat_path <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds"
  seurat_obj <- readRDS(seurat_path)
  refined_annotations <- read.csv("results/cell_type_annotations/refined_classifications.csv")
  
  # Create output directory
  dir.create("results/cell_type_annotations/refined_plots", recursive = TRUE, showWarnings = FALSE)
  
  # Create UMAP plots
  cat("Creating UMAP visualizations...\n")
  umap_plots <- create_refined_umaps(seurat_obj, refined_annotations)
  
  # Save UMAP plots
  for (plot_name in names(umap_plots)) {
    ggsave(paste0("results/cell_type_annotations/refined_plots/umap_", plot_name, ".pdf"),
           umap_plots[[plot_name]],
           width = 10, height = 8)
  }
  
  # Create marker heatmap
  cat("Creating marker expression heatmap...\n")
  key_markers <- c(
    # DA markers
    "TH", "DDC", "SLC6A3", "SLC18A2", "NR4A2", "FOXA2",
    # Vulnerability markers
    "SOX6", "ALDH1A1", "KCNJ6", "MEIS1", "MEIS2",
    # Resilience markers
    "CALB1", "FOXP2", "SORCS3", "RELN", "TMEFF2",
    # Regional markers
    "OTX2", "GBX2", "EN1", "EN2", "HOXD10", "HOXD11",
    # Floor plate
    "LMX1A", "LMX1B", "SHH", "CORIN"
  )
  
  # Filter to existing markers
  key_markers <- key_markers[key_markers %in% rownames(seurat_obj)]
  
  ht <- create_marker_heatmap(seurat_obj, key_markers)
  
  pdf("results/cell_type_annotations/refined_plots/marker_heatmap.pdf", 
      width = 12, height = 10)
  draw(ht)
  dev.off()
  
  # Create summary plots
  cat("Creating summary visualizations...\n")
  summary_plots <- create_summary_plots(refined_annotations)
  
  for (plot_name in names(summary_plots)) {
    ggsave(paste0("results/cell_type_annotations/refined_plots/summary_", plot_name, ".pdf"),
           summary_plots[[plot_name]],
           width = 8, height = 6)
  }
  
  # Create feature plots for key markers
  cat("Creating feature plots for vulnerability markers...\n")
  vulnerability_markers <- c("SOX6", "ALDH1A1", "KCNJ6", "MEIS1", "MEIS2", "CALB1")
  vulnerability_markers <- vulnerability_markers[vulnerability_markers %in% rownames(seurat_obj)]
  
  if (length(vulnerability_markers) > 0) {
    p <- FeaturePlot(seurat_obj, 
                    features = vulnerability_markers,
                    ncol = 3,
                    pt.size = 0.1)
    
    ggsave("results/cell_type_annotations/refined_plots/vulnerability_markers_features.pdf",
           p, width = 15, height = 10)
  }
  
  cat("\n\nVisualizations complete! Saved to:\n")
  cat("  results/cell_type_annotations/refined_plots/\n")
  
  return(list(
    umap_plots = umap_plots,
    summary_plots = summary_plots
  ))
}

# Run if not interactive
if (!interactive()) {
  plots <- main()
}