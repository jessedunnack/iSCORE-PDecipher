#!/usr/bin/env Rscript

# Add Cell Type Annotations to Seurat Object
# Creates new metadata columns with fine cluster cell type assignments

library(Seurat)
library(dplyr)

# Function to add cell type annotations
add_celltype_annotations <- function(seurat_obj, annotations_file = NULL) {
  
  # Load annotations if file provided, otherwise use default
  if (is.null(annotations_file)) {
    annotations_file <- "results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv"
  }
  
  if (!file.exists(annotations_file)) {
    stop("Annotations file not found. Please run comprehensive_cell_type_annotation.R first")
  }
  
  cat("Loading cell type annotations from:", annotations_file, "\n")
  annotations <- read.csv(annotations_file, stringsAsFactors = FALSE)
  
  # Check what clustering column exists
  metadata_cols <- colnames(seurat_obj@meta.data)
  cat("\nAvailable metadata columns:\n")
  print(metadata_cols)
  
  # Find fine cluster column (you may need to adjust this)
  cluster_cols <- grep("cluster|res\\.", metadata_cols, value = TRUE)
  cat("\nPotential cluster columns found:\n")
  print(cluster_cols)
  
  # Look for the fine cluster column by name first
  if ("seurat_clusters_fine" %in% metadata_cols) {
    fine_cluster_col <- "seurat_clusters_fine"
    cat("\nFound 'seurat_clusters_fine' column\n")
  } else {
    # If not found, determine which column has 36 clusters (fine resolution)
    n_clusters <- sapply(cluster_cols, function(col) {
      length(unique(seurat_obj@meta.data[[col]]))
    })
    
    # Find columns with 36 clusters
    fine_cluster_candidates <- names(n_clusters)[which(n_clusters == 36)]
    
    if (length(fine_cluster_candidates) == 0) {
      cat("\nWarning: Could not identify fine cluster column with 36 clusters\n")
      cat("Number of clusters per column:\n")
      print(n_clusters)
      stop("Could not identify fine cluster column")
    } else if (length(fine_cluster_candidates) > 1) {
      cat("\nMultiple columns with 36 clusters found:\n")
      print(fine_cluster_candidates)
      
      # Try to pick the most likely one
      if ("seurat_clusters" %in% fine_cluster_candidates) {
        fine_cluster_col <- "seurat_clusters"
        cat("Using 'seurat_clusters' as fine cluster column\n")
      } else if (any(grepl("SCT_snn_res\\.0\\.8", fine_cluster_candidates))) {
        fine_cluster_col <- grep("SCT_snn_res\\.0\\.8", fine_cluster_candidates, value = TRUE)[1]
        cat("Using", fine_cluster_col, "as fine cluster column\n")
      } else {
        fine_cluster_col <- fine_cluster_candidates[1]
        cat("Using first candidate:", fine_cluster_col, "\n")
      }
    } else {
      fine_cluster_col <- fine_cluster_candidates[1]
    }
  }
  
  cat("\nUsing fine cluster column:", fine_cluster_col, "\n")
  fine_clusters <- seurat_obj@meta.data[[fine_cluster_col]]
  
  # Create lookup tables - ensure keys are character strings
  celltype_lookup <- setNames(annotations$cell_type, as.character(annotations$fine_cluster))
  subtype_lookup <- setNames(annotations$subtype, as.character(annotations$fine_cluster))
  confidence_lookup <- setNames(annotations$confidence, as.character(annotations$fine_cluster))
  markers_lookup <- setNames(annotations$key_markers, as.character(annotations$fine_cluster))
  
  # Add new metadata columns
  cat("\nAdding cell type annotations to metadata...\n")
  
  # Convert fine clusters to character for lookup
  fine_clusters_char <- as.character(fine_clusters)
  
  # Create new metadata vectors with cell names
  celltypes_fine <- celltype_lookup[fine_clusters_char]
  celltypes_fine_subtype <- subtype_lookup[fine_clusters_char]
  celltypes_fine_confidence <- confidence_lookup[fine_clusters_char]
  celltypes_fine_markers <- markers_lookup[fine_clusters_char]
  
  # Set names to match cell barcodes
  names(celltypes_fine) <- rownames(seurat_obj@meta.data)
  names(celltypes_fine_subtype) <- rownames(seurat_obj@meta.data)
  names(celltypes_fine_confidence) <- rownames(seurat_obj@meta.data)
  names(celltypes_fine_markers) <- rownames(seurat_obj@meta.data)
  
  # Add to Seurat object
  seurat_obj[["celltypes_fine"]] <- celltypes_fine
  seurat_obj[["celltypes_fine_subtype"]] <- celltypes_fine_subtype
  seurat_obj[["celltypes_fine_confidence"]] <- celltypes_fine_confidence
  seurat_obj[["celltypes_fine_markers"]] <- celltypes_fine_markers
  
  # Handle any NA values (clusters not in annotation)
  na_count <- sum(is.na(seurat_obj$celltypes_fine))
  if (na_count > 0) {
    cat("\nWarning:", na_count, "cells have no annotation (setting to 'Unknown')\n")
    seurat_obj$celltypes_fine[is.na(seurat_obj$celltypes_fine)] <- "Unknown"
    seurat_obj$celltypes_fine_subtype[is.na(seurat_obj$celltypes_fine_subtype)] <- ""
    seurat_obj$celltypes_fine_confidence[is.na(seurat_obj$celltypes_fine_confidence)] <- "Low"
    seurat_obj$celltypes_fine_markers[is.na(seurat_obj$celltypes_fine_markers)] <- ""
  }
  
  # Create simplified version (major cell types only)
  seurat_obj$celltypes_fine_major <- case_when(
    grepl("Dopaminergic", seurat_obj$celltypes_fine) ~ "Dopaminergic Neurons",
    grepl("Floor Plate|Floorplate", seurat_obj$celltypes_fine) ~ "Floor Plate Progenitors",
    grepl("Neural Progenitor|Neuroblast", seurat_obj$celltypes_fine) ~ "Neural Progenitors",
    grepl("Choroid Plexus", seurat_obj$celltypes_fine) ~ "Choroid Plexus",
    grepl("Leptomeningeal", seurat_obj$celltypes_fine) ~ "Leptomeningeal",
    grepl("Vascular|Pericyte", seurat_obj$celltypes_fine) ~ "Vascular Cells",
    grepl("Mesenchymal|Fibroblast", seurat_obj$celltypes_fine) ~ "Mesenchymal/Fibroblast",
    grepl("Oligodendrocyte", seurat_obj$celltypes_fine) ~ "Oligodendrocytes",
    grepl("Stressed|Dying", seurat_obj$celltypes_fine) ~ "Stressed/Dying Cells",
    grepl("Proliferating", seurat_obj$celltypes_fine) ~ "Proliferating Cells",
    grepl("Hypothalamic", seurat_obj$celltypes_fine) ~ "Hypothalamic Neurons",
    grepl("Caudal|Spinal", seurat_obj$celltypes_fine) ~ "Caudal Neural",
    TRUE ~ "Other/Unknown"
  )
  
  # Print summary
  cat("\n=== Cell Type Annotation Summary ===\n")
  print(table(seurat_obj$celltypes_fine_major))
  
  cat("\n=== Confidence Level Summary ===\n")
  print(table(seurat_obj$celltypes_fine_confidence))
  
  return(seurat_obj)
}

# Function to visualize annotations
visualize_celltype_annotations <- function(seurat_obj) {
  library(ggplot2)
  
  # Create UMAP with cell types
  p1 <- DimPlot(seurat_obj, 
                group.by = "celltypes_fine_major", 
                label = TRUE, 
                repel = TRUE,
                pt.size = 0.5) +
    ggtitle("Cell Types (Major Categories)") +
    theme(legend.position = "bottom")
  
  # Create UMAP with confidence levels
  p2 <- DimPlot(seurat_obj, 
                group.by = "celltypes_fine_confidence",
                cols = c("High" = "darkgreen", "Medium" = "orange", "Low" = "red"),
                pt.size = 0.5) +
    ggtitle("Annotation Confidence")
  
  # Save plots
  dir.create("results/cell_type_annotations/plots", recursive = TRUE, showWarnings = FALSE)
  
  ggsave("results/cell_type_annotations/plots/umap_celltypes_major.pdf", p1, width = 12, height = 10)
  ggsave("results/cell_type_annotations/plots/umap_confidence.pdf", p2, width = 10, height = 8)
  
  return(list(p1, p2))
}

# Main function
main <- function(seurat_path = NULL) {
  cat("Adding Cell Type Annotations to Seurat Object\n")
  cat("============================================\n\n")
  
  # Default path if not provided
  if (is.null(seurat_path)) {
    seurat_path <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"
  }
  
  # Check if file exists
  if (!file.exists(seurat_path)) {
    stop("Seurat object not found at: ", seurat_path)
  }
  
  # Load Seurat object
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_path)
  cat("Loaded object with", ncol(seurat_obj), "cells\n")
  
  # Add annotations
  seurat_obj <- add_celltype_annotations(seurat_obj)
  
  # Save annotated object
  output_path <- gsub("\\.rds$", "_annotated.rds", seurat_path)
  cat("\nSaving annotated Seurat object to:", output_path, "\n")
  saveRDS(seurat_obj, output_path)
  
  # Create visualizations
  cat("\nCreating visualization plots...\n")
  plots <- visualize_celltype_annotations(seurat_obj)
  
  cat("\n=== Success! ===\n")
  cat("New metadata columns added:\n")
  cat("- celltypes_fine: Main cell type assignment\n")
  cat("- celltypes_fine_subtype: Detailed subtype\n")
  cat("- celltypes_fine_confidence: Annotation confidence (High/Medium/Low)\n")
  cat("- celltypes_fine_markers: Key marker genes\n")
  cat("- celltypes_fine_major: Simplified major categories\n")
  
  return(seurat_obj)
}

# Example usage function
example_usage <- function() {
  cat("\n=== Example Usage ===\n\n")
  cat("# Basic usage with default path:\n")
  cat("source('scripts/add_celltype_annotations_to_seurat.R')\n")
  cat("seurat_annotated <- main()\n\n")
  
  cat("# With custom path:\n")
  cat("seurat_annotated <- main('path/to/your/seurat.rds')\n\n")
  
  cat("# Access the new metadata:\n")
  cat("table(seurat_annotated$celltypes_fine_major)\n")
  cat("DimPlot(seurat_annotated, group.by = 'celltypes_fine_major')\n\n")
  
  cat("# Filter for specific cell types:\n")
  cat("da_neurons <- subset(seurat_annotated, celltypes_fine_major == 'Dopaminergic Neurons')\n")
}

# Run example if executed with --example flag
args <- commandArgs(trailingOnly = TRUE)
if ("--example" %in% args) {
  example_usage()
} else if (!interactive()) {
  # Run main function if executed as script
  seurat_obj <- main()
}