#!/usr/bin/env Rscript

# Quick diagnostic to check for common metadata issues

library(Seurat)
library(dplyr)

cat("=================================\n")
cat("QUICK METADATA DIAGNOSTIC CHECK\n")
cat("=================================\n\n")

# Load object
x <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_validated.rds")

# Check available metadata columns
cat("Available metadata columns:\n")
cat(paste(colnames(x@meta.data), collapse=", "), "\n\n")

# 1. Check cluster 5 specifically (should be Dopaminergic)
cat("DIAGNOSTIC 1: Checking cluster 5 (expected to be Dopaminergic)\n")
cluster5_cells <- WhichCells(x, expression = seurat_clusters_fine == 5)
cat("Number of cells in cluster 5:", length(cluster5_cells), "\n")

# Check what cell types are assigned to cluster 5
# First check which celltype columns exist
celltype_cols <- grep("celltype|validated", colnames(x@meta.data), value = TRUE)
cat("Available celltype columns:", paste(celltype_cols, collapse=", "), "\n")

if ("validated_celltype_major" %in% colnames(x@meta.data)) {
  cluster5_types <- table(x$validated_celltype_major[x$seurat_clusters_fine == 5])
  cat("\nCell types in cluster 5:\n")
  print(cluster5_types)
} else {
  cat("\nWARNING: validated_celltype_major column not found!\n")
  cat("Checking other celltype columns...\n")
  for (col in celltype_cols) {
    cat("\n", col, "in cluster 5:\n")
    print(table(x@meta.data[x$seurat_clusters_fine == 5, col]))
  }
}

# Check TH expression in cluster 5
if ("TH" %in% rownames(x)) {
  # Check which assay to use
  if ("RNA" %in% names(x@assays)) {
    th_expr <- AverageExpression(x, features = "TH", group.by = "seurat_clusters_fine", assays = "RNA")$RNA
  } else if ("SCT" %in% names(x@assays)) {
    th_expr <- AverageExpression(x, features = "TH", group.by = "seurat_clusters_fine", assays = "SCT")$SCT
  } else {
    cat("\nNo RNA or SCT assay found\n")
    th_expr <- NULL
  }
  
  if (!is.null(th_expr) && "5" %in% colnames(th_expr)) {
    cat("\nTH expression in cluster 5:", th_expr["TH", "5"], "\n")
    cat("TH expression across all clusters:\n")
    th_vals <- th_expr["TH",]
    print(head(sort(th_vals, decreasing = TRUE), 10))
  } else {
    cat("\nCould not calculate TH expression\n")
  }
}

# 2. Check cluster 4 (should be Choroid Plexus)
cat("\n\nDIAGNOSTIC 2: Checking cluster 4 (expected to be Choroid Plexus)\n")
if ("validated_celltype_major" %in% colnames(x@meta.data)) {
  cluster4_types <- table(x$validated_celltype_major[x$seurat_clusters_fine == 4])
  cat("Cell types in cluster 4:\n")
  print(cluster4_types)
} else {
  cat("validated_celltype_major column not found\n")
}

if ("TTR" %in% rownames(x)) {
  ttr_expr <- AverageExpression(x, features = "TTR", group.by = "seurat_clusters_fine", assays = "RNA")$RNA
  if ("4" %in% colnames(ttr_expr)) {
    cat("\nTTR expression in cluster 4:", ttr_expr["TTR", "4"], "\n")
  }
}

# 3. Check for any obvious misalignments
cat("\n\nDIAGNOSTIC 3: Checking for obvious misalignments\n")

# Sample some cells and show their metadata
set.seed(123)
sample_cells <- sample(colnames(x), 20)

# Get available columns
available_cols <- c("seurat_clusters_fine")
if ("celltypes_fine" %in% colnames(x@meta.data)) available_cols <- c(available_cols, "celltypes_fine")
if ("validated_celltype" %in% colnames(x@meta.data)) available_cols <- c(available_cols, "validated_celltype")
if ("validated_celltype_major" %in% colnames(x@meta.data)) available_cols <- c(available_cols, "validated_celltype_major")

sample_data <- x@meta.data[sample_cells, available_cols]

cat("\nSample of 20 cells:\n")
print(sample_data)

# 4. Check if any clusters have mixed cell types
cat("\n\nDIAGNOSTIC 4: Clusters with mixed cell type assignments\n")

if ("validated_celltype_major" %in% colnames(x@meta.data)) {
  mixed_clusters <- x@meta.data %>%
    group_by(seurat_clusters_fine) %>%
    summarise(
      n_celltypes = n_distinct(validated_celltype_major),
      celltypes = paste(unique(validated_celltype_major), collapse = " | ")
    ) %>%
    filter(n_celltypes > 1)
} else {
  cat("Cannot check mixed clusters - validated_celltype_major not found\n")
  mixed_clusters <- data.frame()
}

if (nrow(mixed_clusters) > 0) {
  cat("WARNING: Some clusters have multiple cell types assigned!\n")
  print(mixed_clusters)
} else {
  cat("Good: Each cluster has a single cell type assigned.\n")
}

# 5. Quick marker check
cat("\n\nDIAGNOSTIC 5: Quick marker expression check\n")
key_markers <- c("TH", "TTR", "TUBB3", "SOX2", "FOXA2")
available <- key_markers[key_markers %in% rownames(x)]

for (marker in available) {
  expr <- AverageExpression(x, features = marker, group.by = "seurat_clusters_fine", assays = "RNA")$RNA
  top_cluster <- names(which.max(expr[marker,]))
  
  if ("validated_celltype_major" %in% colnames(x@meta.data)) {
    celltype <- unique(x$validated_celltype_major[x$seurat_clusters_fine == top_cluster])[1]
    cat("\n", marker, "highest in cluster", top_cluster, "(", celltype, ")")
    
    # Flag potential issues
    if (marker == "TH" && !grepl("Dopaminergic", celltype)) {
      cat(" <- WARNING: TH high but not labeled as Dopaminergic!")
    }
    if (marker == "TTR" && !grepl("Choroid", celltype)) {
      cat(" <- WARNING: TTR high but not labeled as Choroid Plexus!")
    }
  } else {
    cat("\n", marker, "highest in cluster", top_cluster)
  }
}