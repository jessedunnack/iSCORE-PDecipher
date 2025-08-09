#!/usr/bin/env Rscript

# Quick check of marker file structure

cat("Checking marker file structure...\n\n")

# Check cluster markers
markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")

cat("Structure of cluster_markers.rds:\n")
cat("Class:", class(markers), "\n")
cat("Names:", names(markers)[1:5], "...\n")
cat("Length:", length(markers), "\n")

# Check first element
cat("\nFirst element structure:\n")
str(markers[[1]], max.level = 2)

# Check if it has coarse and fine
if ("coarse" %in% names(markers)) {
  cat("\nHas 'coarse' markers\n")
}
if ("fine" %in% names(markers)) {
  cat("\nHas 'fine' markers\n")
}

# Check all markers expression if exists
if (file.exists("inst/extdata/umap_data/iSCORE_PD_CRISPRi_all_markers_expr.rds")) {
  all_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_all_markers_expr.rds")
  cat("\n\nStructure of all_markers_expr.rds:\n")
  str(all_markers, max.level = 1)
}