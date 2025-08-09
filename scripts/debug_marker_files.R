#!/usr/bin/env Rscript

# Debug marker file structure

cat("Checking marker file structure...\n\n")

# Check coarse markers
cat("1. Loading coarse marker file...\n")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")

cat("Structure:\n")
str(coarse_markers, max.level = 2)

cat("\nNames of elements:\n")
print(names(coarse_markers))

# Check if it's a data frame or list
if (is.data.frame(coarse_markers)) {
  cat("\nIt's a data frame. Columns:\n")
  print(colnames(coarse_markers))
  
  cat("\nUnique clusters:\n")
  if ("cluster" %in% colnames(coarse_markers)) {
    print(unique(coarse_markers$cluster))
  }
  
  cat("\nFirst few rows:\n")
  print(head(coarse_markers))
} else if (is.list(coarse_markers)) {
  cat("\nIt's a list with", length(coarse_markers), "elements\n")
  
  # Check first element
  cat("\nFirst element:\n")
  if (length(coarse_markers) > 0) {
    print(names(coarse_markers)[1])
    str(coarse_markers[[1]], max.level = 1)
  }
}