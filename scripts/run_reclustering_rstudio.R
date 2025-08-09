#!/usr/bin/env Rscript

# RUN IN RSTUDIO
# This script re-runs clustering with proper resolutions

cat("=================================================================\n")
cat("RECLUSTERING SCRIPT FOR RSTUDIO\n")
cat("=================================================================\n\n")

cat("This script will:\n")
cat("1. Re-run FindClusters with resolution 0.2 for coarse clusters\n")
cat("2. Re-run FindClusters with default resolution for fine clusters\n")
cat("3. Perform quick investigation of the new clusters\n")
cat("4. Save the reclustered object\n\n")

cat("Expected runtime: 5-10 minutes\n\n")

# Source the reclustering script
source("scripts/recluster_and_investigate.R")

cat("\n\nReclustering complete!\n")
cat("Check results in:\n")
cat("- results/seurat_obj_reclustered.rds (updated Seurat object)\n")
cat("- Console output for cluster summaries\n")