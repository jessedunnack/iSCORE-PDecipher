#!/usr/bin/env Rscript

# RUN IN RSTUDIO - FIXED VERSION
# This script re-runs clustering and OVERWRITES existing columns

cat("=================================================================\n")
cat("RECLUSTERING SCRIPT FOR RSTUDIO - OVERWRITES EXISTING COLUMNS\n")
cat("=================================================================\n\n")

cat("WARNING: This will OVERWRITE:\n")
cat("- seurat_clusters_coarse (with resolution 0.2)\n")
cat("- seurat_clusters_fine (with default resolution)\n\n")

cat("Press Enter to continue or Ctrl+C to cancel...\n")
# readline()  # Uncomment this if you want a pause

# Source the reclustering script
source("scripts/recluster_and_investigate_fixed.R")

cat("\n\nReclustering complete!\n")
cat("The existing cluster columns have been overwritten.\n")