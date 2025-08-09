#!/usr/bin/env Rscript

# RUN_COARSE_INVESTIGATION_RSTUDIO.R
# Simplified version for running in RStudio

cat("=================================================================\n")
cat("COARSE CLUSTER INVESTIGATION FOR RSTUDIO\n")
cat("Run this script in RStudio by sourcing it\n")
cat("=================================================================\n\n")

# Source the main script
cat("Sourcing comprehensive coarse cluster investigation script...\n")
cat("This may take several minutes due to the large Seurat object.\n\n")

# Set working directory if needed
# setwd("E:/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher")

# Run the comprehensive analysis
source("scripts/comprehensive_coarse_cluster_investigation.R")

cat("\n\nAnalysis complete! Check the results in:\n")
cat("- results/comprehensive_coarse_cluster_results.rds\n")
cat("- results/coarse_cluster_identity_summary.csv\n")
cat("- results/comprehensive_coarse_investigation_log.txt (if you saved console output)\n")