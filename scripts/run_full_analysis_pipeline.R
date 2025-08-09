#!/usr/bin/env Rscript

# RUN FULL ANALYSIS PIPELINE
# Complete workflow: recluster -> analyze expression -> map hierarchies

cat("=================================================================\n")
cat("FULL CLUSTER ANALYSIS PIPELINE\n")
cat("=================================================================\n\n")

cat("This pipeline will:\n")
cat("1. Re-cluster with correct resolutions (if not already done)\n")
cat("2. Analyze expression in all coarse clusters\n")
cat("3. Create coarse-to-fine mapping\n")
cat("4. Assign cluster identities\n\n")

# Step 1: Check if reclustering is needed
if (file.exists("results/seurat_obj_reclustered.rds")) {
  cat("Found reclustered object. Skipping reclustering step.\n\n")
} else {
  cat("Running reclustering...\n")
  source("scripts/recluster_and_investigate_fixed.R")
  cat("\nReclustering complete.\n\n")
}

# Step 2: Run comprehensive expression analysis
cat("\nRunning comprehensive expression analysis...\n")
cat("This will check all marker sets across all clusters.\n\n")
source("scripts/analyze_reclustered_expression.R")

cat("\n\n=== PIPELINE COMPLETE ===\n")
cat("Check results in:\n")
cat("- results/reclustered_analysis/coarse_cluster_identities.csv\n")
cat("- results/reclustered_analysis/fine_to_coarse_mapping.csv\n")
cat("- Console output for detailed marker expression\n")