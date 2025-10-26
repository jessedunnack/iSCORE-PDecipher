#!/usr/bin/env Rscript
# PHASE 1: Inspect pooled MixScale cluster RDS structure

cat("\n=== INSPECTING POOLED MIXSCALE CLUSTER RDS ===\n\n")

# Load one pooled cluster file
pooled_file <- "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0/all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds"

cat("Loading file:", pooled_file, "\n")
cat("File size:", file.size(pooled_file) / 1024^2, "MB\n\n")

cluster_data <- readRDS(pooled_file)

# Check structure
cat("=== STRUCTURE ===\n")
cat("Class:", class(cluster_data), "\n")
cat("Number of perturbations:", length(cluster_data), "\n")
cat("First 5 perturbations:", paste(head(names(cluster_data), 5), collapse = ", "), "\n\n")

# Check first perturbation
cat("=== FIRST PERTURBATION STRUCTURE ===\n")
first_pert <- cluster_data[[1]]
cat("Perturbation name:", names(cluster_data)[1], "\n")
cat("Class:", class(first_pert), "\n")
cat("Columns:", paste(colnames(first_pert), collapse = ", "), "\n")
cat("Rows:", nrow(first_pert), "\n\n")

cat("=== FIRST FEW ROWS ===\n")
print(head(first_pert, 5))

cat("\n=== P-VALUE COLUMNS CHECK ===\n")
cat("Has p_weight:", "p_weight" %in% colnames(first_pert), "\n")
cat("Has p_weight_BH:", "p_weight_BH" %in% colnames(first_pert), "\n")
cat("Has p_weight_bonferroni:", "p_weight_bonferroni" %in% colnames(first_pert), "\n")

cat("\n=== COMPARISON WITH REQUIRED FORMAT ===\n")
cat("CURRENT columns: gene_ID, log2FC, p_weight, p_weight_BH, p_weight_bonferroni, ...\n")
cat("REQUIRED columns: (gene names as rownames), avg_log2FC, p_val_adj, pct.1, pct.2\n")
cat("\nCONVERSION NEEDED:\n")
cat("1. Move gene_ID to rownames\n")
cat("2. Rename log2FC -> avg_log2FC\n")
cat("3. Rename p_weight_BH -> p_val_adj (or chosen correction)\n")
cat("4. Add pct.1 and pct.2 columns (can be NA)\n")

cat("\n=== INSPECTION COMPLETE ===\n")
