#!/usr/bin/env Rscript

#' Test Enhanced Direction Analysis with Real Data (v0.2.6)
#' 
#' This script tests the enhanced direction-aware analysis using actual
#' differential expression data from MAST and CRISPRi experiments

# Load required functions
source("R/enhanced_direction_analysis.R")
source("R/gene_harmonization.R")
source("R/signature_analysis.R")
source("R/experiment_weighting.R")

cat("=== REAL DATA DIRECTION ANALYSIS TESTING ===\n\n")

# Load real DE data
de_data_path <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds"
cat("Loading DE data from:", de_data_path, "\n")
full_de <- readRDS(de_data_path)

# Load experiment weights
cat("\nLoading experiment weights...\n")
cell_counts <- load_crispri_cell_counts()
experiment_weights <- calculate_experiment_weights(cell_counts)

# Test specific genes and clusters
test_cases <- list(
  list(gene = "LRRK2", cluster = "cluster_0", expected = "opposing"),
  list(gene = "LRRK2", cluster = "cluster_5", expected = "opposing"),
  list(gene = "SNCA", cluster = "cluster_0", expected = "same"),
  list(gene = "SNCA", cluster = "cluster_5", expected = "same"),
  list(gene = "PINK1", cluster = "cluster_0", expected = "same"),
  list(gene = "PARK7", cluster = "cluster_0", expected = "same")
)

cat("\nTesting with real differential expression data:\n")
cat("==============================================\n")

results_summary <- list()

for (test_case in test_cases) {
  gene <- test_case$gene
  cluster <- test_case$cluster
  expected <- test_case$expected
  
  cat("\n\nAnalyzing", gene, "in", cluster, "(expecting", expected, "pattern):\n")
  cat(paste(rep("-", 60), collapse=""), "\n")
  
  # Get MAST data
  mast_gene <- gene
  if (gene == "SNCA") {
    # Use SNCA_A30P as representative for SNCA
    mast_gene <- "SNCA_A30P"
  } else if (gene == "PARK2") {
    mast_gene <- "PRKN"
  }
  
  # Check if data exists
  if (!mast_gene %in% names(full_de$iSCORE_PD_MAST)) {
    cat("  WARNING: No MAST data for", mast_gene, "\n")
    next
  }
  
  if (!cluster %in% names(full_de$iSCORE_PD_MAST[[mast_gene]])) {
    cat("  WARNING: No MAST data for", mast_gene, "in", cluster, "\n")
    next
  }
  
  mast_data <- full_de$iSCORE_PD_MAST[[mast_gene]][[cluster]]$results
  
  # Get CRISPRi data (need to handle multiple experiments)
  crispri_gene <- gene
  if (gene == "PRKN") {
    crispri_gene <- "PARK2"
  }
  
  if (!crispri_gene %in% names(full_de$CRISPRi_Mixscale)) {
    cat("  WARNING: No CRISPRi data for", crispri_gene, "\n")
    next
  }
  
  if (!cluster %in% names(full_de$CRISPRi_Mixscale[[crispri_gene]])) {
    cat("  WARNING: No CRISPRi data for", crispri_gene, "in", cluster, "\n")
    next
  }
  
  crispri_data <- full_de$CRISPRi_Mixscale[[crispri_gene]][[cluster]]$results
  
  # Extract background genes (all genes tested)
  mast_background <- rownames(mast_data)
  crispri_background <- if (!is.null(rownames(crispri_data))) {
    rownames(crispri_data)
  } else if ("gene" %in% names(crispri_data)) {
    crispri_data$gene
  } else {
    crispri_data[[1]]
  }
  
  background_genes <- intersect(mast_background, crispri_background)
  
  cat("  Data summary:\n")
  cat("    MAST genes tested:", length(mast_background), "\n")
  cat("    CRISPRi genes tested:", length(crispri_background), "\n")
  cat("    Background intersection:", length(background_genes), "\n")
  
  # Run enhanced direction analysis
  direction_result <- enhanced_direction_analysis(
    mast_data = mast_data,
    crispri_data = crispri_data,
    gene_name = gene,
    background_genes = background_genes,
    lfc_threshold = 0.25,
    p_threshold = 0.05
  )
  
  # Store results
  result_key <- paste(gene, cluster, sep = "_")
  results_summary[[result_key]] <- list(
    gene = gene,
    cluster = cluster,
    expected = expected,
    observed = direction_result$primary_pattern,
    biological_expectation = direction_result$biological_expectation,
    same_direction_count = direction_result$same_direction$overlap_count,
    opposite_direction_count = direction_result$opposite_direction$overlap_count,
    same_p = direction_result$same_direction$fisher_p,
    opposite_p = direction_result$opposite_direction$fisher_p,
    combined_p = direction_result$combined_analysis$combined_p
  )
  
  # Print detailed results
  cat("\n  Direction analysis results:\n")
  cat("    Biological expectation:", direction_result$biological_expectation, "\n")
  cat("    Primary pattern detected:", direction_result$primary_pattern, "\n")
  cat("    Same-direction overlap:", direction_result$same_direction$overlap_count, "genes\n")
  cat("    Opposite-direction overlap:", direction_result$opposite_direction$overlap_count, "genes\n")
  cat("    Same-direction p-value:", format(direction_result$same_direction$fisher_p, scientific = TRUE), "\n")
  cat("    Opposite-direction p-value:", format(direction_result$opposite_direction$fisher_p, scientific = TRUE), "\n")
  cat("    Combined p-value:", format(direction_result$combined_analysis$combined_p, scientific = TRUE), "\n")
  
  # Check if result matches expectation
  if (direction_result$biological_expectation == expected) {
    cat("    ✓ PASS: Biological expectation matches expected pattern\n")
  } else {
    cat("    ✗ FAIL: Biological expectation mismatch\n")
  }
  
  # Show example genes if overlaps exist
  if (direction_result$same_direction$overlap_count > 0) {
    example_same <- head(direction_result$same_direction$all_same_direction_genes, 5)
    cat("    Example same-direction genes:", paste(example_same, collapse = ", "), "\n")
  }
  
  if (direction_result$opposite_direction$overlap_count > 0) {
    example_opposite <- head(direction_result$opposite_direction$all_opposite_direction_genes, 5)
    cat("    Example opposite-direction genes:", paste(example_opposite, collapse = ", "), "\n")
  }
}

# Summary report
cat("\n\n=== SUMMARY REPORT ===\n")
cat("======================\n")

pass_count <- 0
fail_count <- 0

for (result in results_summary) {
  status <- if (result$biological_expectation == result$expected) "PASS" else "FAIL"
  if (status == "PASS") pass_count <- pass_count + 1 else fail_count <- fail_count + 1
  
  cat(sprintf("%-8s %-10s: Expected %-8s, Got %-8s [%s]\n",
              result$gene, result$cluster, result$expected, 
              result$biological_expectation, status))
  cat(sprintf("  → Same: %d genes (p=%.2e), Opposite: %d genes (p=%.2e)\n",
              result$same_direction_count, result$same_p,
              result$opposite_direction_count, result$opposite_p))
}

cat("\nOverall: ", pass_count, "PASS,", fail_count, "FAIL\n")

# Test with multiple experiments using enhanced overlap significance
cat("\n\n=== TESTING ENHANCED OVERLAP WITH EXPERIMENT WEIGHTING ===\n")
cat("==========================================================\n")

# Prepare CRISPRi data in expected format (list of experiments)
# Extract data for each experiment separately
crispri_data_raw <- full_de$CRISPRi_Mixscale[["LRRK2"]][["cluster_0"]]$results
experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
crispri_experiments_list <- list()

for (exp in experiments) {
  # Create experiment-specific data frame
  log2fc_col <- paste0("log2FC_", exp)
  pval_col <- paste0("p_cell_type", exp, ":weight")
  
  if (log2fc_col %in% names(crispri_data_raw) && pval_col %in% names(crispri_data_raw)) {
    exp_data <- data.frame(
      log2FC = crispri_data_raw[[log2fc_col]],
      p_val_adj = crispri_data_raw[[pval_col]],
      row.names = rownames(crispri_data_raw),
      stringsAsFactors = FALSE
    )
    crispri_experiments_list[[exp]] <- exp_data
  }
}

# Test LRRK2 with enhanced analysis
test_enhanced <- calculate_enhanced_overlap_significance(
  mast_data = full_de$iSCORE_PD_MAST[["LRRK2"]][["cluster_0"]]$results,
  crispri_experiments_data = crispri_experiments_list,
  gene_name = "LRRK2",
  cluster_id = "cluster_0", 
  experiment_weights = experiment_weights,
  use_enhanced_analysis = TRUE
)

cat("\nLRRK2 Enhanced Analysis Results:\n")
cat("  Weighted meta p-value:", format(test_enhanced$weighted_meta_p, scientific = TRUE), "\n")
cat("  Direction analysis included:", !is.null(test_enhanced$enhanced_statistics), "\n")
if (!is.null(test_enhanced$enhanced_statistics)) {
  cat("  Primary direction pattern:", test_enhanced$enhanced_statistics$primary_pattern, "\n")
  cat("  Biological expectation:", test_enhanced$enhanced_statistics$biological_expectation, "\n")
}

cat("\n=== TESTING COMPLETE ===\n")