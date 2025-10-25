# Test script for pooled MixScale import functions
# Created: October 24, 2025
# Purpose: Validate new import functions with actual FDR-corrected data

# Source the new import functions
source("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/R/import_pooled_mixscale_functions.R")

# Initialize test results
test_results <- list()

cat("===============================================\n")
cat("TESTING POOLED MIXSCALE IMPORT FUNCTIONS\n")
cat("===============================================\n\n")

# Define paths
fpd_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/"
crispri_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/"

# ========================================
# TEST 1: FPD with BH correction (recommended)
# ========================================
cat("TEST 1: Loading FPD data with p_weight_BH...\n")
test_results$fpd_bh <- tryCatch({
  fpd_bh <- import_pooled_mixscale_data(
    fpd_dir,
    pval_column = "p_weight_BH",
    dataset_type = "FPD"
  )

  cat("✓ Successfully loaded FPD with BH correction\n")
  cat("  - Perturbations:", length(fpd_bh), "\n")
  cat("  - Clusters:", length(unique(unlist(lapply(fpd_bh, names)))), "\n")
  cat("  - First perturbation:", names(fpd_bh)[1], "\n")

  # Validate structure
  first_pert <- fpd_bh[[1]]
  first_cluster <- first_pert[[1]]

  cat("  - Metadata keys:", paste(names(first_cluster$metadata), collapse=", "), "\n")
  cat("  - P-value column used:", first_cluster$metadata$pval_column_used, "\n")
  cat("  - Available p-value columns:", paste(first_cluster$available_pval_columns, collapse=", "), "\n")

  list(status = "PASS", data = fpd_bh)
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# TEST 2: FPD with original p-values
# ========================================
cat("TEST 2: Loading FPD data with p_weight (uncorrected)...\n")
test_results$fpd_uncorrected <- tryCatch({
  fpd_uncorrected <- import_pooled_mixscale_data(
    fpd_dir,
    pval_column = "p_weight",
    dataset_type = "FPD"
  )

  cat("✓ Successfully loaded FPD with uncorrected p-values\n")
  cat("  - Perturbations:", length(fpd_uncorrected), "\n")

  list(status = "PASS", data = fpd_uncorrected)
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# TEST 3: FPD with Bonferroni correction
# ========================================
cat("TEST 3: Loading FPD data with p_weight_bonferroni...\n")
test_results$fpd_bonferroni <- tryCatch({
  fpd_bonf <- import_pooled_mixscale_data(
    fpd_dir,
    pval_column = "p_weight_bonferroni",
    dataset_type = "FPD"
  )

  cat("✓ Successfully loaded FPD with Bonferroni correction\n")
  cat("  - Perturbations:", length(fpd_bonf), "\n")

  list(status = "PASS", data = fpd_bonf)
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# TEST 4: CRISPRi with BH correction
# ========================================
cat("TEST 4: Loading CRISPRi data with p_weight_BH...\n")
test_results$crispri_bh <- tryCatch({
  crispri_bh <- import_pooled_mixscale_data(
    crispri_dir,
    pval_column = "p_weight_BH",
    dataset_type = "CRISPRi"
  )

  cat("✓ Successfully loaded CRISPRi with BH correction\n")
  cat("  - Perturbations:", length(crispri_bh), "\n")
  cat("  - Clusters:", length(unique(unlist(lapply(crispri_bh, names)))), "\n")

  list(status = "PASS", data = crispri_bh)
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# TEST 5: Detect format function
# ========================================
cat("TEST 5: Testing detect_mixscale_format()...\n")
test_results$format_detection <- tryCatch({
  # Load a test file directly
  test_file <- file.path(fpd_dir,
    "all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0",
    "all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds")

  test_data <- readRDS(test_file)
  format <- detect_mixscale_format(test_data)

  cat("✓ Format detection successful\n")
  cat("  - Detected format:", format, "\n")

  if (format == "pooled") {
    cat("  - Correctly identified as pooled\n")
    list(status = "PASS", format = format)
  } else {
    cat("  - ERROR: Should be pooled but detected as:", format, "\n")
    list(status = "FAIL", format = format)
  }
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# TEST 6: Data integrity checks
# ========================================
cat("TEST 6: Validating data integrity...\n")
test_results$data_integrity <- tryCatch({
  # Use FPD BH data for validation
  if (test_results$fpd_bh$status == "PASS") {
    fpd_data <- test_results$fpd_bh$data

    # Check a specific perturbation-cluster combination
    first_pert_name <- names(fpd_data)[1]
    first_pert <- fpd_data[[first_pert_name]]
    first_cluster_name <- names(first_pert)[1]
    first_cluster <- first_pert[[first_cluster_name]]

    # Validate structure
    has_results <- !is.null(first_cluster$results)
    has_metadata <- !is.null(first_cluster$metadata)
    has_background <- !is.null(first_cluster$background_genes)

    # Validate columns in results
    results_df <- first_cluster$results
    required_cols <- c("gene_ID", "log2FC", "p_weight_BH")
    has_all_cols <- all(required_cols %in% colnames(results_df))

    # Check for NAs
    has_na_geneid <- any(is.na(results_df$gene_ID))
    has_na_log2fc <- any(is.na(results_df$log2FC))

    cat("✓ Data integrity checks:\n")
    cat("  - Has results:", has_results, "\n")
    cat("  - Has metadata:", has_metadata, "\n")
    cat("  - Has background genes:", has_background, "\n")
    cat("  - Has all required columns:", has_all_cols, "\n")
    cat("  - Number of background genes:", length(first_cluster$background_genes), "\n")
    cat("  - NA in gene_ID:", has_na_geneid, "\n")
    cat("  - NA in log2FC:", has_na_log2fc, "\n")

    if (has_results && has_metadata && has_background && has_all_cols) {
      list(status = "PASS")
    } else {
      list(status = "FAIL", reason = "Missing required components")
    }
  } else {
    cat("  Skipping - FPD BH data not loaded\n")
    list(status = "SKIP")
  }
}, error = function(e) {
  cat("✗ FAILED:", e$message, "\n")
  list(status = "FAIL", error = e$message)
})

cat("\n")

# ========================================
# SUMMARY
# ========================================
cat("===============================================\n")
cat("TEST SUMMARY\n")
cat("===============================================\n")

passed <- sum(sapply(test_results, function(x) x$status == "PASS"))
failed <- sum(sapply(test_results, function(x) x$status == "FAIL"))
skipped <- sum(sapply(test_results, function(x) x$status == "SKIP"))

cat("Passed: ", passed, "/", length(test_results), "\n")
cat("Failed: ", failed, "/", length(test_results), "\n")
cat("Skipped:", skipped, "/", length(test_results), "\n")

if (failed == 0) {
  cat("\n✓ ALL TESTS PASSED!\n")
} else {
  cat("\n✗ SOME TESTS FAILED - Review output above\n")
}

cat("\n")

# Save test results
saveRDS(test_results,
  "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/tests/pooled_import_test_results.rds")
cat("Test results saved to: tests/pooled_import_test_results.rds\n")
