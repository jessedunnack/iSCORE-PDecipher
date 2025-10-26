#!/usr/bin/env Rscript
# Convert All Pooled MixScale Datasets
# Date: October 26, 2025
# Purpose: Create all 6 datasets (FPD_BH already done, creating remaining 5)

# Source the conversion functions
source("convert_pooled_to_full_de.R")

cat("\n=== BATCH CONVERSION: ALL POOLED MIXSCALE DATASETS ===\n\n")

# Define all datasets to create
datasets <- list(
  # FPD datasets
  list(
    name = "FPD_uncorrected",
    source_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit",
    output_dir = "/mnt/e/THESIS/scRNASeq/mixscale/FPD_uncorrected_dataset",
    dataset_name = "FPD",
    pval_column = "p_weight",
    enrichment_source = "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_FPD_p_weight"
  ),
  list(
    name = "FPD_bonferroni",
    source_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit",
    output_dir = "/mnt/e/THESIS/scRNASeq/mixscale/FPD_bonferroni_dataset",
    dataset_name = "FPD",
    pval_column = "p_weight_bonferroni",
    enrichment_source = "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_FPD_p_weight_bonferroni"
  ),

  # CRISPRi datasets
  list(
    name = "CRISPRi_BH",
    source_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit",
    output_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_BH_dataset",
    dataset_name = "CRISPRi",
    pval_column = "p_weight_BH",
    enrichment_source = "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_CRISPRi_p_weight_BH"
  ),
  list(
    name = "CRISPRi_uncorrected",
    source_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit",
    output_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_uncorrected_dataset",
    dataset_name = "CRISPRi",
    pval_column = "p_weight",
    enrichment_source = "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_CRISPRi_p_weight"
  ),
  list(
    name = "CRISPRi_bonferroni",
    source_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit",
    output_dir = "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_bonferroni_dataset",
    dataset_name = "CRISPRi",
    pval_column = "p_weight_bonferroni",
    enrichment_source = "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_CRISPRi_p_weight_bonferroni"
  )
)

# Track success/failure
results <- list()
successful <- 0
failed <- 0

# Process each dataset
for (i in seq_along(datasets)) {
  ds <- datasets[[i]]

  cat("\n")
  cat("=" , rep("=", 70), "=\n", sep = "")
  cat("DATASET ", i, "/", length(datasets), ": ", ds$name, "\n", sep = "")
  cat("=" , rep("=", 70), "=\n", sep = "")
  cat("\n")

  tryCatch({
    # Convert
    convert_pooled_to_full_de(
      source_dir = ds$source_dir,
      output_dir = ds$output_dir,
      dataset_name = ds$dataset_name,
      pval_column = ds$pval_column,
      enrichment_source = ds$enrichment_source
    )

    # Validate
    validate_converted_dataset(ds$output_dir)

    # Record success
    results[[ds$name]] <- "SUCCESS"
    successful <- successful + 1

    cat("\n✓ ", ds$name, " COMPLETE\n", sep = "")

  }, error = function(e) {
    # Record failure
    results[[ds$name]] <- paste("FAILED:", e$message)
    failed <- failed + 1

    cat("\n✗ ", ds$name, " FAILED: ", e$message, "\n", sep = "")
  })
}

# Final summary
cat("\n")
cat("=" , rep("=", 70), "=\n", sep = "")
cat("BATCH CONVERSION SUMMARY\n")
cat("=" , rep("=", 70), "=\n", sep = "")
cat("\n")

cat("Total datasets processed:", length(datasets), "\n")
cat("Successful:", successful, "\n")
cat("Failed:", failed, "\n\n")

cat("Results:\n")
for (name in names(results)) {
  status <- results[[name]]
  symbol <- if (status == "SUCCESS") "✓" else "✗"
  cat("  ", symbol, " ", name, ": ", status, "\n", sep = "")
}

if (failed == 0) {
  cat("\n🎉 ALL DATASETS CONVERTED SUCCESSFULLY! 🎉\n")
  cat("\nNext steps:\n")
  cat("1. Update R/dataset_validator.R → get_dataset_options()\n")
  cat("2. Test launch_app() with each dataset\n")
  cat("3. Final testing and commit\n\n")
} else {
  cat("\n⚠ SOME DATASETS FAILED - Please review errors above\n\n")
}
