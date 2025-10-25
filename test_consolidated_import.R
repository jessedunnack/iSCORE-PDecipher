#!/usr/bin/env Rscript
# Test Consolidated Enrichment Import
# Tests that the updated import_enrichment_with_correction() function works correctly

cat("\n========================================================\n")
cat("TESTING CONSOLIDATED ENRICHMENT IMPORT\n")
cat("========================================================\n\n")

# Source the import functions
source("R/import_pooled_mixscale_functions.R")

# Test 1: Load FPD with BH correction (recommended)
cat("Test 1: Loading FPD with BH correction...\n")
fpd_bh <- tryCatch({
  import_enrichment_with_correction(
    dataset = "FPD",
    pval_correction = "BH"
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(fpd_bh) && is.data.frame(fpd_bh)) {
  cat("✓ SUCCESS: Loaded", nrow(fpd_bh), "enrichment terms\n")
  cat("  - Perturbations:", length(unique(fpd_bh$mutation_perturbation)), "\n")
  cat("  - Clusters:", paste(sort(unique(fpd_bh$cluster)), collapse=", "), "\n")
  cat("  - Enrichment types:", paste(unique(fpd_bh$enrichment_type), collapse=", "), "\n")
  cat("  - Consolidated:", attr(fpd_bh, "is_consolidated"), "\n")
} else {
  cat("✗ FAILED: Could not load FPD BH data\n")
}

cat("\n")

# Test 2: Load CRISPRi with uncorrected p-values
cat("Test 2: Loading CRISPRi with uncorrected p-values...\n")
crispri_none <- tryCatch({
  import_enrichment_with_correction(
    dataset = "CRISPRi",
    pval_correction = "none"
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(crispri_none) && is.data.frame(crispri_none)) {
  cat("✓ SUCCESS: Loaded", nrow(crispri_none), "enrichment terms\n")
  cat("  - Perturbations:", length(unique(crispri_none$mutation_perturbation)), "\n")
  cat("  - Clusters:", paste(sort(unique(crispri_none$cluster)), collapse=", "), "\n")
  cat("  - Consolidated:", attr(crispri_none, "is_consolidated"), "\n")
} else {
  cat("✗ FAILED: Could not load CRISPRi uncorrected data\n")
}

cat("\n")

# Test 3: Load CRISPRi with Bonferroni correction
cat("Test 3: Loading CRISPRi with Bonferroni correction...\n")
crispri_bonf <- tryCatch({
  import_enrichment_with_correction(
    dataset = "CRISPRi",
    pval_correction = "bonferroni"
  )
}, error = function(e) {
  cat("ERROR:", e$message, "\n")
  return(NULL)
})

if (!is.null(crispri_bonf) && is.data.frame(crispri_bonf)) {
  cat("✓ SUCCESS: Loaded", nrow(crispri_bonf), "enrichment terms\n")
  cat("  - Perturbations:", length(unique(crispri_bonf$mutation_perturbation)), "\n")
  cat("  - Clusters:", paste(sort(unique(crispri_bonf$cluster)), collapse=", "), "\n")
  cat("  - Consolidated:", attr(crispri_bonf, "is_consolidated"), "\n")
} else {
  cat("✗ FAILED: Could not load CRISPRi Bonferroni data\n")
}

cat("\n========================================================\n")
cat("TESTING COMPLETE\n")
cat("========================================================\n\n")
