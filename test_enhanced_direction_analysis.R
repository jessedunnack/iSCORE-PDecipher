#!/usr/bin/env Rscript

#' Test Enhanced Direction Analysis v0.2.6
#' 
#' This script tests the enhanced direction-aware analysis to verify:
#' 1. LRRK2 correctly detected as "opposing" direction expectation
#' 2. SNCA variants correctly detected as "same" direction expectation
#' 3. Direction-aware Fisher's tests work properly
#' 4. Biological context weighting functions correctly

# Load required functions
library(dplyr)
source("R/enhanced_direction_analysis.R")
source("R/gene_harmonization.R")
source("R/signature_analysis.R")

cat("=== ENHANCED DIRECTION ANALYSIS TESTING ===\n\n")

# Test 1: Biological Direction Expectations
cat("TEST 1: Biological Direction Expectations\n")
cat("=========================================\n")

test_genes <- c("LRRK2", "SNCA_A30P", "SNCA_A53T", "SNCA", "PINK1", "PARK2", "PRKN", "VPS13C_A444P", "ATP13A2")

for (gene in test_genes) {
  expectation <- get_biological_direction_expectation(gene)
  cat(sprintf("%-12s: %s\n", gene, expectation))
}

cat("\n")

# Test 2: Test with Mock Data
cat("TEST 2: Mock Data Direction Analysis\n")
cat("====================================\n")

# Create mock MAST data (LRRK2 upregulates some genes)
mock_mast_lrrk2 <- data.frame(
  avg_log2FC = c(0.5, -0.3, 0.8, -0.4, 0.6, -0.2, 0.3, -0.5),
  p_val_adj = c(0.01, 0.02, 0.001, 0.03, 0.005, 0.04, 0.02, 0.01),
  row.names = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8"),
  stringsAsFactors = FALSE
)

# Create mock CRISPRi data (LRRK2 knockdown should have opposite effects)
mock_crispri_lrrk2 <- data.frame(
  log2FC_C12_FPD_24 = c(-0.4, 0.4, -0.7, 0.3, -0.5, 0.3, -0.2, 0.6),  # Opposite to MAST
  p_val_adj = c(0.02, 0.01, 0.002, 0.04, 0.008, 0.03, 0.04, 0.015),
  gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8"),
  stringsAsFactors = FALSE
)

# Create mock SNCA data (should show same-direction effects)
mock_mast_snca <- data.frame(
  avg_log2FC = c(0.4, -0.3, 0.6, -0.5, 0.3, -0.4, 0.5, -0.2),
  p_val_adj = c(0.01, 0.02, 0.003, 0.01, 0.04, 0.02, 0.01, 0.045),
  row.names = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8"),
  stringsAsFactors = FALSE
)

mock_crispri_snca <- data.frame(
  log2FC_C12_FPD_24 = c(0.3, -0.4, 0.5, -0.4, 0.4, -0.3, 0.4, -0.3),  # Same direction as MAST
  p_val_adj = c(0.02, 0.01, 0.005, 0.02, 0.03, 0.03, 0.015, 0.04),
  gene = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8"),
  stringsAsFactors = FALSE
)

# Test LRRK2 analysis (expecting opposite-direction to be significant)
cat("Testing LRRK2 (should detect opposing effects):\n")
lrrk2_result <- enhanced_direction_analysis(
  mast_data = mock_mast_lrrk2,
  crispri_data = mock_crispri_lrrk2,
  gene_name = "LRRK2",
  background_genes = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8")
)

cat("\nTesting SNCA (should detect same-direction effects):\n")
snca_result <- enhanced_direction_analysis(
  mast_data = mock_mast_snca,
  crispri_data = mock_crispri_snca,
  gene_name = "SNCA_A30P",
  background_genes = c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8")
)

# Test 3: Detailed Results Analysis
cat("\n\nTEST 3: Detailed Results Analysis\n")
cat("==================================\n")

cat("LRRK2 Results:\n")
cat("  Biological expectation:", lrrk2_result$biological_expectation, "\n")
cat("  Primary pattern:", lrrk2_result$primary_pattern, "\n")
cat("  Same-direction overlap count:", lrrk2_result$same_direction$overlap_count, "\n")
cat("  Opposite-direction overlap count:", lrrk2_result$opposite_direction$overlap_count, "\n")
cat("  Combined p-value:", lrrk2_result$combined_analysis$combined_p, "\n")

if (lrrk2_result$biological_expectation == "opposing" && lrrk2_result$primary_pattern == "opposite") {
  cat("  ✓ PASS: LRRK2 correctly identified as opposing pattern\n")
} else {
  cat("  ✗ FAIL: LRRK2 should be opposing pattern\n")
}

cat("\nSNCA Results:\n")
cat("  Biological expectation:", snca_result$biological_expectation, "\n")
cat("  Primary pattern:", snca_result$primary_pattern, "\n")
cat("  Same-direction overlap count:", snca_result$same_direction$overlap_count, "\n")
cat("  Opposite-direction overlap count:", snca_result$opposite_direction$overlap_count, "\n")
cat("  Combined p-value:", snca_result$combined_analysis$combined_p, "\n")

if (snca_result$biological_expectation == "same") {
  cat("  ✓ PASS: SNCA correctly identified as same-direction expectation\n")
} else {
  cat("  ✗ FAIL: SNCA should be same-direction expectation\n")
}

# Test 4: P-value Combination
cat("\n\nTEST 4: P-value Combination Testing\n")
cat("===================================\n")

# Test weighted combination for LRRK2 (should prioritize opposite-direction)
lrrk2_combination <- combine_direction_pvalues(
  same_direction_p = 0.1,
  opposite_direction_p = 0.01,
  gene_name = "LRRK2"
)

cat("LRRK2 p-value combination:\n")
cat("  Primary test:", lrrk2_combination$primary_test, "\n")
cat("  Primary p-value:", lrrk2_combination$primary_p, "\n")
cat("  Secondary p-value:", lrrk2_combination$secondary_p, "\n")
cat("  Combined p-value:", lrrk2_combination$combined_p, "\n")

if (lrrk2_combination$primary_test == "opposite") {
  cat("  ✓ PASS: LRRK2 correctly prioritizes opposite-direction test\n")
} else {
  cat("  ✗ FAIL: LRRK2 should prioritize opposite-direction test\n")
}

# Test weighted combination for SNCA (should prioritize same-direction)
snca_combination <- combine_direction_pvalues(
  same_direction_p = 0.01,
  opposite_direction_p = 0.1,
  gene_name = "SNCA_A30P"
)

cat("\nSNCA p-value combination:\n")
cat("  Primary test:", snca_combination$primary_test, "\n")
cat("  Primary p-value:", snca_combination$primary_p, "\n")
cat("  Secondary p-value:", snca_combination$secondary_p, "\n")
cat("  Combined p-value:", snca_combination$combined_p, "\n")

if (snca_combination$primary_test == "same") {
  cat("  ✓ PASS: SNCA correctly prioritizes same-direction test\n")
} else {
  cat("  ✗ FAIL: SNCA should prioritize same-direction test\n")
}

cat("\n=== TESTING COMPLETE ===\n")