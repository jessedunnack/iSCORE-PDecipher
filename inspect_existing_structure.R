#!/usr/bin/env Rscript
# PHASE 1: Inspect existing full_DE_results.rds structure
# This tells us EXACTLY what format the app expects

cat("\n=== INSPECTING EXISTING full_DE_results.rds ===\n\n")

# Load the existing file
existing_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds"

cat("Loading file:", existing_file, "\n")
cat("File size:", file.size(existing_file) / 1024^2, "MB\n\n")

full_de <- readRDS(existing_file)

# Check top-level structure
cat("=== TOP-LEVEL STRUCTURE ===\n")
cat("Class:", class(full_de), "\n")
cat("Names:", paste(names(full_de), collapse = ", "), "\n\n")

# Check iSCORE_PD_MAST structure
if ("iSCORE_PD_MAST" %in% names(full_de)) {
  cat("=== iSCORE_PD_MAST STRUCTURE ===\n")
  mast_data <- full_de$iSCORE_PD_MAST
  cat("Is NULL:", is.null(mast_data), "\n")
  if (!is.null(mast_data)) {
    cat("Number of mutations:", length(mast_data), "\n")
    cat("First 3 mutations:", paste(head(names(mast_data), 3), collapse = ", "), "\n")

    # Check first mutation structure
    first_mut <- mast_data[[1]]
    cat("\nFirst mutation clusters:", paste(names(first_mut), collapse = ", "), "\n")

    # Check first cluster structure
    first_cluster <- first_mut[[1]]
    cat("First cluster elements:", paste(names(first_cluster), collapse = ", "), "\n")

    # Check results dataframe
    if ("results" %in% names(first_cluster)) {
      results_df <- first_cluster$results
      cat("Results columns:", paste(colnames(results_df), collapse = ", "), "\n")
      cat("Results rows:", nrow(results_df), "\n")
    }
  }
  cat("\n")
}

# Check CRISPRi_Mixscale structure
if ("CRISPRi_Mixscale" %in% names(full_de)) {
  cat("=== CRISPRi_Mixscale STRUCTURE ===\n")
  mixscale_data <- full_de$CRISPRi_Mixscale
  cat("Is NULL:", is.null(mixscale_data), "\n")
  if (!is.null(mixscale_data)) {
    cat("Number of perturbations:", length(mixscale_data), "\n")
    cat("First 5 perturbations:", paste(head(names(mixscale_data), 5), collapse = ", "), "\n")

    # Check first perturbation structure
    first_pert <- mixscale_data[[1]]
    cat("\nFirst perturbation clusters:", paste(names(first_pert), collapse = ", "), "\n")

    # Check first cluster structure
    first_cluster <- first_pert[[1]]
    cat("First cluster elements:", paste(names(first_cluster), collapse = ", "), "\n")

    # Check results dataframe
    if ("results" %in% names(first_cluster)) {
      results_df <- first_cluster$results
      cat("Results class:", class(results_df), "\n")
      cat("Results columns:", paste(colnames(results_df), collapse = ", "), "\n")
      cat("Results rows:", nrow(results_df), "\n")
      cat("\nFirst few rows:\n")
      print(head(results_df, 3))
    }
  }
  cat("\n")
}

cat("\n=== STRUCTURE INSPECTION COMPLETE ===\n")
cat("This is the EXACT format we need to match!\n")
