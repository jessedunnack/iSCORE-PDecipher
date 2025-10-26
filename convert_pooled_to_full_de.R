#!/usr/bin/env Rscript
# Convert Pooled MixScale Data to launch_app() Format
# Date: October 26, 2025
# Purpose: Transform cluster-organized pooled MixScale data into the
#          perturbation-organized structure expected by launch_app()

cat("\n=== POOLED MIXSCALE TO FULL_DE_RESULTS CONVERTER ===\n\n")

#' Convert pooled MixScale cluster data to full_DE_results format
#'
#' @param source_dir Directory containing cluster subdirectories with *_mixscale_DEGs.rds files
#' @param output_dir Directory where converted dataset will be created
#' @param dataset_name "FPD" or "CRISPRi"
#' @param pval_column Which p-value column to use: "p_weight", "p_weight_BH", or "p_weight_bonferroni"
#' @param enrichment_source Directory containing enrichment results (will be copied)
#' @return Invisible TRUE if successful
convert_pooled_to_full_de <- function(
  source_dir,
  output_dir,
  dataset_name,
  pval_column = "p_weight_BH",
  enrichment_source = NULL
) {

  # Validate inputs
  if (!dir.exists(source_dir)) {
    stop("Source directory does not exist: ", source_dir)
  }

  valid_pval_cols <- c("p_weight", "p_weight_BH", "p_weight_bonferroni")
  if (!pval_column %in% valid_pval_cols) {
    stop("pval_column must be one of: ", paste(valid_pval_cols, collapse = ", "))
  }

  cat("Configuration:\n")
  cat("  Source:", source_dir, "\n")
  cat("  Output:", output_dir, "\n")
  cat("  Dataset:", dataset_name, "\n")
  cat("  P-value column:", pval_column, "\n\n")

  # Find all cluster RDS files
  cat("Searching for cluster RDS files...\n")
  cluster_files <- list.files(
    source_dir,
    pattern = "_mixscale_DEGs\\.rds$",
    recursive = TRUE,
    full.names = TRUE
  )

  if (length(cluster_files) == 0) {
    stop("No *_mixscale_DEGs.rds files found in ", source_dir)
  }

  cat("Found", length(cluster_files), "cluster files\n\n")

  # Initialize output structure
  full_DE_results <- list(
    iSCORE_PD_MAST = NULL,
    CRISPRi_Mixscale = list()
  )

  # Process each cluster file
  for (cluster_file in cluster_files) {
    # Extract cluster number from filename
    cluster_num <- gsub(".*_clust_(\\d+)_mixscale_DEGs\\.rds$", "\\1", basename(cluster_file))
    cluster_name <- paste0("cluster_", cluster_num)

    cat("Processing", cluster_name, "...\n")

    # Load cluster data
    cluster_data <- readRDS(cluster_file)

    cat("  Found", length(cluster_data), "perturbations\n")

    # Process each perturbation in this cluster
    for (pert_name in names(cluster_data)) {
      pert_df <- cluster_data[[pert_name]]

      # Validate required columns exist
      required_cols <- c("gene_ID", "log2FC", pval_column)
      missing_cols <- setdiff(required_cols, colnames(pert_df))
      if (length(missing_cols) > 0) {
        stop("Missing required columns in ", pert_name, ": ",
             paste(missing_cols, collapse = ", "))
      }

      # Transform dataframe
      results <- pert_df

      # Move gene_ID to rownames
      rownames(results) <- results$gene_ID

      # Rename columns to match expected format
      results$avg_log2FC <- results$log2FC
      results$p_val_adj <- results[[pval_column]]
      results$p_val <- results$p_weight  # Keep original for compatibility

      # Add dummy pct columns (required by app, but not available in pooled data)
      results$pct.1 <- NA
      results$pct.2 <- NA

      # Initialize perturbation list if needed
      if (!pert_name %in% names(full_DE_results$CRISPRi_Mixscale)) {
        full_DE_results$CRISPRi_Mixscale[[pert_name]] <- list()
      }

      # Store in nested structure: perturbation → cluster → results
      full_DE_results$CRISPRi_Mixscale[[pert_name]][[cluster_name]] <- list(
        results = results,
        metadata = list(
          source = "pooled_mixscale",
          dataset = dataset_name,
          pval_correction = pval_column,
          cluster = cluster_num
        ),
        background_genes = character()
      )
    }
  }

  # Report summary
  n_perturbations <- length(full_DE_results$CRISPRi_Mixscale)
  n_clusters <- length(unique(gsub(".*_clust_(\\d+)_.*", "\\1", basename(cluster_files))))

  cat("\n=== CONVERSION SUMMARY ===\n")
  cat("Total perturbations:", n_perturbations, "\n")
  cat("Total clusters:", n_clusters, "\n")
  cat("First 5 perturbations:", paste(head(names(full_DE_results$CRISPRi_Mixscale), 5), collapse = ", "), "\n")

  # Create output directory
  if (!dir.exists(output_dir)) {
    cat("\nCreating output directory:", output_dir, "\n")
    dir.create(output_dir, recursive = TRUE)
  }

  # Save converted data
  output_file <- file.path(output_dir, "full_DE_results.rds")
  cat("Saving to:", output_file, "\n")
  saveRDS(full_DE_results, output_file, compress = TRUE)

  cat("File size:", round(file.size(output_file) / 1024^2, 2), "MB\n")

  # Copy enrichment file if provided
  if (!is.null(enrichment_source) && dir.exists(enrichment_source)) {
    cat("\n=== COPYING ENRICHMENT FILE ===\n")

    enrich_file <- file.path(enrichment_source, "all_enrichment_padj005_complete_with_direction.rds")

    if (file.exists(enrich_file)) {
      output_enrich <- file.path(output_dir, "all_enrichment_padj005_complete_with_direction.rds")
      cat("Copying enrichment from:", enrich_file, "\n")
      cat("To:", output_enrich, "\n")
      file.copy(enrich_file, output_enrich, overwrite = TRUE)
      cat("Enrichment file size:", round(file.size(output_enrich) / 1024^2, 2), "MB\n")
    } else {
      warning("Enrichment file not found: ", enrich_file)
    }
  }

  cat("\n=== CONVERSION COMPLETE ===\n")

  invisible(TRUE)
}

#' Validate converted dataset
#'
#' @param data_dir Directory containing full_DE_results.rds
#' @return TRUE if valid, otherwise stops with error
validate_converted_dataset <- function(data_dir) {
  cat("\n=== VALIDATING CONVERTED DATASET ===\n\n")

  # Check file exists
  de_file <- file.path(data_dir, "full_DE_results.rds")
  enrich_file <- file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds")

  cat("Checking files exist...\n")
  if (!file.exists(de_file)) {
    stop("Missing: full_DE_results.rds")
  }
  cat("  ✓ full_DE_results.rds exists\n")

  if (!file.exists(enrich_file)) {
    warning("  ⚠ all_enrichment_padj005_complete_with_direction.rds not found")
  } else {
    cat("  ✓ all_enrichment_padj005_complete_with_direction.rds exists\n")
  }

  # Load and validate structure
  cat("\nValidating structure...\n")
  full_de <- readRDS(de_file)

  # Check top-level structure
  if (!is.list(full_de)) {
    stop("full_DE_results is not a list")
  }

  if (!"iSCORE_PD_MAST" %in% names(full_de)) {
    stop("Missing iSCORE_PD_MAST element")
  }
  cat("  ✓ Has iSCORE_PD_MAST element\n")

  if (!"CRISPRi_Mixscale" %in% names(full_de)) {
    stop("Missing CRISPRi_Mixscale element")
  }
  cat("  ✓ Has CRISPRi_Mixscale element\n")

  # Check perturbations exist
  if (length(full_de$CRISPRi_Mixscale) == 0) {
    stop("No perturbations found in CRISPRi_Mixscale")
  }
  cat("  ✓ Found", length(full_de$CRISPRi_Mixscale), "perturbations\n")

  # Check first perturbation structure
  first_pert <- full_de$CRISPRi_Mixscale[[1]]
  if (length(first_pert) == 0) {
    stop("First perturbation has no clusters")
  }
  cat("  ✓ First perturbation has", length(first_pert), "clusters\n")

  # Check first cluster structure
  first_cluster <- first_pert[[1]]
  if (!"results" %in% names(first_cluster)) {
    stop("First cluster missing 'results' element")
  }
  cat("  ✓ First cluster has 'results' element\n")

  # Check results dataframe
  results <- first_cluster$results
  if (!is.data.frame(results)) {
    stop("results is not a data.frame")
  }
  cat("  ✓ results is a data.frame\n")

  # Check rownames
  if (is.null(rownames(results)) || all(rownames(results) == as.character(1:nrow(results)))) {
    stop("results missing gene symbol rownames")
  }
  cat("  ✓ results has gene symbol rownames\n")

  # Check required columns
  required_cols <- c("avg_log2FC", "p_val_adj")
  missing_cols <- setdiff(required_cols, colnames(results))
  if (length(missing_cols) > 0) {
    stop("results missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  cat("  ✓ results has required columns: avg_log2FC, p_val_adj\n")

  cat("\n=== VALIDATION PASSED ===\n")
  cat("Dataset ready for use with launch_app()\n\n")

  invisible(TRUE)
}

# ============================================================================
# MAIN EXECUTION: Convert FPD with BH correction (TEST DATASET)
# ============================================================================

cat("Starting conversion for TEST DATASET: FPD with BH correction\n\n")

# Configuration
SOURCE_DIR <- "/mnt/e/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit"
OUTPUT_DIR <- "/mnt/e/THESIS/scRNASeq/mixscale/FPD_BH_dataset"
DATASET_NAME <- "FPD"
PVAL_COLUMN <- "p_weight_BH"
ENRICHMENT_SOURCE <- "/mnt/e/THESIS/scRNASeq/mixscale/enrichment_results_FPD_p_weight_BH"

# Convert
tryCatch({
  convert_pooled_to_full_de(
    source_dir = SOURCE_DIR,
    output_dir = OUTPUT_DIR,
    dataset_name = DATASET_NAME,
    pval_column = PVAL_COLUMN,
    enrichment_source = ENRICHMENT_SOURCE
  )

  # Validate
  validate_converted_dataset(OUTPUT_DIR)

  cat("\n✓ TEST DATASET CONVERSION SUCCESSFUL\n")
  cat("Ready to test with: launch_app(data_dir = '", OUTPUT_DIR, "')\n\n", sep = "")

}, error = function(e) {
  cat("\n✗ CONVERSION FAILED\n")
  cat("Error:", conditionMessage(e), "\n\n")
  stop(e)
})
