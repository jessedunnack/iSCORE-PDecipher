#!/usr/bin/env Rscript
# Consolidate Pooled MixScale Enrichment Results
# Creates all_enrichment_padj005_complete_with_direction.rds files
# for FPD and CRISPRi datasets with different p-value corrections
#
# Created: October 24, 2025

library(dplyr)

#' Consolidate enrichment results from directory structure into single RDS
#'
#' @param enrichment_dir Directory containing cluster subdirectories with enrichment results
#' @param output_file Path to output RDS file
#' @param padj_threshold P-value threshold for filtering (default 0.05)
#' @return Data frame with all consolidated enrichment results
consolidate_enrichment_results <- function(
  enrichment_dir,
  output_file,
  padj_threshold = 0.05
) {

  message("\n===============================================")
  message("CONSOLIDATING ENRICHMENT RESULTS")
  message("===============================================\n")
  message("Input directory: ", enrichment_dir)
  message("Output file: ", output_file)
  message("P-value threshold: ", padj_threshold)

  # Find all cluster directories
  cluster_dirs <- list.dirs(enrichment_dir, recursive = FALSE, full.names = TRUE)
  cluster_dirs <- cluster_dirs[!grepl("parallel_joblog", cluster_dirs)]

  message(sprintf("\nFound %d cluster directories", length(cluster_dirs)))

  all_results <- list()
  result_count <- 0

  for (cluster_dir in cluster_dirs) {
    cluster_name <- basename(cluster_dir)
    # Extract cluster number (e.g., "cluster_0" from "all_FPD_no_multiplets_noExptSplit_clust_0")
    cluster_num <- gsub(".*_clust_", "cluster_", cluster_name)

    message(sprintf("\nProcessing %s -> %s", cluster_name, cluster_num))

    # Find all perturbation directories
    pert_dirs <- list.dirs(cluster_dir, recursive = FALSE, full.names = TRUE)

    message(sprintf("  Found %d perturbations", length(pert_dirs)))

    for (pert_dir in pert_dirs) {
      pert_name <- basename(pert_dir)

      # Find all enrichment method directories
      method_dirs <- list.dirs(pert_dir, recursive = FALSE, full.names = TRUE)
      method_dirs <- method_dirs[!grepl("diagnostics", method_dirs)]

      for (method_dir in method_dirs) {
        method_name <- basename(method_dir)

        # Find all RDS files (ALL, UP, DOWN)
        rds_files <- list.files(method_dir, pattern = "\\.rds$", full.names = TRUE)

        for (rds_file in rds_files) {
          # Extract direction from filename
          direction <- gsub(paste0(method_name, "_"), "", gsub("\\.rds", "", basename(rds_file)))

          # Load RDS file
          enrich_obj <- tryCatch({
            readRDS(rds_file)
          }, error = function(e) {
            message("    WARNING: Could not load ", basename(rds_file), ": ", e$message)
            return(NULL)
          })

          if (is.null(enrich_obj)) next

          # Extract results dataframe
          results_df <- tryCatch({
            if (inherits(enrich_obj, "enrichResult")) {
              enrich_obj@result
            } else if (is.list(enrich_obj) && "result" %in% names(enrich_obj)) {
              enrich_obj$result@result
            } else if (is.data.frame(enrich_obj)) {
              enrich_obj
            } else {
              NULL
            }
          }, error = function(e) {
            message("    WARNING: Could not extract results from ", basename(rds_file))
            return(NULL)
          })

          if (is.null(results_df) || nrow(results_df) == 0) next

          # Filter by p-value threshold
          if ("p.adjust" %in% colnames(results_df)) {
            results_df <- results_df[results_df$p.adjust <= padj_threshold, ]
          } else if ("padj" %in% colnames(results_df)) {
            results_df <- results_df[results_df$padj <= padj_threshold, ]
          } else if ("fdr" %in% colnames(results_df)) {
            results_df <- results_df[results_df$fdr <= padj_threshold, ]
          }

          if (nrow(results_df) == 0) next

          # Add metadata columns
          results_df$mutation_perturbation <- pert_name
          results_df$cluster <- cluster_num
          results_df$enrichment_type <- method_name
          results_df$direction <- direction

          # Standardize column names
          if ("ID" %in% colnames(results_df)) {
            results_df$term_id <- results_df$ID
          } else if ("pathway" %in% colnames(results_df)) {
            results_df$term_id <- results_df$pathway
          } else if ("term" %in% colnames(results_df)) {
            results_df$term_id <- results_df$term
          }

          if ("Description" %in% colnames(results_df)) {
            results_df$term_description <- results_df$Description
          } else if ("description" %in% colnames(results_df)) {
            results_df$term_description <- results_df$description
          } else if ("pathway" %in% colnames(results_df)) {
            results_df$term_description <- results_df$pathway
          }

          # Standardize p.adjust column
          if ("padj" %in% colnames(results_df) && !"p.adjust" %in% colnames(results_df)) {
            results_df$p.adjust <- results_df$padj
          } else if ("fdr" %in% colnames(results_df) && !"p.adjust" %in% colnames(results_df)) {
            results_df$p.adjust <- results_df$fdr
          }

          # Add other useful columns if they exist
          if (!"GeneRatio" %in% colnames(results_df) && "gene_ratio" %in% colnames(results_df)) {
            results_df$GeneRatio <- results_df$gene_ratio
          }

          if (!"Count" %in% colnames(results_df) && "count" %in% colnames(results_df)) {
            results_df$Count <- results_df$count
          }

          if (!"geneID" %in% colnames(results_df) && "genes" %in% colnames(results_df)) {
            results_df$geneID <- results_df$genes
          }

          # Store in list
          result_count <- result_count + 1
          all_results[[result_count]] <- results_df
        }
      }
    }
  }

  message(sprintf("\n\nCollected %d result dataframes", length(all_results)))

  if (length(all_results) == 0) {
    stop("ERROR: No enrichment results found!")
  }

  # Combine all results
  message("Combining all results...")
  combined_results <- bind_rows(all_results)

  message(sprintf("Total enrichment terms: %d", nrow(combined_results)))
  message(sprintf("Unique perturbations: %d", length(unique(combined_results$mutation_perturbation))))
  message(sprintf("Unique clusters: %d", length(unique(combined_results$cluster))))
  message(sprintf("Enrichment types: %s", paste(unique(combined_results$enrichment_type), collapse=", ")))
  message(sprintf("Directions: %s", paste(unique(combined_results$direction), collapse=", ")))

  # Save to RDS
  message("\nSaving to: ", output_file)
  saveRDS(combined_results, output_file)

  message("\n✓ CONSOLIDATION COMPLETE!\n")

  return(combined_results)
}

# ============================================================================
# MAIN EXECUTION
# ============================================================================

# Base directory
base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper"

# Define all enrichment directories and their output files
enrichment_configs <- list(
  list(
    dir = file.path(base_dir, "enrichment_results_FPD_p_weight"),
    output = file.path(base_dir, "enrichment_results_FPD_p_weight", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "FPD (uncorrected p-values)"
  ),
  list(
    dir = file.path(base_dir, "enrichment_results_FPD_p_weight_BH"),
    output = file.path(base_dir, "enrichment_results_FPD_p_weight_BH", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "FPD (BH corrected)"
  ),
  list(
    dir = file.path(base_dir, "enrichment_results_FPD_p_weight_bonferroni"),
    output = file.path(base_dir, "enrichment_results_FPD_p_weight_bonferroni", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "FPD (Bonferroni corrected)"
  ),
  list(
    dir = file.path(base_dir, "enrichment_results_CRISPRi_p_weight"),
    output = file.path(base_dir, "enrichment_results_CRISPRi_p_weight", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "CRISPRi (uncorrected p-values)"
  ),
  list(
    dir = file.path(base_dir, "enrichment_results_CRISPRi_p_weight_BH"),
    output = file.path(base_dir, "enrichment_results_CRISPRi_p_weight_BH", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "CRISPRi (BH corrected)"
  ),
  list(
    dir = file.path(base_dir, "enrichment_results_CRISPRi_p_weight_bonferroni"),
    output = file.path(base_dir, "enrichment_results_CRISPRi_p_weight_bonferroni", "all_enrichment_padj005_complete_with_direction.rds"),
    name = "CRISPRi (Bonferroni corrected)"
  )
)

# Process all configurations
message("\n")
message("========================================================================")
message("CONSOLIDATING ALL POOLED MIXSCALE ENRICHMENT RESULTS")
message("========================================================================\n")
message(sprintf("Total configurations to process: %d\n", length(enrichment_configs)))

for (i in seq_along(enrichment_configs)) {
  config <- enrichment_configs[[i]]

  message("\n")
  message("------------------------------------------------------------------------")
  message(sprintf("[%d/%d] Processing: %s", i, length(enrichment_configs), config$name))
  message("------------------------------------------------------------------------")

  if (!dir.exists(config$dir)) {
    message("WARNING: Directory does not exist, skipping: ", config$dir)
    next
  }

  tryCatch({
    results <- consolidate_enrichment_results(
      enrichment_dir = config$dir,
      output_file = config$output,
      padj_threshold = 0.05
    )

    message(sprintf("✓ SUCCESS: %s", config$name))

  }, error = function(e) {
    message(sprintf("✗ ERROR processing %s: %s", config$name, e$message))
  })
}

message("\n")
message("========================================================================")
message("ALL CONSOLIDATIONS COMPLETE")
message("========================================================================\n")
message("Output files created:")
for (config in enrichment_configs) {
  if (file.exists(config$output)) {
    file_size <- file.size(config$output) / 1024 / 1024  # MB
    message(sprintf("  ✓ %s (%.1f MB)", basename(dirname(config$output)), file_size))
  }
}
message("\n")
