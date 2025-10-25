# New functions for importing FDR-corrected pooled MixScale data
# Created: October 24, 2025
# Purpose: Support new Perturb-seq-only datasets with FDR corrections

#' Detect MixScale Data Format
#'
#' Determines if MixScale data uses experiment-split or pooled structure
#'
#' @param de_results Loaded MixScale results (list of perturbations)
#' @return Character: "experiment_split" or "pooled"
#' @export
detect_mixscale_format <- function(de_results) {
  if (!is.list(de_results) || length(de_results) == 0) {
    stop("de_results must be a non-empty list")
  }

  # Check first perturbation's column names
  first_pert <- de_results[[1]]

  if (!is.data.frame(first_pert)) {
    stop("First element of de_results should be a data.frame")
  }

  col_names <- names(first_pert)

  # Experiment-split has patterns like "log2FC_C12_FPD-24"
  if (any(grepl("log2FC_C\\d+_", col_names))) {
    return("experiment_split")
  }

  # Pooled has simple "log2FC" column with p_weight variants
  if ("log2FC" %in% col_names && "p_weight" %in% col_names) {
    # Additional check for FDR columns
    has_fdr <- any(c("p_weight_BH", "p_weight_bonferroni") %in% col_names)
    if (has_fdr) {
      return("pooled")
    } else {
      # Older pooled format without FDR
      return("pooled")
    }
  }

  stop("Unable to detect data format. Expected either experiment-split or pooled structure.")
}


#' Import Pooled MixScale Data with FDR Corrections
#'
#' Imports MixScale differential expression results from pooled analysis with FDR-corrected p-values.
#' Compatible with the new Perturb-seq-only datasets that have simple column structure.
#'
#' @param mixscale_dir Directory containing cluster subdirectories with *_mixscale_DEGs.rds files
#' @param pval_column Which p-value column to use: "p_weight" (uncorrected),
#'                    "p_weight_BH" (Benjamini-Hochberg FDR), or "p_weight_bonferroni"
#' @param dataset_type Optional: "FPD" or "CRISPRi" to help with metadata. Auto-detected if NULL.
#' @return List structure compatible with existing app modules, organized as:
#'         perturbation -> cluster -> list(results, metadata, background_genes, ...)
#' @export
#'
#' @examples
#' \dontrun{
#' # Load FPD data with BH correction (recommended)
#' fpd_data <- import_pooled_mixscale_data(
#'   "../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/",
#'   pval_column = "p_weight_BH"
#' )
#'
#' # Load CRISPRi data with original p-values
#' crispri_data <- import_pooled_mixscale_data(
#'   "../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/",
#'   pval_column = "p_weight"
#' )
#' }
import_pooled_mixscale_data <- function(
  mixscale_dir,
  pval_column = "p_weight_BH",
  dataset_type = NULL
) {

  # Validate pval_column
  valid_pval_cols <- c("p_weight", "p_weight_BH", "p_weight_bonferroni")
  if (!pval_column %in% valid_pval_cols) {
    stop("pval_column must be one of: ", paste(valid_pval_cols, collapse = ", "))
  }

  # Find all RDS files recursively
  rds_files <- list.files(
    mixscale_dir,
    pattern = "_mixscale_DEGs\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(rds_files) == 0) {
    stop("No *_mixscale_DEGs.rds files found in ", mixscale_dir)
  }

  message("Found ", length(rds_files), " RDS files to process")

  # Auto-detect dataset type from path if not specified
  if (is.null(dataset_type)) {
    if (grepl("FPD", mixscale_dir)) {
      dataset_type <- "FPD"
    } else if (grepl("CRISPRi", mixscale_dir)) {
      dataset_type <- "CRISPRi"
    } else {
      dataset_type <- "Unknown"
    }
  }

  # Initialize results list
  results <- list()

  # Process each file
  for (file_path in rds_files) {
    message("Processing: ", basename(file_path))

    # Load data safely
    data <- tryCatch(
      readRDS(file_path),
      error = function(e) {
        warning("Failed to load ", file_path, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(data)) next

    # Extract cluster ID from filename or path
    cluster_id <- extract_cluster_id(file_path)

    # Verify this is pooled format
    format_type <- tryCatch(
      detect_mixscale_format(data),
      error = function(e) {
        warning("Could not detect format for ", file_path, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(format_type)) next

    if (format_type != "pooled") {
      warning("File ", basename(file_path), " appears to be ", format_type,
              " format, not pooled. Skipping.")
      next
    }

    # Process each perturbation in this cluster
    for (pert_name in names(data)) {
      # Skip metadata if present
      if (pert_name == "metadata") next

      pert_data <- data[[pert_name]]

      # Validate it's a data frame
      if (!is.data.frame(pert_data)) {
        warning("Perturbation ", pert_name, " in cluster ", cluster_id,
                " is not a data.frame. Skipping.")
        next
      }

      # Check for required columns
      required_cols <- c("gene_ID", "log2FC", pval_column)
      missing_cols <- setdiff(required_cols, colnames(pert_data))

      if (length(missing_cols) > 0) {
        warning("Perturbation ", pert_name, " in cluster ", cluster_id,
                " is missing columns: ", paste(missing_cols, collapse = ", "), ". Skipping.")
        next
      }

      # Initialize perturbation entry if needed
      if (is.null(results[[pert_name]])) {
        results[[pert_name]] <- list()
      }

      # Get background genes (all genes tested)
      background_genes <- pert_data$gene_ID

      # Create metadata for this perturbation-cluster combination
      cluster_metadata <- list(
        gene = pert_name,
        cluster = cluster_id,
        dataset_type = dataset_type,
        modality = "Perturb-seq",
        is_pooled = TRUE,
        pval_column_used = pval_column,
        file_path = file_path,
        n_genes_tested = nrow(pert_data)
      )

      # Store results in structure compatible with existing app
      # This matches the structure from import_mixscale_data()
      results[[pert_name]][[cluster_id]] <- list(
        results = pert_data,
        metadata = cluster_metadata,
        background_genes = background_genes,
        gene_column = "gene_ID",
        has_gene_rownames = FALSE,
        log2fc_columns = "log2FC",  # Simple column name for pooled
        p_value_columns = pval_column,  # Selected p-value column
        # Store all available p-value columns for reference
        available_pval_columns = intersect(
          valid_pval_cols,
          colnames(pert_data)
        )
      )
    }
  }

  if (length(results) == 0) {
    stop("No valid perturbation data was imported")
  }

  message("Successfully imported ", length(results), " perturbations across ",
          length(unique(sapply(results, names))), " clusters")

  return(results)
}


#' Extract cluster ID from file path
#'
#' Helper function to extract cluster identifier from file paths or filenames
#'
#' @param file_path The full path to a results file
#' @return String: the extracted cluster ID (e.g., "cluster_0")
#' @export
extract_cluster_id <- function(file_path) {
  # Extract from patterns like "clust_0" in the filename
  filename <- basename(file_path)
  if (grepl("clust_([0-9]+)", filename)) {
    cluster_num <- regmatches(filename, regexpr("clust_([0-9]+)", filename))
    return(paste0("cluster_", gsub("clust_", "", cluster_num)))
  }

  # Extract from patterns like "Cluster0" in directory path
  if (grepl("Cluster([0-9]+)", file_path)) {
    cluster_match <- regmatches(file_path, regexpr("Cluster[0-9]+", file_path))
    cluster_num <- gsub("Cluster", "", cluster_match)
    return(paste0("cluster_", cluster_num))
  }

  # Fallback: use a numbered default
  warning("Could not extract cluster ID from ", file_path)
  return("cluster_unknown")
}


#' Import Enrichment Results from Specific P-Value Correction
#'
#' Loads enrichment analysis results that were generated using a specific p-value correction method.
#' The enrichment pipeline generates separate output directories for each correction type.
#'
#' @param base_dir Base directory containing enrichment results folders
#'                 (default: "/mnt/e/ASAP/scRNASeq/PerturbSeq/final")
#' @param dataset "FPD" or "CRISPRi" - which dataset's enrichment to load
#' @param pval_correction "none", "BH", or "bonferroni" - which p-value correction was used
#' @return List of enrichment results organized by perturbation and cluster
#' @export
#'
#' @examples
#' \dontrun{
#' # Load FPD enrichment with BH correction
#' fpd_enrich_bh <- import_enrichment_with_correction(
#'   dataset = "FPD",
#'   pval_correction = "BH"
#' )
#'
#' # Load CRISPRi enrichment with original p-values
#' crispri_enrich <- import_enrichment_with_correction(
#'   dataset = "CRISPRi",
#'   pval_correction = "none"
#' )
#' }
import_enrichment_with_correction <- function(
  base_dir = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper",
  dataset,
  pval_correction = "BH"
) {

  # Validate inputs
  if (!dataset %in% c("FPD", "CRISPRi")) {
    stop("dataset must be 'FPD' or 'CRISPRi'")
  }

  if (!pval_correction %in% c("none", "BH", "bonferroni")) {
    stop("pval_correction must be 'none', 'BH', or 'bonferroni'")
  }

  # Map correction method to directory suffix
  dir_suffix <- switch(pval_correction,
    "none" = "_p_weight",
    "BH" = "_p_weight_BH",
    "bonferroni" = "_p_weight_bonferroni"
  )

  # Construct path to enrichment directory
  enrich_dir <- file.path(
    base_dir,
    paste0("enrichment_results_", dataset, dir_suffix)
  )

  # Check if directory exists
  if (!dir.exists(enrich_dir)) {
    warning("Enrichment directory does not exist: ", enrich_dir)
    warning("Enrichment results may not have been generated yet. Check HPC job status.")
    return(list(
      status = "not_found",
      expected_path = enrich_dir,
      message = "Enrichment data not available. Run HPC enrichment pipeline first."
    ))
  }

  # NEW: Check for consolidated enrichment file first (MUCH FASTER!)
  consolidated_file <- file.path(enrich_dir, "all_enrichment_padj005_complete_with_direction.rds")

  if (file.exists(consolidated_file)) {
    message("Loading consolidated enrichment file: ", basename(consolidated_file))

    enrichment_data <- tryCatch({
      readRDS(consolidated_file)
    }, error = function(e) {
      warning("Failed to load consolidated file: ", e$message)
      warning("Falling back to individual file loading...")
      return(NULL)
    })

    if (!is.null(enrichment_data) && is.data.frame(enrichment_data)) {
      # Validate expected columns
      required_cols <- c("mutation_perturbation", "cluster", "enrichment_type", "direction", "p.adjust")
      missing_cols <- setdiff(required_cols, colnames(enrichment_data))

      if (length(missing_cols) == 0) {
        # Add metadata
        attr(enrichment_data, "dataset") <- dataset
        attr(enrichment_data, "pval_correction") <- pval_correction
        attr(enrichment_data, "import_date") <- Sys.Date()
        attr(enrichment_data, "source_file") <- consolidated_file
        attr(enrichment_data, "is_consolidated") <- TRUE

        message("Successfully loaded ", nrow(enrichment_data), " enrichment terms")
        message("  Perturbations: ", length(unique(enrichment_data$mutation_perturbation)))
        message("  Clusters: ", paste(sort(unique(enrichment_data$cluster)), collapse = ", "))
        message("  Enrichment types: ", paste(unique(enrichment_data$enrichment_type), collapse = ", "))

        return(enrichment_data)
      } else {
        warning("Consolidated file missing required columns: ", paste(missing_cols, collapse = ", "))
        warning("Falling back to individual file loading...")
      }
    }
  } else {
    message("Consolidated file not found: ", basename(consolidated_file))
    message("Falling back to individual file loading...")
  }

  # FALLBACK: Load individual enrichment files (backward compatibility)
  message("Loading individual enrichment files from: ", enrich_dir)

  # Find all enrichment RDS files
  enrich_files <- list.files(
    enrich_dir,
    pattern = "\\.rds$",
    full.names = TRUE,
    recursive = TRUE
  )

  if (length(enrich_files) == 0) {
    warning("No RDS files found in ", enrich_dir)
    return(list(
      status = "empty",
      directory = enrich_dir,
      message = "Directory exists but contains no RDS files"
    ))
  }

  message("Found ", length(enrich_files), " enrichment files")

  # Initialize results structure
  enrichment_results <- list()

  # Load each file
  for (file_path in enrich_files) {
    # Skip the consolidated file if present
    if (basename(file_path) == "all_enrichment_padj005_complete_with_direction.rds") {
      next
    }

    # Extract perturbation and cluster info from path
    # Expected structure: enrichment_results_*/cluster_X/perturbation_enrichment.rds

    path_parts <- strsplit(file_path, "/")[[1]]

    # Find cluster directory
    cluster_dir <- path_parts[grepl("cluster_", path_parts)]
    if (length(cluster_dir) == 0) {
      cluster_dir <- "unknown_cluster"
    } else {
      cluster_dir <- cluster_dir[1]
    }

    # Extract perturbation from filename
    filename <- basename(file_path)
    perturbation <- gsub("_enrichment\\.rds$", "", filename)

    # Load enrichment data
    enrich_data <- tryCatch(
      readRDS(file_path),
      error = function(e) {
        warning("Failed to load ", file_path, ": ", e$message)
        return(NULL)
      }
    )

    if (is.null(enrich_data)) next

    # Store in nested structure
    if (is.null(enrichment_results[[perturbation]])) {
      enrichment_results[[perturbation]] <- list()
    }

    enrichment_results[[perturbation]][[cluster_dir]] <- enrich_data
  }

  if (length(enrichment_results) == 0) {
    warning("No enrichment data was successfully loaded")
    return(list(
      status = "load_failed",
      directory = enrich_dir,
      message = "Files found but could not be loaded"
    ))
  }

  message("Successfully loaded enrichment for ", length(enrichment_results),
          " perturbations")

  # Add metadata about this import
  attr(enrichment_results, "dataset") <- dataset
  attr(enrichment_results, "pval_correction") <- pval_correction
  attr(enrichment_results, "import_date") <- Sys.Date()
  attr(enrichment_results, "source_directory") <- enrich_dir
  attr(enrichment_results, "is_consolidated") <- FALSE

  return(enrichment_results)
}
