# Functions for importing MAST and MixScale data with consistent structure

#' Extract cluster ID from file path
#'
#' @param file_path The full path to a results file
#' @return String: the extracted cluster ID
extract_cluster_id <- function(file_path) {
  # Extract from patterns like "clust_0" in the filename
  filename <- basename(file_path)
  if (grepl("clust_([0-9]+)", filename)) {
    cluster_num <- regmatches(filename, regexpr("clust_([0-9]+)", filename))
    return(paste0("cluster_", gsub("clust_", "", cluster_num)))
  }
  
  # Extract from patterns like "Cluster0" in directory path
  if (grepl("Cluster([0-9]+)", file_path)) {
    return(regmatches(file_path, regexpr("Cluster[0-9]+", file_path)))
  }
  
  # Fallback: use a numbered default
  return(paste0("cluster_unknown"))
}

#' Import MAST differential expression results with consistent structure
#'
#' @param input_dir Directory containing MAST result files
#' @return List of structured MAST results
import_mast_data <- function(input_dir) {
  # Get list of RDS files
  files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  
  # Initialize results list
  results <- list()
  
  # Process each file
  for (file_path in files) {
    # Extract mutation name from filename - updated pattern to handle both old and new naming
    filename <- basename(file_path)
    
    # Try new pattern first (mutation_NAME_results.rds)
    if (grepl("^mutation_(.+)_results\\.rds$", filename)) {
      mutation_name <- gsub("^mutation_(.+)_results\\.rds$", "\\1", filename)
    } 
    # Fall back to old pattern (mutation_NAME_results_RNA_batchspecific.rds)
    else if (grepl("^mutation_(.+)_results_RNA_batchspecific.*\\.rds$", filename)) {
      mutation_name <- gsub("^mutation_(.+)_results_RNA_batchspecific.*\\.rds$", "\\1", filename)
    }
    else {
      # Skip files that don't match expected patterns
      next
    }
    
    # Load RDS file
    data <- readRDS(file_path)
    
    # Skip if no metadata or empty results
    if (is.null(data$metadata)) next
    
    # Initialize mutation entry if needed
    if (is.null(results[[mutation_name]])) {
      results[[mutation_name]] <- list()
    }
    
    # Process clusters
    for (cluster_name in names(data)) {
      if (cluster_name == "metadata") next
      
      # Skip if no results for this cluster or results are error messages
      if (!is.data.frame(data[[cluster_name]]) || 
          (nrow(data[[cluster_name]]) == 1 && "error" %in% colnames(data[[cluster_name]]))) {
        next
      }
      
      # Get background genes - all genes tested in this analysis (rownames of the DEG table)
      background_genes <- rownames(data[[cluster_name]])
      
      # Create standardized structure with background_genes at same level as results and metadata
      results[[mutation_name]][[cluster_name]] <- list(
        results = data[[cluster_name]],
        metadata = list(
          mutation = mutation_name,
          cluster = cluster_name,
          test = "MAST",
          control = data$metadata$control,
          batches_used = data$metadata$batches_used,
          latent_vars = data$metadata$latent_vars,
          date = data$metadata$date
        ),
        background_genes = background_genes
      )
    }
    
    # Add global metadata for this mutation
    results[[mutation_name]]$metadata <- data$metadata
  }
  
  return(results)
}

#' Import MixScale differential expression results with consistent structure
#'
#' @param input_dir Directory containing MixScale result files
#' @param modality Optional: specify "CRISPRi" or "CRISPRa" to import only one modality
#' @return List of structured MixScale results
import_mixscale_data <- function(input_dir, modality = NULL) {
  # Find all DEG results files
  rds_files <- list.files(input_dir, pattern = "\\DEGs.rds$", full.names = TRUE, recursive = TRUE)
  if (length(rds_files) == 0) stop("No RDS files found in ", input_dir)
  
  # Filter by modality if specified
  if (!is.null(modality)) {
    if (!modality %in% c("CRISPRi", "CRISPRa")) {
      stop("modality must be 'CRISPRi' or 'CRISPRa'")
    }
    # Keep only files from the specified modality subdirectory
    rds_files <- rds_files[grepl(paste0("/", modality, "/"), rds_files)]
    if (length(rds_files) == 0) {
      stop(sprintf("No RDS files found for modality '%s' in %s", modality, input_dir))
    }
  }
  
  # Initialize results list
  results <- list()
  
  # Process each file
  for (file_path in rds_files) {
    # Load data safely
    data <- tryCatch(readRDS(file_path), error = function(e) NULL)
    if (is.null(data)) next
    
    # Extract cluster ID using the extract_cluster_id helper function
    cluster_id <- extract_cluster_id(file_path)
    
    # Extract modality from path if not specified
    file_modality <- if (grepl("/CRISPRi/", file_path)) "CRISPRi" else if (grepl("/CRISPRa/", file_path)) "CRISPRa" else NA
    
    # Process each perturbation
    for (pert in names(data)) {
      if (pert == "metadata") next
      
      # Get perturbation data and structure it properly
      pert_data <- data[[pert]]
      
      # Structure results consistently
      if (is.list(pert_data) && !is.data.frame(pert_data) && !is.null(pert_data$results)) {
        pert_results <- pert_data$results
        pert_metadata <- pert_data$metadata
      } else if (is.data.frame(pert_data)) {
        pert_results <- pert_data
        pert_metadata <- NULL
      } else {
        next  # Skip if invalid data
      }
      
      # Skip if no valid results
      if (!is.data.frame(pert_results) || nrow(pert_results) == 0) next
      
      # Initialize perturbation entry if needed
      if (is.null(results[[pert]])) results[[pert]] <- list()
      
      # Extract experiment information
      experiment_ids <- character(0)
      weighted_experiments <- character(0)
      
      # Get gene column name - try common patterns
      gene_col <- intersect(c("gene_ID", "gene_id", "geneid", "gene"), colnames(pert_results))[1]
      if (is.na(gene_col) && !is.null(rownames(pert_results)) && 
          !all(grepl("^[0-9]+$", rownames(pert_results)))) {
        # If no obvious gene column but rownames look like gene IDs, use those
        has_gene_rownames <- TRUE
      } else {
        has_gene_rownames <- FALSE
      }
      
      # Try to get experiment IDs from column names
      log2fc_cols <- grep("^log2FC_", colnames(pert_results), value = TRUE)
      p_value_cols <- grep("^p_cell_type", colnames(pert_results), value = TRUE)
      p_value_weight_cols <- grep(":weight$", p_value_cols, value = TRUE)
      
      # Check for CRISPRa single-experiment structure (generic column names)
      if (length(log2fc_cols) == 0 && "log2FC" %in% colnames(pert_results)) {
        # CRISPRa structure detected - single experiment with generic columns
        # Extract experiment ID from the mixscale object filename if available
        exp_id <- "A15_FPD-24"  # Default for CRISPRa based on observed data
        
        # Try to extract from file path pattern like "A15_FPD-24"
        if (grepl("([A-Z][0-9]+_[A-Z]+-[0-9]+)", file_path)) {
          extracted_id <- regmatches(file_path, regexpr("([A-Z][0-9]+_[A-Z]+-[0-9]+)", file_path))
          if (length(extracted_id) > 0) exp_id <- extracted_id[1]
        }
        
        experiment_ids <- exp_id
        log2fc_cols <- "log2FC"
        
        # Map generic p-value columns to experiment-specific format for consistency
        if ("p_weight" %in% colnames(pert_results)) {
          p_value_cols <- "p_weight"
          p_value_weight_cols <- "p_weight"
          weighted_experiments <- exp_id
        }
        
        # Add mapped columns to maintain consistency with CRISPRi structure
        pert_results[[paste0("log2FC_", exp_id)]] <- pert_results$log2FC
        if ("p_weight" %in% colnames(pert_results)) {
          pert_results[[paste0("p_cell_type", exp_id, ":weight")]] <- pert_results$p_weight
        }
        
      } else if (length(log2fc_cols) > 0) {
        # Original CRISPRi multi-experiment structure
        experiment_ids <- gsub("^log2FC_", "", log2fc_cols)
        
        # Extract experiment IDs from p-value column names for better matching
        if (length(p_value_weight_cols) > 0) {
          experiment_ids_from_p <- unique(gsub("^p_cell_type([^:]+).*$", "\\1", p_value_weight_cols))
          # Only keep experiment IDs that have both log2FC and p-value columns
          experiment_ids <- intersect(experiment_ids, experiment_ids_from_p)
        }
        
        # Check for weighted p-values
        for (exp_id in experiment_ids) {
          if (any(grepl(paste0("p_cell_type", exp_id, ".*:weight"), colnames(pert_results)))) {
            weighted_experiments <- c(weighted_experiments, exp_id)
          }
        }
      } else if ("p_weight" %in% colnames(pert_results)) {
        experiment_ids <- paste0("exp_", cluster_id)
      }
      
      # Create experiment-specific background genes map
      experiment_background_genes <- list()
      
      # Global background - all genes in the dataset
      if (!is.na(gene_col)) {
        global_background_genes <- unique(pert_results[[gene_col]])
      } else if (has_gene_rownames) {
        global_background_genes <- rownames(pert_results)
      } else {
        global_background_genes <- character(0)
      }
      
      # Generate experiment-specific background gene lists (non-NA genes in each experiment)
      for (exp_id in experiment_ids) {
        lfc_col <- paste0("log2FC_", exp_id)
        
        # For CRISPRa, also check the generic column if mapped column doesn't exist
        if (!lfc_col %in% colnames(pert_results)) {
          if ("log2FC" %in% colnames(pert_results) && exp_id %in% experiment_ids) {
            lfc_col <- "log2FC"  # Use generic column for CRISPRa
          } else {
            # Skip if no matching log2FC column
            next
          }
        }
        
        # Background genes for this experiment are those without NA in log2FC column
        if (!is.na(gene_col)) {
          exp_bg_genes <- pert_results[[gene_col]][!is.na(pert_results[[lfc_col]])]
        } else if (has_gene_rownames) {
          exp_bg_genes <- rownames(pert_results)[!is.na(pert_results[[lfc_col]])]
        } else {
          exp_bg_genes <- character(0)
        }
        
        # Store experiment-specific background genes
        experiment_background_genes[[exp_id]] <- exp_bg_genes
      }
      
      # Create metadata
      cluster_metadata <- list(
        gene = pert,
        cluster = cluster_id,
        modality = file_modality,  # Add modality info
        experiments = experiment_ids,
        is_weighted = length(weighted_experiments) > 0,
        weighted_experiments = weighted_experiments,
        file_path = file_path
      )
      
      # Override with original metadata if available
      if (!is.null(pert_metadata)) {
        cluster_metadata <- pert_metadata
        cluster_metadata$gene <- pert
        cluster_metadata$cluster <- cluster_id
        cluster_metadata$file_path <- file_path
        if (is.null(cluster_metadata$modality) && !is.na(file_modality)) {
          cluster_metadata$modality <- file_modality
        }
        if (is.null(cluster_metadata$experiments)) {
          cluster_metadata$experiments <- experiment_ids
        }
        if (is.null(cluster_metadata$is_weighted)) {
          cluster_metadata$is_weighted <- length(weighted_experiments) > 0
        }
      }
      
      # Update column lists to include mapped columns for CRISPRa
      if ("log2FC" %in% colnames(pert_results) && length(log2fc_cols) == 1 && log2fc_cols[1] == "log2FC") {
        # For CRISPRa, include the mapped columns we created
        log2fc_cols <- c(log2fc_cols, grep("^log2FC_", colnames(pert_results), value = TRUE))
        p_value_cols <- c(p_value_cols, grep("^p_cell_type.*:weight$", colnames(pert_results), value = TRUE))
        p_value_weight_cols <- grep(":weight$", p_value_cols, value = TRUE)
      }
      
      # Store results with enhanced background gene tracking
      results[[pert]][[cluster_id]] <- list(
        results = pert_results,
        metadata = cluster_metadata,
        background_genes = global_background_genes,
        experiment_background_genes = experiment_background_genes,
        gene_column = gene_col,
        has_gene_rownames = has_gene_rownames,
        log2fc_columns = log2fc_cols,
        p_value_columns = p_value_cols,
        weighted_p_value_columns = p_value_weight_cols
      )
    }
  }
  
  return(results)
}


# The following code is now handled in the vignette properly with absolute paths
# and is commented out to prevent accidental execution when sourcing this file

# Example of how to use the import functions with new directory structures:
# full <- list()
# full$iSCORE_PD_MAST <- import_mast_data("/path/to/iSCORE-PD_MAST_analysis/")
# full$CRISPRi_Mixscale <- import_mixscale_data("/path/to/PerturbSeq_MixScale_analysis_full_dataset/", modality = "CRISPRi")
# full$CRISPRa_Mixscale <- import_mixscale_data("/path/to/PerturbSeq_MixScale_analysis_full_dataset/", modality = "CRISPRa")
# # Or import both modalities at once:
# full$PerturbSeq_Mixscale <- import_mixscale_data("/path/to/PerturbSeq_MixScale_analysis_full_dataset/")
# saveRDS(full, "/path/to/output/full_DE_results.rds")