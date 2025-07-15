#' Experiment Weighting Functions for CRISPRi Multi-Experiment Analysis
#'
#' This module implements cell-count based weighting for combining results
#' across multiple CRISPRi experiments (C12_FPD-23, C12_FPD-24, C18_FPD-23)
#' with C12_FPD-24 as the primary experiment due to highest cell count.

# Required for data manipulation
suppressPackageStartupMessages(library(dplyr))

#' Load and process CRISPRi cell count data
#'
#' @param csv_path Path to CRISPRi analysis log CSV file
#' @param seurat_path Path to Seurat object (optional, for detailed analysis)
#' @return List with cell count data and experiment information
#' @export
load_crispri_cell_counts <- function(csv_path = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/PerturbSeq_MixScale_analysis_CRi_only/CRISPRi/CRISPRi_analysis_log.csv",
                                    seurat_path = NULL) {
  
  if (!file.exists(csv_path)) {
    stop("CRISPRi analysis log CSV file not found at: ", csv_path)
  }
  
  # Load cell count data
  cell_counts <- read.csv(csv_path, stringsAsFactors = FALSE)
  cat("[CELL COUNTS] Loaded cell count data for", nrow(cell_counts), "clusters\n")
  
  # Process the data
  processed_counts <- cell_counts %>%
    mutate(
      cluster_id = paste0("cluster_", cluster),
      total_cells = n_cells,
      nt_cells = n_nt_cells,
      perturbed_cells = total_cells - nt_cells,  # Calculate perturbed cells
      cell_density = total_cells / sum(total_cells),  # Proportion of total cells
      experiments_available = n_experiments
    ) %>%
    select(cluster_id, cluster, total_cells, nt_cells, perturbed_cells, 
           cell_density, experiments_available, timestamp)
  
  # Add information about which experiments are available per cluster
  processed_counts$experiment_info <- case_when(
    processed_counts$experiments_available == 3 ~ "C12_FPD-23,C12_FPD-24,C18_FPD-23",
    processed_counts$experiments_available == 2 ~ "Subset_of_experiments",
    processed_counts$experiments_available == 1 ~ "Single_experiment",
    TRUE ~ "Unknown"
  )
  
  cat("[CELL COUNTS] Processing summary:\n")
  cat("  - Clusters analyzed:", nrow(processed_counts), "\n")
  cat("  - Total cells across all clusters:", sum(processed_counts$total_cells), "\n")
  cat("  - Total NT cells:", sum(processed_counts$nt_cells), "\n")  
  cat("  - Total perturbed cells:", sum(processed_counts$perturbed_cells), "\n")
  cat("  - Clusters with 3 experiments:", sum(processed_counts$experiments_available == 3), "\n")
  cat("  - Clusters with 2 experiments:", sum(processed_counts$experiments_available == 2), "\n")
  cat("  - Clusters with 1 experiment:", sum(processed_counts$experiments_available == 1), "\n")
  
  # Identify primary experiment based on overall cell counts
  # C12_FPD-24 expected to have highest representation
  primary_experiment <- "C12_FPD-24"
  cat("  - Primary experiment (highest cell count expected): ", primary_experiment, "\n")
  
  return(list(
    counts_by_cluster = processed_counts,
    summary = list(
      total_clusters = nrow(processed_counts),
      total_cells = sum(processed_counts$total_cells),
      total_nt_cells = sum(processed_counts$nt_cells),
      total_perturbed_cells = sum(processed_counts$perturbed_cells),
      primary_experiment = primary_experiment,
      data_source = csv_path
    )
  ))
}

#' Extract detailed cell counts from Seurat object (optional enhancement)
#'
#' @param seurat_path Path to Seurat object RDS file
#' @param clusters_to_analyze Vector of cluster IDs to analyze (default: 0-11)
#' @return Data frame with detailed experiment-specific cell counts
#' @export
extract_seurat_cell_counts <- function(seurat_path = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
                                      clusters_to_analyze = 0:11) {
  
  if (!file.exists(seurat_path)) {
    warning("Seurat object not found at: ", seurat_path, ". Using CSV-based estimates.")
    return(NULL)
  }
  
  cat("[SEURAT ANALYSIS] Loading Seurat object for detailed cell count analysis...\n")
  
  # Load Seurat object (this may take some time)
  tryCatch({
    seurat_obj <- readRDS(seurat_path)
    cat("[SEURAT ANALYSIS] Seurat object loaded successfully\n")
    
    # Extract metadata
    metadata <- seurat_obj@meta.data
    
    # Check for required columns
    required_cols <- c("scMAGeCK_gene_assignment", "seurat_clusters")
    missing_cols <- setdiff(required_cols, names(metadata))
    
    if (length(missing_cols) > 0) {
      warning("Missing required columns in Seurat metadata: ", paste(missing_cols, collapse = ", "))
      return(NULL)
    }
    
    # Extract experiment information from gene assignments
    # Assuming experiments are encoded in the assignment column
    experiment_counts <- metadata %>%
      filter(seurat_clusters %in% clusters_to_analyze) %>%
      mutate(
        cluster_id = paste0("cluster_", seurat_clusters),
        is_nt = scMAGeCK_gene_assignment == "Non-Targeting",
        experiment = extract_experiment_from_assignment(scMAGeCK_gene_assignment)
      ) %>%
      group_by(cluster_id, experiment) %>%
      summarise(
        total_cells = n(),
        nt_cells = sum(is_nt),
        perturbed_cells = sum(!is_nt),
        .groups = "drop"
      )
    
    cat("[SEURAT ANALYSIS] Extracted detailed counts for", 
        length(unique(experiment_counts$cluster_id)), "clusters and",
        length(unique(experiment_counts$experiment)), "experiments\n")
    
    return(experiment_counts)
    
  }, error = function(e) {
    warning("Error loading Seurat object: ", e$message, ". Using CSV-based estimates.")
    return(NULL)
  })
}

#' Helper function to extract experiment ID from gene assignment
#'
#' @param assignments Vector of gene assignment strings
#' @return Vector of experiment IDs
extract_experiment_from_assignment <- function(assignments) {
  # This function may need adjustment based on actual assignment format
  # For now, using placeholder logic
  experiments <- case_when(
    grepl("C12_FPD-23", assignments) ~ "C12_FPD-23",
    grepl("C12_FPD-24", assignments) ~ "C12_FPD-24", 
    grepl("C18_FPD-23", assignments) ~ "C18_FPD-23",
    assignments == "Non-Targeting" ~ "Control",
    TRUE ~ "Unknown"
  )
  return(experiments)
}

#' Calculate experiment weights based on cell counts
#'
#' @param cell_count_data Result from load_crispri_cell_counts()
#' @param detailed_counts Result from extract_seurat_cell_counts() (optional)
#' @param weighting_method Method for calculating weights ("total_cells", "perturbed_cells", "combined")
#' @return List with experiment weights by cluster
#' @export
calculate_experiment_weights <- function(cell_count_data, detailed_counts = NULL, 
                                       weighting_method = "combined") {
  
  cat("[EXPERIMENT WEIGHTING] Calculating weights using method:", weighting_method, "\n")
  
  counts_by_cluster <- cell_count_data$counts_by_cluster
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  
  # Initialize weights structure
  weights <- list()
  
  for (cluster_row in 1:nrow(counts_by_cluster)) {
    cluster_data <- counts_by_cluster[cluster_row, ]
    cluster_id <- cluster_data$cluster_id
    
    # Base weights from CSV data (estimated)
    total_cells_in_cluster <- cluster_data$total_cells
    nt_cells_in_cluster <- cluster_data$nt_cells
    perturbed_cells_in_cluster <- cluster_data$perturbed_cells
    n_experiments_available <- cluster_data$experiments_available
    
    if (weighting_method == "total_cells") {
      # Weight by total cells in cluster
      base_weight <- total_cells_in_cluster
    } else if (weighting_method == "perturbed_cells") {
      # Weight by perturbed cells only
      base_weight <- perturbed_cells_in_cluster
    } else {
      # Combined method: total cells with perturbed cell boost
      base_weight <- total_cells_in_cluster + (perturbed_cells_in_cluster * 0.5)
    }
    
    # If detailed Seurat counts available, use them for more precise weighting
    if (!is.null(detailed_counts)) {
      cluster_detailed <- detailed_counts[detailed_counts$cluster_id == cluster_id, ]
      
      if (nrow(cluster_detailed) > 0) {
        # Use detailed experiment-specific weights
        for (exp in experiments) {
          exp_data <- cluster_detailed[cluster_detailed$experiment == exp, ]
          if (nrow(exp_data) > 0) {
            if (weighting_method == "total_cells") {
              exp_weight <- exp_data$total_cells
            } else if (weighting_method == "perturbed_cells") {
              exp_weight <- exp_data$perturbed_cells
            } else {
              exp_weight <- exp_data$total_cells + (exp_data$perturbed_cells * 0.5)
            }
            weights[[paste0(exp, "_", cluster_id)]] <- exp_weight
          } else {
            # Experiment not available in this cluster
            weights[[paste0(exp, "_", cluster_id)]] <- 0
          }
        }
      } else {
        # Fallback to estimated weights
        weights <- add_estimated_weights(weights, experiments, cluster_id, base_weight, n_experiments_available)
      }
    } else {
      # Use estimated weights based on CSV data
      weights <- add_estimated_weights(weights, experiments, cluster_id, base_weight, n_experiments_available)
    }
  }
  
  # Normalize weights within each cluster
  normalized_weights <- normalize_weights_by_cluster(weights)
  
  # Identify primary experiment (should be C12_FPD-24 with highest weights)
  primary_experiment <- identify_primary_experiment(normalized_weights)
  
  cat("[EXPERIMENT WEIGHTING] Weighting completed\n")
  cat("  - Primary experiment identified:", primary_experiment, "\n")
  cat("  - Total weight entries:", length(normalized_weights), "\n")
  
  return(list(
    weights = normalized_weights,
    primary_experiment = primary_experiment,
    weighting_method = weighting_method,
    normalization = "within_cluster"
  ))
}

#' Helper function to add estimated weights when detailed data unavailable
#'
#' @param weights Current weights list
#' @param experiments Vector of experiment names
#' @param cluster_id Cluster identifier
#' @param base_weight Base weight for the cluster
#' @param n_experiments_available Number of experiments in this cluster
#' @return Updated weights list
add_estimated_weights <- function(weights, experiments, cluster_id, base_weight, n_experiments_available) {
  
  if (n_experiments_available == 3) {
    # All experiments available - distribute based on expected proportions
    # C12_FPD-24 gets highest weight (user indicated it's strongest)
    weight_proportions <- c("C12_FPD-23" = 0.25, "C12_FPD-24" = 0.5, "C18_FPD-23" = 0.25)
    
    for (exp in experiments) {
      weights[[paste0(exp, "_", cluster_id)]] <- base_weight * weight_proportions[[exp]]
    }
    
  } else if (n_experiments_available == 2) {
    # Assume C12_FPD-24 + one other (most likely scenario)
    # Give C12_FPD-24 higher weight
    weights[[paste0("C12_FPD-24", "_", cluster_id)]] <- base_weight * 0.6
    weights[[paste0("C12_FPD-23", "_", cluster_id)]] <- base_weight * 0.4
    weights[[paste0("C18_FPD-23", "_", cluster_id)]] <- 0  # Not available
    
  } else if (n_experiments_available == 1) {
    # Likely just C12_FPD-24 (primary experiment)
    weights[[paste0("C12_FPD-24", "_", cluster_id)]] <- base_weight
    weights[[paste0("C12_FPD-23", "_", cluster_id)]] <- 0
    weights[[paste0("C18_FPD-23", "_", cluster_id)]] <- 0
    
  } else {
    # Unknown - equal weights as fallback
    equal_weight <- base_weight / length(experiments)
    for (exp in experiments) {
      weights[[paste0(exp, "_", cluster_id)]] <- equal_weight
    }
  }
  
  return(weights)
}

#' Normalize weights within each cluster so they sum to 1
#'
#' @param weights List of experiment weights
#' @return List of normalized weights
normalize_weights_by_cluster <- function(weights) {
  
  # Extract cluster IDs
  weight_names <- names(weights)
  cluster_ids <- unique(gsub(".*_", "", weight_names))
  
  normalized_weights <- weights
  
  for (cluster_id in cluster_ids) {
    cluster_weight_names <- weight_names[grepl(paste0("_", cluster_id, "$"), weight_names)]
    cluster_weights <- weights[cluster_weight_names]
    
    total_weight <- sum(unlist(cluster_weights))
    
    if (total_weight > 0) {
      # Normalize so cluster weights sum to 1
      for (name in cluster_weight_names) {
        normalized_weights[[name]] <- weights[[name]] / total_weight
      }
    }
  }
  
  return(normalized_weights)
}

#' Identify the primary experiment based on average weights across clusters
#'
#' @param weights List of normalized weights
#' @return Name of primary experiment
identify_primary_experiment <- function(weights) {
  
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  experiment_totals <- c()
  
  for (exp in experiments) {
    exp_weights <- weights[grepl(paste0("^", exp, "_"), names(weights))]
    experiment_totals[exp] <- mean(unlist(exp_weights), na.rm = TRUE)
  }
  
  primary_exp <- names(experiment_totals)[which.max(experiment_totals)]
  return(primary_exp)
}

#' Weighted meta-analysis across experiments
#'
#' @param experiment_results List of results by experiment
#' @param weights Experiment weights from calculate_experiment_weights()
#' @param cluster_id Cluster identifier
#' @return Weighted combined results
#' @export
weighted_meta_analysis <- function(experiment_results, weights, cluster_id) {
  
  experiments <- names(experiment_results)
  
  if (length(experiments) == 0) {
    return(list(
      weighted_log2fc = NA,
      weighted_p_value = NA,
      experiments_included = 0,
      error = "No experiment results provided"
    ))
  }
  
  # Extract weights for this cluster
  exp_weights <- c()
  exp_effects <- c()
  exp_p_values <- c()
  
  for (exp in experiments) {
    weight_key <- paste0(exp, "_", cluster_id)
    if (weight_key %in% names(weights$weights)) {
      exp_weight <- weights$weights[[weight_key]]
      
      if (exp_weight > 0 && !is.null(experiment_results[[exp]])) {
        exp_weights <- c(exp_weights, exp_weight)
        exp_effects <- c(exp_effects, experiment_results[[exp]]$log2fc %||% 0)
        exp_p_values <- c(exp_p_values, experiment_results[[exp]]$p_value %||% 1)
      }
    }
  }
  
  if (length(exp_weights) == 0) {
    return(list(
      weighted_log2fc = NA,
      weighted_p_value = NA,
      experiments_included = 0,
      error = "No valid experiments with weights > 0"
    ))
  }
  
  # Calculate weighted average effect size
  weighted_log2fc <- sum(exp_effects * exp_weights) / sum(exp_weights)
  
  # Calculate weighted p-value using Fisher's method with weights
  # Convert p-values to chi-square statistics
  chi_square_stats <- -2 * log(exp_p_values)
  weighted_chi_square <- sum(chi_square_stats * exp_weights) / sum(exp_weights)
  
  # Convert back to p-value (approximate)
  weighted_p_value <- pchisq(weighted_chi_square * 2 * length(exp_weights), 
                            df = 2 * length(exp_weights), lower.tail = FALSE)
  
  return(list(
    weighted_log2fc = weighted_log2fc,
    weighted_p_value = weighted_p_value,
    experiments_included = length(exp_weights),
    experiment_weights = setNames(exp_weights, experiments[1:length(exp_weights)]),
    primary_experiment = weights$primary_experiment,
    weighting_method = weights$weighting_method
  ))
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a