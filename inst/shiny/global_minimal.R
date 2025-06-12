# Global configuration and functions for PerturbSeq Enrichment Analysis App
# This file contains core data loading and filtering functions

# Load required libraries
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(dplyr)
library(ggplot2)
library(plotly)
library(glue)

# Global configuration
options(shiny.maxRequestSize = 500*1024^2)  # 500MB upload limit

# Get data paths from environment variables
data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
de_file <- Sys.getenv("ISCORE_DE_FILE", "")
enrichment_file <- Sys.getenv("ISCORE_ENRICHMENT_FILE", "")
enrichment_dir <- Sys.getenv("ISCORE_ENRICHMENT_DIR", "")

# ============================================================================
# CORE DATA FUNCTIONS
# ============================================================================

#' Load consolidated enrichment data
#' @param dataset_name Name of dataset to load (or path to custom file)
#' @return data frame with enrichment results
load_enrichment_data <- function(dataset_name = "all_modalities") {
  
  cat("Loading dataset:", dataset_name, "\n")
  
  # First check if we have an environment variable set
  if (dataset_name == "all_modalities" && enrichment_file != "" && file.exists(enrichment_file)) {
    cat("Loading from environment variable path:", enrichment_file, "\n")
    dataset_file <- enrichment_file
  } else if (dataset_name %in% c("all_modalities", "mast_only", "mast_crispi", "crispi_only", "crispa_only")) {
    # Try datasets directory first
    dataset_file <- file.path("datasets", paste0(dataset_name, ".rds"))
    
    # Fallback to main directory for all_modalities
    if (!file.exists(dataset_file) && dataset_name == "all_modalities") {
      dataset_file <- "../all_enrichment_padj005_complete_with_direction.rds"
    }
  } else {
    # Custom file path
    dataset_file <- dataset_name
  }
  
  if (!file.exists(dataset_file)) {
    stop("Dataset file not found: ", dataset_file)
  }
  
  # Load the data
  data <- readRDS(dataset_file)
  cat("Loaded", nrow(data), "enrichment terms\n")
  
  # Ensure required columns exist
  required_cols <- c("mutation_perturbation", "cluster", "enrichment_type", "direction", "p.adjust", "Description")
  
  # Handle alternative column names
  if (!"mutation_perturbation" %in% names(data) && "gene" %in% names(data)) {
    data$mutation_perturbation <- data$gene
  }
  
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    warning("Missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Add modality column if missing (for backward compatibility)
  if (!"modality" %in% names(data) && "method" %in% names(data)) {
    data$modality <- case_when(
      data$method == "MixScale_CRISPRa" ~ "CRISPRa",
      data$method == "MixScale" ~ "CRISPRi", 
      data$method == "MAST" ~ "MAST",
      TRUE ~ NA_character_
    )
  }
  
  # Add analysis_type column if missing
  if (!"analysis_type" %in% names(data) && "method" %in% names(data)) {
    data$analysis_type <- case_when(
      grepl("MixScale", data$method) ~ "MixScale",
      TRUE ~ data$method
    )
  }
  
  # Filter to significant results only
  data <- data[!is.na(data$p.adjust) & data$p.adjust <= 0.05, ]
  
  cat("Filtered to", nrow(data), "significant terms (p.adjust <= 0.05)\n")
  
  return(data)
}

#' Get available datasets
#' @return list of available dataset options
get_available_datasets <- function() {
  
  # Check what datasets are available
  datasets_dir <- "datasets"
  available_datasets <- list()
  
  # Predefined datasets
  dataset_files <- c(
    "all_modalities" = "datasets/all_modalities.rds",
    "mast_only" = "datasets/mast_only.rds", 
    "mast_crispi" = "datasets/mast_crispi.rds",
    "crispi_only" = "datasets/crispi_only.rds",
    "crispa_only" = "datasets/crispa_only.rds"
  )
  
  # Check which files exist
  for (name in names(dataset_files)) {
    file_path <- dataset_files[[name]]
    
    # Special case for all_modalities - check main directory too
    if (name == "all_modalities" && !file.exists(file_path)) {
      fallback_path <- "../all_enrichment_padj005_complete_with_direction.rds"
      if (file.exists(fallback_path)) {
        available_datasets[[name]] <- list(
          name = "All Modalities",
          description = "Complete dataset (MAST + CRISPRi + CRISPRa)",
          file = fallback_path,
          size = round(file.info(fallback_path)$size / (1024^2), 1)
        )
      }
    } else if (file.exists(file_path)) {
      available_datasets[[name]] <- list(
        name = switch(name,
                     "all_modalities" = "All Modalities",
                     "mast_only" = "MAST Only", 
                     "mast_crispi" = "MAST + CRISPRi",
                     "crispi_only" = "CRISPRi Only",
                     "crispa_only" = "CRISPRa Only"),
        description = switch(name,
                           "all_modalities" = "Complete dataset (MAST + CRISPRi + CRISPRa)",
                           "mast_only" = "Genetic mutations only",
                           "mast_crispi" = "Mutations + knockdowns", 
                           "crispi_only" = "CRISPR knockdowns only",
                           "crispa_only" = "CRISPR activation only"),
        file = file_path,
        size = round(file.info(file_path)$size / (1024^2), 1)
      )
    }
  }
  
  return(available_datasets)
}

#' Filter enrichment data based on selections
#' @param data Input data frame
#' @param gene Gene/mutation to filter by
#' @param cluster Cluster to filter by  
#' @param enrichment_type Enrichment type to filter by
#' @param direction Direction to filter by (UP/DOWN/ALL)
#' @param modality Modality to filter by (optional)
#' @return Filtered data frame
get_significant_terms_from_consolidated <- function(data, gene = NULL, cluster = NULL, 
                                                   enrichment_type = NULL, direction = "ALL",
                                                   modality = NULL, analysis_type = NULL,
                                                   experiment = NULL, pval_threshold = 0.05) {
  
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame())
  }
  
  # Start with all data
  filtered_data <- data
  
  # Filter by gene/mutation
  if (!is.null(gene) && gene != "All" && gene != "") {
    # Check which column to use
    if ("gene" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$gene == gene, ]
    } else if ("mutation_perturbation" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$mutation_perturbation == gene, ]
    }
  }
  
  # Filter by cluster
  if (!is.null(cluster) && cluster != "All") {
    filtered_data <- filtered_data[filtered_data$cluster == cluster, ]
  }
  
  # Filter by enrichment type
  if (!is.null(enrichment_type) && enrichment_type != "All") {
    filtered_data <- filtered_data[filtered_data$enrichment_type == enrichment_type, ]
  }
  
  # Filter by direction (fixed bug from original)
  if (!is.null(direction) && direction != "ALL") {
    filtered_data <- filtered_data[filtered_data$direction == direction, ]
  }
  
  # Filter by modality if specified
  if (!is.null(modality) && modality != "All" && "modality" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$modality == modality, ]
  }
  
  # Filter by analysis type (method) if specified
  if (!is.null(analysis_type) && analysis_type != "All" && "method" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$method == analysis_type, ]
  }
  
  # Filter by experiment if specified
  if (!is.null(experiment) && experiment != "All" && experiment != "default" && "experiment" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$experiment == experiment, ]
  }
  
  # Filter by p-value threshold
  if (!is.null(pval_threshold) && "p.adjust" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$p.adjust <= pval_threshold, ]
  }
  
  # Sort by p-value and limit for performance
  if (nrow(filtered_data) > 0) {
    filtered_data <- filtered_data[order(filtered_data$p.adjust), ]
    
    # Limit to top 1000 for performance
    if (nrow(filtered_data) > 1000) {
      filtered_data <- filtered_data[1:1000, ]
    }
  }
  
  return(filtered_data)
}

#' Get dataset summary statistics
#' @param data Input data frame
#' @return List with summary statistics
get_dataset_summary <- function(data) {
  
  if (is.null(data) || nrow(data) == 0) {
    return(list(
      total_terms = 0,
      genes = 0,
      clusters = 0,
      enrichment_types = 0,
      modalities = 0
    ))
  }
  
  summary_stats <- list(
    total_terms = nrow(data),
    genes = length(unique(data$mutation_perturbation)),
    clusters = length(unique(data$cluster)),
    enrichment_types = length(unique(data$enrichment_type))
  )
  
  # Add modality info if available
  if ("modality" %in% names(data)) {
    summary_stats$modalities = length(unique(data$modality[!is.na(data$modality)]))
    summary_stats$modality_breakdown = table(data$modality, useNA = "ifany")
  }
  
  # Add method info if available  
  if ("method" %in% names(data)) {
    summary_stats$method_breakdown = table(data$method, useNA = "ifany")
  }
  
  return(summary_stats)
}

#' Format dataset info for display
#' @param summary_stats Summary statistics from get_dataset_summary
#' @return Formatted text string
format_dataset_info <- function(summary_stats) {
  
  info_text <- paste(
    paste("Terms:", format(summary_stats$total_terms, big.mark = ",")),
    paste("Genes:", summary_stats$genes),
    paste("Clusters:", summary_stats$clusters),
    paste("Databases:", summary_stats$enrichment_types),
    sep = "\n"
  )
  
  if (!is.null(summary_stats$modalities) && summary_stats$modalities > 0) {
    info_text <- paste(info_text, paste("Modalities:", summary_stats$modalities), sep = "\n")
  }
  
  return(info_text)
}

# ============================================================================
# VISUALIZATION HELPERS  
# ============================================================================

#' Create a summary table for enrichment types
#' @param data Input data frame
#' @return Data frame for display
create_enrichment_summary_table <- function(data) {
  
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame(Enrichment_Type = character(), Count = integer()))
  }
  
  summary_data <- as.data.frame(table(data$enrichment_type))
  names(summary_data) <- c("Enrichment_Type", "Count")
  summary_data <- summary_data[order(summary_data$Count, decreasing = TRUE), ]
  
  return(summary_data)
}

#' Create a summary table for genes
#' @param data Input data frame  
#' @return Data frame for display
create_gene_summary_table <- function(data) {
  
  if (is.null(data) || nrow(data) == 0) {
    return(data.frame(Gene = character(), Total_Terms = integer()))
  }
  
  gene_summary <- as.data.frame(table(data$mutation_perturbation))
  names(gene_summary) <- c("Gene", "Total_Terms") 
  gene_summary <- gene_summary[order(gene_summary$Total_Terms, decreasing = TRUE), ]
  
  return(gene_summary)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Clean column names for display
#' @param df Input data frame
#' @return Data frame with cleaned column names
clean_display_columns <- function(df) {
  
  # Select key columns for display
  display_cols <- c("Description", "p.adjust", "direction")
  
  # Add optional columns if they exist
  optional_cols <- c("count", "fold_enrichment", "geneID", "BgRatio", "GeneRatio")
  for (col in optional_cols) {
    if (col %in% names(df)) {
      display_cols <- c(display_cols, col)
    }
  }
  
  # Return subset with available columns
  available_cols <- intersect(display_cols, names(df))
  return(df[, available_cols, drop = FALSE])
}

#' Validate uploaded dataset
#' @param file_path Path to uploaded file
#' @return List with validation results
validate_dataset <- function(file_path) {
  
  tryCatch({
    
    # Load the data
    data <- readRDS(file_path)
    
    # Check if it's a data frame
    if (!is.data.frame(data)) {
      return(list(valid = FALSE, message = "File must contain a data frame"))
    }
    
    # Check minimum size
    if (nrow(data) < 10) {
      return(list(valid = FALSE, message = "Dataset too small (minimum 10 rows)"))
    }
    
    # Check required columns
    required_cols <- c("mutation_perturbation", "cluster", "enrichment_type", "direction", "p.adjust", "Description")
    
    # Handle alternative column names
    if (!"mutation_perturbation" %in% names(data) && "gene" %in% names(data)) {
      names(data)[names(data) == "gene"] <- "mutation_perturbation"
    }
    
    missing_cols <- setdiff(required_cols, names(data))
    if (length(missing_cols) > 0) {
      return(list(valid = FALSE, message = paste("Missing columns:", paste(missing_cols, collapse = ", "))))
    }
    
    # Check p.adjust values
    if (any(is.na(data$p.adjust)) || any(data$p.adjust < 0) || any(data$p.adjust > 1)) {
      return(list(valid = FALSE, message = "Invalid p.adjust values (must be 0-1)"))
    }
    
    return(list(valid = TRUE, data = data, message = "Dataset validated successfully"))
    
  }, error = function(e) {
    return(list(valid = FALSE, message = paste("Error loading file:", e$message)))
  })
}

#' Load differential expression results
#' @return Data frame with DE results or NULL
load_de_results <- function() {
  # First check environment variable
  if (de_file != "" && file.exists(de_file)) {
    cat("Loading DE results from environment variable path:", de_file, "\n")
    return(readRDS(de_file))
  }
  
  # Try default location
  default_path <- file.path(dirname(getwd()), "full_DE_results.rds")
  if (file.exists(default_path)) {
    cat("Loading DE results from default path:", default_path, "\n")
    return(readRDS(default_path))
  }
  
  cat("Warning: Could not find full_DE_results.rds\n")
  return(NULL)
}

cat("Global functions loaded successfully\n")