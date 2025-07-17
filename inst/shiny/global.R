# Consolidated Global Configuration for iSCORE-PDecipher Shiny App
# Single point for all library loading and configuration

# =============================================================================
# LIBRARY LOADING - All required packages in one place
# =============================================================================

# Core Shiny libraries (order matters!)
library(shiny)
library(shinyjs)  # Load shinyjs early before UI definition
library(shinyWidgets)
library(shinycssloaders)
library(DT)

# Data manipulation and visualization
library(dplyr)
library(ggplot2)
library(plotly)
library(glue)
library(tibble)
library(tidyr)
library(stringr)

# Color and styling
library(viridis)
library(RColorBrewer)
library(colourpicker)

# Advanced visualization (lazy load these - they're heavy)
# library(ComplexHeatmap)  # Load on demand
# library(circlize)        # Load on demand

# =============================================================================
# GLOBAL CONFIGURATION
# =============================================================================

# Shiny app settings
options(shiny.maxRequestSize = 500*1024^2)  # 500MB upload limit

# Environment variables
data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
de_file <- Sys.getenv("ISCORE_DE_FILE", "")
enrichment_file <- Sys.getenv("ISCORE_ENRICHMENT_FILE", "")
enrichment_dir <- Sys.getenv("ISCORE_ENRICHMENT_DIR", "")

# Application configuration
APP_CONFIG <- list(
  # Data paths
  enrichment_base_dir = enrichment_dir,
  
  # Analysis parameters
  default_pval_threshold = 0.05,
  gsea_pval_threshold = 0.25,
  max_terms_display = 50,
  max_heatmap_rows = 500,
  
  # Visualization settings
  default_colors = c(
    "UP" = "#fee8c8",
    "DOWN" = "#deebf7",
    "ALL" = "#d9f0a3"
  ),
  
  # Available options
  analysis_types = c("MAST", "MixScale"),
  enrichment_types = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
  directions = c("ALL", "UP", "DOWN"),
  
  # Method comparison options
  comparison_methods = c(
    "All" = "all",
    "MAST only" = "MAST",
    "MixScale only" = "MixScale",
    "Union" = "union",
    "Intersection" = "intersection"
  )
)

# =============================================================================
# UTILITY FUNCTIONS (Essential only - moved complex ones to modules)
# =============================================================================

#' Get significant terms with filtering
#' Simplified version for core filtering needs
get_significant_terms_simple <- function(data, 
                                       analysis_type = NULL,
                                       gene = NULL, 
                                       cluster = NULL,
                                       enrichment_type = NULL,
                                       direction = NULL,
                                       modality = NULL) {
  
  if (is.null(data) || nrow(data) == 0) return(data.frame())
  
  # Apply filters
  filtered_data <- data
  
  if (!is.null(analysis_type) && analysis_type != "ALL") {
    if ("analysis_type" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$analysis_type == analysis_type, ]
    } else if ("method" %in% names(filtered_data)) {
      # Map method to analysis_type
      method_filter <- switch(analysis_type,
                             "MAST" = "MAST",
                             "MixScale" = c("MixScale", "MixScale_CRISPRa"),
                             analysis_type)
      filtered_data <- filtered_data[filtered_data$method %in% method_filter, ]
    }
  }
  
  if (!is.null(gene) && gene != "ALL") {
    filtered_data <- filtered_data[filtered_data$mutation_perturbation == gene, ]
  }
  
  if (!is.null(cluster) && cluster != "ALL") {
    filtered_data <- filtered_data[filtered_data$cluster == cluster, ]
  }
  
  if (!is.null(enrichment_type) && enrichment_type != "ALL") {
    filtered_data <- filtered_data[filtered_data$enrichment_type == enrichment_type, ]
  }
  
  if (!is.null(direction) && direction != "ALL") {
    filtered_data <- filtered_data[filtered_data$direction == direction, ]
  }
  
  if (!is.null(modality) && modality != "ALL") {
    if ("modality" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$modality == modality, ]
    }
  }
  
  return(filtered_data)
}

#' Load heavy packages on demand
#' @param package_name Name of package to load
load_package_on_demand <- function(package_name) {
  if (!requireNamespace(package_name, quietly = TRUE)) {
    warning(paste("Package", package_name, "not available"))
    return(FALSE)
  }
  
  if (!package_name %in% (.packages())) {
    library(package_name, character.only = TRUE)
  }
  return(TRUE)
}

# =============================================================================
# FILTERING FUNCTIONS
# =============================================================================

#' Filter consolidated enrichment data based on selection criteria
#' @param data Consolidated enrichment data frame
#' @param gene Gene/mutation to filter by
#' @param cluster Cell cluster to filter by 
#' @param enrichment_type Enrichment database to filter by
#' @param direction Gene regulation direction (ALL, UP, DOWN)
#' @param modality Analysis modality (if applicable)
#' @param analysis_type Analysis method (MAST, MixScale)
#' @param experiment Experiment ID to filter by
#' @param pval_threshold P-value threshold for significance
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
  
  # Filter by direction (FIXED: always filter to prevent inflation)
  if (!is.null(direction)) {
    filtered_data <- filtered_data[filtered_data$direction == direction, ]
  }
  
  # Filter by modality if specified (FIXED: handle CRISPRi/CRISPRa properly)
  if (!is.null(modality) && modality != "All") {
    if ("modality" %in% names(filtered_data)) {
      # Use modality column if it exists
      filtered_data <- filtered_data[filtered_data$modality == modality, ]
    } else {
      # Fallback: filter by analysis_type or method patterns for CRISPRi/CRISPRa
      if (modality == "CRISPRi") {
        if ("analysis_type" %in% names(filtered_data)) {
          filtered_data <- filtered_data[grepl("MixScale.*CRISPRi|CRISPRi", filtered_data$analysis_type, ignore.case = TRUE), ]
        } else if ("method" %in% names(filtered_data)) {
          filtered_data <- filtered_data[grepl("MixScale.*CRISPRi|CRISPRi", filtered_data$method, ignore.case = TRUE), ]
        }
      } else if (modality == "CRISPRa") {
        if ("analysis_type" %in% names(filtered_data)) {
          filtered_data <- filtered_data[grepl("MixScale.*CRISPRa|CRISPRa", filtered_data$analysis_type, ignore.case = TRUE), ]
        } else if ("method" %in% names(filtered_data)) {
          filtered_data <- filtered_data[grepl("MixScale.*CRISPRa|CRISPRa", filtered_data$method, ignore.case = TRUE), ]
        }
      }
    }
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

# =============================================================================

# =============================================================================
# GENE ASSOCIATION LOOKUP
# =============================================================================

# Load gene association functions if available in package
tryCatch({
  # Check if we're in package mode (installed) or development mode
  if (requireNamespace('iSCORE.PDecipher', quietly = TRUE)) {
    # Package mode - functions should be exported
    message('Loading gene associations from package...')
    success <- iSCORE.PDecipher::load_gene_associations()
    if (success) {
      message('Gene associations loaded successfully from package')
    } else {
      message('Failed to load gene associations from package')
    }
  } else {
    # Development mode - source the file
    gene_lookup_file <- '../../R/gene_association_lookup.R'
    if (file.exists(gene_lookup_file)) {
      message('Loading gene associations from development file...')
      source(gene_lookup_file)
      success <- load_gene_associations()
      if (success) {
        message('Gene associations loaded successfully from development file')
      } else {
        message('Failed to load gene associations from development file')
      }
    } else {
      # Try alternative path for when running from project root
      alt_gene_lookup_file <- 'R/gene_association_lookup.R'
      if (file.exists(alt_gene_lookup_file)) {
        message('Loading gene associations from alternative development path...')
        source(alt_gene_lookup_file)
        success <- load_gene_associations()
        if (success) {
          message('Gene associations loaded successfully from alternative path')
        } else {
          message('Failed to load gene associations from alternative path')
        }
      } else {
        message('Gene association file not found - gene display will be unavailable')
      }
    }
  }
}, error = function(e) {
  message('Could not load gene associations: ', e$message)
  message('Gene display functionality will be limited')
})

# INITIALIZATION
# =============================================================================

cat("Global configuration loaded successfully\n")
cat("Environment variables configured\n")
cat("Ready for data loading\n")
