# Consolidated Global Configuration for iSCORE-PDecipher Shiny App
# Single point for all library loading and configuration

# =============================================================================
# LIBRARY LOADING - All required packages in one place
# =============================================================================

# Core Shiny libraries
library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
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
# INITIALIZATION
# =============================================================================

cat("Global configuration loaded successfully\n")
cat("Environment variables configured\n")
cat("Ready for data loading\n")