# Global configuration and functions for PD Enrichment Explorer

# Load core analysis functions
source("data_import_functions.R")
source("enrichment_analysis.R")
source("file_access.R")
source("term_extraction_fixed.R")
source("visualization.R")
source("plot_enrichment_results_v7.51.R")
source("unified_enrichment_heatmaps.R")

# Additional required libraries
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(RColorBrewer)
library(tidyr)
library(stringr)
library(tibble)

# Global configuration
APP_CONFIG <- list(
  # Data paths
  enrichment_base_dir = "/Users/hockemeyer/Desktop/Functional Enrichment/enrichment_results",
  
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
    "Intersection" = "intersection",
    "Union" = "union"
  ),
  
  # Heatmap metric options
  heatmap_metrics = c(
    "P-value (-log10)" = "p.adjust",
    "Fold Enrichment" = "FoldEnrichment",
    "Z-score" = "zScore",
    "NES (GSEA)" = "NES"
  )
)

# Utility functions

#' Get available genes from enrichment results directory
get_available_genes <- function(base_dir = APP_CONFIG$enrichment_base_dir, 
                               analysis_type = NULL) {
  if (!dir.exists(base_dir)) {
    return(character(0))
  }
  
  if (!is.null(analysis_type)) {
    type_dir <- file.path(base_dir, analysis_type)
    if (dir.exists(type_dir)) {
      genes <- list.dirs(type_dir, full.names = FALSE, recursive = FALSE)
      return(sort(unique(genes[genes != ""])))
    }
  } else {
    # Get genes from all analysis types
    all_genes <- character(0)
    for (type in APP_CONFIG$analysis_types) {
      type_dir <- file.path(base_dir, type)
      if (dir.exists(type_dir)) {
        genes <- list.dirs(type_dir, full.names = FALSE, recursive = FALSE)
        all_genes <- c(all_genes, genes[genes != ""])
      }
    }
    return(sort(unique(all_genes)))
  }
}

#' Get available clusters for a gene
get_available_clusters <- function(base_dir = APP_CONFIG$enrichment_base_dir,
                                 analysis_type, gene) {
  gene_dir <- file.path(base_dir, analysis_type, gene)
  if (!dir.exists(gene_dir)) {
    return(character(0))
  }
  
  clusters <- list.dirs(gene_dir, full.names = FALSE, recursive = FALSE)
  clusters <- clusters[grepl("^cluster_", clusters)]
  
  # Sort numerically
  cluster_nums <- as.numeric(gsub("cluster_", "", clusters))
  clusters <- clusters[order(cluster_nums)]
  
  return(clusters)
}

#' Get available experiments for a gene/cluster
get_available_experiments <- function(base_dir = APP_CONFIG$enrichment_base_dir,
                                    analysis_type, gene, cluster) {
  cluster_dir <- file.path(base_dir, analysis_type, gene, cluster)
  if (!dir.exists(cluster_dir)) {
    return(character(0))
  }
  
  experiments <- list.dirs(cluster_dir, full.names = FALSE, recursive = FALSE)
  experiments <- experiments[experiments != ""]
  
  # For MAST, should primarily be "default"
  # For MixScale, various experiment codes
  return(sort(experiments))
}

#' Load enrichment result safely
load_enrichment_safe <- function(file_path) {
  tryCatch({
    if (file.exists(file_path)) {
      result <- readRDS(file_path)
      return(result)
    } else {
      return(NULL)
    }
  }, error = function(e) {
    message("Error loading file: ", file_path)
    message("Error: ", e$message)
    return(NULL)
  })
}

#' Create a summary of enrichment results
summarize_enrichment <- function(enrichment_data, 
                               enrichment_type,
                               pval_threshold = 0.05) {
  if (is.null(enrichment_data)) {
    return(data.frame(
      total_terms = 0,
      significant_terms = 0,
      top_term = NA,
      top_pvalue = NA
    ))
  }
  
  # Extract the result based on enrichment type
  if (enrichment_type == "STRING") {
    if (is.list(enrichment_data) && "enrichment" %in% names(enrichment_data)) {
      result_df <- enrichment_data$enrichment
      pval_col <- "fdr"
      desc_col <- "description"
    } else {
      return(data.frame(
        total_terms = 0,
        significant_terms = 0,
        top_term = NA,
        top_pvalue = NA
      ))
    }
  } else if (enrichment_type == "GSEA") {
    # Handle GSEA format
    pval_col <- "padj"
    desc_col <- "pathway"
    result_df <- enrichment_data
  } else {
    # Standard enrichResult S4 object
    if (methods::is(enrichment_data, "enrichResult")) {
      result_df <- enrichment_data@result
      pval_col <- "p.adjust"
      desc_col <- "Description"
    } else {
      return(data.frame(
        total_terms = 0,
        significant_terms = 0,
        top_term = NA,
        top_pvalue = NA
      ))
    }
  }
  
  # Calculate summary statistics
  total_terms <- nrow(result_df)
  
  if (pval_col %in% names(result_df)) {
    significant_terms <- sum(result_df[[pval_col]] < pval_threshold, na.rm = TRUE)
    
    if (significant_terms > 0) {
      top_idx <- which.min(result_df[[pval_col]])
      top_term <- result_df[[desc_col]][top_idx]
      top_pvalue <- result_df[[pval_col]][top_idx]
    } else {
      top_term <- NA
      top_pvalue <- NA
    }
  } else {
    significant_terms <- 0
    top_term <- NA
    top_pvalue <- NA
  }
  
  return(data.frame(
    total_terms = total_terms,
    significant_terms = significant_terms,
    top_term = top_term,
    top_pvalue = top_pvalue
  ))
}

#' Format p-value for display
format_pvalue <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return(formatC(p, format = "e", digits = 2))
  return(formatC(p, format = "f", digits = 4))
}

#' Create a color palette for heatmaps
get_heatmap_colors <- function(metric_type) {
  if (metric_type == "NES") {
    # Divergent colors for NES
    colorRamp2(c(-3, 0, 3), c("blue", "white", "red"))
  } else {
    # Sequential colors for other metrics
    colorRamp2(c(0, 2, 5, 10, 20), 
              c("white", "#fee5d9", "#fcae91", "#fb6a4a", "#cb181d"))
  }
}

#' Validate gene list input
validate_gene_list <- function(genes, species = "human") {
  # Remove empty strings and whitespace
  genes <- trimws(genes)
  genes <- genes[genes != ""]
  
  # Remove duplicates
  genes <- unique(genes)
  
  # Basic validation
  valid_genes <- genes[nchar(genes) > 1]
  
  return(list(
    valid_genes = valid_genes,
    n_input = length(genes),
    n_valid = length(valid_genes),
    n_invalid = length(genes) - length(valid_genes)
  ))
}

#' Create safe file name from user input
safe_filename <- function(name) {
  # Remove special characters and spaces
  name <- gsub("[^[:alnum:]_-]", "_", name)
  # Remove multiple underscores
  name <- gsub("_+", "_", name)
  # Add timestamp
  paste0(name, "_", format(Sys.time(), "%Y%m%d_%H%M%S"))
}

# Set options for the app
options(
  shiny.maxRequestSize = 50*1024^2,  # 50MB file upload limit
  scipen = 999  # Avoid scientific notation
)

# Create necessary directories if they don't exist
if (!dir.exists("www")) dir.create("www")
if (!dir.exists("data")) dir.create("data")
if (!dir.exists("modules")) dir.create("modules")

print("Global configuration loaded successfully")