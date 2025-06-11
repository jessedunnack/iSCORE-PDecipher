# Minimal global configuration for PD Enrichment Explorer
# This version avoids loading complex dependencies that aren't immediately needed

# Core libraries only
library(tibble)
library(tidyr)

# Check for optional packages used in gene symbol conversion
if (requireNamespace("clusterProfiler", quietly = TRUE)) {
  suppressPackageStartupMessages(library(clusterProfiler))
}
if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  suppressPackageStartupMessages(library(org.Hs.eg.db))
}

# Check for optional packages used in pathview visualization
if (requireNamespace("pathview", quietly = TRUE)) {
  suppressPackageStartupMessages(library(pathview))
}
if (requireNamespace("base64enc", quietly = TRUE)) {
  suppressPackageStartupMessages(library(base64enc))
}

# Global configuration
APP_CONFIG <- list(
  # Data paths
  enrichment_results_path = "/Users/hockemeyer/Desktop/Functional Enrichment/enrichment_results",
  consolidated_data_path = "/Users/hockemeyer/Desktop/Functional Enrichment/all_enrichment_padj005_complete_with_direction.rds",
  
  # Analysis types
  analysis_types = c("MAST" = "MAST", "MixScale" = "MixScale"),
  
  # Enrichment types
  enrichment_types = c(
    "GO Biological Process" = "GO_BP",
    "GO Cellular Component" = "GO_CC", 
    "GO Molecular Function" = "GO_MF",
    "KEGG Pathways" = "KEGG",
    "Reactome Pathways" = "Reactome",
    "WikiPathways" = "WikiPathways",
    "STRING PPI" = "STRING",
    "GSEA" = "GSEA"
  ),
  
  # Direction types
  direction_types = c("All" = "ALL", "Up-regulated" = "UP", "Down-regulated" = "DOWN"),
  
  # Analysis parameters
  default_pval_threshold = 0.05,
  gsea_pval_threshold = 0.25,
  max_terms_display = 50,
  max_heatmap_rows = 500,
  
  # Performance settings
  enable_caching = TRUE,
  cache_dir = "cache",
  max_cache_size_mb = 500,
  
  # UI settings
  default_theme = "blue",
  show_debug = FALSE,
  
  # File handling
  max_file_size_mb = 50,
  supported_formats = c(".rds", ".csv", ".txt", ".xlsx"),
  
  # Heatmap configuration
  heatmap_metrics = c("P-value" = "p.adjust", "Fold Enrichment" = "fold_enrichment", "Gene Count" = "count"),
  heatmap_scope = c("Terms across all clusters (this mutation)" = "across_clusters", 
                    "Terms across all mutations (this cluster)" = "across_mutations")
)

# Utility functions
get_available_genes <- function() {
  base_dir <- APP_CONFIG$enrichment_results_path
  genes <- c()
  
  # Try to get available genes from directory structure
  if (dir.exists(file.path(base_dir, "MAST"))) {
    mast_genes <- list.dirs(file.path(base_dir, "MAST"), full.names = FALSE, recursive = FALSE)
    genes <- c(genes, mast_genes)
  }
  
  if (dir.exists(file.path(base_dir, "MixScale"))) {
    mixscale_genes <- list.dirs(file.path(base_dir, "MixScale"), full.names = FALSE, recursive = FALSE)
    genes <- c(genes, mixscale_genes)
  }
  
  unique(genes[genes != ""])
}

get_available_clusters <- function(gene = NULL) {
  if (is.null(gene)) return(character())
  
  base_dir <- APP_CONFIG$enrichment_results_path
  clusters <- c()
  
  # Check both MAST and MixScale
  for (analysis_type in c("MAST", "MixScale")) {
    gene_dir <- file.path(base_dir, analysis_type, gene)
    if (dir.exists(gene_dir)) {
      gene_clusters <- list.dirs(gene_dir, full.names = FALSE, recursive = FALSE)
      clusters <- c(clusters, gene_clusters)
    }
  }
  
  unique(clusters[clusters != ""])
}

get_available_experiments <- function(gene = NULL, cluster = NULL, analysis_type = NULL) {
  if (is.null(gene) || is.null(cluster)) return(character())
  
  base_dir <- APP_CONFIG$enrichment_results_path
  experiments <- c()
  
  # If analysis_type is provided, only check that type
  types_to_check <- if (!is.null(analysis_type)) analysis_type else c("MAST", "MixScale")
  
  for (type in types_to_check) {
    cluster_dir <- file.path(base_dir, type, gene, cluster)
    if (dir.exists(cluster_dir)) {
      exp_dirs <- list.dirs(cluster_dir, full.names = FALSE, recursive = FALSE)
      # Filter out empty strings and diagnostic directories
      exp_dirs <- exp_dirs[!exp_dirs %in% c("", ".", "..", "diagnostics")]
      
      # For MixScale, sort to put "combined" first if it exists
      if (type == "MixScale" && "combined" %in% exp_dirs) {
        exp_dirs <- c("combined", setdiff(exp_dirs, "combined"))
      }
      
      experiments <- c(experiments, exp_dirs)
    }
  }
  
  unique(experiments[experiments != ""])
}

# Summarize enrichment results
summarize_enrichment <- function(result, enrichment_type, pval_threshold = 0.05) {
  if (is.null(result)) {
    return(list(
      total_terms = 0,
      significant_terms = 0,
      top_term = "N/A",
      top_pvalue = 1
    ))
  }
  
  # Extract data based on enrichment type
  if (enrichment_type == "STRING") {
    df <- result$enrichment
    total <- nrow(df)
    sig <- sum(df$fdr < pval_threshold, na.rm = TRUE)
    top_idx <- which.min(df$fdr)
    top_term <- if (length(top_idx) > 0) df$description[top_idx] else "N/A"
    top_pval <- if (length(top_idx) > 0) df$fdr[top_idx] else 1
  } else if (methods::is(result, "enrichResult")) {
    df <- result@result
    total <- nrow(df)
    sig <- sum(df$p.adjust < pval_threshold, na.rm = TRUE)
    top_idx <- which.min(df$p.adjust)
    top_term <- if (length(top_idx) > 0) df$Description[top_idx] else "N/A"
    top_pval <- if (length(top_idx) > 0) df$p.adjust[top_idx] else 1
  } else {
    return(list(
      total_terms = 0,
      significant_terms = 0,
      top_term = "N/A",
      top_pvalue = 1
    ))
  }
  
  list(
    total_terms = total,
    significant_terms = sig,
    top_term = top_term,
    top_pvalue = top_pval
  )
}

# Format p-value for display
format_pvalue <- function(pval) {
  if (is.na(pval) || is.null(pval)) return("N/A")
  if (pval < 0.001) {
    return(sprintf("%.2e", pval))
  } else {
    return(sprintf("%.4f", pval))
  }
}

# Load enrichment result safely
load_enrichment_safe <- function(file_path) {
  # Check if file_path is NULL or empty
  if (is.null(file_path) || length(file_path) == 0 || file_path == "") {
    message("Error: file_path is NULL or empty")
    return(NULL)
  }
  
  tryCatch({
    if (file.exists(file_path)) {
      result <- readRDS(file_path)
      return(result)
    } else {
      message("File does not exist: ", file_path)
      return(NULL)
    }
  }, error = function(e) {
    message("Error loading file: ", file_path)
    message("Error: ", e$message)
    return(NULL)
  })
}

# Simple file existence check
check_enrichment_file <- function(analysis_type, gene, cluster, experiment = "default", 
                                  enrichment_type, direction) {
  base_dir <- APP_CONFIG$enrichment_results_path
  file_path <- file.path(base_dir, analysis_type, gene, cluster, experiment, enrichment_type,
                         paste0(enrichment_type, "_", direction, ".rds"))
  file.exists(file_path)
}

# Simple data loader
load_enrichment_result <- function(analysis_type, gene, cluster, experiment = "default", 
                                   enrichment_type, direction) {
  base_dir <- APP_CONFIG$enrichment_results_path
  file_path <- file.path(base_dir, analysis_type, gene, cluster, experiment, enrichment_type,
                         paste0(enrichment_type, "_", direction, ".rds"))
  
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    return(NULL)
  }
}

# Load consolidated enrichment data
load_consolidated_data <- function() {
  path <- APP_CONFIG$consolidated_data_path
  if (file.exists(path)) {
    tryCatch({
      data <- readRDS(path)
      message("Loaded consolidated data with ", nrow(data), " rows")
      return(data)
    }, error = function(e) {
      message("Error loading consolidated data: ", e$message)
      return(NULL)
    })
  } else {
    message("Consolidated data file not found: ", path)
    return(NULL)
  }
}

# Get significant terms from consolidated data
get_significant_terms_from_consolidated <- function(consolidated_data, 
                                                    analysis_type = NULL,
                                                    gene = NULL,
                                                    cluster = NULL,
                                                    experiment = NULL,
                                                    enrichment_type = NULL,
                                                    direction = NULL,
                                                    pval_threshold = 0.05) {
  
  if (is.null(consolidated_data)) return(NULL)
  
  # Filter the consolidated data based on parameters
  filtered_data <- consolidated_data
  
  if (!is.null(analysis_type)) {
    filtered_data <- filtered_data[filtered_data$method == analysis_type, ]
  }
  
  if (!is.null(gene)) {
    filtered_data <- filtered_data[filtered_data$gene == gene, ]
  }
  
  if (!is.null(cluster)) {
    filtered_data <- filtered_data[filtered_data$cluster == cluster, ]
  }
  
  if (!is.null(experiment)) {
    filtered_data <- filtered_data[filtered_data$experiment == experiment, ]
  }
  
  if (!is.null(enrichment_type)) {
    filtered_data <- filtered_data[filtered_data$enrichment_type == enrichment_type, ]
  }
  
  if (!is.null(direction)) {
    filtered_data <- filtered_data[filtered_data$direction == direction, ]
  }
  
  # Filter by p-value threshold
  # Handle different p-value column names
  if ("p.adjust" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$p.adjust < pval_threshold, ]
  } else if ("fdr" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$fdr < pval_threshold, ]
  } else if ("padj" %in% names(filtered_data)) {
    filtered_data <- filtered_data[filtered_data$padj < pval_threshold, ]
  }
  
  return(filtered_data)
}

# Initialize app data
initialize_app_data <- function() {
  # Load consolidated data on app initialization
  consolidated_data <- load_consolidated_data()
  
  list(
    available_genes = get_available_genes(),
    config = APP_CONFIG,
    session_id = format(Sys.time(), "%Y%m%d_%H%M%S"),
    consolidated_data = consolidated_data
  )
}