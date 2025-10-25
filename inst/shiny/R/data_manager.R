# Centralized Data Manager
# Single point for data loading and caching to prevent multiple loads

# Global data cache
.data_cache <- new.env()

#' Initialize Data Manager
#' Sets up the data caching system
init_data_manager <- function() {
  .data_cache$enrichment_data <- NULL
  .data_cache$de_data <- NULL
  .data_cache$pooled_mixscale_data <- NULL
  .data_cache$pval_correction_type <- "p_weight_BH"  # Default to BH correction
  .data_cache$dataset_type <- NULL  # FPD or CRISPRi
  .data_cache$loading_status <- "not_loaded"
  .data_cache$load_time <- NULL
}

#' Get Enrichment Data (Cached)
#' Loads data once and caches it for all modules
#' @param force_reload Force reload even if cached
#' @return Enrichment data frame
get_enrichment_data <- function(force_reload = FALSE) {
  
  # Return cached data if available
  if (!force_reload && !is.null(.data_cache$enrichment_data)) {
    return(.data_cache$enrichment_data)
  }
  
  # Load data for the first time
  if (.data_cache$loading_status != "loading") {
    .data_cache$loading_status <- "loading"
    
    tryCatch({
      # Use environment variable path
      enrichment_file <- Sys.getenv("ISCORE_ENRICHMENT_FILE", "")
      
      if (enrichment_file == "" || !file.exists(enrichment_file)) {
        stop("Enrichment file not found. Please check ISCORE_ENRICHMENT_FILE environment variable.")
      }
      
      cat("Loading enrichment data from:", enrichment_file, "\n")
      data <- readRDS(enrichment_file)
      
      # Process and clean data
      required_cols <- c("mutation_perturbation", "cluster", "enrichment_type", "direction", "p.adjust", "Description")
      
      # Handle alternative column names
      if (!"mutation_perturbation" %in% names(data) && "gene" %in% names(data)) {
        data$mutation_perturbation <- data$gene
      }
      
      # Add method/modality columns if missing
      if (!"method" %in% names(data) && "analysis_type" %in% names(data)) {
        data$method <- data$analysis_type
      }
      
      # Filter to significant results
      data <- data[!is.na(data$p.adjust) & data$p.adjust <= 0.05, ]
      
      # Cache the data
      .data_cache$enrichment_data <- data
      .data_cache$loading_status <- "loaded"
      .data_cache$load_time <- Sys.time()
      
      cat("Successfully cached", nrow(data), "enrichment terms\n")
      
      return(data)
      
    }, error = function(e) {
      .data_cache$loading_status <- "error"
      cat("Error loading enrichment data:", e$message, "\n")
      return(NULL)
    })
  }
  
  # Return NULL if loading failed
  return(.data_cache$enrichment_data)
}

#' Get Data Loading Status
#' @return Status of data loading
get_loading_status <- function() {
  return(.data_cache$loading_status)
}

#' Get Available Genes
#' @return Vector of available genes/mutations
get_available_genes <- function() {
  data <- get_enrichment_data()
  if (is.null(data)) return(character())
  return(unique(data$mutation_perturbation))
}

#' Get Available Clusters  
#' @return Vector of available clusters
#' Natural Sort Clusters
#' @description Sort cluster names naturally so cluster_10 comes after cluster_9
#' @param clusters Character vector of cluster names
#' @return Naturally sorted cluster names
natural_sort_clusters <- function(clusters) {
  if (length(clusters) == 0) return(clusters)
  
  # Extract numeric parts from cluster names
  numeric_parts <- as.numeric(gsub(".*_", "", clusters))
  
  # If all clusters follow the cluster_X pattern, sort by numeric value
  if (!any(is.na(numeric_parts))) {
    return(clusters[order(numeric_parts)])
  }
  
  # Otherwise, fall back to regular sort
  return(sort(clusters))
}

get_available_clusters <- function() {
  data <- get_enrichment_data()
  if (is.null(data)) return(character())
  return(natural_sort_clusters(unique(data$cluster)))
}

#' Get Data Summary
#' @return List with summary statistics
get_data_summary <- function() {
  data <- get_enrichment_data()
  if (is.null(data)) {
    return(list(
      total_terms = 0,
      genes = 0,
      clusters = 0,
      experiments = 0,
      enrichment_types = 0,
      load_time = NULL
    ))
  }
  
  return(list(
    total_terms = nrow(data),
    genes = length(unique(data$mutation_perturbation)),
    clusters = length(unique(data$cluster)),
    experiments = if("experiment" %in% names(data)) length(unique(data$experiment)) else 1,
    enrichment_types = length(unique(data$enrichment_type)),
    load_time = .data_cache$load_time
  ))
}

#' ==============================================================================
#' POOLED MIXSCALE DATA FUNCTIONS (Added October 24, 2025)
#' ==============================================================================

#' Load pooled MixScale import functions
load_pooled_mixscale_functions <- function() {
  # Source the pooled MixScale import functions
  import_file <- "../../R/import_pooled_mixscale_functions.R"
  if (file.exists(import_file)) {
    source(import_file)
    return(TRUE)
  } else {
    # Try alternative path
    alt_import_file <- "R/import_pooled_mixscale_functions.R"
    if (file.exists(alt_import_file)) {
      source(alt_import_file)
      return(TRUE)
    }
  }
  warning("Could not load pooled MixScale import functions")
  return(FALSE)
}

#' Get Pooled MixScale Data with FDR Correction
#' @param mixscale_dir Directory containing pooled MixScale RDS files
#' @param pval_column Which p-value column: "p_weight", "p_weight_BH", or "p_weight_bonferroni"
#' @param dataset_type "FPD" or "CRISPRi"
#' @param force_reload Force reload even if cached
#' @return Pooled MixScale data
get_pooled_mixscale_data <- function(
  mixscale_dir = NULL,
  pval_column = "p_weight_BH",
  dataset_type = NULL,
  force_reload = FALSE
) {

  # Check if we should use cached data
  cache_key <- paste(mixscale_dir, pval_column, dataset_type, sep = "_")

  if (!force_reload &&
      !is.null(.data_cache$pooled_mixscale_data) &&
      identical(.data_cache$pval_correction_type, pval_column) &&
      identical(.data_cache$dataset_type, dataset_type)) {
    return(.data_cache$pooled_mixscale_data)
  }

  # Load import functions if not already loaded
  if (!exists("import_pooled_mixscale_data")) {
    if (!load_pooled_mixscale_functions()) {
      stop("Could not load pooled MixScale import functions")
    }
  }

  # Use environment variable if directory not specified
  if (is.null(mixscale_dir)) {
    mixscale_dir <- Sys.getenv("ISCORE_POOLED_MIXSCALE_DIR", "")
    if (mixscale_dir == "") {
      stop("mixscale_dir not provided and ISCORE_POOLED_MIXSCALE_DIR not set")
    }
  }

  # Auto-detect dataset type if not specified
  if (is.null(dataset_type)) {
    dataset_type <- Sys.getenv("ISCORE_DATASET_TYPE", "")
    if (dataset_type == "") {
      # Try to detect from directory name
      if (grepl("FPD", mixscale_dir)) {
        dataset_type <- "FPD"
      } else if (grepl("CRISPRi", mixscale_dir)) {
        dataset_type <- "CRISPRi"
      } else {
        warning("Could not detect dataset type, defaulting to FPD")
        dataset_type <- "FPD"
      }
    }
  }

  cat("Loading pooled MixScale data...\n")
  cat("  Directory:", mixscale_dir, "\n")
  cat("  P-value column:", pval_column, "\n")
  cat("  Dataset type:", dataset_type, "\n")

  tryCatch({
    data <- import_pooled_mixscale_data(
      mixscale_dir = mixscale_dir,
      pval_column = pval_column,
      dataset_type = dataset_type
    )

    # Cache the data
    .data_cache$pooled_mixscale_data <- data
    .data_cache$pval_correction_type <- pval_column
    .data_cache$dataset_type <- dataset_type
    .data_cache$load_time <- Sys.time()

    cat("Successfully loaded and cached pooled MixScale data\n")

    return(data)

  }, error = function(e) {
    cat("Error loading pooled MixScale data:", e$message, "\n")
    return(NULL)
  })
}

#' Get Pooled Enrichment Data with FDR Correction
#' @param dataset "FPD" or "CRISPRi"
#' @param pval_correction "none", "BH", or "bonferroni"
#' @param force_reload Force reload even if cached
#' @return Enrichment data
get_pooled_enrichment_data <- function(
  dataset = NULL,
  pval_correction = "BH",
  force_reload = FALSE
) {

  # Load import functions if not already loaded
  if (!exists("import_enrichment_with_correction")) {
    if (!load_pooled_mixscale_functions()) {
      stop("Could not load pooled MixScale import functions")
    }
  }

  # Use cached dataset type if not specified
  if (is.null(dataset)) {
    dataset <- .data_cache$dataset_type
    if (is.null(dataset)) {
      dataset <- Sys.getenv("ISCORE_DATASET_TYPE", "FPD")
    }
  }

  cat("Loading pooled enrichment data...\n")
  cat("  Dataset:", dataset, "\n")
  cat("  P-value correction:", pval_correction, "\n")

  tryCatch({
    data <- import_enrichment_with_correction(
      dataset = dataset,
      pval_correction = pval_correction
    )

    return(data)

  }, error = function(e) {
    cat("Error loading pooled enrichment data:", e$message, "\n")
    return(NULL)
  })
}

#' Set P-value Correction Type
#' @param pval_type "p_weight", "p_weight_BH", or "p_weight_bonferroni"
set_pval_correction <- function(pval_type) {
  valid_types <- c("p_weight", "p_weight_BH", "p_weight_bonferroni")
  if (!pval_type %in% valid_types) {
    stop("Invalid p-value type. Must be one of: ", paste(valid_types, collapse = ", "))
  }
  .data_cache$pval_correction_type <- pval_type
  cat("P-value correction type set to:", pval_type, "\n")
}

#' Get Current P-value Correction Type
#' @return Current p-value correction type
get_pval_correction <- function() {
  if (is.null(.data_cache$pval_correction_type)) {
    return("p_weight_BH")  # Default
  }
  return(.data_cache$pval_correction_type)
}

#' Set Dataset Type
#' @param dataset_type "FPD" or "CRISPRi"
set_dataset_type <- function(dataset_type) {
  if (!dataset_type %in% c("FPD", "CRISPRi")) {
    stop("Invalid dataset type. Must be 'FPD' or 'CRISPRi'")
  }
  .data_cache$dataset_type <- dataset_type
  cat("Dataset type set to:", dataset_type, "\n")
}

#' Get Current Dataset Type
#' @return Current dataset type
get_dataset_type <- function() {
  if (is.null(.data_cache$dataset_type)) {
    return(Sys.getenv("ISCORE_DATASET_TYPE", "FPD"))  # Default
  }
  return(.data_cache$dataset_type)
}

# Initialize the data manager
init_data_manager()