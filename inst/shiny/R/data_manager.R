# Centralized Data Manager
# Single point for data loading and caching to prevent multiple loads

# Global data cache
.data_cache <- new.env()

#' Initialize Data Manager
#' Sets up the data caching system
init_data_manager <- function() {
  .data_cache$enrichment_data <- NULL
  .data_cache$de_data <- NULL
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
get_available_clusters <- function() {
  data <- get_enrichment_data()
  if (is.null(data)) return(character())
  return(sort(unique(data$cluster)))
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

# Initialize the data manager
init_data_manager()