# Startup Configuration for iSCORE-PDecipher Shiny App
# Manages dataset selection and initial app configuration

#' Detect cross-platform path to datasets
#' @return Base path string appropriate for the current platform
detect_cross_platform_path <- function() {
  # Try to get from user configuration first
  if (requireNamespace("iSCORE.PDecipher", quietly = TRUE)) {
    config_dir <- tryCatch({
      iSCORE.PDecipher::get_parent_data_dir()
    }, error = function(e) NULL)
    
    if (!is.null(config_dir) && dir.exists(config_dir)) {
      return(config_dir)
    }
  }
  
  # Platform-specific fallback paths
  if (.Platform$OS.type == "windows") {
    # Windows paths
    candidates <- c(
      "E:/ASAP/scRNASeq/PerturbSeq/final",
      "C:/ASAP/scRNASeq/PerturbSeq/final",
      "D:/ASAP/scRNASeq/PerturbSeq/final"
    )
  } else {
    # Unix-like paths (Linux, macOS, WSL)
    candidates <- c(
      "/mnt/e/ASAP/scRNASeq/PerturbSeq/final",
      "/mnt/c/ASAP/scRNASeq/PerturbSeq/final", 
      "/mnt/d/ASAP/scRNASeq/PerturbSeq/final",
      "~/ASAP/scRNASeq/PerturbSeq/final",
      "/opt/ASAP/scRNASeq/PerturbSeq/final"
    )
  }
  
  # Find first existing path
  for (path in candidates) {
    if (dir.exists(path)) {
      return(path)
    }
  }
  
  # If nothing found, return the default for the platform
  if (.Platform$OS.type == "windows") {
    return("E:/ASAP/scRNASeq/PerturbSeq/final")
  } else {
    return("/mnt/e/ASAP/scRNASeq/PerturbSeq/final")
  }
}

#' Create startup configuration
#' @return List with startup configuration
create_startup_config <- function() {
  
  # Detect platform and set appropriate base directory
  base_dir <- detect_cross_platform_path()
  
  # Check which datasets are available
  datasets <- list()
  
  # iSCORE-PD (mutations only)
  iscore_pd_dir <- file.path(base_dir, "iSCORE-PD")
  if (dir.exists(iscore_pd_dir)) {
    has_de <- file.exists(file.path(iscore_pd_dir, "full_DE_results.rds"))
    has_enrichment <- file.exists(file.path(iscore_pd_dir, "all_enrichment_padj005_complete_with_direction.rds"))
    has_seurat <- file.exists(file.path(iscore_pd_dir, "iSCORE-PD_final.rds"))
    
    if (has_de && has_enrichment) {
      datasets$iscore_pd <- list(
        name = "iSCORE-PD (Mutations Only)",
        path = iscore_pd_dir,
        description = "13 PD mutations analyzed across ~210K iPSC-derived neurons",
        stats = list(
          n_cells = 210000,
          n_genes = 36601,
          n_clusters = 16,
          n_mutations = 13,
          has_crispr = FALSE
        ),
        files = list(
          de = file.path(iscore_pd_dir, "full_DE_results.rds"),
          enrichment = file.path(iscore_pd_dir, "all_enrichment_padj005_complete_with_direction.rds"),
          seurat = if(has_seurat) file.path(iscore_pd_dir, "iSCORE-PD_final.rds") else NULL,
          enrichment_dir = file.path(iscore_pd_dir, "enrichment_results")
        ),
        status = "ready"
      )
    }
  }
  
  # iSCORE-PD + CRISPRi
  iscore_crispri_dir <- file.path(base_dir, "iSCORE-PD_plus_CRISPRi")
  if (dir.exists(iscore_crispri_dir)) {
    has_de <- file.exists(file.path(iscore_crispri_dir, "full_DE_results.rds"))
    has_enrichment <- file.exists(file.path(iscore_crispri_dir, "all_enrichment_padj005_complete_with_direction.rds"))
    has_seurat <- file.exists(file.path(iscore_crispri_dir, "iSCORE-PD_plus_CRISPRi_final.rds"))
    
    if (has_de && has_enrichment) {
      datasets$iscore_pd_plus_crispri <- list(
        name = "iSCORE-PD + CRISPRi",
        path = iscore_crispri_dir,
        description = "Combined dataset with 13 mutations + 10 CRISPRi perturbations (~230K cells)",
        stats = list(
          n_cells = 231874,
          n_genes = 36601,
          n_clusters = 16,
          n_mutations = 13,
          n_perturbations = 10,
          has_crispr = TRUE
        ),
        files = list(
          de = file.path(iscore_crispri_dir, "full_DE_results.rds"),
          enrichment = file.path(iscore_crispri_dir, "all_enrichment_padj005_complete_with_direction.rds"),
          seurat = if(has_seurat) file.path(iscore_crispri_dir, "iSCORE-PD_plus_CRISPRi_final.rds") else NULL,
          enrichment_dir = file.path(iscore_crispri_dir, "enrichment_results")
        ),
        status = "ready"
      )
    }
  }
  
  # Determine default dataset (prefer CRISPRi if available)
  default_dataset <- if ("iscore_pd_plus_crispri" %in% names(datasets)) {
    "iscore_pd_plus_crispri"
  } else if ("iscore_pd" %in% names(datasets)) {
    "iscore_pd"
  } else {
    NULL
  }
  
  # Check if we have at least one dataset
  app_ready <- length(datasets) > 0
  
  config <- list(
    app_ready = app_ready,
    available_datasets = datasets,
    default_dataset = default_dataset,
    base_dir = base_dir,
    timestamp = Sys.time()
  )
  
  return(config)
}

#' Load dataset into environment
#' @param dataset_key Key of the dataset to load
#' @param config Startup configuration
#' @return TRUE if successful
load_dataset <- function(dataset_key, config) {
  
  if (!dataset_key %in% names(config$available_datasets)) {
    stop("Dataset not found: ", dataset_key)
  }
  
  dataset <- config$available_datasets[[dataset_key]]
  
  # Set environment variables for global.R
  Sys.setenv(ISCORE_DATA_DIR = dataset$path)
  Sys.setenv(ISCORE_DE_FILE = dataset$files$de)
  Sys.setenv(ISCORE_ENRICHMENT_FILE = dataset$files$enrichment)
  Sys.setenv(ISCORE_ENRICHMENT_DIR = dataset$files$enrichment_dir)
  
  # Store dataset info globally
  options(iscore.current_dataset = dataset_key)
  options(iscore.dataset_info = dataset)
  
  cat("Loaded dataset:", dataset$name, "\n")
  cat("  Path:", dataset$path, "\n")
  cat("  Cells:", format(dataset$stats$n_cells, big.mark = ","), "\n")
  cat("  Clusters:", dataset$stats$n_clusters, "\n")
  
  if (dataset$stats$has_crispr) {
    cat("  CRISPRi perturbations:", dataset$stats$n_perturbations, "\n")
  }
  
  return(TRUE)
}

#' Get current dataset info
#' @return Current dataset information
get_current_dataset <- function() {
  dataset_key <- getOption("iscore.current_dataset")
  dataset_info <- getOption("iscore.dataset_info")
  
  if (is.null(dataset_key) || is.null(dataset_info)) {
    return(NULL)
  }
  
  return(list(
    key = dataset_key,
    info = dataset_info
  ))
}

cat("Startup configuration module loaded successfully\n")