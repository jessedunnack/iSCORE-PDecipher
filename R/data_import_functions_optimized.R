# Optimized functions for importing MAST and MixScale data with memory efficiency for 230K+ cells

#' Memory-efficient cluster ID extraction from file path
#'
#' @param file_path The full path to a results file
#' @return String: the extracted cluster ID
extract_cluster_id_fast <- function(file_path) {
  # More efficient regex extraction
  filename <- basename(file_path)
  
  # Pre-compiled patterns for better performance
  clust_pattern <- "clust_([0-9]+)"
  cluster_pattern <- "Cluster([0-9]+)"
  
  if (regexpr(clust_pattern, filename) > 0) {
    cluster_num <- regmatches(filename, regexpr(clust_pattern, filename))
    return(paste0("cluster_", gsub("clust_", "", cluster_num)))
  }
  
  if (regexpr(cluster_pattern, file_path) > 0) {
    return(regmatches(file_path, regexpr(cluster_pattern, file_path)))
  }
  
  return("cluster_unknown")
}

#' Optimized import of MAST differential expression results with memory efficiency
#'
#' @param input_dir Directory containing MAST result files
#' @param lazy_loading Logical, whether to use lazy loading strategies
#' @param memory_efficient Logical, use memory-efficient processing
#' @param parallel_loading Logical, use parallel loading for multiple files
#' @param cache_results Logical, cache processed results for faster re-loading
#' @return List of structured MAST results
import_mast_data_optimized <- function(input_dir, 
                                     lazy_loading = TRUE,
                                     memory_efficient = TRUE,
                                     parallel_loading = TRUE,
                                     cache_results = TRUE) {
  
  # Memory management setup
  if (memory_efficient) {
    options(future.globals.maxSize = 4 * 10^9)  # 4GB limit for imports
    gc(verbose = FALSE)
  }
  
  # Check for cached results first
  cache_file <- file.path(input_dir, "mast_data_cache.rds")
  if (cache_results && file.exists(cache_file)) {
    cache_info <- file.info(cache_file)
    cat(sprintf("Loading MAST data from cache (%s, %.1f MB)...\n", 
               format(cache_info$mtime, "%Y-%m-%d %H:%M"), 
               cache_info$size / 1e6))
    return(readRDS(cache_file))
  }
  
  # Get list of RDS files with optimized file discovery
  files <- list.files(input_dir, pattern = "\\.rds$", full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No RDS files found in ", input_dir)
    return(list())
  }
  
  cat(sprintf("Processing %d MAST result files...\n", length(files)))
  
  # Initialize results list
  results <- list()
  
  # Define single file processing function
  process_mast_file <- function(file_path) {
    tryCatch({
      # Extract mutation name from filename - optimized pattern matching
      filename <- basename(file_path)
      
      # Pre-compile regex patterns for better performance
      new_pattern <- "^mutation_(.+)_results\\.rds$"
      old_pattern <- "^mutation_(.+)_results_RNA_batchspecific.*\\.rds$"
      
      mutation_name <- NULL
      if (grepl(new_pattern, filename)) {
        mutation_name <- gsub(new_pattern, "\\1", filename)
      } else if (grepl(old_pattern, filename)) {
        mutation_name <- gsub(old_pattern, "\\1", filename)
      }
      
      if (is.null(mutation_name)) {
        return(NULL)  # Skip files that don't match expected patterns
      }
      
      # Load RDS file with error handling
      data <- readRDS(file_path)
      
      # Skip if no metadata or empty results
      if (is.null(data$metadata)) return(NULL)
      
      # Initialize file results
      file_results <- list()
      
      # Process clusters with memory efficiency
      for (cluster_name in names(data)) {
        if (cluster_name == "metadata") next
        
        # Skip if no results for this cluster or results are error messages
        if (!is.data.frame(data[[cluster_name]]) || 
            (nrow(data[[cluster_name]]) == 1 && "error" %in% colnames(data[[cluster_name]]))) {
          next
        }
        
        # Get background genes - all genes tested in this analysis (rownames of the DEG table)
        background_genes <- rownames(data[[cluster_name]])
        
        # Create optimized structure with reduced memory footprint
        if (lazy_loading) {
          # Store minimal metadata and file path for lazy loading
          file_results[[cluster_name]] <- list(
            lazy_data = list(
              file_path = file_path,
              cluster_name = cluster_name,
              n_genes = length(background_genes),
              top_genes = head(background_genes, 10)  # Preview
            ),
            metadata = list(
              mutation = mutation_name,
              cluster = cluster_name,
              test = "MAST",
              control = data$metadata$control,
              batches_used = data$metadata$batches_used,
              latent_vars = data$metadata$latent_vars,
              date = data$metadata$date,
              optimization_mode = "lazy_loading"
            ),
            background_genes = background_genes
          )
        } else {
          # Standard loading with optimized structure
          file_results[[cluster_name]] <- list(
            results = data[[cluster_name]],
            metadata = list(
              mutation = mutation_name,
              cluster = cluster_name,
              test = "MAST",
              control = data$metadata$control,
              batches_used = data$metadata$batches_used,
              latent_vars = data$metadata$latent_vars,
              date = data$metadata$date,
              optimization_mode = "standard"
            ),
            background_genes = background_genes
          )
        }
      }
      
      # Add global metadata for this mutation
      if (length(file_results) > 0) {
        file_results$metadata <- data$metadata
        return(list(mutation = mutation_name, data = file_results))
      }
      
      return(NULL)
      
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", basename(file_path), e$message))
      return(NULL)
    })
  }
  
  # Process files with optional parallelization
  if (parallel_loading && length(files) > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("future.apply not available, using sequential processing")
      parallel_loading <- FALSE
    } else {
      # Setup parallel processing
      old_plan <- future::plan()
      on.exit(future::plan(old_plan))
      future::plan(future::multisession, workers = min(4, length(files)))
      
      cat(sprintf("Using parallel loading with %d workers...\n", min(4, length(files))))
      file_results <- future.apply::future_lapply(files, process_mast_file)
    }
  }
  
  if (!parallel_loading) {
    # Sequential processing with progress
    file_results <- list()
    for (i in seq_along(files)) {
      if (i %% max(1, length(files) %/% 10) == 0) {
        cat(sprintf("Processing file %d of %d...\n", i, length(files)))
      }
      file_results[[i]] <- process_mast_file(files[i])
    }
  }
  
  # Organize results by mutation
  for (file_result in file_results) {
    if (!is.null(file_result)) {
      mutation_name <- file_result$mutation
      if (is.null(results[[mutation_name]])) {
        results[[mutation_name]] <- list()
      }
      
      # Merge cluster results
      for (cluster_name in names(file_result$data)) {
        if (cluster_name != "metadata") {
          results[[mutation_name]][[cluster_name]] <- file_result$data[[cluster_name]]
        }
      }
      
      # Add metadata if not present
      if (is.null(results[[mutation_name]]$metadata)) {
        results[[mutation_name]]$metadata <- file_result$data$metadata
      }
    }
  }
  
  # Memory cleanup
  if (memory_efficient) {
    gc(verbose = FALSE)
  }
  
  # Cache results if requested
  if (cache_results && length(results) > 0) {
    saveRDS(results, cache_file)
    cat(sprintf("Cached MAST data to %s\n", cache_file))
  }
  
  cat(sprintf("Successfully imported data for %d mutations\n", length(results)))
  return(results)
}

#' Optimized import of MixScale differential expression results with memory efficiency
#'
#' @param input_dir Directory containing MixScale result files
#' @param modality Optional: specify "CRISPRi" or "CRISPRa" to import only one modality
#' @param lazy_loading Logical, whether to use lazy loading strategies
#' @param memory_efficient Logical, use memory-efficient processing
#' @param parallel_loading Logical, use parallel loading for multiple files
#' @param cache_results Logical, cache processed results for faster re-loading
#' @return List of structured MixScale results
import_mixscale_data_optimized <- function(input_dir, 
                                         modality = NULL,
                                         lazy_loading = TRUE,
                                         memory_efficient = TRUE,
                                         parallel_loading = TRUE,
                                         cache_results = TRUE) {
  
  # Memory management setup
  if (memory_efficient) {
    options(future.globals.maxSize = 4 * 10^9)  # 4GB limit
    gc(verbose = FALSE)
  }
  
  # Check for cached results first
  cache_suffix <- ifelse(is.null(modality), "all", modality)
  cache_file <- file.path(input_dir, sprintf("mixscale_data_%s_cache.rds", cache_suffix))
  if (cache_results && file.exists(cache_file)) {
    cache_info <- file.info(cache_file)
    cat(sprintf("Loading MixScale data from cache (%s, %.1f MB)...\n", 
               format(cache_info$mtime, "%Y-%m-%d %H:%M"), 
               cache_info$size / 1e6))
    return(readRDS(cache_file))
  }
  
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
  
  cat(sprintf("Processing %d MixScale result files...\n", length(rds_files)))
  
  # Initialize results list
  results <- list()
  
  # Define single file processing function
  process_mixscale_file <- function(file_path) {
    tryCatch({
      # Load data safely
      data <- readRDS(file_path)
      if (is.null(data)) return(NULL)
      
      # Extract cluster ID using the optimized helper function
      cluster_id <- extract_cluster_id_fast(file_path)
      
      # Extract modality from path if not specified
      file_modality <- if (grepl("/CRISPRi/", file_path)) "CRISPRi" else if (grepl("/CRISPRa/", file_path)) "CRISPRa" else NA
      
      file_results <- list()
      
      # Process each perturbation with memory efficiency
      for (pert in names(data)) {
        if (pert == "metadata") next
        
        # Get perturbation data and structure it efficiently
        pert_data <- data[[pert]]
        
        # Structure results consistently with memory optimization
        if (is.list(pert_data) && !is.data.frame(pert_data) && !is.null(pert_data$results)) {
          pert_results <- pert_data$results
          pert_metadata <- pert_data$metadata
        } else if (is.data.frame(pert_data)) {
          pert_results <- pert_data
          pert_metadata <- NULL
        } else {
          next  # Skip invalid entries
        }
        
        # Skip empty results
        if (is.null(pert_results) || (is.data.frame(pert_results) && nrow(pert_results) == 0)) {
          next
        }
        
        # Get background genes for enrichment analysis
        if (is.data.frame(pert_results) && "gene_ID" %in% colnames(pert_results)) {
          background_genes <- unique(pert_results$gene_ID)
        } else if (is.data.frame(pert_results)) {
          background_genes <- rownames(pert_results)
        } else {
          background_genes <- character(0)
        }
        
        # Create optimized structure
        if (lazy_loading && is.data.frame(pert_results) && nrow(pert_results) > 1000) {
          # For large results, store summary and file path for lazy loading
          file_results[[pert]] <- list(
            lazy_data = list(
              file_path = file_path,
              perturbation = pert,
              cluster = cluster_id,
              n_genes = length(background_genes),
              top_genes = head(background_genes, 10),
              has_weighted_pvals = any(grepl(":weight", colnames(pert_results)))
            ),
            metadata = list(
              perturbation = pert,
              cluster = cluster_id,
              modality = file_modality,
              test = "MixScale",
              n_genes_tested = length(background_genes),
              optimization_mode = "lazy_loading",
              original_metadata = pert_metadata
            ),
            background_genes = background_genes
          )
        } else {
          # Standard loading for smaller datasets
          file_results[[pert]] <- list(
            results = pert_results,
            metadata = list(
              perturbation = pert,
              cluster = cluster_id,
              modality = file_modality,
              test = "MixScale",
              n_genes_tested = length(background_genes),
              optimization_mode = "standard",
              original_metadata = pert_metadata
            ),
            background_genes = background_genes
          )
        }
      }
      
      return(list(
        cluster_id = cluster_id,
        modality = file_modality,
        data = file_results
      ))
      
    }, error = function(e) {
      warning(sprintf("Error processing file %s: %s", basename(file_path), e$message))
      return(NULL)
    })
  }
  
  # Process files with optional parallelization
  if (parallel_loading && length(rds_files) > 1) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      warning("future.apply not available, using sequential processing")
      parallel_loading <- FALSE
    } else {
      # Setup parallel processing
      old_plan <- future::plan()
      on.exit(future::plan(old_plan))
      future::plan(future::multisession, workers = min(4, length(rds_files)))
      
      cat(sprintf("Using parallel loading with %d workers...\n", min(4, length(rds_files))))
      file_results <- future.apply::future_lapply(rds_files, process_mixscale_file)
    }
  }
  
  if (!parallel_loading) {
    # Sequential processing with progress
    file_results <- list()
    for (i in seq_along(rds_files)) {
      if (i %% max(1, length(rds_files) %/% 10) == 0) {
        cat(sprintf("Processing file %d of %d...\n", i, length(rds_files)))
      }
      file_results[[i]] <- process_mixscale_file(rds_files[i])
    }
  }
  
  # Organize results by modality and perturbation
  for (file_result in file_results) {
    if (!is.null(file_result)) {
      modality_key <- file_result$modality
      cluster_id <- file_result$cluster_id
      
      if (is.na(modality_key)) modality_key <- "Unknown"
      
      # Initialize modality structure
      if (is.null(results[[modality_key]])) {
        results[[modality_key]] <- list()
      }
      
      # Add perturbations
      for (pert in names(file_result$data)) {
        if (is.null(results[[modality_key]][[pert]])) {
          results[[modality_key]][[pert]] <- list()
        }
        results[[modality_key]][[pert]][[cluster_id]] <- file_result$data[[pert]]
      }
    }
  }
  
  # Memory cleanup
  if (memory_efficient) {
    gc(verbose = FALSE)
  }
  
  # Cache results if requested
  if (cache_results && length(results) > 0) {
    saveRDS(results, cache_file)
    cat(sprintf("Cached MixScale data to %s\n", cache_file))
  }
  
  cat(sprintf("Successfully imported MixScale data for %d modalities\n", length(results)))
  return(results)
}

#' Load lazily-stored data on demand
#'
#' @param lazy_data_entry Entry with lazy loading information
#' @return The actual data loaded from file
load_lazy_data <- function(lazy_data_entry) {
  if (is.null(lazy_data_entry$lazy_data)) {
    warning("Not a lazy-loaded entry")
    return(lazy_data_entry)
  }
  
  lazy_info <- lazy_data_entry$lazy_data
  
  # Load the original file
  full_data <- readRDS(lazy_info$file_path)
  
  # Extract the specific cluster/perturbation data
  if ("cluster_name" %in% names(lazy_info)) {
    # MAST data
    target_data <- full_data[[lazy_info$cluster_name]]
  } else if ("perturbation" %in% names(lazy_info)) {
    # MixScale data
    target_data <- full_data[[lazy_info$perturbation]]
  } else {
    stop("Unknown lazy data format")
  }
  
  # Replace lazy_data with actual results
  lazy_data_entry$results <- target_data
  lazy_data_entry$lazy_data <- NULL
  lazy_data_entry$metadata$optimization_mode <- "lazy_loaded"
  
  return(lazy_data_entry)
}

#' Validate optimized data import results
#'
#' @param original_data Data from standard import functions
#' @param optimized_data Data from optimized import functions  
#' @param check_lazy_loading Whether to test lazy loading functionality
#'
#' @return Logical indicating validation success
validate_optimized_import <- function(original_data, optimized_data, check_lazy_loading = TRUE) {
  
  validation_passed <- TRUE
  
  # Check basic structure
  if (length(original_data) != length(optimized_data)) {
    warning("Different number of mutations/modalities found")
    validation_passed <- FALSE
  }
  
  # Check common entries
  common_names <- intersect(names(original_data), names(optimized_data))
  
  for (name in common_names) {
    orig_entry <- original_data[[name]]
    opt_entry <- optimized_data[[name]]
    
    # Check cluster/perturbation structure
    if (name != "metadata") {
      orig_clusters <- names(orig_entry)[names(orig_entry) != "metadata"]
      opt_clusters <- names(opt_entry)[names(opt_entry) != "metadata"]
      
      if (length(setdiff(orig_clusters, opt_clusters)) > 0) {
        warning(sprintf("Missing clusters in optimized data for %s", name))
        validation_passed <- FALSE
      }
      
      # Test lazy loading if enabled
      if (check_lazy_loading) {
        for (cluster in intersect(orig_clusters, opt_clusters)) {
          opt_cluster_data <- opt_entry[[cluster]]
          
          if (!is.null(opt_cluster_data$lazy_data)) {
            # Test lazy loading
            loaded_data <- load_lazy_data(opt_cluster_data)
            
            if (is.null(loaded_data$results)) {
              warning(sprintf("Lazy loading failed for %s %s", name, cluster))
              validation_passed <- FALSE
            }
          }
        }
      }
    }
  }
  
  cat(sprintf("Import validation %s\n", ifelse(validation_passed, "PASSED", "FAILED")))
  return(validation_passed)
}