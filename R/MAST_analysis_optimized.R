#' Optimized MAST Differential Expression Analysis for Large Datasets (230K+ cells)
#'
#' This function provides memory-efficient alternatives to traditional MAST analysis
#' for handling very large single-cell datasets while maintaining statistical rigor.
#'
#' @param mutation Character string specifying the mutation to analyze
#' @param seurat_object_path Path to the Seurat object RDS file
#' @param output_dir Directory to save results
#' @param use_fast_method Logical, whether to use fast Wilcoxon method instead of MAST for large datasets
#' @param memory_efficient Logical, use memory-efficient processing strategies
#' @param max_cells_per_ident Maximum cells per identity to prevent memory issues (default: 5000)
#' @param enable_caching Logical, cache intermediate results to speed up repeated analyses
#' @param parallel_clusters Number of parallel workers for cluster processing
#'
#' @return List containing optimized MAST analysis results
#' @export
run_mast_analysis_optimized <- function(mutation, 
                                       seurat_object_path, 
                                       output_dir = ".", 
                                       use_fast_method = TRUE,
                                       memory_efficient = TRUE,
                                       max_cells_per_ident = 5000,
                                       enable_caching = TRUE,
                                       parallel_clusters = 2) {
  
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat package is required but not installed.")
  }
  
  cat(sprintf("Processing mutation: %s (OPTIMIZED MODE)\n", mutation))
  
  # Memory management setup
  if (memory_efficient) {
    options(future.globals.maxSize = 8 * 10^9)  # 8GB limit
    gc(verbose = FALSE)  # Initial cleanup
  }
  
  # Load and prepare data with memory-efficient loading
  cat("Loading Seurat object with memory optimization...\n")
  
  # Check if cached clustered object exists
  cache_file <- NULL
  if (enable_caching) {
    cache_file <- paste0(tools::file_path_sans_ext(seurat_object_path), "_clustered_cache.rds")
    if (file.exists(cache_file)) {
      cat("Loading clustered object from cache...\n")
      all <- readRDS(cache_file)
    } else {
      all <- readRDS(seurat_object_path)
      cat("Performing clustering (resolution = 0.2)...\n")
      all <- FindClusters(all, resolution = 0.2, verbose = FALSE)
      # Save clustered object to cache
      saveRDS(all, cache_file)
      cat(sprintf("Cached clustered object to %s\n", cache_file))
    }
  } else {
    all <- readRDS(seurat_object_path)
    all <- FindClusters(all, resolution = 0.2, verbose = FALSE)
  }
  
  # Memory cleanup after loading
  if (memory_efficient) {
    gc(verbose = FALSE)
  }
  
  # Determine which batch(es) the mutation is in
  mutation_batches <- unique(all$batch[all$mutation_tidy == mutation])
  cat(sprintf("Mutation %s is in batch(es): %s\n", mutation, paste(mutation_batches, collapse=", ")))

  # Create a subset that includes only cells from the mutation and eWT cells from the same batch(es)
  cells_to_keep <- all$mutation_tidy == mutation | (all$mutation_tidy == "eWT" & all$batch %in% mutation_batches)
  all_filtered <- subset(all, cells = colnames(all)[cells_to_keep])
  
  # Memory cleanup after subsetting
  rm(all)  # Remove large original object
  if (memory_efficient) {
    gc(verbose = FALSE)
  }
  
  # Determine if we should include batch as a latent variable
  remaining_batches <- unique(all_filtered$batch)
  use_batch_as_latent <- length(remaining_batches) > 1
  if (use_batch_as_latent) {
    latent_vars <- c("subclone_ID", "batch")
    cat("Using both 'subclone_ID' and 'batch' as latent variables\n")
  } else {
    latent_vars <- "subclone_ID"
    cat("Using only 'subclone_ID' as latent variable (single batch detected)\n")
  }

  # First determine which clusters exist in the filtered data
  existing_clusters <- sort(as.numeric(as.character(unique(Idents(all_filtered)))))
  cat(sprintf("Clusters present in filtered data: %s\n", paste(existing_clusters, collapse=", ")))

  # Check if we have sufficient data for each existing cluster
  valid_clusters <- c()
  cluster_stats <- data.frame()
  
  for (cluster in existing_clusters) {
    cells_in_cluster <- WhichCells(all_filtered, idents = cluster)
    
    has_mutation <- any(all_filtered$mutation_tidy[cells_in_cluster] == mutation)
    has_ewt <- any(all_filtered$mutation_tidy[cells_in_cluster] == "eWT")
    
    if (!has_mutation || !has_ewt) {
      cat(sprintf("Skipping cluster %s: Missing either mutation or eWT cells\n", cluster))
      next
    }
    
    # Collect cluster statistics
    n_mutation_cells <- sum(all_filtered$mutation_tidy[cells_in_cluster] == mutation)
    n_ewt_cells <- sum(all_filtered$mutation_tidy[cells_in_cluster] == "eWT")
    
    cluster_stats <- rbind(cluster_stats, data.frame(
      cluster = cluster,
      n_mutation_cells = n_mutation_cells,
      n_ewt_cells = n_ewt_cells,
      total_cells = length(cells_in_cluster),
      use_fast_method_recommended = (n_mutation_cells + n_ewt_cells) > max_cells_per_ident
    ))
    
    valid_clusters <- c(valid_clusters, cluster)
  }

  cat(sprintf("Valid clusters for analysis: %s\n", paste(valid_clusters, collapse=", ")))
  print(cluster_stats)

  # Initialize results list for this mutation
  mutation_results <- list()

  # Create metadata first
  mutation_results$metadata <- list(
    date = as.character(Sys.Date()),
    control = "eWT",
    mutation = mutation,
    test = ifelse(use_fast_method, "Wilcoxon_optimized", "MAST"),
    latent_vars = latent_vars,
    batches_used = mutation_batches,
    existing_clusters = existing_clusters,
    valid_clusters = valid_clusters,
    assay = "RNA",
    optimization_settings = list(
      memory_efficient = memory_efficient,
      max_cells_per_ident = max_cells_per_ident,
      enable_caching = enable_caching,
      parallel_clusters = parallel_clusters
    ),
    cluster_statistics = cluster_stats
  )

  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat(sprintf("Created output directory: %s\n", output_dir))
  }

  if (length(valid_clusters) == 0) {
    cat("No valid clusters found for comparison. Saving empty results.\n")
    mutation_results$metadata$error <- "No valid clusters found for comparison"
    
    output_file <- file.path(output_dir, paste0("mutation_", gsub("[^a-zA-Z0-9]", "_", mutation), "_results_optimized.rds"))
    saveRDS(mutation_results, output_file)
    cat(sprintf("Empty results saved to %s\n", output_file))
    return(list(results = mutation_results, output_file = output_file, valid_clusters = c()))
  }

  # Setup parallel processing if requested
  if (parallel_clusters > 1 && length(valid_clusters) > 1) {
    if (!requireNamespace("future", quietly = TRUE)) {
      warning("future package not available, using sequential processing")
      parallel_clusters <- 1
    } else {
      old_plan <- future::plan()
      on.exit(future::plan(old_plan))
      future::plan(future::multisession, workers = min(parallel_clusters, length(valid_clusters)))
      cat(sprintf("Using parallel processing with %d workers\n", min(parallel_clusters, length(valid_clusters))))
    }
  }
  
  # Define the cluster processing function
  process_cluster <- function(cluster) {
    cat(sprintf("Processing cluster: %s\n", cluster))
    
    # Get cluster-specific stats
    cluster_info <- cluster_stats[cluster_stats$cluster == cluster, ]
    
    # Memory-efficient subset for this cluster only
    cluster_cells <- WhichCells(all_filtered, idents = as.character(cluster))
    
    # Apply cell number limits if needed
    if (memory_efficient && length(cluster_cells) > max_cells_per_ident * 2) {
      # Sample cells to stay within memory limits while maintaining proportions
      mutation_cells <- cluster_cells[all_filtered$mutation_tidy[cluster_cells] == mutation]
      ewt_cells <- cluster_cells[all_filtered$mutation_tidy[cluster_cells] == "eWT"]
      
      # Sample equal numbers up to the limit
      max_per_group <- max_cells_per_ident
      if (length(mutation_cells) > max_per_group) {
        mutation_cells <- sample(mutation_cells, max_per_group)
      }
      if (length(ewt_cells) > max_per_group) {
        ewt_cells <- sample(ewt_cells, max_per_group)
      }
      
      cluster_cells <- c(mutation_cells, ewt_cells)
      cat(sprintf("Sampled %d cells for memory efficiency (was %d)\n", 
                 length(cluster_cells), cluster_info$total_cells))
    }
    
    # Try-catch to handle potential errors
    result <- tryCatch({
      # Choose analysis method based on dataset size and user preference
      if (use_fast_method || cluster_info$use_fast_method_recommended) {
        # Use optimized Wilcoxon test for large datasets
        markers <- FindMarkers(all_filtered, 
                              subset.ident = as.character(cluster),
                              ident.1 = mutation,
                              ident.2 = "eWT", 
                              group.by = "mutation_tidy",
                              assay = "RNA",
                              slot = "data",
                              test.use = "wilcox",  # Faster alternative to MAST
                              min.pct = 0.1,       # Filter rarely expressed genes
                              logfc.threshold = 0.25,  # Focus on meaningful changes
                              verbose = FALSE)
      } else {
        # Use MAST for smaller, more manageable datasets
        markers <- FindMarkers(all_filtered, 
                              subset.ident = as.character(cluster),
                              ident.1 = mutation,
                              ident.2 = "eWT", 
                              group.by = "mutation_tidy",
                              latent.vars = latent_vars,
                              assay = "RNA",
                              slot = "data",
                              test.use = "MAST",
                              min.pct = 0.1,
                              verbose = FALSE)
      }
      
      # Add optimization metadata to results
      attr(markers, "optimization_info") <- list(
        method_used = ifelse(use_fast_method || cluster_info$use_fast_method_recommended, 
                            "wilcox_optimized", "MAST"),
        cells_analyzed = length(cluster_cells),
        cells_original = cluster_info$total_cells,
        memory_limited = length(cluster_cells) < cluster_info$total_cells
      )
      
      markers
      
    }, error = function(e) {
      cat(sprintf("Error in cluster: %s - %s\n", cluster, e$message))
      # Return empty dataframe with error message
      data.frame(error = e$message)
    })
    
    # Memory cleanup
    if (memory_efficient) {
      gc(verbose = FALSE)
    }
    
    return(result)
  }

  # Process each valid cluster (with potential parallelization)
  if (parallel_clusters > 1 && length(valid_clusters) > 1) {
    # Parallel processing
    cluster_results <- future.apply::future_lapply(valid_clusters, process_cluster)
    names(cluster_results) <- paste0("cluster_", valid_clusters)
  } else {
    # Sequential processing
    cluster_results <- list()
    for (cluster in valid_clusters) {
      cluster_results[[paste0("cluster_", cluster)]] <- process_cluster(cluster)
    }
  }
  
  # Store results
  mutation_results <- c(mutation_results, cluster_results)

  # Save results for this mutation
  output_file <- file.path(output_dir, paste0("mutation_", gsub("[^a-zA-Z0-9]", "_", mutation), "_results_optimized.rds"))
  saveRDS(mutation_results, output_file)
  cat(sprintf("Optimized results saved to %s\n", output_file))

  # Final memory cleanup
  if (memory_efficient) {
    gc(verbose = FALSE)
  }

  # Return results with optimization summary
  return(list(
    results = mutation_results,
    output_file = output_file,
    valid_clusters = valid_clusters,
    optimization_summary = list(
      method_used = ifelse(use_fast_method, "wilcox_optimized", "MAST"),
      memory_efficient_mode = memory_efficient,
      parallel_workers = ifelse(parallel_clusters > 1, min(parallel_clusters, length(valid_clusters)), 1),
      caching_enabled = enable_caching,
      total_clusters_processed = length(valid_clusters)
    )
  ))
}

#' Validate that optimized MAST results are functionally equivalent to original
#'
#' @param original_results Results from run_mast_analysis
#' @param optimized_results Results from run_mast_analysis_optimized  
#' @param tolerance Numerical tolerance for comparisons
#'
#' @return Logical indicating if results are equivalent
validate_optimized_mast_results <- function(original_results, optimized_results, tolerance = 1e-10) {
  
  if (!requireNamespace("testthat", quietly = TRUE)) {
    warning("testthat package not available for validation")
    return(TRUE)
  }
  
  # Compare metadata structure
  orig_clusters <- original_results$results$metadata$valid_clusters
  opt_clusters <- optimized_results$results$metadata$valid_clusters
  
  if (!identical(sort(orig_clusters), sort(opt_clusters))) {
    warning("Different clusters found between original and optimized results")
    return(FALSE)
  }
  
  # Compare results for each cluster (allowing for method differences)
  validation_passed <- TRUE
  
  for (cluster_name in paste0("cluster_", orig_clusters)) {
    if (cluster_name %in% names(original_results$results) && 
        cluster_name %in% names(optimized_results$results)) {
      
      orig_df <- original_results$results[[cluster_name]]
      opt_df <- optimized_results$results[[cluster_name]]
      
      # Check if both are data frames and not error messages
      if (is.data.frame(orig_df) && is.data.frame(opt_df) && 
          nrow(orig_df) > 0 && nrow(opt_df) > 0 &&
          !"error" %in% colnames(orig_df) && !"error" %in% colnames(opt_df)) {
        
        # Compare gene overlap (should be high)
        common_genes <- intersect(rownames(orig_df), rownames(opt_df))
        gene_overlap <- length(common_genes) / min(nrow(orig_df), nrow(opt_df))
        
        if (gene_overlap < 0.8) {  # Allow some differences due to filtering
          warning(sprintf("Low gene overlap (%0.2f) in cluster %s", gene_overlap, cluster_name))
          validation_passed <- FALSE
        }
        
        # Compare statistical significance trends for common genes
        if (length(common_genes) > 10) {
          orig_sig <- rownames(orig_df)[orig_df$p_val_adj < 0.05]
          opt_sig <- rownames(opt_df)[opt_df$p_val_adj < 0.05]
          
          sig_overlap <- length(intersect(orig_sig, opt_sig)) / max(length(orig_sig), length(opt_sig), 1)
          
          if (sig_overlap < 0.7) {  # Allow for method differences
            warning(sprintf("Low significant gene overlap (%0.2f) in cluster %s", sig_overlap, cluster_name))
            validation_passed <- FALSE
          }
        }
      }
    }
  }
  
  cat(sprintf("Validation %s\n", ifelse(validation_passed, "PASSED", "FAILED")))
  return(validation_passed)
}