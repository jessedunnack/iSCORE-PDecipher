#' Data Sampling Utilities for Large Datasets
#' 
#' Functions to create sampled subsets of large single-cell datasets
#' for preview mode and performance optimization.
#'
#' @author Claude
#' @name data_sampling

#' Sample cells from a Seurat object
#'
#' @param seurat_obj Seurat object to sample from
#' @param n_cells Number of cells to sample (default: 50000)
#' @param seed Random seed for reproducibility (default: 42)
#' @param preserve_proportions Preserve cluster proportions (default: TRUE)
#' @param min_cells_per_cluster Minimum cells per cluster (default: 100)
#' 
#' @return Sampled Seurat object
#' @export
sample_seurat_cells <- function(seurat_obj, 
                               n_cells = 50000,
                               seed = 42,
                               preserve_proportions = TRUE,
                               min_cells_per_cluster = 100) {
  
  set.seed(seed)
  
  total_cells <- ncol(seurat_obj)
  
  # If dataset is already smaller than n_cells, return as is
  if (total_cells <= n_cells) {
    message(sprintf("Dataset has %d cells, no sampling needed", total_cells))
    return(seurat_obj)
  }
  
  message(sprintf("Sampling %d cells from %d total cells (%.1f%%)", 
                  n_cells, total_cells, 100 * n_cells / total_cells))
  
  if (preserve_proportions && !is.null(seurat_obj$seurat_clusters)) {
    # Sample proportionally from each cluster
    clusters <- as.character(seurat_obj$seurat_clusters)
    cluster_table <- table(clusters)
    
    # Calculate cells per cluster
    cluster_proportions <- cluster_table / total_cells
    cells_per_cluster <- round(cluster_proportions * n_cells)
    
    # Ensure minimum cells per cluster
    cells_per_cluster[cells_per_cluster < min_cells_per_cluster] <- min_cells_per_cluster
    
    # Adjust if we exceed n_cells
    if (sum(cells_per_cluster) > n_cells) {
      # Scale down proportionally
      scale_factor <- n_cells / sum(cells_per_cluster)
      cells_per_cluster <- round(cells_per_cluster * scale_factor)
      cells_per_cluster[cells_per_cluster < min_cells_per_cluster] <- min_cells_per_cluster
    }
    
    # Sample from each cluster
    sampled_cells <- c()
    
    for (cluster_id in names(cluster_table)) {
      cluster_cells <- which(clusters == cluster_id)
      n_sample <- min(cells_per_cluster[cluster_id], length(cluster_cells))
      
      sampled_cells <- c(sampled_cells, 
                        sample(cluster_cells, n_sample, replace = FALSE))
    }
    
    # If we need more cells, sample randomly from remaining
    if (length(sampled_cells) < n_cells) {
      remaining_cells <- setdiff(1:total_cells, sampled_cells)
      n_additional <- n_cells - length(sampled_cells)
      sampled_cells <- c(sampled_cells,
                        sample(remaining_cells, n_additional, replace = FALSE))
    }
    
  } else {
    # Simple random sampling
    sampled_cells <- sample(1:total_cells, n_cells, replace = FALSE)
  }
  
  # Subset the Seurat object
  sampled_obj <- seurat_obj[, sampled_cells]
  
  # Add sampling metadata
  sampled_obj@misc$sampling_info <- list(
    original_n_cells = total_cells,
    sampled_n_cells = length(sampled_cells),
    sampling_fraction = length(sampled_cells) / total_cells,
    seed = seed,
    timestamp = Sys.time(),
    preserve_proportions = preserve_proportions
  )
  
  message(sprintf("Sampling complete: %d cells selected", ncol(sampled_obj)))
  
  return(sampled_obj)
}

#' Create preview dataset for Shiny app
#'
#' @param seurat_obj Full Seurat object
#' @param preview_cells Number of cells for preview (default: 50000)
#' @param cache_dir Directory to cache preview data (default: "cache/")
#' @param force_recreate Force recreation even if cache exists (default: FALSE)
#' 
#' @return List with full and preview Seurat objects
#' @export
create_preview_dataset <- function(seurat_obj, 
                                  preview_cells = 50000,
                                  cache_dir = "cache/",
                                  force_recreate = FALSE) {
  
  # Create cache directory if needed
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
  }
  
  # Generate cache filename based on dataset characteristics
  dataset_hash <- digest::digest(list(
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    preview_cells = preview_cells
  ))
  
  cache_file <- file.path(cache_dir, sprintf("preview_%s.rds", dataset_hash))
  
  # Check if cached preview exists
  if (!force_recreate && file.exists(cache_file)) {
    message("Loading cached preview dataset...")
    preview_obj <- readRDS(cache_file)
    
    return(list(
      full = seurat_obj,
      preview = preview_obj,
      is_preview = TRUE,
      cache_file = cache_file
    ))
  }
  
  # Create preview
  message("Creating preview dataset...")
  preview_obj <- sample_seurat_cells(seurat_obj, n_cells = preview_cells)
  
  # Save to cache
  message(sprintf("Saving preview to cache: %s", cache_file))
  saveRDS(preview_obj, cache_file)
  
  return(list(
    full = seurat_obj,
    preview = preview_obj,
    is_preview = TRUE,
    cache_file = cache_file
  ))
}

#' Extract UMAP data for fast plotting
#'
#' @param seurat_obj Seurat object
#' @param sample_n Optional number of cells to sample
#' @param metadata_cols Additional metadata columns to include
#' 
#' @return Data frame with UMAP coordinates and metadata
#' @export
extract_umap_data <- function(seurat_obj, 
                             sample_n = NULL,
                             metadata_cols = c("seurat_clusters", "orig.ident")) {
  
  # Sample if requested
  if (!is.null(sample_n) && sample_n < ncol(seurat_obj)) {
    cell_indices <- sample(1:ncol(seurat_obj), sample_n)
  } else {
    cell_indices <- 1:ncol(seurat_obj)
  }
  
  # Extract UMAP coordinates
  if ("umap" %in% names(seurat_obj@reductions)) {
    umap_coords <- seurat_obj@reductions$umap@cell.embeddings[cell_indices, ]
  } else if ("UMAP" %in% names(seurat_obj@reductions)) {
    umap_coords <- seurat_obj@reductions$UMAP@cell.embeddings[cell_indices, ]
  } else {
    stop("No UMAP reduction found in Seurat object")
  }
  
  # Create base data frame
  umap_data <- data.frame(
    cell = rownames(umap_coords),
    UMAP1 = umap_coords[, 1],
    UMAP2 = umap_coords[, 2],
    stringsAsFactors = FALSE
  )
  
  # Add metadata columns
  for (col in metadata_cols) {
    if (col %in% colnames(seurat_obj@meta.data)) {
      umap_data[[col]] <- seurat_obj@meta.data[[col]][cell_indices]
    }
  }
  
  return(umap_data)
}

#' Progressive loading strategy for UMAP plots
#'
#' @param seurat_obj Seurat object
#' @param stages Vector of cell counts for progressive loading
#' 
#' @return List of UMAP data frames at different resolutions
#' @export
create_progressive_umap <- function(seurat_obj, 
                                   stages = c(1000, 5000, 20000, 50000, Inf)) {
  
  total_cells <- ncol(seurat_obj)
  
  # Adjust stages to dataset size
  stages <- stages[stages < total_cells]
  stages <- c(stages, total_cells)
  
  progressive_data <- list()
  
  for (i in seq_along(stages)) {
    n_cells <- min(stages[i], total_cells)
    
    message(sprintf("Creating stage %d: %d cells", i, n_cells))
    
    progressive_data[[paste0("stage_", i)]] <- list(
      n_cells = n_cells,
      data = extract_umap_data(seurat_obj, sample_n = n_cells)
    )
  }
  
  return(progressive_data)
}

#' Memory-efficient data loader for Shiny
#'
#' @param data_path Path to Seurat object RDS file
#' @param preview_mode Start in preview mode (default: TRUE)
#' @param preview_cells Number of cells for preview (default: 50000)
#' 
#' @return Reactive dataset object for Shiny
#' @export
create_reactive_dataset <- function(data_path, 
                                   preview_mode = TRUE,
                                   preview_cells = 50000) {
  
  # This would be used in Shiny server function
  # Returns a function that can be called to get the appropriate dataset
  
  dataset_loader <- function(use_preview = preview_mode) {
    
    if (!file.exists(data_path)) {
      stop(sprintf("Data file not found: %s", data_path))
    }
    
    # Load full dataset
    message("Loading dataset...")
    full_data <- readRDS(data_path)
    
    if (use_preview && ncol(full_data) > preview_cells) {
      # Create or load preview
      result <- create_preview_dataset(
        full_data, 
        preview_cells = preview_cells
      )
      
      return(list(
        current = result$preview,
        full_available = TRUE,
        is_preview = TRUE,
        n_cells_current = ncol(result$preview),
        n_cells_full = ncol(full_data)
      ))
      
    } else {
      # Use full dataset
      return(list(
        current = full_data,
        full_available = FALSE,
        is_preview = FALSE,
        n_cells_current = ncol(full_data),
        n_cells_full = ncol(full_data)
      ))
    }
  }
  
  return(dataset_loader)
}

#' Estimate memory usage for dataset
#'
#' @param seurat_obj Seurat object
#' @param include_assays Include assay data in calculation (default: TRUE)
#' 
#' @return List with memory usage statistics
#' @export
estimate_memory_usage <- function(seurat_obj, include_assays = TRUE) {
  
  memory_stats <- list()
  
  # Basic dimensions
  memory_stats$n_cells <- ncol(seurat_obj)
  memory_stats$n_genes <- nrow(seurat_obj)
  
  # Estimate memory for different components
  if (include_assays) {
    # RNA assay (sparse matrix)
    assay_size <- object.size(seurat_obj@assays$RNA)
    memory_stats$assay_mb <- as.numeric(assay_size) / 1024^2
  }
  
  # Metadata
  metadata_size <- object.size(seurat_obj@meta.data)
  memory_stats$metadata_mb <- as.numeric(metadata_size) / 1024^2
  
  # Reductions
  reduction_sizes <- list()
  for (red_name in names(seurat_obj@reductions)) {
    red_size <- object.size(seurat_obj@reductions[[red_name]])
    reduction_sizes[[red_name]] <- as.numeric(red_size) / 1024^2
  }
  memory_stats$reductions_mb <- reduction_sizes
  
  # Total
  memory_stats$total_mb <- sum(unlist(memory_stats[grep("_mb$", names(memory_stats))]))
  
  # Recommendations
  memory_stats$recommended_ram_gb <- ceiling(memory_stats$total_mb * 2 / 1024)
  
  return(memory_stats)
}