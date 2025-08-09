# Tests for Optimized Functions in iSCORE-PDecipher
# Created: August 2025
#
# These tests validate the optimized functions maintain functional equivalence
# while providing performance improvements for 230k+ cell datasets

test_that("optimized MAST analysis maintains functional equivalence", {
  
  skip_if_not_installed("Seurat")
  
  # Create test dataset
  test_data <- create_mast_test_data("LRRK2_G2019S")
  seurat_obj <- test_data$seurat_object
  
  # Scale up for realistic testing
  if (ncol(seurat_obj) < 800) {
    larger_obj <- create_test_seurat(n_cells = 800, n_genes = 100)
    larger_obj$mutation_tidy[1:200] <- "LRRK2_G2019S"
    larger_obj$mutation_tidy[201:400] <- "eWT"
    seurat_obj <- larger_obj
  }
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(seurat_obj, temp_file)
  
  # Run original MAST analysis
  original_results <- run_mast_analysis(
    mutation = "LRRK2_G2019S",
    seurat_object_path = temp_file,
    output_dir = tempdir()
  )
  
  # Run optimized MAST analysis
  optimized_results <- run_mast_analysis_optimized(
    mutation = "LRRK2_G2019S",
    seurat_object_path = temp_file,
    output_dir = tempdir(),
    use_fast_method = FALSE,  # Use MAST for direct comparison
    memory_efficient = TRUE,
    parallel_clusters = 1     # Single threaded for deterministic comparison
  )
  
  # Validate functional equivalence
  expect_true(validate_optimized_mast_results(original_results, optimized_results))
  
  # Validate structure preservation
  expect_true(validate_mast_results(optimized_results))
  expect_equal(length(original_results$valid_clusters), length(optimized_results$valid_clusters))
  
  message("Optimized MAST analysis maintains functional equivalence")
})

test_that("optimized MAST analysis provides performance improvements", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pryr")
  
  # Create larger test dataset for performance testing
  test_obj <- create_test_seurat(n_cells = 1500, n_genes = 150)
  test_obj$mutation_tidy[1:500] <- "SNCA_A53T"
  test_obj$mutation_tidy[501:1000] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  # Measure original performance
  gc()
  mem_before_orig <- pryr::mem_used()
  time_orig <- system.time({
    original_results <- run_mast_analysis(
      mutation = "SNCA_A53T",
      seurat_object_path = temp_file,
      output_dir = tempdir()
    )
  })
  mem_after_orig <- pryr::mem_used()
  
  # Measure optimized performance
  gc()
  mem_before_opt <- pryr::mem_used()
  time_opt <- system.time({
    optimized_results <- run_mast_analysis_optimized(
      mutation = "SNCA_A53T",
      seurat_object_path = temp_file,
      output_dir = tempdir(),
      use_fast_method = TRUE,  # Use fast method for speed
      memory_efficient = TRUE,
      parallel_clusters = 2
    )
  })
  mem_after_opt <- pryr::mem_used()
  
  # Calculate improvements
  speedup <- time_orig[["elapsed"]] / time_opt[["elapsed"]]
  memory_orig <- as.numeric(mem_after_orig - mem_before_orig)
  memory_opt <- as.numeric(mem_after_opt - mem_before_opt)
  memory_improvement <- 1 - (memory_opt / memory_orig)
  
  # Validate improvements
  expect_true(speedup >= 1.2, info = sprintf("Expected speedup >=1.2x, got %.2fx", speedup))
  expect_true(memory_improvement >= 0.1, info = sprintf("Expected memory reduction >=10%%, got %.1f%%", memory_improvement * 100))
  
  # Validate both analyses succeeded
  expect_true(validate_mast_results(original_results))
  expect_true(validate_mast_results(optimized_results))
  
  message(sprintf("Performance improvements: %.2fx speedup, %.1f%% memory reduction", 
                 speedup, memory_improvement * 100))
})

test_that("optimized data import maintains functional equivalence", {
  
  # Create test data directory with mock MAST results
  test_dir <- tempdir()
  import_test_dir <- file.path(test_dir, "import_test")
  dir.create(import_test_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(import_test_dir, recursive = TRUE))
  
  # Create mock result files
  for (i in 1:3) {
    mutation <- paste0("TEST_MUTATION_", i)
    
    mock_results <- list(
      metadata = list(
        date = Sys.Date(),
        control = "eWT",
        mutation = mutation,
        test = "MAST",
        batches_used = c("batch1", "batch2"),
        latent_vars = c("subclone_ID", "batch")
      )
    )
    
    # Add cluster results
    for (cluster in 0:2) {
      n_genes <- 100
      mock_results[[paste0("cluster_", cluster)]] <- data.frame(
        p_val = runif(n_genes, 0, 0.1),
        avg_log2FC = rnorm(n_genes, 0, 1),
        p_val_adj = runif(n_genes, 0, 0.2),
        pct.1 = runif(n_genes, 0, 1),
        pct.2 = runif(n_genes, 0, 1),
        row.names = paste0("GENE_", 1:n_genes, "_", cluster)
      )
    }
    
    saveRDS(mock_results, file.path(import_test_dir, paste0("mutation_", mutation, "_results.rds")))
  }
  
  # Import with original function
  original_data <- import_mast_data(import_test_dir)
  
  # Import with optimized function
  optimized_data <- import_mast_data_optimized(
    import_test_dir,
    lazy_loading = FALSE,  # Disable lazy loading for direct comparison
    memory_efficient = TRUE,
    parallel_loading = TRUE,
    cache_results = FALSE   # Disable caching for clean comparison
  )
  
  # Validate functional equivalence
  expect_true(validate_optimized_import(original_data, optimized_data, check_lazy_loading = FALSE))
  
  # Validate basic structure
  expect_equal(length(original_data), length(optimized_data))
  expect_equal(sort(names(original_data)), sort(names(optimized_data)))
  
  # Check specific mutations
  for (mutation in names(original_data)) {
    if (mutation != "metadata") {
      expect_true(mutation %in% names(optimized_data))
      
      # Check cluster structure
      orig_clusters <- names(original_data[[mutation]])[names(original_data[[mutation]]) != "metadata"]
      opt_clusters <- names(optimized_data[[mutation]])[names(optimized_data[[mutation]]) != "metadata"]
      expect_equal(sort(orig_clusters), sort(opt_clusters))
      
      # Check data integrity for first cluster
      if (length(orig_clusters) > 0) {
        first_cluster <- orig_clusters[1]
        orig_genes <- original_data[[mutation]][[first_cluster]]$background_genes
        opt_genes <- optimized_data[[mutation]][[first_cluster]]$background_genes
        expect_equal(sort(orig_genes), sort(opt_genes))
      }
    }
  }
  
  message("Optimized data import maintains functional equivalence")
})

test_that("optimized data import provides performance benefits", {
  
  skip_if_not_installed("pryr")
  
  # Create larger test dataset
  test_dir <- tempdir()
  import_test_dir <- file.path(test_dir, "perf_import_test")
  dir.create(import_test_dir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(import_test_dir, recursive = TRUE))
  
  # Create multiple larger result files
  for (i in 1:8) {  # More files for realistic testing
    mutation <- paste0("PERF_TEST_MUT_", i)
    
    mock_results <- list(
      metadata = list(
        date = Sys.Date(),
        control = "eWT",
        mutation = mutation,
        test = "MAST"
      )
    )
    
    # Add more clusters with more genes
    for (cluster in 0:4) {
      n_genes <- 800  # Larger gene sets
      mock_results[[paste0("cluster_", cluster)]] <- data.frame(
        p_val = runif(n_genes),
        avg_log2FC = rnorm(n_genes),
        p_val_adj = runif(n_genes),
        pct.1 = runif(n_genes),
        pct.2 = runif(n_genes),
        row.names = paste0("GENE_", 1:n_genes, "_C", cluster, "_M", i)
      )
    }
    
    saveRDS(mock_results, file.path(import_test_dir, paste0("mutation_", mutation, "_results.rds")))
  }
  
  # Benchmark original import
  gc()
  mem_before_orig <- pryr::mem_used()
  time_orig <- system.time({
    original_data <- import_mast_data(import_test_dir)
  })
  mem_after_orig <- pryr::mem_used()
  
  # Benchmark optimized import
  gc()
  mem_before_opt <- pryr::mem_used()
  time_opt <- system.time({
    optimized_data <- import_mast_data_optimized(
      import_test_dir,
      lazy_loading = TRUE,
      memory_efficient = TRUE,
      parallel_loading = TRUE,
      cache_results = TRUE
    )
  })
  mem_after_opt <- pryr::mem_used()
  
  # Calculate improvements
  speedup <- time_orig[["elapsed"]] / time_opt[["elapsed"]]
  memory_orig <- as.numeric(mem_after_orig - mem_before_orig)
  memory_opt <- as.numeric(mem_after_opt - mem_before_opt)
  memory_reduction <- 1 - (memory_opt / memory_orig)
  
  # Validate improvements
  expect_true(speedup >= 1.1, info = sprintf("Expected speedup >=1.1x, got %.2fx", speedup))
  expect_true(memory_reduction >= 0.05, info = sprintf("Expected memory reduction >=5%%, got %.1f%%", memory_reduction * 100))
  
  # Test lazy loading functionality
  if (length(optimized_data) > 0) {
    first_mutation <- names(optimized_data)[1]
    if (first_mutation != "metadata") {
      first_cluster <- names(optimized_data[[first_mutation]])[1]
      cluster_data <- optimized_data[[first_mutation]][[first_cluster]]
      
      if (!is.null(cluster_data$lazy_data)) {
        # Test lazy loading
        loaded_data <- load_lazy_data(cluster_data)
        expect_true(!is.null(loaded_data$results))
        expect_null(loaded_data$lazy_data)
        expect_equal(loaded_data$metadata$optimization_mode, "lazy_loaded")
      }
    }
  }
  
  message(sprintf("Import performance improvements: %.2fx speedup, %.1f%% memory reduction", 
                 speedup, memory_reduction * 100))
})

test_that("caching mechanisms work correctly", {
  
  # Create test dataset for caching test
  test_obj <- create_test_seurat(n_cells = 500, n_genes = 80)
  test_obj$mutation_tidy[1:150] <- "CACHE_TEST_MUT"
  test_obj$mutation_tidy[151:300] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  # First run should create cache
  result1 <- run_mast_analysis_optimized(
    mutation = "CACHE_TEST_MUT",
    seurat_object_path = temp_file,
    output_dir = tempdir(),
    enable_caching = TRUE,
    memory_efficient = TRUE
  )
  
  # Check if cache file was created
  cache_file <- paste0(tools::file_path_sans_ext(temp_file), "_clustered_cache.rds")
  expect_true(file.exists(cache_file))
  on.exit(unlink(cache_file), add = TRUE)
  
  # Second run should use cache (should be faster)
  time_with_cache <- system.time({
    result2 <- run_mast_analysis_optimized(
      mutation = "CACHE_TEST_MUT",
      seurat_object_path = temp_file,
      output_dir = tempdir(),
      enable_caching = TRUE,
      memory_efficient = TRUE
    )
  })
  
  # Validate both results are valid and similar
  expect_true(validate_mast_results(result1))
  expect_true(validate_mast_results(result2))
  expect_equal(length(result1$valid_clusters), length(result2$valid_clusters))
  
  # Cache usage should be faster than typical analysis
  expect_true(time_with_cache[["elapsed"]] < 60)  # Should complete quickly with cache
  
  message("Caching mechanisms validated successfully")
})

test_that("parallel processing works correctly", {
  
  skip_if_not_installed("future.apply")
  
  # Create test dataset with multiple clusters
  test_obj <- create_test_seurat(n_cells = 1000, n_genes = 100)
  test_obj$mutation_tidy[1:300] <- "PARALLEL_TEST_MUT"
  test_obj$mutation_tidy[301:600] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  # Test single-threaded
  time_single <- system.time({
    result_single <- run_mast_analysis_optimized(
      mutation = "PARALLEL_TEST_MUT",
      seurat_object_path = temp_file,
      output_dir = tempdir(),
      parallel_clusters = 1,
      memory_efficient = TRUE
    )
  })
  
  # Test multi-threaded
  time_parallel <- system.time({
    result_parallel <- run_mast_analysis_optimized(
      mutation = "PARALLEL_TEST_MUT",
      seurat_object_path = temp_file,
      output_dir = tempdir(),
      parallel_clusters = 2,
      memory_efficient = TRUE
    )
  })
  
  # Validate both results
  expect_true(validate_mast_results(result_single))
  expect_true(validate_mast_results(result_parallel))
  expect_equal(length(result_single$valid_clusters), length(result_parallel$valid_clusters))
  
  # Parallel should be at least as fast (may not always be faster due to overhead)
  speedup_ratio <- time_single[["elapsed"]] / time_parallel[["elapsed"]]
  expect_true(speedup_ratio >= 0.8)  # Allow some overhead
  
  message(sprintf("Parallel processing validation: %.2fx speedup achieved", speedup_ratio))
})

test_that("memory efficiency optimizations work under stress", {
  
  skip_if_not_installed("pryr")
  
  # Create larger test dataset to stress memory management
  test_obj <- create_test_seurat(n_cells = 2000, n_genes = 200)
  test_obj$mutation_tidy[1:600] <- "MEMORY_STRESS_MUT"  
  test_obj$mutation_tidy[601:1200] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  # Monitor memory usage during analysis
  gc()
  mem_start <- pryr::mem_used()
  
  result <- run_mast_analysis_optimized(
    mutation = "MEMORY_STRESS_MUT",
    seurat_object_path = temp_file,
    output_dir = tempdir(),
    memory_efficient = TRUE,
    max_cells_per_ident = 800,  # Limit cells to control memory
    use_fast_method = TRUE      # Use fast method for large dataset
  )
  
  mem_peak <- pryr::mem_used()
  gc()
  mem_end <- pryr::mem_used()
  
  # Validate analysis succeeded
  expect_true(validate_mast_results(result))
  
  # Memory should return close to starting level (within 20%)
  memory_cleanup_ratio <- as.numeric(mem_end - mem_start) / as.numeric(mem_peak - mem_start)
  expect_true(memory_cleanup_ratio < 0.3, 
             info = sprintf("Memory cleanup ratio: %.2f (should be <0.3)", memory_cleanup_ratio))
  
  # Peak memory increase should be reasonable
  peak_increase_mb <- as.numeric(mem_peak - mem_start) / 1e6
  expect_true(peak_increase_mb < 2000,  # Less than 2GB increase
             info = sprintf("Peak memory increase: %.1f MB", peak_increase_mb))
  
  message(sprintf("Memory stress test passed: %.1f MB peak increase, %.1f%% cleanup", 
                 peak_increase_mb, (1 - memory_cleanup_ratio) * 100))
})

test_that("error handling and robustness improvements work", {
  
  # Test with problematic dataset (very few cells)
  small_obj <- create_test_seurat(n_cells = 50, n_genes = 20)
  small_obj$mutation_tidy[1:10] <- "ROBUST_TEST_MUT"
  small_obj$mutation_tidy[11:20] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(small_obj, temp_file)
  
  # This should handle the small dataset gracefully
  result <- run_mast_analysis_optimized(
    mutation = "ROBUST_TEST_MUT",
    seurat_object_path = temp_file,
    output_dir = tempdir(),
    memory_efficient = TRUE,
    use_fast_method = TRUE
  )
  
  # Should complete without errors and provide informative metadata
  expect_true(!is.null(result))
  expect_true(!is.null(result$optimization_summary))
  expect_true(result$optimization_summary$method_used %in% c("wilcox_optimized", "MAST"))
  
  # Test with non-existent file
  expect_error(
    run_mast_analysis_optimized(
      mutation = "NONEXISTENT",
      seurat_object_path = "nonexistent_file.rds",
      output_dir = tempdir()
    ),
    class = "error"
  )
  
  message("Error handling and robustness validated")
})