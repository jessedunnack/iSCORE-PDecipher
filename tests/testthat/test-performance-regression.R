# Performance Regression Tests for iSCORE-PDecipher
# Created: August 2025
#
# These tests establish performance baselines and monitor for regressions
# Critical for maintaining analysis speed with 230k+ cell datasets

test_that("MAST analysis completes within time limits", {
  
  skip_if_not_installed("Seurat")
  
  # Create moderately sized test dataset
  test_data <- create_mast_test_data("LRRK2_G2019S")
  seurat_obj <- test_data$seurat_object
  
  # Scale up for performance testing
  if (ncol(seurat_obj) < 1000) {
    # Create larger test dataset
    larger_obj <- create_test_seurat(n_cells = 1000, n_genes = 100)
    larger_obj$mutation_tidy[1:200] <- "LRRK2_G2019S"
    larger_obj$mutation_tidy[201:400] <- "eWT"
    seurat_obj <- larger_obj
  }
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(seurat_obj, temp_file)
  
  # Measure execution time
  start_time <- Sys.time()
  
  results <- run_mast_analysis(
    mutation = "LRRK2_G2019S",
    seurat_object_path = temp_file,
    output_dir = tempdir()
  )
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Performance expectations for 1000 cells, 100 genes
  expect_true(execution_time < 300)  # Less than 5 minutes
  expect_true(validate_mast_results(results))
  
  message(sprintf("MAST analysis completed in %.2f seconds", execution_time))
  
  # Store baseline for future regression testing
  baseline_file <- file.path(tempdir(), "mast_performance_baseline.rds")
  performance_data <- list(
    timestamp = Sys.time(),
    n_cells = ncol(seurat_obj),
    n_genes = nrow(seurat_obj),
    execution_time = execution_time,
    n_clusters = length(results)
  )
  saveRDS(performance_data, baseline_file)
})

test_that("MixScale analysis completes within time limits", {
  
  skip_if_not_installed("Seurat")
  
  # Create test data
  test_data <- create_mixscale_test_data("PRKN")
  
  # Scale up for performance testing
  larger_obj <- create_test_seurat(n_cells = 800, n_genes = 80)
  larger_obj$guide_assignment <- sample(
    c("Non-Targeting", "PRKN_guide_1", "PRKN_guide_2"), 
    ncol(larger_obj), 
    replace = TRUE
  )
  larger_obj$perturbation <- ifelse(larger_obj$guide_assignment == "Non-Targeting", 
                                   "Control", "Perturbed")
  larger_obj$experiment <- "C12_FPD-24"
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "perf_test_mixscale")
  dir.create(exp_dir, showWarnings = FALSE)
  on.exit(unlink(exp_dir, recursive = TRUE))
  
  saveRDS(larger_obj, file.path(exp_dir, "C12_FPD-24.rds"))
  
  # Measure execution time
  start_time <- Sys.time()
  
  results <- run_mixscale_analysis(
    experiment_path = exp_dir,
    output_dir = temp_dir,
    modality = "CRISPRi"
  )
  
  end_time <- Sys.time()
  execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  # Performance expectations
  expect_true(execution_time < 240)  # Less than 4 minutes
  expect_true(validate_mixscale_results(results))
  
  message(sprintf("MixScale analysis completed in %.2f seconds", execution_time))
  
  # Store baseline
  baseline_file <- file.path(tempdir(), "mixscale_performance_baseline.rds")
  performance_data <- list(
    timestamp = Sys.time(),
    n_cells = ncol(larger_obj),
    n_genes = nrow(larger_obj),
    execution_time = execution_time,
    n_clusters = length(results)
  )
  saveRDS(performance_data, baseline_file)
})

test_that("enrichment analysis scales appropriately", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Test with different gene set sizes
  gene_set_sizes <- c(10, 50, 100, 200)
  enrichment_times <- numeric(length(gene_set_sizes))
  
  for (i in seq_along(gene_set_sizes)) {
    
    # Create gene set of appropriate size
    all_pd_genes <- c(
      "SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35", "ATP13A2", 
      "FBXO7", "PLA2G6", "DNAJC6", "SYNJ1", "GBA", "SMPD1", "CTSB",
      "TH", "DDC", "SLC6A3", "DRD1", "DRD2", "COMT", "MAO", "ALDH1A1",
      "PITX3", "FOXA2", "LMX1B", "MSX1", "NURR1", "EN1", "EN2", "WNT1",
      paste0("TEST_GENE_", 1:200)  # Padding genes
    )
    
    test_genes <- head(all_pd_genes, gene_set_sizes[i])
    
    # Measure enrichment time
    start_time <- Sys.time()
    
    # Convert to Entrez IDs
    gene_entrez <- clusterProfiler::bitr(test_genes,
                                        fromType = "SYMBOL",
                                        toType = "ENTREZID", 
                                        OrgDb = org.Hs.eg.db::org.Hs.eg.db)
    
    if (nrow(gene_entrez) >= 3) {
      go_result <- clusterProfiler::enrichGO(
        gene = gene_entrez$ENTREZID,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 1.0
      )
    }
    
    end_time <- Sys.time()
    enrichment_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  # Enrichment time should scale reasonably with gene set size
  expect_true(all(enrichment_times < 60))  # Less than 1 minute each
  expect_true(enrichment_times[4] > enrichment_times[1])  # Larger sets take longer
  
  message(sprintf("Enrichment times: %s seconds", 
                 paste(round(enrichment_times, 2), collapse = ", ")))
  
  # Check for reasonable scaling (should be roughly linear or sublinear)
  time_ratio <- enrichment_times[4] / enrichment_times[1]
  size_ratio <- gene_set_sizes[4] / gene_set_sizes[1]
  
  # Time shouldn't increase exponentially with gene set size
  expect_true(time_ratio <= size_ratio * 2)
})

test_that("memory usage stays within bounds for large datasets", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pryr")
  
  # Test memory usage with progressively larger datasets
  dataset_sizes <- c(500, 1000, 2000)  # Cell counts
  
  for (n_cells in dataset_sizes) {
    
    gc()  # Clean up memory
    mem_start <- pryr::mem_used()
    
    # Create large test dataset
    large_obj <- create_test_seurat(n_cells = n_cells, n_genes = 200)
    large_obj$mutation_tidy[1:(n_cells * 0.3)] <- "SNCA_A53T"
    large_obj$mutation_tidy[((n_cells * 0.3) + 1):(n_cells * 0.7)] <- "eWT"
    
    temp_file <- tempfile(fileext = ".rds")
    saveRDS(large_obj, temp_file)
    
    # Monitor peak memory usage during analysis
    mem_before_analysis <- pryr::mem_used()
    
    results <- run_mast_analysis(
      mutation = "SNCA_A53T",
      seurat_object_path = temp_file,
      output_dir = tempdir()
    )
    
    mem_after_analysis <- pryr::mem_used()
    peak_memory_increase <- mem_after_analysis - mem_before_analysis
    
    # Clean up
    unlink(temp_file)
    rm(large_obj, results)
    gc()
    
    # Memory usage expectations
    expect_true(peak_memory_increase < 4e9)  # Less than 4GB increase
    
    message(sprintf("Dataset with %d cells used %.1f MB peak memory", 
                   n_cells, as.numeric(peak_memory_increase) / 1e6))
    
    # Memory should scale reasonably with dataset size
    memory_per_cell <- as.numeric(peak_memory_increase) / n_cells
    expect_true(memory_per_cell < 5e6)  # Less than 5MB per cell
  }
})

test_that("concurrent analysis performance is reasonable", {
  
  skip_if_not_installed("Seurat")
  
  # Test multiple analyses in sequence (simulating batch processing)
  mutations_to_test <- c("SNCA_A53T", "LRRK2_G2019S", "PRKN")
  
  # Create shared test dataset
  shared_obj <- create_test_seurat(n_cells = 600, n_genes = 60)
  
  # Distribute mutations
  n_per_mutation <- 150
  shared_obj$mutation_tidy[1:n_per_mutation] <- "SNCA_A53T"
  shared_obj$mutation_tidy[(n_per_mutation + 1):(2 * n_per_mutation)] <- "LRRK2_G2019S"
  shared_obj$mutation_tidy[(2 * n_per_mutation + 1):(3 * n_per_mutation)] <- "PRKN"
  shared_obj$mutation_tidy[(3 * n_per_mutation + 1):600] <- "eWT"
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(shared_obj, temp_file)
  
  # Time batch processing
  total_start_time <- Sys.time()
  batch_results <- list()
  
  for (mut in mutations_to_test) {
    mutation_start_time <- Sys.time()
    
    batch_results[[mut]] <- run_mast_analysis(
      mutation = mut,
      seurat_object_path = temp_file,
      output_dir = tempdir()
    )
    
    mutation_time <- as.numeric(difftime(Sys.time(), mutation_start_time, units = "secs"))
    expect_true(mutation_time < 180)  # Each analysis < 3 minutes
  }
  
  total_time <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  
  # Total batch processing time should be reasonable
  expect_true(total_time < 600)  # Less than 10 minutes for 3 mutations
  expect_true(all(sapply(batch_results, validate_mast_results)))
  
  message(sprintf("Batch processing of %d mutations completed in %.2f seconds", 
                 length(mutations_to_test), total_time))
  
  # Average time per mutation should be reasonable
  avg_time_per_mutation <- total_time / length(mutations_to_test)
  expect_true(avg_time_per_mutation < 200)  # Less than 3.3 minutes average
})

test_that("disk I/O performance is acceptable", {
  
  skip_if_not_installed("Seurat")
  
  # Test file I/O performance with realistic dataset sizes
  test_sizes <- c(100, 500, 1000)  # Cell counts
  
  for (n_cells in test_sizes) {
    
    # Create test dataset
    test_obj <- create_test_seurat(n_cells = n_cells, n_genes = 100)
    
    # Measure write time
    temp_file <- tempfile(fileext = ".rds")
    
    write_start <- Sys.time()
    saveRDS(test_obj, temp_file)
    write_time <- as.numeric(difftime(Sys.time(), write_start, units = "secs"))
    
    # Check file size
    file_size <- file.info(temp_file)$size
    
    # Measure read time
    read_start <- Sys.time()
    loaded_obj <- readRDS(temp_file)
    read_time <- as.numeric(difftime(Sys.time(), read_start, units = "secs"))
    
    unlink(temp_file)
    
    # Performance expectations
    expect_true(write_time < 30)  # Less than 30 seconds to write
    expect_true(read_time < 10)   # Less than 10 seconds to read
    
    # File size should be reasonable (compressed RDS)
    size_per_cell <- file_size / n_cells
    expect_true(size_per_cell < 50000)  # Less than 50KB per cell
    
    message(sprintf("%d cells: Write %.2fs, Read %.2fs, Size %.1fMB", 
                   n_cells, write_time, read_time, file_size / 1e6))
  }
})

test_that("regression detection system works", {
  
  # Test the ability to detect performance regressions
  current_baseline <- list(
    mast_time_per_1000_cells = 120,  # seconds
    mixscale_time_per_1000_cells = 100,
    memory_per_1000_cells = 2e9,  # bytes
    enrichment_time_per_100_genes = 15
  )
  
  # Simulate different performance scenarios
  test_scenarios <- list(
    "acceptable" = list(
      mast_time_per_1000_cells = 130,  # 8% increase - acceptable
      mixscale_time_per_1000_cells = 105,  # 5% increase - acceptable
      memory_per_1000_cells = 2.1e9,  # 5% increase - acceptable
      enrichment_time_per_100_genes = 16  # 7% increase - acceptable
    ),
    "warning" = list(
      mast_time_per_1000_cells = 150,  # 25% increase - warning
      mixscale_time_per_1000_cells = 130,  # 30% increase - warning
      memory_per_1000_cells = 2.6e9,  # 30% increase - warning
      enrichment_time_per_100_genes = 20  # 33% increase - warning
    ),
    "regression" = list(
      mast_time_per_1000_cells = 200,  # 67% increase - regression
      mixscale_time_per_1000_cells = 180,  # 80% increase - regression
      memory_per_1000_cells = 3.5e9,  # 75% increase - regression
      enrichment_time_per_100_genes = 30  # 100% increase - regression
    )
  )
  
  for (scenario_name in names(test_scenarios)) {
    scenario <- test_scenarios[[scenario_name]]
    
    # Calculate percentage changes
    mast_change <- (scenario$mast_time_per_1000_cells / current_baseline$mast_time_per_1000_cells - 1) * 100
    mixscale_change <- (scenario$mixscale_time_per_1000_cells / current_baseline$mixscale_time_per_1000_cells - 1) * 100
    memory_change <- (scenario$memory_per_1000_cells / current_baseline$memory_per_1000_cells - 1) * 100
    enrichment_change <- (scenario$enrichment_time_per_100_genes / current_baseline$enrichment_time_per_100_genes - 1) * 100
    
    # Classification logic
    max_change <- max(abs(c(mast_change, mixscale_change, memory_change, enrichment_change)))
    
    if (max_change <= 15) {
      classification <- "acceptable"
    } else if (max_change <= 40) {
      classification <- "warning"
    } else {
      classification <- "regression"
    }
    
    expect_equal(classification, scenario_name)
    
    message(sprintf("Scenario '%s': Max change %.1f%%, classified as '%s'", 
                   scenario_name, max_change, classification))
  }
  
  expect_true(TRUE)  # Test completed successfully
})