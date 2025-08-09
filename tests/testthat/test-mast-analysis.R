# Unit Tests for MAST Analysis Functions
# Created: August 2025
#
# These tests validate the core MAST differential expression analysis
# functions used for iSCORE-PD mutation analysis

test_that("run_mast_analysis creates valid output structure", {
  
  # Skip if Seurat not available
  skip_if_not_installed("Seurat")
  
  # Create test data
  test_data <- create_mast_test_data("SNCA_A53T")
  
  # Create temporary file for test Seurat object
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_data$seurat_object, temp_file)
  
  # Create temporary output directory
  temp_dir <- tempdir()
  
  # Run MAST analysis (should not crash)
  expect_no_error({
    results <- run_mast_analysis(
      mutation = test_data$mutation,
      seurat_object_path = temp_file,
      output_dir = temp_dir
    )
  })
  
  # Validate result structure
  expect_true(is.list(results))
  expect_true(validate_mast_results(results))
})

test_that("run_mast_analysis handles batch effects correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Create test data with multiple batches
  seurat_obj <- create_test_seurat(n_cells = 300)
  
  # Set up multi-batch scenario
  seurat_obj$mutation_tidy[1:100] <- "LRRK2_G2019S"
  seurat_obj$mutation_tidy[101:200] <- "eWT"
  seurat_obj$batch[1:150] <- "batch1"
  seurat_obj$batch[151:300] <- "batch2"
  
  # Ensure LRRK2 and eWT are in both batches (special case in code)
  seurat_obj$batch[seurat_obj$mutation_tidy == "LRRK2_G2019S"] <- c(rep("batch1", 50), rep("batch2", 50))
  seurat_obj$batch[seurat_obj$mutation_tidy == "eWT"] <- c(rep("batch1", 50), rep("batch2", 50))
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(seurat_obj, temp_file)
  
  # LRRK2 should allow cross-batch analysis (special case in original code)
  expect_no_error({
    results <- run_mast_analysis(
      mutation = "LRRK2_G2019S",
      seurat_object_path = temp_file,
      output_dir = tempdir()
    )
  })
  
  expect_true(is.list(results))
})

test_that("run_mast_analysis validates input parameters correctly", {
  
  # Test NULL mutation parameter
  expect_error(
    run_mast_analysis(NULL, "dummy_path.rds"),
    class = "error"
  )
  
  # Test non-existent file path
  expect_error(
    run_mast_analysis("SNCA_A53T", "non_existent_file.rds"),
    class = "error"
  )
  
  # Test invalid output directory (if we can make one)
  if (.Platform$OS.type != "windows") {  # Skip on Windows due to permission issues
    expect_error(
      run_mast_analysis("SNCA_A53T", tempfile(fileext = ".rds"), "/root/invalid_dir"),
      class = "error"
    )
  }
})

test_that("MAST analysis produces statistically valid results", {
  
  skip_if_not_installed("Seurat")
  
  # Create test data with known differential expression pattern
  test_data <- create_mast_test_data("SNCA_A53T")
  seurat_obj <- test_data$seurat_object
  
  # Artificially create differential expression for SNCA gene
  if ("SNCA" %in% rownames(seurat_obj)) {
    # Increase SNCA expression in mutation cells
    snca_cells <- which(seurat_obj$mutation_tidy == "SNCA_A53T")
    if (length(snca_cells) > 0) {
      seurat_obj@assays$RNA@data["SNCA", snca_cells] <- 
        seurat_obj@assays$RNA@data["SNCA", snca_cells] + 2
    }
  }
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(seurat_obj, temp_file)
  
  results <- run_mast_analysis(
    mutation = "SNCA_A53T",
    seurat_object_path = temp_file,
    output_dir = tempdir()
  )
  
  # Check statistical validity
  for (cluster_result in results) {
    if (is.data.frame(cluster_result) && nrow(cluster_result) > 0) {
      
      # P-values should be between 0 and 1
      expect_true(all(cluster_result$p_val >= 0 & cluster_result$p_val <= 1, na.rm = TRUE))
      expect_true(all(cluster_result$p_val_adj >= 0 & cluster_result$p_val_adj <= 1, na.rm = TRUE))
      
      # Adjusted p-values should be >= raw p-values (Benjamini-Hochberg)
      expect_true(all(cluster_result$p_val_adj >= cluster_result$p_val, na.rm = TRUE))
      
      # Percentages should be between 0 and 1
      expect_true(all(cluster_result$pct.1 >= 0 & cluster_result$pct.1 <= 1, na.rm = TRUE))
      expect_true(all(cluster_result$pct.2 >= 0 & cluster_result$pct.2 <= 1, na.rm = TRUE))
      
      # Log fold changes should be numeric and finite
      expect_true(all(is.finite(cluster_result$avg_log2FC) | is.na(cluster_result$avg_log2FC)))
    }
  }
})

test_that("MAST analysis handles edge cases gracefully", {
  
  skip_if_not_installed("Seurat")
  
  # Test with very few cells
  small_obj <- create_edge_case_data("single_cell")
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(small_obj, temp_file)
  
  # Should handle gracefully (may warn but shouldn't crash)
  expect_no_error({
    results <- run_mast_analysis(
      mutation = "eWT",  # Use a mutation that exists
      seurat_object_path = temp_file,
      output_dir = tempdir()
    )
  })
  
  # Test with empty cluster scenario
  empty_cluster_obj <- create_edge_case_data("empty_cluster")
  temp_file2 <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file2), add = TRUE)
  saveRDS(empty_cluster_obj, temp_file2)
  
  expect_no_error({
    results <- run_mast_analysis(
      mutation = "eWT",
      seurat_object_path = temp_file2,
      output_dir = tempdir()
    )
  })
})

test_that("MAST analysis clustering resolution is correct", {
  
  skip_if_not_installed("Seurat")
  
  test_data <- create_mast_test_data("PARK7")
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_data$seurat_object, temp_file)
  
  results <- run_mast_analysis(
    mutation = "PARK7",
    seurat_object_path = temp_file,
    output_dir = tempdir()
  )
  
  # Should use resolution 0.2 as specified in the original code
  # Check that we get reasonable number of clusters (not too many or too few)
  n_clusters <- length(results)
  expect_true(n_clusters >= 1)
  expect_true(n_clusters <= 20)  # Reasonable upper bound for 500 test cells
})

test_that("MAST analysis mutation detection works correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Create test data with specific mutations
  seurat_obj <- create_test_seurat(n_cells = 200)
  
  # Set specific mutation pattern
  mutations <- c("SNCA_A53T", "LRRK2_G2019S", "eWT", "PRKN", "PARK7")
  seurat_obj$mutation_tidy <- sample(mutations, ncol(seurat_obj), replace = TRUE,
                                   prob = c(0.3, 0.2, 0.3, 0.1, 0.1))
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(seurat_obj, temp_file)
  
  # Test that each mutation can be analyzed
  for (mut in c("SNCA_A53T", "LRRK2_G2019S", "PRKN")) {
    expect_no_error({
      results <- run_mast_analysis(
        mutation = mut,
        seurat_object_path = temp_file,
        output_dir = tempdir()
      )
    })
    
    expect_true(is.list(results))
  }
})

test_that("MAST analysis memory usage is reasonable", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pryr")
  
  test_data <- create_mast_test_data("SNCA_A53T")
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_data$seurat_object, temp_file)
  
  # Monitor memory usage
  memory_test <- check_memory_usage(
    run_mast_analysis,
    mutation = "SNCA_A53T",
    seurat_object_path = temp_file,
    output_dir = tempdir()
  )
  
  # Memory increase should be reasonable for test data
  expect_true(memory_test$memory_used < 1e9)  # Less than 1GB for test data
  expect_true(validate_mast_results(memory_test$result))
})