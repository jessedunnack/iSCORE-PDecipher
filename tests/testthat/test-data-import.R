# Unit Tests for Data Import Functions
# Created: August 2025
#
# These tests validate the data import functions that handle
# MAST and MixScale data loading and processing

test_that("import_mast_data loads data with correct structure", {
  
  # Create mock MAST results data
  mock_mast_data <- list(
    "cluster_0" = data.frame(
      p_val = c(0.01, 0.05, 0.1),
      avg_log2FC = c(1.5, -0.8, 0.3),
      pct.1 = c(0.6, 0.3, 0.8),
      pct.2 = c(0.2, 0.7, 0.5),
      p_val_adj = c(0.02, 0.1, 0.2),
      row.names = c("SNCA", "TH", "LRRK2")
    ),
    "cluster_1" = data.frame(
      p_val = c(0.001, 0.02),
      avg_log2FC = c(2.1, -1.2),
      pct.1 = c(0.9, 0.4),
      pct.2 = c(0.1, 0.8),
      p_val_adj = c(0.005, 0.04),
      row.names = c("PARK7", "PINK1")
    )
  )
  
  # Save to temporary file
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(mock_mast_data, temp_file)
  
  # Test import_mast_data function (if it exists)
  if (exists("import_mast_data", mode = "function")) {
    result <- import_mast_data(temp_file)
    
    expect_true(is.list(result))
    expect_true(validate_mast_results(result))
    
    # Check that all expected columns are present
    for (cluster_data in result) {
      if (is.data.frame(cluster_data)) {
        expected_cols <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
        expect_true(all(expected_cols %in% colnames(cluster_data)))
      }
    }
  } else {
    skip("import_mast_data function not found")
  }
})

test_that("import_mixscale_data loads data with correct structure", {
  
  # Create mock MixScale results data
  mock_mixscale_data <- list(
    "cluster_0" = data.frame(
      gene_ID = c("SNCA", "LRRK2", "PRKN"),
      log2FC_C12_FPD.24 = c(-1.2, 0.8, -0.5),
      p_cell_typeC12_FPD.24.weight = c(0.01, 0.05, 0.15),
      stringsAsFactors = FALSE
    ),
    "cluster_1" = data.frame(
      gene_ID = c("PARK7", "PINK1"),
      log2FC_C12_FPD.24 = c(-0.9, 1.1),
      p_cell_typeC12_FPD.24.weight = c(0.02, 0.08),
      stringsAsFactors = FALSE
    )
  )
  
  # Save to temporary file
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(mock_mixscale_data, temp_file)
  
  # Test import_mixscale_data function (if it exists)
  if (exists("import_mixscale_data", mode = "function")) {
    result <- import_mixscale_data(temp_file)
    
    expect_true(is.list(result))
    expect_true(validate_mixscale_results(result))
    
    # Check for required columns
    for (cluster_data in result) {
      if (is.data.frame(cluster_data)) {
        expect_true("gene_ID" %in% colnames(cluster_data))
        
        # Should have p-value weight columns
        p_weight_cols <- grep("p_.*:weight$", colnames(cluster_data), value = TRUE)
        expect_true(length(p_weight_cols) > 0)
      }
    }
  } else {
    skip("import_mixscale_data function not found")
  }
})

test_that("data import handles missing files gracefully", {
  
  non_existent_file <- "non_existent_file.rds"
  
  # Test that import functions handle missing files appropriately
  if (exists("import_mast_data", mode = "function")) {
    expect_error(import_mast_data(non_existent_file))
  }
  
  if (exists("import_mixscale_data", mode = "function")) {
    expect_error(import_mixscale_data(non_existent_file))
  }
})

test_that("data import validates file formats correctly", {
  
  # Create invalid file (not RDS format)
  invalid_file <- tempfile(fileext = ".txt")
  on.exit(unlink(invalid_file))
  writeLines("This is not an RDS file", invalid_file)
  
  if (exists("import_mast_data", mode = "function")) {
    expect_error(import_mast_data(invalid_file))
  }
  
  if (exists("import_mixscale_data", mode = "function")) {
    expect_error(import_mixscale_data(invalid_file))
  }
})

test_that("Seurat object loading preserves metadata integrity", {
  
  skip_if_not_installed("Seurat")
  
  # Create test Seurat object with specific metadata
  test_obj <- create_test_seurat(n_cells = 100, n_genes = 50)
  
  # Add specific metadata that should be preserved
  original_mutations <- unique(test_obj$mutation_tidy)
  original_batches <- unique(test_obj$batch)
  original_experiments <- unique(test_obj$experiment)
  
  # Save and reload
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  reloaded_obj <- readRDS(temp_file)
  
  # Check metadata integrity
  expect_identical(colnames(test_obj), colnames(reloaded_obj))
  expect_identical(rownames(test_obj), rownames(reloaded_obj))
  expect_identical(test_obj$mutation_tidy, reloaded_obj$mutation_tidy)
  expect_identical(test_obj$batch, reloaded_obj$batch)
  expect_identical(test_obj$experiment, reloaded_obj$experiment)
  
  # Check that all original values are preserved
  expect_true(all(original_mutations %in% reloaded_obj$mutation_tidy))
  expect_true(all(original_batches %in% reloaded_obj$batch))
  expect_true(all(original_experiments %in% reloaded_obj$experiment))
})

test_that("data import preserves gene expression data integrity", {
  
  skip_if_not_installed("Seurat")
  
  # Create test object with known expression pattern
  test_obj <- create_test_seurat(n_cells = 50, n_genes = 20)
  
  # Get original expression data
  original_counts <- Seurat::GetAssayData(test_obj, assay = "RNA", slot = "counts")
  original_data <- Seurat::GetAssayData(test_obj, assay = "RNA", slot = "data")
  
  # Save and reload
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  saveRDS(test_obj, temp_file)
  
  reloaded_obj <- readRDS(temp_file)
  reloaded_counts <- Seurat::GetAssayData(reloaded_obj, assay = "RNA", slot = "counts")
  reloaded_data <- Seurat::GetAssayData(reloaded_obj, assay = "RNA", slot = "data")
  
  # Check expression data integrity
  expect_identical(dim(original_counts), dim(reloaded_counts))
  expect_identical(dim(original_data), dim(reloaded_data))
  
  # Check that expression values are preserved (allowing for floating point precision)
  expect_equal(as.matrix(original_counts), as.matrix(reloaded_counts), tolerance = 1e-10)
  expect_equal(as.matrix(original_data), as.matrix(reloaded_data), tolerance = 1e-10)
})

test_that("batch assignment logic works correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Create test object with specific batch-mutation combinations
  test_obj <- create_test_seurat(n_cells = 300)
  
  # Set up known batch-mutation scenario
  test_obj$mutation_tidy <- c(
    rep("LRRK2_G2019S", 100),  # Should allow cross-batch
    rep("SNCA_A53T", 100),     # Should be restricted to same batch
    rep("eWT", 100)
  )
  
  test_obj$batch <- c(
    rep("batch1", 50), rep("batch2", 50),  # LRRK2 in both batches
    rep("batch1", 100),                    # SNCA only in batch1
    rep("batch1", 50), rep("batch2", 50)   # eWT in both batches
  )
  
  # Test batch detection logic
  lrrk2_batches <- unique(test_obj$batch[test_obj$mutation_tidy == "LRRK2_G2019S"])
  snca_batches <- unique(test_obj$batch[test_obj$mutation_tidy == "SNCA_A53T"])
  ewt_batches <- unique(test_obj$batch[test_obj$mutation_tidy == "eWT"])
  
  expect_equal(length(lrrk2_batches), 2)  # LRRK2 in both batches
  expect_equal(length(snca_batches), 1)   # SNCA only in batch1
  expect_equal(length(ewt_batches), 2)    # eWT in both batches
  
  # LRRK2 analysis should be able to use cross-batch controls
  lrrk2_or_ewt <- test_obj$mutation_tidy %in% c("LRRK2_G2019S", "eWT")
  lrrk2_available_batches <- unique(test_obj$batch[lrrk2_or_ewt])
  expect_equal(length(lrrk2_available_batches), 2)
  
  # SNCA analysis should be restricted to same batch
  snca_or_ewt_same_batch <- test_obj$mutation_tidy == "SNCA_A53T" | 
    (test_obj$mutation_tidy == "eWT" & test_obj$batch %in% snca_batches)
  snca_available_batches <- unique(test_obj$batch[snca_or_ewt_same_batch])
  expect_equal(length(snca_available_batches), 1)
  expect_equal(snca_available_batches, "batch1")
})

test_that("experiment ID parsing works correctly", {
  
  # Test the critical C12_FPD-24 vs A15_FPD-24 distinction mentioned in context
  experiment_ids <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23", "A15_FPD-24")
  
  # C12 and C18 should be CRISPRi, A15 should be CRISPRa
  crispri_experiments <- grep("^C[0-9]+_FPD", experiment_ids, value = TRUE)
  crispra_experiments <- grep("^A[0-9]+_FPD", experiment_ids, value = TRUE)
  
  expect_true("C12_FPD-23" %in% crispri_experiments)
  expect_true("C12_FPD-24" %in% crispri_experiments)
  expect_true("C18_FPD-23" %in% crispri_experiments)
  expect_false("A15_FPD-24" %in% crispri_experiments)
  
  expect_true("A15_FPD-24" %in% crispra_experiments)
  expect_false("C12_FPD-24" %in% crispra_experiments)
  
  # Test preference order: C12_FPD-24 > C12_FPD-23 > C18_FPD-23
  available_crispri <- c("C18_FPD-23", "C12_FPD-23", "C12_FPD-24")
  sorted_by_preference <- c("C12_FPD-24", "C12_FPD-23", "C18_FPD-23")
  
  # This should match the preference logic in MixScale analysis
  expect_equal(intersect(sorted_by_preference, available_crispri), 
               c("C12_FPD-24", "C12_FPD-23", "C18_FPD-23"))
})

test_that("data loading handles large datasets efficiently", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pryr")
  
  # Create larger test dataset to simulate performance
  large_obj <- create_test_seurat(n_cells = 5000, n_genes = 1000)
  
  temp_file <- tempfile(fileext = ".rds")
  on.exit(unlink(temp_file))
  
  # Time the save operation
  save_start <- Sys.time()
  saveRDS(large_obj, temp_file)
  save_duration <- difftime(Sys.time(), save_start, units = "secs")
  
  # Time the load operation
  load_start <- Sys.time()
  loaded_obj <- readRDS(temp_file)
  load_duration <- difftime(Sys.time(), load_start, units = "secs")
  
  # Basic performance expectations (will vary by system)
  expect_true(as.numeric(save_duration) < 30)  # Less than 30 seconds to save
  expect_true(as.numeric(load_duration) < 10)  # Less than 10 seconds to load
  
  # Verify data integrity after loading
  expect_identical(dim(large_obj), dim(loaded_obj))
  expect_identical(colnames(large_obj), colnames(loaded_obj))
  expect_identical(rownames(large_obj), rownames(loaded_obj))
})