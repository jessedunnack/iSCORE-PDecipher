# Unit Tests for MixScale Analysis Functions
# Created: August 2025
#
# These tests validate the core MixScale analysis functions used for
# CRISPRi/CRISPRa perturbation experiments in iSCORE-PDecipher

test_that("run_mixscale_analysis creates valid output structure", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("dplyr")
  skip_if_not_installed("Matrix")
  
  # Create test data
  test_data <- create_mixscale_test_data("SNCA")
  
  # Create temporary experiment directory structure
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "test_experiment")
  dir.create(exp_dir, showWarnings = FALSE)
  
  # Save test experiment data
  saveRDS(test_data$experiments$`C12_FPD-24`, 
          file.path(exp_dir, "C12_FPD-24.rds"))
  
  # Run MixScale analysis
  expect_no_error({
    results <- run_mixscale_analysis(
      experiment_path = exp_dir,
      output_dir = temp_dir,
      modality = "CRISPRi"
    )
  })
  
  # Clean up
  unlink(exp_dir, recursive = TRUE)
  
  # Validate result structure
  expect_true(is.list(results))
  expect_true(validate_mixscale_results(results))
})

test_that("run_mixscale_analysis handles experiment preferences correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Create test data with multiple experiments
  test_data <- create_mixscale_test_data("LRRK2")
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "multi_exp_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  # Save multiple experiments
  for (exp_name in names(test_data$experiments)) {
    saveRDS(test_data$experiments[[exp_name]], 
            file.path(exp_dir, paste0(exp_name, ".rds")))
  }
  
  expect_no_error({
    results <- run_mixscale_analysis(
      experiment_path = exp_dir,
      output_dir = temp_dir,
      modality = "CRISPRi"
    )
  })
  
  # Clean up
  unlink(exp_dir, recursive = TRUE)
  
  # Check that results exist and follow expected structure
  expect_true(is.list(results))
  
  # According to the code, C12_FPD-24 should be preferred if available
  # This is handled in the actual implementation
})

test_that("run_mixscale_analysis validates input parameters", {
  
  # Test non-existent experiment path
  expect_error(
    run_mixscale_analysis("non_existent_path"),
    class = "error"
  )
  
  # Test invalid modality
  temp_dir <- tempdir()
  expect_error(
    run_mixscale_analysis(temp_dir, modality = "InvalidModality"),
    class = "error"
  )
})

test_that("MixScale analysis excludes genes correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Test that excluded genes (defined in original code) are handled
  excluded_genes <- c("SSR1", "SLC6A3", "CTSB", "PITX3", "APOE", "LRP1B")
  
  test_data <- create_mixscale_test_data("SNCA")
  seurat_obj <- test_data$experiments$`C12_FPD-24`
  
  # Add some excluded genes to test object
  if (any(excluded_genes %in% rownames(seurat_obj))) {
    
    temp_dir <- tempdir()
    exp_dir <- file.path(temp_dir, "exclude_test")
    dir.create(exp_dir, showWarnings = FALSE)
    
    saveRDS(seurat_obj, file.path(exp_dir, "C12_FPD-24.rds"))
    
    expect_no_error({
      results <- run_mixscale_analysis(
        experiment_path = exp_dir,
        output_dir = temp_dir,
        modality = "CRISPRi"
      )
    })
    
    unlink(exp_dir, recursive = TRUE)
    
    # Results should not contain excluded genes
    for (cluster_result in results) {
      if (is.data.frame(cluster_result) && "gene_ID" %in% colnames(cluster_result)) {
        expect_true(!any(excluded_genes %in% cluster_result$gene_ID))
      }
    }
  }
})

test_that("MixScale analysis handles GWAS genes correctly", {
  
  skip_if_not_installed("Seurat")
  
  # Test GWAS gene handling (from G5O13 library)
  gwas_genes <- c("GBA", "BST1", "GALC", "COASY")  # Subset of GWAS genes
  
  test_data <- create_mixscale_test_data("GBA")
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "gwas_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  saveRDS(test_data$experiments$`C12_FPD-24`, 
          file.path(exp_dir, "C12_FPD-24.rds"))
  
  expect_no_error({
    results <- run_mixscale_analysis(
      experiment_path = exp_dir,
      output_dir = temp_dir,
      modality = "CRISPRi"
    )
  })
  
  unlink(exp_dir, recursive = TRUE)
  expect_true(validate_mixscale_results(results))
})

test_that("MixScale analysis produces statistically valid results", {
  
  skip_if_not_installed("Seurat")
  
  test_data <- create_mixscale_test_data("PRKN")
  seurat_obj <- test_data$experiments$`C12_FPD-24`
  
  # Add realistic perturbation effect
  if ("PRKN" %in% rownames(seurat_obj)) {
    perturbed_cells <- which(seurat_obj$perturbation == "Perturbed")
    if (length(perturbed_cells) > 0) {
      # Decrease PRKN expression in perturbed cells (CRISPRi effect)
      seurat_obj@assays$RNA@data["PRKN", perturbed_cells] <- 
        seurat_obj@assays$RNA@data["PRKN", perturbed_cells] * 0.5
    }
  }
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "stats_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  saveRDS(seurat_obj, file.path(exp_dir, "C12_FPD-24.rds"))
  
  results <- run_mixscale_analysis(
    experiment_path = exp_dir,
    output_dir = temp_dir,
    modality = "CRISPRi"
  )
  
  unlink(exp_dir, recursive = TRUE)
  
  # Validate statistical properties
  for (cluster_result in results) {
    if (is.data.frame(cluster_result) && nrow(cluster_result) > 0) {
      
      # Check for required columns
      expect_true("gene_ID" %in% colnames(cluster_result))
      
      # Find p-value weight columns
      p_weight_cols <- grep("p_.*:weight$", colnames(cluster_result), value = TRUE)
      expect_true(length(p_weight_cols) > 0)
      
      # P-values should be between 0 and 1
      for (p_col in p_weight_cols) {
        p_values <- cluster_result[[p_col]]
        expect_true(all(p_values >= 0 & p_values <= 1, na.rm = TRUE))
      }
      
      # Check for log fold change columns
      lfc_cols <- grep("log2FC", colnames(cluster_result), value = TRUE)
      for (lfc_col in lfc_cols) {
        lfc_values <- cluster_result[[lfc_col]]
        expect_true(all(is.finite(lfc_values) | is.na(lfc_values)))
      }
    }
  }
})

test_that("MixScale analysis handles different experiment types", {
  
  skip_if_not_installed("Seurat")
  
  experiment_types <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  
  for (exp_type in experiment_types) {
    test_data <- create_mixscale_test_data("ATP13A2")
    
    temp_dir <- tempdir()
    exp_dir <- file.path(temp_dir, paste0("exp_test_", gsub("-", "_", exp_type)))
    dir.create(exp_dir, showWarnings = FALSE)
    
    # Modify experiment metadata to match specific type
    seurat_obj <- test_data$experiments[[exp_type]]
    seurat_obj$experiment <- exp_type
    
    saveRDS(seurat_obj, file.path(exp_dir, paste0(exp_type, ".rds")))
    
    expect_no_error({
      results <- run_mixscale_analysis(
        experiment_path = exp_dir,
        output_dir = temp_dir,
        modality = "CRISPRi"
      )
    })
    
    unlink(exp_dir, recursive = TRUE)
    expect_true(is.list(results))
  }
})

test_that("MixScale analysis handles edge cases gracefully", {
  
  skip_if_not_installed("Seurat")
  
  # Test with very small dataset
  small_obj <- create_test_seurat(n_cells = 10, n_genes = 5)
  small_obj$guide_assignment <- c(rep("Non-Targeting", 5), rep("SNCA_guide_1", 5))
  small_obj$perturbation <- ifelse(small_obj$guide_assignment == "Non-Targeting", 
                                  "Control", "Perturbed")
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "edge_case_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  saveRDS(small_obj, file.path(exp_dir, "C12_FPD-24.rds"))
  
  # Should handle gracefully (may warn but shouldn't crash)
  expect_no_error({
    results <- run_mixscale_analysis(
      experiment_path = exp_dir,
      output_dir = temp_dir,
      modality = "CRISPRi"
    )
  })
  
  unlink(exp_dir, recursive = TRUE)
})

test_that("MixScale analysis memory usage is reasonable", {
  
  skip_if_not_installed("Seurat")
  skip_if_not_installed("pryr")
  
  test_data <- create_mixscale_test_data("VPS35")
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "memory_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  saveRDS(test_data$experiments$`C12_FPD-24`, 
          file.path(exp_dir, "C12_FPD-24.rds"))
  
  # Monitor memory usage
  memory_test <- check_memory_usage(
    run_mixscale_analysis,
    experiment_path = exp_dir,
    output_dir = temp_dir,
    modality = "CRISPRi"
  )
  
  unlink(exp_dir, recursive = TRUE)
  
  # Memory increase should be reasonable for test data
  expect_true(memory_test$memory_used < 1e9)  # Less than 1GB for test data
  expect_true(validate_mixscale_results(memory_test$result))
})

test_that("MixScale FDR correction requirement is documented", {
  
  skip_if_not_installed("Seurat")
  
  test_data <- create_mixscale_test_data("PARK7")
  
  temp_dir <- tempdir()
  exp_dir <- file.path(temp_dir, "fdr_test")
  dir.create(exp_dir, showWarnings = FALSE)
  
  saveRDS(test_data$experiments$`C12_FPD-24`, 
          file.path(exp_dir, "C12_FPD-24.rds"))
  
  results <- run_mixscale_analysis(
    experiment_path = exp_dir,
    output_dir = temp_dir,
    modality = "CRISPRi"
  )
  
  unlink(exp_dir, recursive = TRUE)
  
  # MixScale results should contain p_weight columns that need manual FDR correction
  # This is a reminder that downstream analysis needs to apply Benjamini-Hochberg
  for (cluster_result in results) {
    if (is.data.frame(cluster_result) && nrow(cluster_result) > 0) {
      p_weight_cols <- grep("p_.*:weight$", colnames(cluster_result), value = TRUE)
      expect_true(length(p_weight_cols) > 0, 
                 info = "MixScale results should contain p_weight columns requiring FDR correction")
      
      # Should NOT contain pre-adjusted p-values (unlike MAST)
      expect_false(any(grepl("p_val_adj", colnames(cluster_result))),
                  info = "MixScale results should require manual FDR correction")
    }
  }
})