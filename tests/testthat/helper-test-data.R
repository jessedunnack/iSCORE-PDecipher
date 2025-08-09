# Test Helper Functions for iSCORE-PDecipher
# Created: August 2025
#
# This file contains helper functions for creating realistic test data
# that mimics the patterns found in actual iSCORE-PD datasets

#' Create a minimal Seurat object for testing
#'
#' @param n_cells Number of cells to simulate (default: 1000)
#' @param n_genes Number of genes to simulate (default: 100)
#' @param n_clusters Number of clusters to create (default: 5)
#' @return A Seurat object with realistic metadata
create_test_seurat <- function(n_cells = 1000, n_genes = 100, n_clusters = 5) {
  
  # Skip if Seurat not available (for CRAN checks)
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    skip("Seurat not available")
  }
  
  # Create sparse expression matrix with realistic patterns
  set.seed(12345)  # Reproducible results
  
  # Simulate sparse single-cell expression data
  # Most entries are 0, some cells have moderate expression
  expression_matrix <- Matrix::rsparsematrix(
    nrow = n_genes, 
    ncol = n_cells,
    density = 0.1,  # 10% of entries are non-zero (realistic for scRNA-seq)
    nnz = n_genes * n_cells * 0.1
  )
  
  # Make expression values realistic (log-normalized counts)
  expression_matrix@x <- abs(expression_matrix@x) * 3  # Scale to realistic range
  
  # Create gene names matching real iSCORE-PD patterns
  # Ensure no duplicates and no underscores (Seurat requirement)
  base_genes <- c(
    # Key PD genes that should be in tests
    "SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35", 
    "ATP13A2", "FBXO7", "PLA2G6", "DNAJC6", "SYNJ1",
    # Common marker genes
    "TH", "DAT", "VMAT2", "AADC", "PITX3", "FOXA2", "LMX1B",
    # Neuronal markers
    "RBFOX3", "MAP2", "TUBB3", "SYP", "SYN1"
  )
  
  # Generate additional unique gene names without underscores
  additional_genes <- paste0("GENE", sprintf("%04d", 1:n_genes))
  all_genes <- c(base_genes, additional_genes)
  
  # Take unique genes and ensure we have exactly n_genes
  gene_names <- unique(all_genes)[1:n_genes]
  
  # Ensure no duplicates
  if (length(unique(gene_names)) != length(gene_names)) {
    gene_names <- make.unique(gene_names)
  }
  
  rownames(expression_matrix) <- gene_names
  colnames(expression_matrix) <- paste0("CELL_", sprintf("%04d", 1:n_cells))
  
  # Create Seurat object
  seurat_obj <- Seurat::CreateSeuratObject(
    counts = expression_matrix,
    project = "iSCORE_PD_Test"
  )
  
  # Add realistic metadata matching iSCORE-PD structure
  mutations <- c("eWT", "SNCA_A53T", "LRRK2_G2019S", "PRKN", "PARK7")
  batches <- c("batch1", "batch2", "batch3")
  subclones <- paste0("subclone_", 1:10)
  
  # Assign metadata with realistic distributions
  seurat_obj$mutation_tidy <- sample(mutations, n_cells, replace = TRUE, 
                                   prob = c(0.4, 0.2, 0.2, 0.1, 0.1))
  seurat_obj$batch <- sample(batches, n_cells, replace = TRUE)
  seurat_obj$subclone_ID <- sample(subclones, n_cells, replace = TRUE)
  
  # Add experiment-specific metadata
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  seurat_obj$experiment <- sample(experiments, n_cells, replace = TRUE)
  
  # Perform basic processing
  seurat_obj <- Seurat::NormalizeData(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::FindVariableFeatures(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::ScaleData(seurat_obj, verbose = FALSE)
  
  # FIX: Adjust PCA dimensions based on available features
  # PCA can compute at most min(n_genes, n_cells) - 1 components
  n_pcs <- min(50, n_genes - 1, n_cells - 1)
  
  # Only run PCA if we have enough dimensions
  if (n_pcs >= 2) {
    seurat_obj <- Seurat::RunPCA(seurat_obj, npcs = n_pcs, verbose = FALSE)
    
    # Use available PCs for neighbors (max 10 or available)
    n_dims <- min(10, n_pcs)
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, dims = 1:n_dims, verbose = FALSE)
  } else {
    # For very small datasets, skip PCA and use raw features
    warning("Dataset too small for PCA, using raw features for clustering")
    seurat_obj <- Seurat::FindNeighbors(seurat_obj, features = rownames(seurat_obj)[1:min(10, nrow(seurat_obj))], verbose = FALSE)
  }
  seurat_obj <- Seurat::FindClusters(seurat_obj, resolution = 0.2, verbose = FALSE)
  
  return(seurat_obj)
}

#' Create MAST-specific test data
#'
#' @param mutation Character string specifying mutation to test
#' @return List with Seurat object and expected analysis parameters
create_mast_test_data <- function(mutation = "SNCA_A53T") {
  
  # Create base Seurat object
  seurat_obj <- create_test_seurat(n_cells = 500, n_genes = 50)
  
  # Ensure mutation is present in reasonable numbers
  n_mut_cells <- 100
  n_wt_cells <- 300
  
  # Set specific mutation pattern
  cell_assignments <- c(
    rep(mutation, n_mut_cells),
    rep("eWT", n_wt_cells),
    rep("OTHER_MUT", ncol(seurat_obj) - n_mut_cells - n_wt_cells)
  )
  
  seurat_obj$mutation_tidy <- sample(cell_assignments)
  
  # Ensure both mutation and eWT are in same batch for valid comparison
  seurat_obj$batch[seurat_obj$mutation_tidy %in% c(mutation, "eWT")] <- "batch1"
  
  return(list(
    seurat_object = seurat_obj,
    mutation = mutation,
    expected_clusters = length(unique(Seurat::Idents(seurat_obj))),
    expected_latent_vars = c("subclone_ID")
  ))
}

#' Create MixScale-specific test data
#'
#' @param gene Character string specifying gene to test
#' @return List with experiment data structure
create_mixscale_test_data <- function(gene = "SNCA") {
  
  # Create realistic MixScale experimental structure
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  
  test_data <- list()
  
  for (exp in experiments) {
    # Create Seurat object for this experiment
    seurat_obj <- create_test_seurat(n_cells = 200, n_genes = 30)
    
    # Add perturbation-specific metadata
    guide_types <- c("Non-Targeting", paste0(gene, "_guide_", 1:3))
    seurat_obj$guide_assignment <- sample(guide_types, ncol(seurat_obj), replace = TRUE)
    seurat_obj$experiment <- exp
    seurat_obj$perturbation <- ifelse(seurat_obj$guide_assignment == "Non-Targeting", 
                                    "Control", "Perturbed")
    
    test_data[[exp]] <- seurat_obj
  }
  
  return(list(
    experiments = test_data,
    target_gene = gene,
    expected_experiments = experiments,
    modality = "CRISPRi"
  ))
}

#' Create edge case test data
#'
#' @param case_type Type of edge case ("empty_cluster", "single_cell", "missing_genes")
#' @return Appropriate test data for edge case
create_edge_case_data <- function(case_type = "empty_cluster") {
  
  base_obj <- create_test_seurat(n_cells = 100, n_genes = 20)
  
  switch(case_type,
    "empty_cluster" = {
      # Create situation where one cluster has no cells of interest
      base_obj$mutation_tidy[Seurat::Idents(base_obj) == "0"] <- "RARE_MUTATION"
      base_obj
    },
    "single_cell" = {
      # Create cluster with single cell - use minimum viable numbers
      single_cell_obj <- create_test_seurat(n_cells = 10, n_genes = 10)
      single_cell_obj
    },
    "missing_genes" = {
      # Create object missing key genes
      genes_to_remove <- c("SNCA", "TH", "LRRK2")
      remaining_genes <- setdiff(rownames(base_obj), genes_to_remove)
      subset(base_obj, features = remaining_genes)
    },
    {
      base_obj
    }
  )
}

#' Validate MAST analysis results structure
#'
#' @param results MAST analysis results to validate
#' @return Boolean indicating if structure is valid
validate_mast_results <- function(results) {
  
  # Check if results is a list
  if (!is.list(results)) return(FALSE)
  
  # Check for required columns in differential expression results
  required_columns <- c("p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  
  for (cluster_result in results) {
    if (is.data.frame(cluster_result)) {
      if (!all(required_columns %in% colnames(cluster_result))) {
        return(FALSE)
      }
      
      # Check that p-values are between 0 and 1
      if (any(cluster_result$p_val < 0 | cluster_result$p_val > 1, na.rm = TRUE)) {
        return(FALSE)
      }
      
      # Check that percentages are between 0 and 1
      if (any(cluster_result$pct.1 < 0 | cluster_result$pct.1 > 1, na.rm = TRUE)) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}

#' Validate MixScale analysis results structure  
#'
#' @param results MixScale analysis results to validate
#' @return Boolean indicating if structure is valid
validate_mixscale_results <- function(results) {
  
  if (!is.list(results)) return(FALSE)
  
  # MixScale results have different structure - check for p_weight columns
  for (cluster_result in results) {
    if (is.data.frame(cluster_result)) {
      
      # Check for gene_ID column
      if (!"gene_ID" %in% colnames(cluster_result)) {
        return(FALSE)
      }
      
      # Check for p-value weight columns (experiment-specific)
      p_weight_cols <- grep("p_.*:weight$", colnames(cluster_result), value = TRUE)
      if (length(p_weight_cols) == 0) {
        return(FALSE)
      }
    }
  }
  
  return(TRUE)
}

#' Check memory usage during test
#'
#' @param test_function Function to test
#' @param ... Arguments to pass to test_function
#' @return List with result and memory usage
check_memory_usage <- function(test_function, ...) {
  
  gc()  # Clean up before test
  mem_before <- pryr::mem_used()
  
  result <- test_function(...)
  
  mem_after <- pryr::mem_used()
  mem_increase <- mem_after - mem_before
  
  return(list(
    result = result,
    memory_used = mem_increase,
    memory_before = mem_before,
    memory_after = mem_after
  ))
}