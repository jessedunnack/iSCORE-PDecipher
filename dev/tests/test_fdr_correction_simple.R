#!/usr/bin/env Rscript

#' Simple FDR Correction Validation Test (v0.2.6)
#' 
#' This script performs a focused test of our enhanced hierarchical FDR correction
#' to verify Type I error control with a smaller, more manageable simulation.

# Load required functions
source("R/manuscript_signature_discovery.R")

cat("=== SIMPLE FDR CORRECTION VALIDATION ===\n\n")

# Reduced simulation parameters for faster testing
n_simulations <- 100
n_tests <- 100  # Total tests per simulation

# Test different alpha levels
alpha_levels <- c(0.01, 0.05, 0.10)

cat("Simulation Parameters:\n")
cat("  Simulations:", n_simulations, "\n")
cat("  Tests per simulation:", n_tests, "\n")
cat("  Alpha levels:", paste(alpha_levels, collapse=", "), "\n\n")

# Function to create minimal test data
create_test_data <- function(n_tests) {
  
  # Generate uniform p-values under null hypothesis
  raw_p <- runif(n_tests)
  
  # Create minimal signature rankings data frame
  signature_rankings <- data.frame(
    gene_pair = paste0("pair_", 1:n_tests),
    cluster = paste0("cluster_", sample(0:5, n_tests, replace = TRUE)),
    gene_fisher_p = raw_p,
    intersection_fisher_p = raw_p,
    union_fisher_p = raw_p,
    stringsAsFactors = FALSE
  )
  
  return(signature_rankings)
}

# Function to test a specific correction method
test_fdr_method <- function(method_name, correction_func, p_values, alpha) {
  
  # Apply correction
  corrected_p <- correction_func(p_values)
  
  # Count significant results (false discoveries under null)
  n_significant <- sum(corrected_p < alpha, na.rm = TRUE)
  
  # Calculate observed FDR
  total_tests <- length(p_values)
  observed_fdr <- n_significant / total_tests
  
  return(list(
    method = method_name,
    n_significant = n_significant,
    total_tests = total_tests,
    observed_fdr = observed_fdr,
    target_fdr = alpha
  ))
}

# Test different correction methods
cat("Testing FDR correction methods...\n\n")

# Store results
results_summary <- list()

for (alpha in alpha_levels) {
  cat("Testing α =", alpha, "\n")
  
  # Storage for this alpha level
  method_results <- list()
  
  for (sim in 1:n_simulations) {
    
    # Generate test data
    test_data <- create_test_data(n_tests)
    
    # Test 1: Standard Benjamini-Hochberg
    bh_corrected <- p.adjust(test_data$gene_fisher_p, method = "BH")
    bh_result <- test_fdr_method("Benjamini-Hochberg", 
                                function(p) p.adjust(p, method = "BH"),
                                test_data$gene_fisher_p, alpha)
    
    # Test 2: Benjamini-Yekutieli (more conservative)
    by_corrected <- p.adjust(test_data$gene_fisher_p, method = "BY")
    by_result <- test_fdr_method("Benjamini-Yekutieli",
                                function(p) p.adjust(p, method = "BY"),
                                test_data$gene_fisher_p, alpha)
    
    # Test 3: Our enhanced hierarchical method (suppress output)
    enhanced_result <- tryCatch({
      
      # Suppress verbose output
      options(verbose = FALSE)
      capture.output({
        corrected_data <- apply_enhanced_fdr_correction_v026(
          signature_rankings = test_data,
          use_enhanced_method = TRUE
        )
      })
      
      # Extract hierarchical FDR corrected p-values
      if ("gene_fisher_p_fdr_enhanced_hierarchical" %in% colnames(corrected_data)) {
        enhanced_p <- corrected_data$gene_fisher_p_fdr_enhanced_hierarchical
      } else {
        # Try intersection approach
        enhanced_p <- corrected_data$intersection_fisher_p_fdr_enhanced_hierarchical
      }
      
      test_fdr_method("Enhanced Hierarchical",
                     function(p) enhanced_p,
                     test_data$gene_fisher_p, alpha)
      
    }, error = function(e) {
      list(method = "Enhanced Hierarchical", n_significant = NA, 
           total_tests = n_tests, observed_fdr = NA, target_fdr = alpha)
    })
    
    # Store results
    if (sim == 1) {
      method_results[["BH"]] <- list()
      method_results[["BY"]] <- list() 
      method_results[["Enhanced"]] <- list()
    }
    
    method_results[["BH"]][[sim]] <- bh_result
    method_results[["BY"]][[sim]] <- by_result
    method_results[["Enhanced"]][[sim]] <- enhanced_result
  }
  
  results_summary[[paste0("alpha_", alpha)]] <- method_results
}

# Analyze results
cat("\n=== RESULTS SUMMARY ===\n")
cat("========================\n\n")

for (alpha in alpha_levels) {
  cat(sprintf("α = %.2f:\n", alpha))
  cat("--------\n")
  
  alpha_key <- paste0("alpha_", alpha)
  alpha_results <- results_summary[[alpha_key]]
  
  for (method in names(alpha_results)) {
    method_data <- alpha_results[[method]]
    
    # Extract observed FDRs
    observed_fdrs <- sapply(method_data, function(x) x$observed_fdr)
    valid_fdrs <- observed_fdrs[!is.na(observed_fdrs)]
    
    if (length(valid_fdrs) == 0) {
      cat(sprintf("  %s: No valid results\n", method))
      next
    }
    
    # Calculate summary statistics
    mean_fdr <- mean(valid_fdrs)
    median_fdr <- median(valid_fdrs)
    sd_fdr <- sd(valid_fdrs)
    
    # Check if FDR is controlled (with 20% tolerance for small simulations)
    fdr_controlled <- mean_fdr <= alpha * 1.2
    
    cat(sprintf("  %s:\n", method))
    cat(sprintf("    Valid sims: %d/%d\n", length(valid_fdrs), n_simulations))
    cat(sprintf("    Mean FDR: %.4f (target: %.2f)\n", mean_fdr, alpha))
    cat(sprintf("    Median FDR: %.4f\n", median_fdr))
    cat(sprintf("    SD FDR: %.4f\n", sd_fdr))
    
    if (fdr_controlled) {
      cat("    ✓ FDR CONTROLLED\n")
    } else {
      cat("    ✗ FDR INFLATED\n")
    }
    cat("\n")
  }
  cat("\n")
}

# Overall assessment
cat("=== OVERALL ASSESSMENT ===\n")

# Check each method across all alpha levels
methods_controlled <- c()

for (method in c("BH", "BY", "Enhanced")) {
  method_controlled <- TRUE
  
  for (alpha in alpha_levels) {
    alpha_key <- paste0("alpha_", alpha)
    if (alpha_key %in% names(results_summary)) {
      method_data <- results_summary[[alpha_key]][[method]]
      
      if (is.null(method_data)) {
        method_controlled <- FALSE
        next
      }
      
      observed_fdrs <- sapply(method_data, function(x) x$observed_fdr)
      valid_fdrs <- observed_fdrs[!is.na(observed_fdrs)]
      
      if (length(valid_fdrs) == 0) {
        method_controlled <- FALSE
        next
      }
      
      mean_fdr <- mean(valid_fdrs)
      if (mean_fdr > alpha * 1.2) {  # 20% tolerance
        method_controlled <- FALSE
        break
      }
    }
  }
  
  methods_controlled <- c(methods_controlled, method_controlled)
  names(methods_controlled)[length(methods_controlled)] <- method
}

cat("Method Assessment:\n")
for (method in names(methods_controlled)) {
  status <- if (methods_controlled[method]) "✓ CONTROLLED" else "✗ INFLATED"
  cat(sprintf("  %s: %s\n", method, status))
}

if (methods_controlled["Enhanced"]) {
  cat("\n✅ SUCCESS: Enhanced hierarchical FDR correction maintains Type I error control\n")
} else {
  cat("\n❌ FAILURE: Enhanced hierarchical FDR correction shows inflated Type I error\n")
}

cat("\n=== TEST COMPLETE ===\n")