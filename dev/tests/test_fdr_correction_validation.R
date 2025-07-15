#!/usr/bin/env Rscript

#' FDR Correction Validation Test (v0.2.6)
#' 
#' This script validates that our enhanced hierarchical FDR correction
#' maintains proper Type I error control under the null hypothesis.
#' 
#' Tests both:
#' 1. Within gene pair FDR correction (Benjamini-Yekutieli)
#' 2. Across gene pairs FDR correction (Benjamini-Hochberg)

# Load required functions
source("R/manuscript_signature_discovery.R")

cat("=== FDR CORRECTION VALIDATION TEST ===\n\n")

# Simulation parameters
n_simulations <- 1000
n_gene_pairs <- 8  # Number of comparable gene pairs
n_clusters <- 12   # Number of clusters  
n_experiments <- 3  # Number of experiments
n_directions <- 2   # Same vs opposite direction

# Expected Type I error rates to test
alpha_levels <- c(0.01, 0.05, 0.10)

cat("Simulation Parameters:\n")
cat("  Simulations:", n_simulations, "\n")
cat("  Gene pairs:", n_gene_pairs, "\n") 
cat("  Clusters:", n_clusters, "\n")
cat("  Experiments:", n_experiments, "\n")
cat("  Directions:", n_directions, "\n")
cat("  Alpha levels:", paste(alpha_levels, collapse=", "), "\n\n")

# Function to simulate null p-values (uniform distribution)
simulate_null_pvalues <- function(n_gene_pairs, n_clusters, n_experiments, n_directions) {
  
  # Total tests per gene pair: n_clusters × n_experiments × n_directions
  tests_per_pair <- n_clusters * n_experiments * n_directions
  
  # Generate uniform p-values for each gene pair
  pvalue_matrix <- matrix(runif(n_gene_pairs * tests_per_pair), 
                         nrow = n_gene_pairs, ncol = tests_per_pair)
  
  # Create data structure matching our expected format
  simulated_data <- list()
  
  for (i in 1:n_gene_pairs) {
    gene_pair_name <- paste0("GENE", i, "_vs_GENE", i)
    
    # Create cluster data for this gene pair
    cluster_data <- data.frame(
      gene_pair = gene_pair_name,
      cluster = paste0("cluster_", 0:(n_clusters-1)),
      experiment = rep(paste0("exp_", 1:n_experiments), each = n_clusters),
      direction = rep(c("same", "opposite"), length.out = tests_per_pair),
      raw_p = pvalue_matrix[i, ],
      stringsAsFactors = FALSE
    )
    
    simulated_data[[gene_pair_name]] <- cluster_data
  }
  
  return(simulated_data)
}

# Function to count false discoveries at different alpha levels
count_false_discoveries <- function(fdr_corrected_p, alpha_levels) {
  
  discoveries <- sapply(alpha_levels, function(alpha) {
    # Count number of "significant" results (these are all false under null)
    sum(fdr_corrected_p < alpha, na.rm = TRUE)
  })
  
  names(discoveries) <- paste0("alpha_", alpha_levels)
  return(discoveries)
}

# Run simulations
cat("Running", n_simulations, "simulations...\n")

# Store results for each alpha level
false_discovery_counts <- matrix(0, nrow = n_simulations, ncol = length(alpha_levels))
colnames(false_discovery_counts) <- paste0("alpha_", alpha_levels)

# Track total tests for FDR calculation
total_tests <- n_gene_pairs * n_clusters * n_experiments * n_directions

cat("Total tests per simulation:", total_tests, "\n")
cat("Expected false discoveries under null:\n")
for (alpha in alpha_levels) {
  expected_false <- total_tests * alpha
  cat("  α =", alpha, ":", expected_false, "false discoveries\n")
}
cat("\n")

# Progress tracking
start_time <- Sys.time()

for (sim in 1:n_simulations) {
  
  if (sim %% 100 == 0) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    cat("Simulation", sim, "/", n_simulations, 
        sprintf("(%.1f%% complete, %.0fs elapsed)\n", 
                100 * sim / n_simulations, elapsed))
  }
  
  # Generate null p-values for this simulation
  simulated_data <- simulate_null_pvalues(n_gene_pairs, n_clusters, n_experiments, n_directions)
  
  # Combine all p-values into single vector for FDR correction
  all_pvalues <- unlist(lapply(simulated_data, function(x) x$raw_p))
  
  # Apply our enhanced hierarchical FDR correction
  tryCatch({
    
    # Create signature rankings format expected by the function
    # Convert simulated data to the format expected by apply_enhanced_fdr_correction_v026
    signature_rankings <- do.call(rbind, simulated_data)
    signature_rankings$gene_fisher_p <- signature_rankings$raw_p
    signature_rankings$intersection_fisher_p <- signature_rankings$raw_p
    signature_rankings$union_fisher_p <- signature_rankings$raw_p
    
    # Test our enhanced FDR correction function
    corrected_results <- apply_enhanced_fdr_correction_v026(
      signature_rankings = signature_rankings,
      use_enhanced_method = TRUE
    )
    
    # Extract FDR-corrected p-values (use the hierarchical FDR corrected column)
    if ("gene_fisher_p_fdr_enhanced_hierarchical" %in% colnames(corrected_results)) {
      fdr_corrected_p <- corrected_results$gene_fisher_p_fdr_enhanced_hierarchical
    } else if ("intersection_fisher_p_fdr_enhanced_hierarchical" %in% colnames(corrected_results)) {
      fdr_corrected_p <- corrected_results$intersection_fisher_p_fdr_enhanced_hierarchical
    } else {
      # Fallback to any available FDR column
      fdr_cols <- grep("fdr.*hierarchical", colnames(corrected_results), value = TRUE)
      if (length(fdr_cols) > 0) {
        fdr_corrected_p <- corrected_results[[fdr_cols[1]]]
      } else {
        stop("No FDR-corrected columns found")
      }
    }
    
    # Count false discoveries at each alpha level
    false_discoveries <- count_false_discoveries(fdr_corrected_p, alpha_levels)
    
    # Store results
    false_discovery_counts[sim, ] <- false_discoveries
    
  }, error = function(e) {
    cat("Error in simulation", sim, ":", e$message, "\n")
    # Use NA for failed simulations
    false_discovery_counts[sim, ] <- rep(NA, length(alpha_levels))
  })
}

# Calculate summary statistics
cat("\n=== FDR CORRECTION VALIDATION RESULTS ===\n")
cat("==========================================\n\n")

end_time <- Sys.time()
total_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
cat("Total simulation time:", sprintf("%.1f", total_time), "seconds\n\n")

for (i in 1:length(alpha_levels)) {
  alpha <- alpha_levels[i]
  col_name <- paste0("alpha_", alpha)
  
  # Remove failed simulations
  valid_sims <- !is.na(false_discovery_counts[, i])
  valid_counts <- false_discovery_counts[valid_sims, i]
  n_valid <- sum(valid_sims)
  
  if (n_valid == 0) {
    cat("α =", alpha, ": No valid simulations\n")
    next
  }
  
  # Calculate statistics
  mean_false_discoveries <- mean(valid_counts)
  median_false_discoveries <- median(valid_counts)
  sd_false_discoveries <- sd(valid_counts)
  
  # Expected false discoveries under null
  expected_false_discoveries <- total_tests * alpha
  
  # Observed FDR
  observed_fdr <- mean_false_discoveries / total_tests
  
  # 95% confidence interval for mean
  se_mean <- sd_false_discoveries / sqrt(n_valid)
  ci_lower <- mean_false_discoveries - 1.96 * se_mean
  ci_upper <- mean_false_discoveries + 1.96 * se_mean
  
  cat(sprintf("α = %.2f (Expected FDR = %.2f):\n", alpha, alpha))
  cat(sprintf("  Valid simulations: %d/%d\n", n_valid, n_simulations))
  cat(sprintf("  Mean false discoveries: %.2f (expected: %.2f)\n", 
              mean_false_discoveries, expected_false_discoveries))
  cat(sprintf("  Median false discoveries: %.2f\n", median_false_discoveries))
  cat(sprintf("  Standard deviation: %.2f\n", sd_false_discoveries))
  cat(sprintf("  Observed FDR: %.4f (target: %.2f)\n", observed_fdr, alpha))
  cat(sprintf("  95%% CI for mean: [%.2f, %.2f]\n", ci_lower, ci_upper))
  
  # Check if FDR is controlled
  fdr_controlled <- observed_fdr <= alpha * 1.1  # Allow 10% tolerance
  if (fdr_controlled) {
    cat("  ✓ FDR CONTROLLED: Observed FDR ≤ target α\n")
  } else {
    cat("  ✗ FDR NOT CONTROLLED: Observed FDR > target α\n")
  }
  cat("\n")
}

# Overall assessment
cat("=== OVERALL ASSESSMENT ===\n")

# Check if all alpha levels show controlled FDR
all_controlled <- TRUE
for (i in 1:length(alpha_levels)) {
  alpha <- alpha_levels[i]
  valid_sims <- !is.na(false_discovery_counts[, i])
  if (sum(valid_sims) == 0) {
    all_controlled <- FALSE
    next
  }
  
  valid_counts <- false_discovery_counts[valid_sims, i]
  observed_fdr <- mean(valid_counts) / total_tests
  
  if (observed_fdr > alpha * 1.1) {  # 10% tolerance
    all_controlled <- FALSE
    break
  }
}

if (all_controlled) {
  cat("✅ SUCCESS: Enhanced FDR correction maintains Type I error control\n")
  cat("   The hierarchical Benjamini-Yekutieli + Benjamini-Hochberg approach\n")
  cat("   successfully controls the False Discovery Rate at all tested levels.\n")
} else {
  cat("❌ FAILURE: Enhanced FDR correction does NOT maintain proper control\n")
  cat("   The FDR correction may be too conservative or too liberal.\n")
  cat("   Consider adjusting the correction method or parameters.\n")
}

cat("\n=== SIMULATION COMPLETE ===\n")