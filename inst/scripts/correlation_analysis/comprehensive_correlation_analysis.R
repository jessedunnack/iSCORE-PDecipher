#!/usr/bin/env Rscript

#' Comprehensive Correlation Analysis Across All Gene Pairs, Clusters, and Experiments
#' 
#' Tests different gene filtering approaches for correlation analysis:
#' 1. ALL genes (current app approach)
#' 2. Significant DE genes only  
#' 3. Top 100/200/500/1000 genes by |log2FC|
#' 4. Intersection of significant genes

cat("=== COMPREHENSIVE CORRELATION ANALYSIS ===\n")

# Load required data
if (!file.exists("full_DE_results.rds")) {
  stop("ERROR: full_DE_results.rds not found. Run from correct directory.")
}

de_data <- readRDS("full_DE_results.rds")
cat("✓ Loaded DE data\n")

# Define gene pairs to test
gene_pairs <- data.frame(
  mast_gene = c("PARK7", "LRRK2", "PRKN", "VPS13C_W395C", "PINK1", "VPS13C_A444P", 
                "SNCA_A30P", "SYNJ1", "ATP13A2", "DNAJC6", "FBXO7", "SNCA_A53T"),
  crispri_gene = c("PARK7", "LRRK2", "PARK2", "VPS13C", "PINK1", "VPS13C", 
                   "SNCA", "SYNJ1", "ATP13A2", "DNAJC6", "FBXO7", "SNCA"),
  stringsAsFactors = FALSE
)

# Define experiments and clusters to test
experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
test_clusters <- c("cluster_0", "cluster_1", "cluster_2", "cluster_4", "cluster_5")

cat("Testing", nrow(gene_pairs), "gene pairs ×", length(experiments), "experiments ×", length(test_clusters), "clusters\n")
cat("Total combinations:", nrow(gene_pairs) * length(experiments) * length(test_clusters), "\n\n")

# Helper function to calculate correlation with error handling
calc_correlation_safe <- function(mast_df, crispri_df, approach_name, gene_pair, experiment, cluster) {
  tryCatch({
    merged <- merge(mast_df[, c("gene_name", "log2FC")], 
                    crispri_df[, c("gene_name", "log2FC")],
                    by = "gene_name", suffixes = c("_mast", "_crispri"))
    merged <- merged[complete.cases(merged), ]
    
    if (nrow(merged) < 3) {
      return(data.frame(
        gene_pair = gene_pair,
        experiment = experiment, 
        cluster = cluster,
        approach = approach_name,
        n_genes = nrow(merged),
        correlation = NA,
        p_value = NA,
        abs_correlation = NA,
        significant = FALSE,
        note = "Too few overlapping genes",
        stringsAsFactors = FALSE
      ))
    }
    
    cor_test <- cor.test(merged$log2FC_mast, merged$log2FC_crispri)
    
    return(data.frame(
      gene_pair = gene_pair,
      experiment = experiment,
      cluster = cluster, 
      approach = approach_name,
      n_genes = nrow(merged),
      correlation = round(cor_test$estimate, 4),
      p_value = cor_test$p.value,
      abs_correlation = round(abs(cor_test$estimate), 4),
      significant = cor_test$p.value < 0.05,
      note = ifelse(cor_test$p.value < 0.05, "Significant", "Non-significant"),
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    return(data.frame(
      gene_pair = gene_pair,
      experiment = experiment,
      cluster = cluster,
      approach = approach_name, 
      n_genes = 0,
      correlation = NA,
      p_value = NA,
      abs_correlation = NA,
      significant = FALSE,
      note = paste("Error:", e$message),
      stringsAsFactors = FALSE
    ))
  })
}

# Function to apply different filtering approaches
test_gene_pair_combination <- function(gene_pair_row, experiment, cluster) {
  mast_gene <- gene_pair_row$mast_gene
  crispri_gene <- gene_pair_row$crispri_gene
  gene_pair_name <- paste0(mast_gene, "_vs_", crispri_gene)
  
  cat("  Testing:", gene_pair_name, "-", experiment, "-", cluster, "\n")
  
  # Check if data exists
  if (!mast_gene %in% names(de_data$iSCORE_PD_MAST) ||
      !cluster %in% names(de_data$iSCORE_PD_MAST[[mast_gene]]) ||
      !crispri_gene %in% names(de_data$CRISPRi_Mixscale) ||
      !cluster %in% names(de_data$CRISPRi_Mixscale[[crispri_gene]])) {
    return(data.frame(
      gene_pair = gene_pair_name,
      experiment = experiment,
      cluster = cluster,
      approach = "Data not available",
      n_genes = 0, correlation = NA, p_value = NA, abs_correlation = NA,
      significant = FALSE, note = "Missing data", stringsAsFactors = FALSE
    ))
  }
  
  # Extract MAST data
  mast_results <- de_data$iSCORE_PD_MAST[[mast_gene]][[cluster]]$results
  mast_data_all <- data.frame(
    gene_name = rownames(mast_results),
    log2FC = mast_results$avg_log2FC,
    p_val_adj = mast_results$p_val_adj,
    stringsAsFactors = FALSE
  )
  mast_data_all <- mast_data_all[!is.na(mast_data_all$log2FC), ]
  
  # Extract CRISPRi data
  crispri_results <- de_data$CRISPRi_Mixscale[[crispri_gene]][[cluster]]$results
  exp_col <- paste0("log2FC_", experiment)
  p_col <- paste0("p_cell_type", experiment)
  
  if (!exp_col %in% colnames(crispri_results)) {
    return(data.frame(
      gene_pair = gene_pair_name,
      experiment = experiment,
      cluster = cluster,
      approach = "Experiment not available",
      n_genes = 0, correlation = NA, p_value = NA, abs_correlation = NA,
      significant = FALSE, note = paste("Missing", experiment, "data"), stringsAsFactors = FALSE
    ))
  }
  
  crispri_data_all <- data.frame(
    gene_name = crispri_results$gene_ID,
    log2FC = crispri_results[[exp_col]],
    p_val = if(p_col %in% colnames(crispri_results)) crispri_results[[p_col]] else NA,
    stringsAsFactors = FALSE
  )
  crispri_data_all <- crispri_data_all[!is.na(crispri_data_all$log2FC), ]
  
  # Test different approaches
  results <- list()
  
  # Approach 1: ALL genes (current app approach)
  results[[1]] <- calc_correlation_safe(mast_data_all, crispri_data_all, 
                                       "1_ALL_genes", gene_pair_name, experiment, cluster)
  
  # Approach 2: Significant DE genes only
  mast_sig <- mast_data_all[!is.na(mast_data_all$p_val_adj) & 
                           mast_data_all$p_val_adj < 0.05 & 
                           abs(mast_data_all$log2FC) > 0.25, ]
  
  if (!all(is.na(crispri_data_all$p_val))) {
    crispri_sig <- crispri_data_all[!is.na(crispri_data_all$p_val) & 
                                   crispri_data_all$p_val < 0.05 & 
                                   abs(crispri_data_all$log2FC) > 0.25, ]
  } else {
    # Fallback: top 20% by absolute log2FC
    threshold <- quantile(abs(crispri_data_all$log2FC), 0.8, na.rm = TRUE)
    crispri_sig <- crispri_data_all[abs(crispri_data_all$log2FC) > threshold, ]
  }
  
  results[[2]] <- calc_correlation_safe(mast_sig, crispri_sig,
                                       "2_Significant_DE", gene_pair_name, experiment, cluster)
  
  # Approach 3-6: Top N genes by absolute log2FC
  for (i in 1:4) {
    top_n <- c(100, 200, 500, 1000)[i]
    mast_top <- head(mast_data_all[order(abs(mast_data_all$log2FC), decreasing = TRUE), ], top_n)
    crispri_top <- head(crispri_data_all[order(abs(crispri_data_all$log2FC), decreasing = TRUE), ], top_n)
    
    results[[2 + i]] <- calc_correlation_safe(mast_top, crispri_top,
                                             paste0(3 + i - 1, "_Top_", top_n), gene_pair_name, experiment, cluster)
  }
  
  # Approach 7: Intersection of significant genes
  if (nrow(mast_sig) > 0 && nrow(crispri_sig) > 0) {
    common_genes <- intersect(mast_sig$gene_name, crispri_sig$gene_name)
    if (length(common_genes) > 3) {
      mast_intersect <- mast_sig[mast_sig$gene_name %in% common_genes, ]
      crispri_intersect <- crispri_sig[crispri_sig$gene_name %in% common_genes, ]
      results[[7]] <- calc_correlation_safe(mast_intersect, crispri_intersect,
                                           "7_Intersection_Sig", gene_pair_name, experiment, cluster)
    } else {
      results[[7]] <- data.frame(
        gene_pair = gene_pair_name, experiment = experiment, cluster = cluster,
        approach = "7_Intersection_Sig", n_genes = length(common_genes),
        correlation = NA, p_value = NA, abs_correlation = NA,
        significant = FALSE, note = "Too few intersecting genes", stringsAsFactors = FALSE
      )
    }
  } else {
    results[[7]] <- data.frame(
      gene_pair = gene_pair_name, experiment = experiment, cluster = cluster,
      approach = "7_Intersection_Sig", n_genes = 0,
      correlation = NA, p_value = NA, abs_correlation = NA,
      significant = FALSE, note = "No significant genes", stringsAsFactors = FALSE
    )
  }
  
  return(do.call(rbind, results))
}

# Run comprehensive analysis
cat("Starting comprehensive analysis...\n")
all_results <- data.frame()
total_combinations <- nrow(gene_pairs) * length(experiments) * length(test_clusters)
current_combination <- 0

for (i in 1:nrow(gene_pairs)) {
  gene_pair_row <- gene_pairs[i, ]
  cat("Gene pair", i, "of", nrow(gene_pairs), ":", gene_pair_row$mast_gene, "vs", gene_pair_row$crispri_gene, "\n")
  
  for (experiment in experiments) {
    for (cluster in test_clusters) {
      current_combination <- current_combination + 1
      cat("  Progress:", current_combination, "/", total_combinations, 
          "(", round(current_combination/total_combinations*100, 1), "%)\n")
      
      combination_results <- test_gene_pair_combination(gene_pair_row, experiment, cluster)
      all_results <- rbind(all_results, combination_results)
    }
  }
}

# Save detailed results
write.csv(all_results, "comprehensive_correlation_results.csv", row.names = FALSE)
cat("✓ Detailed results saved to comprehensive_correlation_results.csv\n")

# Generate summary analysis
cat("\n=== COMPREHENSIVE ANALYSIS SUMMARY ===\n")

# Filter to valid results only
valid_results <- all_results[!is.na(all_results$correlation) & 
                            !all_results$approach %in% c("Data not available", "Experiment not available"), ]

cat("Total tests run:", nrow(all_results), "\n")
cat("Valid correlation results:", nrow(valid_results), "\n")
cat("Success rate:", round(nrow(valid_results)/nrow(all_results)*100, 1), "%\n\n")

# Compare approaches overall
cat("=== APPROACH COMPARISON (OVERALL) ===\n")
approach_summary <- aggregate(
  cbind(abs_correlation, significant) ~ approach, 
  data = valid_results, 
  FUN = function(x) c(mean = round(mean(x, na.rm = TRUE), 4), 
                      count = sum(!is.na(x)))
)

approach_stats <- data.frame(
  approach = approach_summary$approach,
  mean_abs_corr = approach_summary$abs_correlation[, "mean"],
  success_rate = round(approach_summary$significant[, "mean"] * 100, 1),
  n_tests = approach_summary$abs_correlation[, "count"],
  stringsAsFactors = FALSE
)
approach_stats <- approach_stats[order(approach_stats$mean_abs_corr, decreasing = TRUE), ]

print(approach_stats, row.names = FALSE)

# Best approach per gene pair
cat("\n=== BEST APPROACH PER GENE PAIR ===\n")
gene_pair_best <- aggregate(abs_correlation ~ gene_pair, data = valid_results, 
                           FUN = function(x) ifelse(all(is.na(x)), NA, max(x, na.rm = TRUE)))

for (gp in unique(valid_results$gene_pair)) {
  gp_data <- valid_results[valid_results$gene_pair == gp, ]
  if (nrow(gp_data) > 0 && any(!is.na(gp_data$abs_correlation))) {
    best_idx <- which.max(gp_data$abs_correlation)
    best_result <- gp_data[best_idx, ]
    cat(sprintf("%-20s: %s (r=%.3f, %s, %s, %s)\n", 
               gp, best_result$approach, best_result$correlation,
               best_result$experiment, best_result$cluster, best_result$note))
  }
}

# Experiment comparison
cat("\n=== EXPERIMENT COMPARISON ===\n")
exp_summary <- aggregate(abs_correlation ~ experiment, data = valid_results, 
                        FUN = function(x) round(mean(x, na.rm = TRUE), 4))
exp_summary <- exp_summary[order(exp_summary$abs_correlation, decreasing = TRUE), ]
print(exp_summary, row.names = FALSE)

# Strong correlations (>0.3)
cat("\n=== STRONG CORRELATIONS (|r| > 0.3) ===\n")
strong_corr <- valid_results[!is.na(valid_results$abs_correlation) & 
                            valid_results$abs_correlation > 0.3, ]
if (nrow(strong_corr) > 0) {
  strong_corr <- strong_corr[order(strong_corr$abs_correlation, decreasing = TRUE), ]
  print(strong_corr[, c("gene_pair", "experiment", "cluster", "approach", "correlation", "n_genes")], 
        row.names = FALSE)
} else {
  cat("No correlations > 0.3 found\n")
}

# Recommendations
cat("\n=== RECOMMENDATIONS ===\n")

best_approach <- approach_stats$approach[1]
best_mean_corr <- approach_stats$mean_abs_corr[1]
improvement_vs_all <- best_mean_corr / approach_stats$mean_abs_corr[approach_stats$approach == "1_ALL_genes"]

cat("BEST OVERALL APPROACH:", best_approach, "\n")
cat("Mean |correlation|:", best_mean_corr, "\n")
cat("Improvement over ALL genes:", round(improvement_vs_all, 2), "x\n")

if (best_mean_corr > 0.25) {
  cat("\n✅ STRONG RECOMMENDATION: Implement", best_approach, "filtering\n")
  cat("   - Shows consistent improvement across gene pairs\n")
  cat("   - Correlations reach moderate strength (0.25+)\n")
  cat("   - Worth implementing in UI with interactive features\n")
} else if (best_mean_corr > 0.15) {
  cat("\n⚠️  MODERATE RECOMMENDATION: Consider implementing", best_approach, "filtering\n")
  cat("   - Some improvement over current approach\n")
  cat("   - Set appropriate user expectations about correlation strength\n")
  cat("   - Focus on specific gene pairs that show stronger correlations\n")
} else {
  cat("\n❌ WEAK RECOMMENDATION: Correlations remain poor even with filtering\n")
  cat("   - Consider removing correlation plots entirely\n")
  cat("   - Focus on gene overlap and pathway overlap metrics instead\n")
  cat("   - If keeping, clearly label as 'exploratory' with low expectations\n")
}

# Gene pair recommendations
good_gene_pairs <- unique(strong_corr$gene_pair)
if (length(good_gene_pairs) > 0) {
  cat("\n🎯 GENE PAIRS WORTH HIGHLIGHTING:\n")
  for (gp in good_gene_pairs) {
    cat("  -", gp, "\n")
  }
  cat("\nConsider showing only these gene pairs in the UI for correlation analysis.\n")
}

cat("\nAnalysis complete! Check comprehensive_correlation_results.csv for detailed data.\n")