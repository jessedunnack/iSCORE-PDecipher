#!/usr/bin/env Rscript

#' Test Different Correlation Approaches for Gene Pair Analysis
#' 
#' This script compares correlation strength using different gene filtering approaches:
#' 1. ALL genes (current app approach) 
#' 2. Significant DE genes only (recommended)
#' 3. Top N significant genes
#' 4. Genes in significant enrichment terms

cat("=== CORRELATION APPROACH COMPARISON ===\n")

# Load required data
if (!file.exists("full_DE_results.rds")) {
  stop("ERROR: full_DE_results.rds not found. Run from correct directory.")
}

de_data <- readRDS("full_DE_results.rds")
cat("✓ Loaded DE data\n")

# Test gene pair: LRRK2_vs_LRRK2, cluster_0, experiment C12_FPD-23
mast_gene <- "LRRK2"
crispri_gene <- "LRRK2" 
test_cluster <- "cluster_0"
test_experiment <- "C12_FPD-23"

cat("\nTesting:", mast_gene, "vs", crispri_gene, "in", test_cluster, "for", test_experiment, "\n")

# Extract MAST data
mast_results <- de_data$iSCORE_PD_MAST[[mast_gene]][[test_cluster]]$results
mast_data_all <- data.frame(
  gene_name = rownames(mast_results),
  log2FC = mast_results$avg_log2FC,
  p_val_adj = mast_results$p_val_adj,
  stringsAsFactors = FALSE
)
mast_data_all <- mast_data_all[!is.na(mast_data_all$log2FC), ]

# Extract CRISPRi data  
crispri_results <- de_data$CRISPRi_Mixscale[[crispri_gene]][[test_cluster]]$results
exp_col <- paste0("log2FC_", test_experiment)
crispri_data_all <- data.frame(
  gene_name = crispri_results$gene_ID,
  log2FC = crispri_results[[exp_col]],
  stringsAsFactors = FALSE
)
crispri_data_all <- crispri_data_all[!is.na(crispri_data_all$log2FC), ]

cat("MAST genes (all):", nrow(mast_data_all), "\n")
cat("CRISPRi genes (all):", nrow(crispri_data_all), "\n")

# Helper function to calculate correlation
calc_correlation <- function(mast_df, crispri_df, label) {
  merged <- merge(mast_df[, c("gene_name", "log2FC")], 
                  crispri_df[, c("gene_name", "log2FC")],
                  by = "gene_name", suffixes = c("_mast", "_crispri"))
  merged <- merged[complete.cases(merged), ]
  
  if (nrow(merged) < 10) {
    return(data.frame(approach = label, n_genes = nrow(merged), 
                     correlation = NA, p_value = NA, note = "Too few genes"))
  }
  
  cor_test <- cor.test(merged$log2FC_mast, merged$log2FC_crispri)
  
  return(data.frame(
    approach = label,
    n_genes = nrow(merged), 
    correlation = round(cor_test$estimate, 4),
    p_value = cor_test$p.value,
    note = ifelse(cor_test$p.value < 0.05, "Significant", "Non-significant"),
    stringsAsFactors = FALSE
  ))
}

# Approach 1: ALL genes (current app approach)
result1 <- calc_correlation(mast_data_all, crispri_data_all, "1. ALL genes (current)")

# Approach 2: Significant DE genes only  
mast_sig <- mast_data_all[mast_data_all$p_val_adj < 0.05 & abs(mast_data_all$log2FC) > 0.25, ]
# For CRISPRi, we need to find the p-value column
p_col <- paste0("p_cell_type", test_experiment)
if (p_col %in% colnames(crispri_results)) {
  crispri_sig_idx <- crispri_results[[p_col]] < 0.05 & abs(crispri_results[[exp_col]]) > 0.25
  crispri_sig_idx[is.na(crispri_sig_idx)] <- FALSE
  crispri_sig <- crispri_data_all[crispri_sig_idx, ]
} else {
  # Fallback: use top 20% of genes by absolute log2FC
  crispri_sig <- crispri_data_all[abs(crispri_data_all$log2FC) > quantile(abs(crispri_data_all$log2FC), 0.8, na.rm = TRUE), ]
}

result2 <- calc_correlation(mast_sig, crispri_sig, "2. Significant DE genes only")

# Approach 3: Top 500 genes by absolute log2FC
mast_top500 <- mast_data_all[order(abs(mast_data_all$log2FC), decreasing = TRUE), ]
mast_top500 <- head(mast_top500, 500)
crispri_top500 <- crispri_data_all[order(abs(crispri_data_all$log2FC), decreasing = TRUE), ]
crispri_top500 <- head(crispri_top500, 500)

result3 <- calc_correlation(mast_top500, crispri_top500, "3. Top 500 by |log2FC|")

# Approach 4: Top 1000 genes by absolute log2FC
mast_top1000 <- head(mast_data_all[order(abs(mast_data_all$log2FC), decreasing = TRUE), ], 1000)
crispri_top1000 <- head(crispri_data_all[order(abs(crispri_data_all$log2FC), decreasing = TRUE), ], 1000)

result4 <- calc_correlation(mast_top1000, crispri_top1000, "4. Top 1000 by |log2FC|")

# Approach 5: Intersection of significant genes (both methods must be significant)
common_genes <- intersect(mast_sig$gene_name, crispri_sig$gene_name)
if (length(common_genes) > 10) {
  mast_intersect <- mast_sig[mast_sig$gene_name %in% common_genes, ]
  crispri_intersect <- crispri_sig[crispri_sig$gene_name %in% common_genes, ]
  result5 <- calc_correlation(mast_intersect, crispri_intersect, "5. Intersection significant")
} else {
  result5 <- data.frame(approach = "5. Intersection significant", n_genes = length(common_genes), 
                       correlation = NA, p_value = NA, note = "Too few overlapping genes", stringsAsFactors = FALSE)
}

# Combine results
all_results <- rbind(result1, result2, result3, result4, result5)

cat("\n=== CORRELATION COMPARISON RESULTS ===\n")
print(all_results, row.names = FALSE)

# Detailed analysis for best approaches
cat("\n=== DETAILED ANALYSIS ===\n")

best_idx <- which.max(abs(all_results$correlation))
if (length(best_idx) > 0 && !is.na(all_results$correlation[best_idx])) {
  best_approach <- all_results$approach[best_idx]
  best_cor <- all_results$correlation[best_idx]
  
  cat("BEST APPROACH:", best_approach, "\n")
  cat("Correlation:", best_cor, "\n")
  cat("Strength interpretation:", 
      ifelse(abs(best_cor) < 0.1, "Negligible",
             ifelse(abs(best_cor) < 0.3, "Weak", 
                    ifelse(abs(best_cor) < 0.5, "Moderate",
                           ifelse(abs(best_cor) < 0.7, "Strong", "Very Strong")))), "\n")
}

# Check if significant genes approach gives better results
if (!is.na(result2$correlation) && !is.na(result1$correlation)) {
  improvement <- abs(result2$correlation) / abs(result1$correlation)
  cat("\nSignificant genes vs ALL genes:\n")
  cat("Fold improvement:", round(improvement, 2), "x\n")
  cat("Absolute correlation change:", round(abs(result2$correlation) - abs(result1$correlation), 4), "\n")
}

cat("\n=== RECOMMENDATIONS ===\n")

if (any(!is.na(all_results$correlation) & abs(all_results$correlation) > 0.3)) {
  strong_approaches <- all_results[!is.na(all_results$correlation) & abs(all_results$correlation) > 0.3, ]
  cat("✓ STRONG correlations found with:\n")
  for (i in 1:nrow(strong_approaches)) {
    cat("  -", strong_approaches$approach[i], ":", round(strong_approaches$correlation[i], 3), "\n")
  }
  cat("\n✅ RECOMMENDATION: Use significant DE genes approach for much stronger correlations\n")
} else if (any(!is.na(all_results$correlation) & abs(all_results$correlation) > 0.15)) {
  cat("⚠️  Moderate correlations found - consider keeping feature with improvements\n")
  cat("✅ RECOMMENDATION: Use top N genes or significant genes approach\n")
} else {
  cat("❌ WEAK correlations across all approaches\n") 
  cat("💡 RECOMMENDATION: Consider removing correlation feature or setting lower expectations\n")
}

cat("\n=== IMPLEMENTATION SUGGESTIONS ===\n")
cat("1. ADD FILTERING OPTIONS:\n")
cat("   - 'All genes' (current - for comparison)\n")
cat("   - 'Significant DE genes only' (recommended)\n")
cat("   - 'Top N most changed genes'\n")
cat("\n2. ADD INTERACTIVE FEATURES:\n")
cat("   - Hover to see gene name, both log2FC values\n")
cat("   - Color points by significance status\n")
cat("   - Brushing/zooming for detailed exploration\n")
cat("\n3. SET USER EXPECTATIONS:\n")
cat("   - Label correlation strength (weak/moderate/strong)\n")
cat("   - Explain that cross-method correlations are typically weak\n")
cat("   - Focus on biological interpretation of significant overlaps\n")

cat("\nAnalysis complete!\n")