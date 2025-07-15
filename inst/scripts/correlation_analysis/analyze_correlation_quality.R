#!/usr/bin/env Rscript

#' Comprehensive Correlation Quality Analysis
#' 
#' This script analyzes correlation strength across all gene pairs, experiments, 
#' and clusters to identify which correlations are worth showing.

# Load required libraries (using base R to avoid dependency issues)
# library(dplyr)
# library(readr)

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# Source the correlation calculation function
if (file.exists("R/signature_analysis.R")) {
  source("R/signature_analysis.R")
} else {
  cat("ERROR: Cannot find R/signature_analysis.R - run from package root directory\n")
  quit(save = "no", status = 1)
}

cat("=== CORRELATION QUALITY ANALYSIS ===\n")
cat("Loading data and analyzing correlation patterns...\n\n")

# Load DE data (same as app uses)
if (!file.exists("full_DE_results.rds")) {
  cat("ERROR: full_DE_results.rds not found. Run the main analysis first.\n")
  quit(save = "no", status = 1)
}

de_data <- readRDS("full_DE_results.rds")
cat("✓ Loaded DE data with", length(de_data), "analysis types\n")

# Define gene pairs to analyze (same as app)
gene_pairs <- data.frame(
  mast_gene = c("PARK7", "LRRK2", "PRKN", "VPS13C_W395C", "PINK1", "VPS13C_A444P", 
                "SNCA_A30P", "SYNJ1", "ATP13A2", "DNAJC6", "FBXO7", "SNCA_A53T"),
  crispri_gene = c("PARK7", "LRRK2", "PARK2", "VPS13C", "PINK1", "VPS13C", 
                   "SNCA", "SYNJ1", "ATP13A2", "DNAJC6", "FBXO7", "SNCA"),
  stringsAsFactors = FALSE
)

# Define experiments to test
experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")

# Helper function to calculate correlation for specific conditions
calculate_correlation_stats <- function(mast_gene, crispri_gene, experiment = NULL, cluster = NULL) {
  
  tryCatch({
    # Get MAST data
    if (!mast_gene %in% names(de_data$iSCORE_PD_MAST)) {
      return(data.frame(correlation = NA, p_value = NA, n_genes = 0, error = "MAST gene not found"))
    }
    
    # Get CRISPRi data  
    if (!crispri_gene %in% names(de_data$CRISPRi_Mixscale)) {
      return(data.frame(correlation = NA, p_value = NA, n_genes = 0, error = "CRISPRi gene not found"))
    }
    
    mast_data <- de_data$iSCORE_PD_MAST[[mast_gene]]
    crispri_data <- de_data$CRISPRi_Mixscale[[crispri_gene]]
    
    # Filter by cluster if specified
    if (!is.null(cluster)) {
      if (cluster %in% names(mast_data)) {
        mast_data <- list(mast_data[[cluster]])
        names(mast_data) <- cluster
      } else {
        return(data.frame(correlation = NA, p_value = NA, n_genes = 0, error = "MAST cluster not found"))
      }
      
      if (cluster %in% names(crispri_data)) {
        crispri_data <- list(crispri_data[[cluster]])
        names(crispri_data) <- cluster
      } else {
        return(data.frame(correlation = NA, p_value = NA, n_genes = 0, error = "CRISPRi cluster not found"))
      }
    }
    
    # Combine data across clusters (or use single cluster)
    combined_mast <- data.frame()
    combined_crispri <- data.frame()
    
    for (clust in names(mast_data)) {
      if (clust %in% names(crispri_data)) {
        mast_df <- mast_data[[clust]]$results
        crispri_df <- crispri_data[[clust]]$results
        
        # Filter by experiment if specified (CRISPRi uses gene_ID column)
        if (!is.null(experiment)) {
          exp_col <- paste0("log2FC_", experiment)
          if (exp_col %in% colnames(crispri_df)) {
            crispri_df <- crispri_df[, c("gene_ID", exp_col)]
            colnames(crispri_df) <- c("gene_name", "log2FC")
          } else {
            next
          }
        } else {
          # Use first available experiment column
          lfc_cols <- grep("^log2FC_C", colnames(crispri_df), value = TRUE)
          if (length(lfc_cols) > 0) {
            crispri_df <- crispri_df[, c("gene_ID", lfc_cols[1])]
            colnames(crispri_df) <- c("gene_name", "log2FC")
          } else {
            next
          }
        }
        
        # Prepare MAST data (genes are rownames, not in a column)
        if ("avg_log2FC" %in% colnames(mast_df)) {
          mast_df$gene_name <- rownames(mast_df)
          mast_df <- mast_df[, c("gene_name", "avg_log2FC")]
          colnames(mast_df)[2] <- "log2FC"
        } else {
          next
        }
        
        combined_mast <- rbind(combined_mast, mast_df)
        combined_crispri <- rbind(combined_crispri, crispri_df)
      }
    }
    
    # Find overlapping genes
    overlap_genes <- intersect(combined_mast$gene_name, combined_crispri$gene_name)
    
    if (length(overlap_genes) < 10) {
      return(data.frame(correlation = NA, p_value = NA, n_genes = length(overlap_genes), error = "Too few overlapping genes"))
    }
    
    # Merge data for correlation
    mast_subset <- combined_mast[combined_mast$gene_name %in% overlap_genes, ]
    crispri_subset <- combined_crispri[combined_crispri$gene_name %in% overlap_genes, ]
    
    merged_data <- merge(mast_subset, crispri_subset, by = "gene_name", suffixes = c("_mast", "_crispri"))
    merged_data <- merged_data[complete.cases(merged_data), ]
    
    if (nrow(merged_data) < 10) {
      return(data.frame(correlation = NA, p_value = NA, n_genes = nrow(merged_data), error = "Too few complete cases"))
    }
    
    # Calculate correlation
    cor_test <- cor.test(merged_data$log2FC_mast, merged_data$log2FC_crispri, method = "pearson")
    
    return(data.frame(
      correlation = cor_test$estimate,
      p_value = cor_test$p.value,
      n_genes = nrow(merged_data),
      error = NA
    ))
    
  }, error = function(e) {
    return(data.frame(correlation = NA, p_value = NA, n_genes = 0, error = as.character(e$message)))
  })
}

# Run comprehensive analysis
cat("Analyzing correlations...\n")
results <- data.frame()

for (i in 1:nrow(gene_pairs)) {
  mast_gene <- gene_pairs$mast_gene[i]
  crispri_gene <- gene_pairs$crispri_gene[i]
  gene_pair <- paste0(mast_gene, "_vs_", crispri_gene)
  
  cat("Processing", gene_pair, "...\n")
  
  for (experiment in experiments) {
    # Pan-cluster correlation (current app approach)
    pan_result <- calculate_correlation_stats(mast_gene, crispri_gene, experiment, cluster = NULL)
    results <- rbind(results, data.frame(
      gene_pair = gene_pair,
      mast_gene = mast_gene,
      crispri_gene = crispri_gene,
      experiment = experiment,
      cluster = "Pan-cluster",
      correlation = pan_result$correlation,
      p_value = pan_result$p_value,
      n_genes = pan_result$n_genes,
      error = pan_result$error %||% NA,
      stringsAsFactors = FALSE
    ))
    
    # Cluster-specific correlations
    # Get available clusters for this gene pair
    if (mast_gene %in% names(de_data$iSCORE_PD_MAST) && 
        crispri_gene %in% names(de_data$CRISPRi_Mixscale)) {
      
      mast_clusters <- names(de_data$iSCORE_PD_MAST[[mast_gene]])
      crispri_clusters <- names(de_data$CRISPRi_Mixscale[[crispri_gene]])
      common_clusters <- intersect(mast_clusters, crispri_clusters)
      
      for (cluster in common_clusters) {
        cluster_result <- calculate_correlation_stats(mast_gene, crispri_gene, experiment, cluster)
        results <- rbind(results, data.frame(
          gene_pair = gene_pair,
          mast_gene = mast_gene,
          crispri_gene = crispri_gene,
          experiment = experiment,
          cluster = cluster,
          correlation = cluster_result$correlation,
          p_value = cluster_result$p_value,
          n_genes = cluster_result$n_genes,
          error = cluster_result$error %||% NA,
          stringsAsFactors = FALSE
        ))
      }
    }
  }
}

# Clean and summarize results
results$correlation <- as.numeric(results$correlation)
results$p_value <- as.numeric(results$p_value)
results$n_genes <- as.numeric(results$n_genes)

# Add significance and strength categories
results$significant <- results$p_value < 0.05 & !is.na(results$p_value)
results$correlation_strength <- cut(abs(results$correlation), 
                                   breaks = c(0, 0.1, 0.3, 0.5, 0.7, 1.0),
                                   labels = c("Negligible", "Weak", "Moderate", "Strong", "Very Strong"),
                                   include.lowest = TRUE)

# Save detailed results
write.csv(results, "correlation_quality_analysis.csv", row.names = FALSE)
cat("✓ Detailed results saved to correlation_quality_analysis.csv\n\n")

# Print summary
cat("=== SUMMARY STATISTICS ===\n")

# Best correlations overall
valid_significant <- results[!is.na(results$correlation) & results$significant, ]
if (nrow(valid_significant) > 0) {
  # Order by absolute correlation (descending)
  best_correlations <- valid_significant[order(abs(valid_significant$correlation), decreasing = TRUE), ]
  best_correlations <- head(best_correlations, 20)
} else {
  best_correlations <- data.frame()
}

if (nrow(best_correlations) > 0) {
  cat("TOP 20 SIGNIFICANT CORRELATIONS:\n")
  print(best_correlations[, c("gene_pair", "experiment", "cluster", "correlation", "p_value", "n_genes")], row.names = FALSE)
  cat("\n")
} else {
  cat("⚠️  NO SIGNIFICANT CORRELATIONS FOUND!\n\n")
}

# Compare pan-cluster vs cluster-specific
cat("CORRELATION STRENGTH COMPARISON:\n")
pan_cluster_results <- results[results$cluster == "Pan-cluster" & !is.na(results$correlation), ]
cluster_specific_results <- results[results$cluster != "Pan-cluster" & !is.na(results$correlation), ]

cat("Pan-cluster correlations:\n")
cat("  Mean |r|:", round(mean(abs(pan_cluster_results$correlation), na.rm = TRUE), 3), "\n")
cat("  Significant:", sum(pan_cluster_results$significant, na.rm = TRUE), "/", nrow(pan_cluster_results), "\n")

cat("Cluster-specific correlations:\n")
cat("  Mean |r|:", round(mean(abs(cluster_specific_results$correlation), na.rm = TRUE), 3), "\n")
cat("  Significant:", sum(cluster_specific_results$significant, na.rm = TRUE), "/", nrow(cluster_specific_results), "\n\n")

# Best gene pairs
cat("BEST GENE PAIRS (by max |correlation|):\n")
valid_results <- results[!is.na(results$correlation), ]
gene_pairs_unique <- unique(valid_results$gene_pair)
gene_pair_summary <- data.frame()

for (gp in gene_pairs_unique) {
  gp_data <- valid_results[valid_results$gene_pair == gp, ]
  max_idx <- which.max(abs(gp_data$correlation))
  
  gene_pair_summary <- rbind(gene_pair_summary, data.frame(
    gene_pair = gp,
    max_abs_correlation = max(abs(gp_data$correlation), na.rm = TRUE),
    best_correlation = gp_data$correlation[max_idx],
    best_p_value = gp_data$p_value[max_idx],
    best_experiment = gp_data$experiment[max_idx],
    best_cluster = gp_data$cluster[max_idx],
    significant_count = sum(gp_data$significant, na.rm = TRUE),
    total_tests = nrow(gp_data),
    stringsAsFactors = FALSE
  ))
}

gene_pair_summary <- gene_pair_summary[order(gene_pair_summary$max_abs_correlation, decreasing = TRUE), ]

print(gene_pair_summary, n = 15)

cat("\n=== RECOMMENDATIONS ===\n")
good_correlations <- gene_pair_summary[gene_pair_summary$max_abs_correlation > 0.3, ]
if (nrow(good_correlations) > 0) {
  cat("✓ Gene pairs worth showing (|r| > 0.3):\n")
  for (i in 1:min(5, nrow(good_correlations))) {
    cat("  -", good_correlations$gene_pair[i], 
        sprintf("(r=%.3f, %s, %s)\n", 
                good_correlations$best_correlation[i],
                good_correlations$best_experiment[i],
                good_correlations$best_cluster[i]))
  }
} else {
  cat("⚠️  NO GENE PAIRS show strong correlations (|r| > 0.3)\n")
  cat("   Consider removing correlation plots or lowering expectations.\n")
}

if (mean(abs(cluster_specific_results$correlation), na.rm = TRUE) > 
    mean(abs(pan_cluster_results$correlation), na.rm = TRUE)) {
  cat("✓ Cluster-specific correlations are stronger than pan-cluster\n")
  cat("   Recommend implementing cluster-specific correlation option.\n")
} else {
  cat("→ Pan-cluster correlations are as good as cluster-specific\n")
  cat("   Current approach is reasonable.\n")
}

cat("\nAnalysis complete! Check correlation_quality_analysis.csv for details.\n")