# Analyze most significant Fisher's exact test results
# Focus on DE gene overlap between CRISPRi and iSCORE-PD mutations
# Account for directionality and avoid double-counting

library(dplyr)
library(tidyr)

# Load the data
cat("Loading enrichment data...\n")
data_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds'
enrichment_data <- readRDS(data_file)

# Also load DE results for Fisher's test data
de_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds'
de_data <- readRDS(de_file)

cat("Data loaded successfully\n")
cat("Enrichment data rows:", nrow(enrichment_data), "\n")
cat("DE data structure:", names(de_data), "\n\n")

# Function to run enhanced Fisher's exact test with directionality
run_direction_aware_fisher <- function(mast_gene, crispri_gene, cluster, de_data) {
  tryCatch({
    # Get MAST DE results
    if (!mast_gene %in% names(de_data$iSCORE_PD_MAST) ||
        !cluster %in% names(de_data$iSCORE_PD_MAST[[mast_gene]])) {
      return(NULL)
    }
    
    mast_results <- de_data$iSCORE_PD_MAST[[mast_gene]][[cluster]]$results
    mast_sig <- mast_results[!is.na(mast_results$p_val_adj) & 
                             mast_results$p_val_adj < 0.05, ]
    
    # Get CRISPRi DE results
    if (!crispri_gene %in% names(de_data$CRISPRi_Mixscale) ||
        !cluster %in% names(de_data$CRISPRi_Mixscale[[crispri_gene]])) {
      return(NULL)
    }
    
    crispri_results <- de_data$CRISPRi_Mixscale[[crispri_gene]][[cluster]]$results
    
    # Get experiment columns
    exp_cols <- grep("^log2FC_", colnames(crispri_results), value = TRUE)
    p_cols <- grep("^p_cell_type.*:weight$", colnames(crispri_results), value = TRUE)
    
    # Analyze each experiment
    experiment_results <- list()
    
    for (i in seq_along(exp_cols)) {
      exp_col <- exp_cols[i]
      p_col <- p_cols[i]
      
      if (!p_col %in% colnames(crispri_results)) next
      
      # Get significant genes for this experiment
      crispri_sig <- crispri_results[!is.na(crispri_results[[p_col]]) & 
                                     crispri_results[[p_col]] < 0.05, ]
      
      # Track directionality
      mast_up <- rownames(mast_sig)[mast_sig$avg_log2FC > 0]
      mast_down <- rownames(mast_sig)[mast_sig$avg_log2FC < 0]
      
      crispri_up <- rownames(crispri_sig)[crispri_sig[[exp_col]] > 0]
      crispri_down <- rownames(crispri_sig)[crispri_sig[[exp_col]] < 0]
      
      # Count overlaps by direction
      same_up <- length(intersect(mast_up, crispri_up))
      same_down <- length(intersect(mast_down, crispri_down))
      opposite_up_down <- length(intersect(mast_up, crispri_down))
      opposite_down_up <- length(intersect(mast_down, crispri_up))
      
      same_direction <- same_up + same_down
      opposite_direction <- opposite_up_down + opposite_down_up
      
      # Get all overlapping genes (avoiding double counting)
      all_overlap <- union(
        union(intersect(mast_up, crispri_up), intersect(mast_down, crispri_down)),
        union(intersect(mast_up, crispri_down), intersect(mast_down, crispri_up))
      )
      
      # Union background (all genes tested in either method)
      union_background <- length(union(rownames(mast_results), rownames(crispri_results)))
      
      # Fisher's exact test on union background
      contingency_matrix <- matrix(c(
        length(all_overlap),  # Both methods
        nrow(mast_sig) - length(all_overlap),  # MAST only
        nrow(crispri_sig) - length(all_overlap),  # CRISPRi only
        union_background - nrow(mast_sig) - nrow(crispri_sig) + length(all_overlap)  # Neither
      ), nrow = 2)
      
      fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
      
      experiment_results[[exp_col]] <- list(
        total_overlap = length(all_overlap),
        same_direction = same_direction,
        opposite_direction = opposite_direction,
        fisher_p = fisher_result$p.value,
        odds_ratio = fisher_result$estimate,
        overlapping_genes = all_overlap
      )
    }
    
    return(experiment_results)
    
  }, error = function(e) {
    cat("Error processing", mast_gene, "vs", crispri_gene, "in", cluster, ":", e$message, "\n")
    return(NULL)
  })
}

# Analyze top gene pairs
gene_pairs <- list(
  c("ATP13A2", "ATP13A2"),
  c("PARK7", "PARK7"),
  c("FBXO7", "FBXO7"),
  c("LRRK2", "LRRK2"),
  c("PINK1", "PINK1"),
  c("PRKN", "PARK2"),  # Note: different names for same gene
  c("SNCA", "SNCA"),
  c("SYNJ1", "SYNJ1"),
  c("DNAJC6", "DNAJC6"),
  c("VPS13C", "VPS13C")
)

# Get all clusters
all_clusters <- unique(enrichment_data$cluster)
cat("Analyzing", length(gene_pairs), "gene pairs across", length(all_clusters), "clusters\n\n")

# Collect results
all_results <- list()

for (pair in gene_pairs) {
  mast_gene <- pair[1]
  crispri_gene <- pair[2]
  pair_name <- paste0(mast_gene, "_vs_", crispri_gene)
  
  cat("Analyzing", pair_name, "...\n")
  
  pair_results <- list()
  
  for (cluster in all_clusters) {
    result <- run_direction_aware_fisher(mast_gene, crispri_gene, cluster, de_data)
    if (!is.null(result)) {
      pair_results[[cluster]] <- result
    }
  }
  
  all_results[[pair_name]] <- pair_results
}

# Find most significant results
cat("\n=== MOST SIGNIFICANT FISHER'S EXACT TEST RESULTS ===\n\n")

# Flatten results for easier analysis
flat_results <- data.frame()

for (pair_name in names(all_results)) {
  for (cluster in names(all_results[[pair_name]])) {
    for (exp in names(all_results[[pair_name]][[cluster]])) {
      exp_data <- all_results[[pair_name]][[cluster]][[exp]]
      
      flat_results <- rbind(flat_results, data.frame(
        gene_pair = pair_name,
        cluster = cluster,
        experiment = exp,
        total_overlap = exp_data$total_overlap,
        same_direction = exp_data$same_direction,
        opposite_direction = exp_data$opposite_direction,
        fisher_p = exp_data$fisher_p,
        odds_ratio = exp_data$odds_ratio,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# Sort by p-value
flat_results <- flat_results[order(flat_results$fisher_p), ]

# Show top 20 most significant
cat("Top 20 most significant overlaps (sorted by p-value):\n")
top_results <- head(flat_results, 20)
print(top_results)

# Aggregate by gene pair (avoiding double counting)
cat("\n\n=== AGGREGATED RESULTS BY GENE PAIR ===\n")

gene_pair_summary <- flat_results %>%
  group_by(gene_pair) %>%
  summarise(
    n_significant = sum(fisher_p < 0.05),
    min_p_value = min(fisher_p),
    max_overlap = max(total_overlap),
    total_clusters = n_distinct(cluster),
    avg_same_direction = mean(same_direction),
    avg_opposite_direction = mean(opposite_direction),
    .groups = 'drop'
  ) %>%
  arrange(min_p_value)

print(gene_pair_summary)

# Extract genes from most significant overlaps
cat("\n\n=== GENES IN MOST SIGNIFICANT OVERLAPS ===\n")

# Get top 5 most significant results
top5 <- head(flat_results, 5)

for (i in 1:nrow(top5)) {
  row <- top5[i, ]
  cat("\n", i, ". ", row$gene_pair, " - ", row$cluster, " (", row$experiment, ")\n", sep="")
  cat("   p-value:", format(row$fisher_p, scientific = TRUE), "\n")
  cat("   Overlap:", row$total_overlap, "genes (", row$same_direction, "same direction,", 
      row$opposite_direction, "opposite)\n")
  
  # Get the actual gene list
  pair <- strsplit(row$gene_pair, "_vs_")[[1]]
  result <- all_results[[row$gene_pair]][[row$cluster]][[row$experiment]]
  
  if (!is.null(result$overlapping_genes)) {
    genes <- head(sort(result$overlapping_genes), 10)
    cat("   Top genes:", paste(genes, collapse=", "), 
        ifelse(length(result$overlapping_genes) > 10, "...", ""), "\n")
  }
}

# Save results
saveRDS(all_results, "fisher_test_results_detailed.rds")
write.csv(flat_results, "fisher_test_results_flat.csv", row.names = FALSE)
write.csv(gene_pair_summary, "fisher_test_results_summary.csv", row.names = FALSE)

cat("\n\nAnalysis complete. Results saved to:\n")
cat("- fisher_test_results_detailed.rds\n")
cat("- fisher_test_results_flat.csv\n")
cat("- fisher_test_results_summary.csv\n")