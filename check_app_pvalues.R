# Check what p-values the Shiny app would show
# Compare with my Fisher's test analysis

# Check what p-values the Shiny app would show

# Load the same data the app uses
enrichment_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds'
de_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/full_DE_results.rds'

enrichment_data <- readRDS(enrichment_file)
de_data <- readRDS(de_file)

# Run signature analysis like the app does
cat("Running signature analysis...\n")

# Get unique clusters
clusters <- unique(enrichment_data$cluster)

# Run for a few key gene pairs
gene_pairs <- list(
  c("SYNJ1", "SYNJ1"),
  c("ATP13A2", "ATP13A2"),
  c("FBXO7", "FBXO7"),
  c("LRRK2", "LRRK2"),
  c("PARK7", "PARK7")
)

# Analyze each pair
results_summary <- data.frame()

for (pair in gene_pairs) {
  mast_gene <- pair[1]
  crispri_gene <- pair[2]
  pair_name <- paste0(mast_gene, "_vs_", crispri_gene)
  
  cat("\nAnalyzing", pair_name, ":\n")
  
  # This is what the app does - it runs Fisher's test for each cluster
  # and combines experiments somehow
  
  for (cluster in clusters[1:5]) {  # Just first 5 clusters for speed
    
    # The app's signature_analysis.R function combines all experiments
    # Let's see what it calculates
    
    # Get MAST genes
    mast_results <- NULL
    if (mast_gene %in% names(de_data$iSCORE_PD_MAST) &&
        cluster %in% names(de_data$iSCORE_PD_MAST[[mast_gene]])) {
      mast_results <- de_data$iSCORE_PD_MAST[[mast_gene]][[cluster]]$results
      mast_sig <- rownames(mast_results[!is.na(mast_results$p_val_adj) & 
                                       mast_results$p_val_adj < 0.05, ])
    }
    
    # Get CRISPRi genes (combining all experiments)
    crispri_sig_all <- c()
    if (crispri_gene %in% names(de_data$CRISPRi_Mixscale) &&
        cluster %in% names(de_data$CRISPRi_Mixscale[[crispri_gene]])) {
      crispri_results <- de_data$CRISPRi_Mixscale[[crispri_gene]][[cluster]]$results
      
      # Get all p-value columns
      p_cols <- grep("^p_cell_type.*:weight$", colnames(crispri_results), value = TRUE)
      
      # Combine significant genes from all experiments
      for (p_col in p_cols) {
        sig_genes <- rownames(crispri_results[!is.na(crispri_results[[p_col]]) & 
                                             crispri_results[[p_col]] < 0.05, ])
        crispri_sig_all <- union(crispri_sig_all, sig_genes)
      }
    }
    
    if (!is.null(mast_results) && length(crispri_sig_all) > 0) {
      # Calculate overlaps
      overlap <- length(intersect(mast_sig, crispri_sig_all))
      
      # Union background
      union_background <- length(union(rownames(mast_results), 
                                      rownames(crispri_results)))
      
      # Fisher's test on union background (what app shows by default now)
      contingency_matrix <- matrix(c(
        overlap,
        length(mast_sig) - overlap,
        length(crispri_sig_all) - overlap,
        union_background - length(mast_sig) - length(crispri_sig_all) + overlap
      ), nrow = 2)
      
      fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
      
      results_summary <- rbind(results_summary, data.frame(
        gene_pair = pair_name,
        cluster = cluster,
        mast_sig_genes = length(mast_sig),
        crispri_sig_genes = length(crispri_sig_all),
        overlap = overlap,
        union_background = union_background,
        fisher_p_union = fisher_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

cat("\n\n=== SUMMARY: What the App Shows ===\n")
print(results_summary[order(results_summary$fisher_p_union), ])

cat("\n\nKey differences from my previous analysis:\n")
cat("1. App combines ALL CRISPRi experiments into one p-value per cluster\n")
cat("2. App defaults to 'union' background (all genes tested in either method)\n")
cat("3. App may apply FDR correction (hierarchical) to these p-values\n")
cat("4. My analysis showed individual experiment p-values (C12_FPD-24, etc.)\n")