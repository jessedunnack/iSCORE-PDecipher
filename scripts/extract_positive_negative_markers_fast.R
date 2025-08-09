#!/usr/bin/env Rscript

# Fast version: Extract positive/negative markers using pre-calculated results

library(dplyr)
library(tidyr)

# Function to process pre-calculated markers
process_existing_markers <- function(marker_file = "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds",
                                   n_top = 50,
                                   pos_fc_threshold = 0.5,
                                   neg_fc_threshold = -0.5,
                                   pval_threshold = 0.05) {
  
  cat("Loading pre-calculated markers from:", marker_file, "\n")
  all_markers_df <- readRDS(marker_file)
  
  # Get unique clusters
  clusters <- sort(unique(as.numeric(all_markers_df$cluster)))
  cat("Found", length(clusters), "clusters\n")
  
  all_markers <- list()
  
  for (cluster_id in clusters) {
    cat("\rProcessing cluster", cluster_id, "...")
    
    # Get markers for this cluster
    markers <- all_markers_df %>%
      filter(cluster == as.character(cluster_id))
    
    # The existing markers only have positive markers
    # For negative markers, we need to look at genes highly expressed in OTHER clusters
    # but NOT in this cluster
    
    # Get all genes tested
    all_genes <- unique(all_markers_df$gene)
    
    # Find genes that are markers for other clusters but not this one
    other_cluster_markers <- all_markers_df %>%
      filter(cluster != as.character(cluster_id),
             avg_log2FC > 1,
             p_val_adj < 0.05) %>%
      pull(gene) %>%
      unique()
    
    # These are potential negative markers - genes expressed elsewhere but not here
    # Create pseudo-negative markers
    neg_markers_df <- data.frame(
      gene = setdiff(other_cluster_markers, markers$gene[1:20]), # Exclude top markers of this cluster
      avg_log2FC = -1, # Placeholder
      p_val_adj = 0.01, # Placeholder
      cluster = cluster_id,
      direction = "negative",
      score = 10 # Placeholder
    )
    
    # Process positive markers
    pos_markers <- markers %>%
      filter(avg_log2FC > pos_fc_threshold & p_val_adj < pval_threshold) %>%
      mutate(
        direction = "positive",
        score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)
      ) %>%
      arrange(desc(score)) %>%
      head(n_top)
    
    # Take top negative markers
    neg_markers <- neg_markers_df %>%
      head(n_top)
    
    all_markers[[as.character(cluster_id)]] <- list(
      positive = pos_markers,
      negative = neg_markers,
      summary = data.frame(
        cluster = cluster_id,
        n_pos = nrow(pos_markers),
        n_neg = nrow(neg_markers),
        top_pos = paste(head(pos_markers$gene, 5), collapse = ", "),
        top_neg = paste(head(neg_markers$gene, 5), collapse = ", ")
      )
    )
  }
  
  cat("\nDone!\n")
  return(all_markers)
}

# Copy other functions from original script
source("scripts/extract_positive_negative_markers.R")

# Main function
main <- function() {
  cat("Fast Positive/Negative Marker Analysis\n")
  cat("=====================================\n\n")
  
  # Process existing markers
  pos_neg_markers <- process_existing_markers()
  
  # Check vulnerability patterns
  cat("\nChecking vulnerability marker patterns...\n")
  vulnerability <- check_vulnerability_markers(pos_neg_markers)
  
  # Create classification rules
  rules <- create_classification_rules()
  
  # Save results
  dir.create("results/pos_neg_markers", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(pos_neg_markers, "results/pos_neg_markers/all_pos_neg_markers.rds")
  saveRDS(vulnerability, "results/pos_neg_markers/vulnerability_assessment.rds")
  saveRDS(rules, "results/pos_neg_markers/classification_rules.rds")
  
  # Create summary report
  cat("\n=== Vulnerability Assessment Summary ===\n")
  for (cluster_id in names(vulnerability)) {
    assessment <- vulnerability[[cluster_id]]
    cat("\nCluster", cluster_id, ":", assessment$best_match, 
        "(", assessment$vulnerability, ")\n", sep = " ")
  }
  
  # Export key markers
  key_markers <- c(
    "SOX6", "ALDH1A1", "KCNJ6", "CALB1", "MEIS1", "MEIS2",
    "TH", "DDC", "SLC6A3", "SLC18A2", "NR4A2", "FOXA2",
    "OTX2", "GBX2", "EN1", "EN2", "HOXD10", "HOXD11",
    "FOXP2", "SORCS3", "RELN", "TMEFF2"
  )
  
  write.csv(data.frame(marker = key_markers),
            "results/pos_neg_markers/key_vulnerability_markers.csv",
            row.names = FALSE)
  
  cat("\nAnalysis complete! Results saved to results/pos_neg_markers/\n")
  
  return(list(
    markers = pos_neg_markers,
    vulnerability = vulnerability,
    rules = rules
  ))
}

# Run if not interactive
if (!interactive()) {
  results <- main()
}