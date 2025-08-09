#!/usr/bin/env Rscript

# Analyze marker strength and annotation confidence for each cluster

library(dplyr)
library(ggplot2)

# Function to analyze marker strength
analyze_marker_strength <- function(marker_file, annotation_file) {
  
  # Load data
  markers <- readRDS(marker_file)
  annotations <- read.csv(annotation_file, stringsAsFactors = FALSE)
  
  # Calculate marker statistics per cluster
  marker_stats <- markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    summarise(
      n_sig_genes = n(),
      n_strong_markers = sum(avg_log2FC > 1),
      n_very_strong = sum(avg_log2FC > 2),
      max_log2FC = max(avg_log2FC),
      mean_log2FC_top10 = mean(head(sort(avg_log2FC, decreasing = TRUE), 10)),
      
      # Specificity: how many clusters express each marker
      # (would need full expression matrix for true specificity)
      
      .groups = 'drop'
    ) %>%
    mutate(cluster_num = as.numeric(as.character(cluster))) %>%
    arrange(cluster_num)
  
  # Merge with annotations
  result <- merge(
    marker_stats,
    annotations[, c("fine_cluster", "cell_type", "confidence")],
    by.x = "cluster_num",
    by.y = "fine_cluster",
    all.x = TRUE
  )
  
  return(result)
}

# Function to calculate confidence scores
calculate_confidence_scores <- function(markers, cluster_id) {
  
  cluster_markers <- markers %>%
    filter(cluster == cluster_id, p_val_adj < 0.05) %>%
    arrange(desc(avg_log2FC))
  
  # Score components (0-100 scale)
  scores <- list()
  
  # 1. Number of significant markers (max 30 points)
  n_sig <- nrow(cluster_markers)
  scores$n_markers <- min(30, n_sig / 10)  # 300+ markers = full points
  
  # 2. Strength of top markers (max 30 points)
  if (n_sig >= 10) {
    top10_fc <- mean(head(cluster_markers$avg_log2FC, 10))
    scores$marker_strength <- min(30, top10_fc * 10)  # FC>3 = full points
  } else {
    scores$marker_strength <- 0
  }
  
  # 3. Statistical significance (max 20 points)
  if (n_sig >= 10) {
    top10_pval <- mean(-log10(head(cluster_markers$p_val_adj + 1e-300, 10)))
    scores$significance <- min(20, top10_pval / 15)  # -log10(p) > 300 = full points
  } else {
    scores$significance <- 0
  }
  
  # 4. Marker exclusivity (max 20 points)
  # This is a placeholder - would need full expression matrix
  scores$exclusivity <- 10  # Default middle score
  
  # Total score
  scores$total <- sum(unlist(scores))
  scores$confidence_level <- case_when(
    scores$total >= 70 ~ "High",
    scores$total >= 40 ~ "Medium",
    TRUE ~ "Low"
  )
  
  return(scores)
}

# Function to identify canonical markers
identify_canonical_markers <- function() {
  # Define canonical markers for major cell types
  canonical <- list(
    dopaminergic = c("TH", "DDC", "SLC6A3", "SLC18A2", "LMX1A", "FOXA2", "EN1", "NR4A2"),
    floor_plate = c("FOXA2", "SHH", "CORIN", "LMX1A", "LMX1B", "ARX"),
    neural_prog = c("SOX2", "NES", "VIM", "PAX6", "FABP7"),
    neurons_pan = c("MAP2", "TUBB3", "SYN1", "SNAP25", "RBFOX3"),
    oligodendro = c("OLIG2", "SOX10", "MBP", "PLP1", "MOG", "MAG"),
    astrocyte = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2"),
    choroid = c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1"),
    vascular = c("TAGLN", "ACTA2", "MYL9", "PECAM1", "VWF"),
    proliferating = c("MKI67", "TOP2A", "CDK1", "PCNA", "MCM2"),
    stress = c("DDIT3", "ATF3", "ATF4", "HSPA1A", "FOS", "JUN")
  )
  
  return(canonical)
}

# Main analysis
main <- function() {
  cat("Analyzing Annotation Confidence\n")
  cat("==============================\n\n")
  
  # File paths
  marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
  annotation_file <- "results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv"
  
  # Analyze marker strength
  cat("Analyzing marker strength per cluster...\n")
  marker_strength <- analyze_marker_strength(marker_file, annotation_file)
  
  # Sort by confidence and marker strength
  marker_strength <- marker_strength %>%
    mutate(
      confidence_order = factor(confidence, levels = c("High", "Medium", "Low")),
      marker_score = n_strong_markers + n_very_strong * 2
    ) %>%
    arrange(confidence_order, desc(marker_score))
  
  # Print summary
  cat("\n=== Marker Strength by Confidence Level ===\n")
  
  conf_summary <- marker_strength %>%
    group_by(confidence) %>%
    summarise(
      n_clusters = n(),
      avg_sig_genes = mean(n_sig_genes),
      avg_strong_markers = mean(n_strong_markers),
      avg_max_FC = mean(max_log2FC),
      .groups = 'drop'
    )
  
  print(conf_summary)
  
  # Identify clusters that might need re-evaluation
  cat("\n\n=== Clusters for Potential Re-evaluation ===\n")
  
  # Low confidence with many markers
  low_conf_strong <- marker_strength %>%
    filter(confidence == "Low", n_strong_markers > 10) %>%
    select(cluster_num, cell_type, n_sig_genes, n_strong_markers, max_log2FC)
  
  if (nrow(low_conf_strong) > 0) {
    cat("\nLow confidence clusters with strong markers:\n")
    print(low_conf_strong)
  }
  
  # High confidence with few markers
  high_conf_weak <- marker_strength %>%
    filter(confidence == "High", n_strong_markers < 5) %>%
    select(cluster_num, cell_type, n_sig_genes, n_strong_markers)
  
  if (nrow(high_conf_weak) > 0) {
    cat("\nHigh confidence clusters with few strong markers:\n")
    print(high_conf_weak)
  }
  
  # Calculate detailed confidence scores
  cat("\n\n=== Detailed Confidence Scores ===\n")
  markers <- readRDS(marker_file)
  
  detailed_scores <- data.frame()
  for (i in 0:35) {
    scores <- calculate_confidence_scores(markers, as.character(i))
    detailed_scores <- rbind(detailed_scores, data.frame(
      cluster = i,
      total_score = scores$total,
      n_markers_score = scores$n_markers,
      strength_score = scores$marker_strength,
      significance_score = scores$significance,
      calculated_confidence = scores$confidence_level,
      stringsAsFactors = FALSE
    ))
  }
  
  # Merge with annotations
  final_analysis <- merge(
    marker_strength,
    detailed_scores,
    by.x = "cluster_num",
    by.y = "cluster"
  )
  
  # Compare original vs calculated confidence
  cat("\nConfidence comparison:\n")
  conf_comparison <- final_analysis %>%
    select(cluster_num, cell_type, confidence, calculated_confidence, total_score) %>%
    filter(confidence != calculated_confidence) %>%
    arrange(cluster_num)
  
  if (nrow(conf_comparison) > 0) {
    cat("\nClusters with confidence mismatch:\n")
    print(conf_comparison)
  }
  
  # Save results
  write.csv(final_analysis, 
            "results/cell_type_annotations/annotation_confidence_analysis.csv",
            row.names = FALSE)
  
  # Create visualization
  cat("\n\nCreating confidence visualization...\n")
  
  p1 <- ggplot(final_analysis, aes(x = reorder(cluster_num, -total_score), 
                                    y = total_score, 
                                    fill = confidence)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("High" = "darkgreen", "Medium" = "orange", "Low" = "red")) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Annotation Confidence Scores by Cluster",
         x = "Cluster", y = "Confidence Score (0-100)")
  
  ggsave("results/cell_type_annotations/confidence_scores_barplot.pdf", p1, width = 12, height = 6)
  
  # Return canonical markers for reference
  canonical <- identify_canonical_markers()
  
  cat("\n\n=== Canonical Marker Reference ===\n")
  for (type in names(canonical)) {
    cat(sprintf("%-15s: %s\n", type, paste(canonical[[type]], collapse = ", ")))
  }
  
  return(list(
    analysis = final_analysis,
    canonical_markers = canonical,
    confidence_comparison = conf_comparison
  ))
}

# Run if not interactive
if (!interactive()) {
  results <- main()
}