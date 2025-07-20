# PD Signature Visualization Functions
# Purpose: Create publication-ready visualizations for PD signatures
# Date: 2025-07-19

library(ggplot2)
library(dplyr)
library(tidyr)
library(pheatmap)
library(RColorBrewer)

# Helper function for word wrapping
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) {
    if (is.na(x)) return('')
    stringr::str_wrap(x, width = width)
  })
}


# Set theme
theme_set(theme_minimal(base_size = 12))

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Create a bar plot of top signatures
#' @param data Data frame with signature results
#' @param title Plot title
#' @param n_top Number of top signatures to show
#' @param color_by Column to use for coloring
create_signature_barplot <- function(data, title, n_top = 20, color_by = "enrichment_type") {
  # Ensure Description column exists
  if (!"Description" %in% names(data)) {
    stop("Description column not found in data")
  }
  
  # Prepare data
  plot_data <- data %>%
    head(n_top) %>%
    mutate(
      Description = as.character(Description),
      Description_short = ifelse(nchar(Description) > 50,
                                paste0(wrap_text(Description, width = 50), "..."),
                                Description)
    ) %>%
    arrange(desc(mean_neg_log_p))
  
  # Create plot
  p <- ggplot(plot_data, aes(x = reorder(Description_short, mean_neg_log_p), 
                             y = mean_neg_log_p,
                             fill = .data[[color_by]])) +
    geom_bar(stat = "identity") +
    coord_flip() +
    labs(
      title = title,
      x = "",
      y = "Mean -log10(adjusted p-value)",
      fill = "Enrichment Type"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10)
    ) +
    scale_fill_brewer(palette = "Set2")
  
  return(p)
}

#' Create a heatmap of pathway enrichment across genes
#' @param enrichment_data Full enrichment data
#' @param pathways Vector of pathway descriptions to include
#' @param method "MAST" or "MixScale" or "both"
#' @param value_col Column to use for heatmap values
create_pathway_gene_heatmap <- function(enrichment_data, pathways, method = "both", 
                                       value_col = "neg_log_p") {
  # Filter data
  filtered_data <- enrichment_data %>%
    filter(Description %in% pathways)
  
  if (method != "both") {
    filtered_data <- filtered_data %>%
      filter(method == method)
  }
  
  # Add value column
  if (value_col == "neg_log_p") {
    filtered_data$value <- -log10(filtered_data$p.adjust)
  } else if (value_col == "fold_enrichment") {
    filtered_data$value <- log2(filtered_data$FoldEnrichment + 1)
  } else {
    filtered_data$value <- filtered_data[[value_col]]
  }
  
  # Create matrix
  heatmap_data <- filtered_data %>%
    group_by(Description, gene) %>%
    summarise(value = mean(value, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = gene, values_from = value, values_fill = 0)
  
  # Convert to matrix
  mat <- as.matrix(heatmap_data[, -1])
  rownames(mat) <- heatmap_data$Description
  
  # Truncate row names if too long
  rownames(mat) <- wrap_text(rownames(mat), width = 60)
  
  # Create heatmap
  pheatmap(
    mat,
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    color = colorRampPalette(c("white", "yellow", "orange", "red"))(50),
    border_color = NA,
    cellwidth = 20,
    cellheight = 12,
    fontsize = 10,
    main = paste("Pathway Enrichment Across Genes -", method),
    annotation_legend = TRUE,
    show_colnames = TRUE,
    show_rownames = TRUE
  )
}

#' Create a comparison plot between MAST and MixScale
#' @param convergent_data Data frame with convergent signatures
#' @param n_top Number of top signatures to show
create_method_comparison_plot <- function(convergent_data, n_top = 15) {
  # Prepare data for plotting
  plot_data <- convergent_data %>%
    head(n_top) %>%
    mutate(
      Description_short = ifelse(nchar(Description) > 40,
                                paste0(wrap_text(Description, width = 40), "..."),
                                Description),
      total_genes = n_genes_mast + n_genes_mixscale
    ) %>%
    arrange(desc(total_genes))
  
  # Create data for grouped bar plot
  mast_data <- plot_data %>%
    select(Description_short, n_genes = n_genes_mast) %>%
    mutate(Method = "MAST")
  
  mixscale_data <- plot_data %>%
    select(Description_short, n_genes = n_genes_mixscale) %>%
    mutate(Method = "MixScale")
  
  combined_data <- rbind(mast_data, mixscale_data)
  
  # Create plot
  p <- ggplot(combined_data, aes(x = reorder(Description_short, n_genes), 
                                 y = n_genes, 
                                 fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_flip() +
    labs(
      title = "Convergent Pathways: Gene Coverage by Method",
      x = "",
      y = "Number of Genes",
      fill = "Method"
    ) +
    theme(
      legend.position = "bottom",
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 10)
    ) +
    scale_fill_manual(values = c("MAST" = "#1f77b4", "MixScale" = "#ff7f0e"))
  
  return(p)
}

#' Create a summary visualization for presentation
#' @param mast_top Top MAST signatures
#' @param mixscale_top Top MixScale signatures
#' @param convergent_top Top convergent signatures
create_summary_visualization <- function(mast_top, mixscale_top, convergent_top) {
  # Prepare data for three-panel plot
  mast_summary <- mast_top %>%
    head(10) %>%
    mutate(
      Category = "MAST-only",
      Description_short = ifelse(nchar(Description) > 40,
                                paste0(wrap_text(Description, width = 40), "..."),
                                Description),
      Score = mean_neg_log_p
    )
  
  mixscale_summary <- mixscale_top %>%
    head(10) %>%
    mutate(
      Category = "MixScale-only",
      Description_short = ifelse(nchar(Description) > 40,
                                paste0(wrap_text(Description, width = 40), "..."),
                                Description),
      Score = mean_neg_log_p
    )
  
  convergent_summary <- convergent_top %>%
    head(10) %>%
    mutate(
      Category = "Convergent",
      Description_short = ifelse(nchar(Description) > 40,
                                paste0(wrap_text(Description, width = 40), "..."),
                                Description),
      Score = mean_neg_log_p
    )
  
  # Combine
  all_summary <- rbind(
    mast_summary %>% select(Category, Description_short, Score),
    mixscale_summary %>% select(Category, Description_short, Score),
    convergent_summary %>% select(Category, Description_short, Score)
  )
  
  # Create faceted plot
  p <- ggplot(all_summary, aes(x = reorder(Description_short, Score), 
                               y = Score, 
                               fill = Category)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    facet_wrap(~ Category, scales = "free", ncol = 1) +
    labs(
      title = "Top Parkinson's Disease Pathway Signatures",
      subtitle = "Comparison across MAST mutations and CRISPRi knockdowns",
      x = "",
      y = "Mean -log10(adjusted p-value)"
    ) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 16),
      plot.subtitle = element_text(size = 12),
      strip.text = element_text(face = "bold", size = 12),
      axis.text.y = element_text(size = 9)
    ) +
    scale_fill_manual(values = c(
      "MAST-only" = "#1f77b4",
      "MixScale-only" = "#ff7f0e", 
      "Convergent" = "#2ca02c"
    ))
  
  return(p)
}

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if (sys.nframe() == 0) {  # Only run if script is executed directly
  
  # Set paths
  results_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
  output_dir <- file.path(results_dir, "visualizations")
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # Load results
  cat("Loading signature results...\n")
  mast_top <- read.csv(file.path(results_dir, "mast_top_fast.csv"))
  mixscale_top <- read.csv(file.path(results_dir, "mixscale_top_fast.csv"))
  convergent_top <- read.csv(file.path(results_dir, "convergent_top_fast.csv"))
  
  # Create visualizations
  cat("Creating visualizations...\n")
  
  # 1. Bar plots for each category
  p1 <- create_signature_barplot(mast_top, "Top MAST-only PD Pathways", n_top = 15)
  ggsave(file.path(output_dir, "mast_only_barplot.pdf"), p1, width = 10, height = 8)
  
  p2 <- create_signature_barplot(mixscale_top, "Top MixScale-only PD Pathways", n_top = 15)
  ggsave(file.path(output_dir, "mixscale_only_barplot.pdf"), p2, width = 10, height = 8)
  
  # 2. Method comparison plot
  p3 <- create_method_comparison_plot(convergent_top, n_top = 15)
  ggsave(file.path(output_dir, "convergent_comparison.pdf"), p3, width = 10, height = 8)
  
  # 3. Summary visualization
  p4 <- create_summary_visualization(mast_top, mixscale_top, convergent_top)
  ggsave(file.path(output_dir, "summary_three_panel.pdf"), p4, width = 10, height = 12)
  
  cat("Visualizations saved to:", output_dir, "\n")
  
  # For heatmaps, we need the full enrichment data
  # This would be loaded if running the full analysis
  cat("\nNote: For full heatmap visualizations, run with complete enrichment data.\n")
}
