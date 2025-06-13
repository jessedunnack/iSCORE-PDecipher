# Progressive Visualization System
# Three-tier architecture for enrichment visualizations

library(ggplot2)
library(dplyr)

# Load ComplexHeatmap conditionally
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  library(ComplexHeatmap)
}

#' Create Tier 1 visualizations (instant from consolidated data)
#' 
#' @param filtered_data Filtered consolidated enrichment data
#' @return List of ggplot objects
create_tier1_visualizations <- function(filtered_data) {
  if (is.null(filtered_data) || nrow(filtered_data) == 0) {
    return(list(
      summary_table = NULL,
      term_counts = NULL,
      pvalue_distribution = NULL,
      top_terms_bar = NULL
    ))
  }
  
  # Get p-value column name
  pval_col <- if ("p.adjust" %in% names(filtered_data)) "p.adjust" 
              else if ("fdr" %in% names(filtered_data)) "fdr" 
              else if ("padj" %in% names(filtered_data)) "padj"
              else NULL
  
  # Get description column
  desc_col <- if ("Description" %in% names(filtered_data)) "Description"
              else if ("description" %in% names(filtered_data)) "description"
              else if ("term" %in% names(filtered_data)) "term"
              else NULL
  
  results <- list()
  
  # 1. Summary statistics table
  results$summary_table <- create_summary_table(filtered_data)
  
  # 2. Term count by category
  if ("enrichment_type" %in% names(filtered_data)) {
    results$term_counts <- ggplot(filtered_data, 
                                 aes(x = enrichment_type)) +
      geom_bar(fill = "#3c8dbc") +
      labs(title = "Significant Terms by Database",
           x = "Enrichment Database",
           y = "Number of Terms") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  # 3. P-value distribution
  if (!is.null(pval_col)) {
    results$pvalue_distribution <- ggplot(filtered_data, 
                                         aes_string(x = pval_col)) +
      geom_histogram(bins = 30, fill = "#3c8dbc", alpha = 0.7) +
      scale_x_continuous(trans = "log10") +
      labs(title = "P-value Distribution",
           x = "Adjusted P-value (log10)",
           y = "Count") +
      theme_minimal()
  }
  
  # 4. Top terms bar plot
  if (!is.null(pval_col) && !is.null(desc_col)) {
    top_terms <- filtered_data %>%
      arrange(!!sym(pval_col)) %>%
      head(20)
    
    results$top_terms_bar <- ggplot(top_terms, 
                                   aes_string(x = sprintf("reorder(%s, -%s)", desc_col, pval_col),
                                             y = sprintf("-log10(%s)", pval_col))) +
      geom_bar(stat = "identity", fill = "#3c8dbc") +
      coord_flip() +
      labs(title = "Top 20 Significant Terms",
           x = "",
           y = "-log10(adjusted p-value)") +
      theme_minimal() +
      theme(axis.text.y = element_text(size = 8))
  }
  
  return(results)
}

#' Create summary table from filtered data
#' 
#' @param data Filtered enrichment data
#' @return Summary data frame
create_summary_table <- function(data) {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  summary_stats <- data %>%
    group_by(enrichment_type, direction) %>%
    summarise(
      n_terms = n(),
      min_pval = min(p.adjust, na.rm = TRUE),
      median_pval = median(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(enrichment_type, direction)
  
  return(summary_stats)
}

#' Create Tier 2 visualizations (unified heatmaps from consolidated data)
#' 
#' @param filtered_data Filtered consolidated enrichment data
#' @param params Visualization parameters
#' @return List of ComplexHeatmap objects
create_tier2_visualizations <- function(filtered_data, params) {
  # Source unified heatmap functions if not already loaded
  if (!exists("create_unified_enrichment_heatmap")) {
    source("unified_enrichment_heatmaps.R")
  }
  
  results <- list()
  
  # Default parameters
  default_params <- list(
    cluster_val = params$cluster,
    method_val = params$method_val %||% "all",
    enrichment_types_to_include = params$enrichment_type,
    show_direction_annotation = params$show_directions %||% TRUE,
    max_terms_overall = params$max_terms %||% 25
  )
  
  # P-value heatmap
  results$pvalue_heatmap <- do.call(create_pvalue_heatmap, 
                                   c(list(data_full = filtered_data), 
                                     default_params))
  
  # Fold enrichment heatmap (if column exists)
  if ("FoldEnrichment" %in% names(filtered_data)) {
    results$enrichment_heatmap <- do.call(create_fold_enrichment_heatmap,
                                         c(list(data_full = filtered_data), 
                                           default_params))
  }
  
  # Z-score heatmap (if column exists)
  if ("zScore" %in% names(filtered_data)) {
    results$zscore_heatmap <- do.call(create_zscore_heatmap,
                                     c(list(data_full = filtered_data), 
                                       default_params))
  }
  
  # Direction comparison (if requested)
  if (params$show_direction_comparison %||% FALSE) {
    results$direction_comparison <- create_direction_comparison_heatmap(
      filtered_data, 
      cluster_val = params$cluster,
      enrichment_type = params$enrichment_type
    )
  }
  
  return(results)
}

#' Create direction comparison heatmap
#' 
#' @param data Filtered enrichment data
#' @param cluster_val Cluster to analyze
#' @param enrichment_type Enrichment type to include
#' @return ComplexHeatmap object
create_direction_comparison_heatmap <- function(data, cluster_val, enrichment_type) {
  # Filter for specific cluster and enrichment type
  subset_data <- data %>%
    filter(cluster == cluster_val,
           enrichment_type == enrichment_type)
  
  if (nrow(subset_data) == 0) return(NULL)
  
  # Create comparison configs for UP vs DOWN vs ALL
  configs <- list(
    list(
      heatmap_type = "standard",
      metric_type = "p.adjust",
      direction_filter = "UP",
      custom_title = "UP-regulated"
    ),
    list(
      heatmap_type = "standard", 
      metric_type = "p.adjust",
      direction_filter = "DOWN",
      custom_title = "DOWN-regulated"
    ),
    list(
      heatmap_type = "standard",
      metric_type = "p.adjust",
      direction_filter = "ALL",
      custom_title = "ALL genes"
    )
  )
  
  # Create comparison heatmaps
  heatmaps <- lapply(configs, function(config) {
    dir_data <- subset_data %>%
      filter(direction == config$direction_filter)
    
    if (nrow(dir_data) > 0) {
      create_unified_enrichment_heatmap(
        dir_data,
        heatmap_type = config$heatmap_type,
        metric_type = config$metric_type,
        custom_title = config$custom_title,
        show_direction_annotation = FALSE
      )
    } else {
      NULL
    }
  })
  
  # Remove NULL heatmaps
  heatmaps <- Filter(Negate(is.null), heatmaps)
  
  # Combine if multiple exist
  if (length(heatmaps) > 1) {
    return(Reduce(`+`, heatmaps))
  } else if (length(heatmaps) == 1) {
    return(heatmaps[[1]])
  } else {
    return(NULL)
  }
}

#' Create Tier 3 visualizations (requiring RDS file loading)
#' 
#' @param file_path Path to RDS file
#' @param plot_type Type of enrichPlot visualization
#' @param params Additional parameters
#' @return Plot object or NULL
create_tier3_visualizations <- function(file_path, plot_type, params = list()) {
  # Check if enrichplot is available
  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    message("enrichplot package not available for Tier 3 visualizations")
    return(NULL)
  }
  
  # Load RDS data with caching
  if (exists("load_enrichment_with_cache")) {
    rds_data <- load_enrichment_with_cache(
      analysis_type = params$analysis_type,
      gene = params$gene,
      cluster = params$cluster,
      experiment = params$experiment,
      enrichment_type = params$enrichment_type,
      direction = params$direction
    )
  } else {
    rds_data <- load_enrichment_safe(file_path)
  }
  
  if (is.null(rds_data)) {
    message("Could not load RDS data for Tier 3 visualization")
    return(NULL)
  }
  
  # Create enrichplot visualization based on type
  result <- tryCatch({
    switch(plot_type,
      "dotplot" = enrichplot::dotplot(rds_data, showCategory = params$n_terms %||% 20),
      "barplot" = enrichplot::barplot(rds_data, showCategory = params$n_terms %||% 20),
      "cnetplot" = enrichplot::cnetplot(rds_data, showCategory = params$n_terms %||% 10),
      "emapplot" = enrichplot::emapplot(rds_data, showCategory = params$n_terms %||% 30),
      "heatplot" = enrichplot::heatplot(rds_data, showCategory = params$n_terms %||% 20),
      "upsetplot" = enrichplot::upsetplot(rds_data),
      NULL
    )
  }, error = function(e) {
    message("Error creating enrichplot visualization: ", e$message)
    NULL
  })
  
  return(result)
}

#' Helper function for NULL default values
#' 
#' @param x Value to check
#' @param y Default value if x is NULL
#' @return x if not NULL, otherwise y
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' Create simple bar plot from consolidated data
#' 
#' @param data Enrichment data
#' @param n_terms Number of terms to show
#' @param metric Metric to plot (p.adjust, FoldEnrichment, etc.)
#' @return ggplot object
create_simple_bar_plot <- function(data, n_terms = 20, metric = "p.adjust") {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  # Get appropriate columns
  desc_col <- if ("Description" %in% names(data)) "Description"
              else if ("description" %in% names(data)) "description"
              else if ("term" %in% names(data)) "term"
              else NULL
  
  if (is.null(desc_col) || !metric %in% names(data)) return(NULL)
  
  # Prepare data
  plot_data <- data %>%
    arrange(!!sym(metric)) %>%
    head(n_terms) %>%
    mutate(Description = factor(!!sym(desc_col), 
                               levels = rev(!!sym(desc_col))))
  
  # Create plot based on metric
  if (metric %in% c("p.adjust", "pvalue", "fdr", "padj")) {
    # For p-values, use -log10 transformation
    p <- ggplot(plot_data, aes(x = Description, y = -log10(!!sym(metric)))) +
      geom_bar(stat = "identity", fill = "#3c8dbc") +
      labs(y = paste0("-log10(", metric, ")"))
  } else {
    # For other metrics, use as is
    p <- ggplot(plot_data, aes(x = Description, y = !!sym(metric))) +
      geom_bar(stat = "identity", fill = "#3c8dbc") +
      labs(y = metric)
  }
  
  p + coord_flip() +
    labs(x = "", title = paste("Top", n_terms, "terms by", metric)) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9))
}

#' Create simple dot plot from consolidated data
#' 
#' @param data Enrichment data
#' @param n_terms Number of terms to show
#' @param x_metric Metric for x-axis
#' @param color_metric Metric for color
#' @param size_metric Metric for size
#' @return ggplot object
create_simple_dot_plot <- function(data, n_terms = 20, 
                                  x_metric = "GeneRatio",
                                  color_metric = "p.adjust", 
                                  size_metric = "Count") {
  if (is.null(data) || nrow(data) == 0) return(NULL)
  
  # Get description column
  desc_col <- if ("Description" %in% names(data)) "Description"
              else if ("description" %in% names(data)) "description"
              else if ("term" %in% names(data)) "term"
              else NULL
  
  if (is.null(desc_col)) return(NULL)
  
  # Handle GeneRatio if it needs to be calculated
  if (x_metric == "GeneRatio" && !x_metric %in% names(data)) {
    if (all(c("Count", "BgRatio") %in% names(data))) {
      # Parse BgRatio and calculate GeneRatio
      data <- data %>%
        mutate(
          BgSize = as.numeric(sub(".*/", "", BgRatio)),
          GeneRatio = Count / BgSize
        )
    } else {
      # Fallback to Count
      x_metric <- "Count"
    }
  }
  
  # Filter to top terms
  plot_data <- data %>%
    arrange(!!sym(color_metric)) %>%
    head(n_terms) %>%
    mutate(Description = factor(!!sym(desc_col), 
                               levels = rev(!!sym(desc_col))))
  
  # Create dot plot
  p <- ggplot(plot_data, aes_string(x = x_metric, 
                                    y = "Description",
                                    color = color_metric,
                                    size = size_metric)) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(x = x_metric,
         y = "",
         color = color_metric,
         size = size_metric,
         title = paste("Top", n_terms, "enriched terms")) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 9))
  
  return(p)
}