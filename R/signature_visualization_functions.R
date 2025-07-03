#' Comprehensive Signature Visualization Functions
#'
#' This module provides enhanced visualizations for signature nomination results,
#' focusing on clear and straightforward presentation of cross-method comparisons.

# Dependencies handled via DESCRIPTION file imports
# library(ggplot2)
# library(plotly) 
# library(dplyr)

#' Create Gene vs Pathway P-value Scatter Plot
#'
#' @param signature_data Data frame with signature results
#' @param interactive Logical, whether to return interactive plotly plot
#' @return ggplot2 or plotly object
#' @export
create_gene_pathway_pvalue_scatter <- function(signature_data, interactive = TRUE) {
  
  # Validate input data
  validate_signature_data(signature_data)
  
  if (nrow(signature_data) == 0) {
    return(ggplot() + 
           geom_text(aes(x = 0.5, y = 0.5, label = "No signature data available"), 
                     size = 6, color = "gray50") +
           theme_void())
  }
  
  # Prepare data with safe column access
  plot_data <- signature_data %>%
    mutate(
      gene_pval = ifelse("gene_fisher_p" %in% colnames(.), gene_fisher_p, NA),
      pathway_pval = ifelse("pathway_fisher_p" %in% colnames(.), pathway_fisher_p, NA),
      signature_strength = get_signature_strength(.),
      cluster_info = get_cluster_info(.)
    ) %>%
    filter(!is.na(gene_pval) | !is.na(pathway_pval))
  
  if (nrow(plot_data) == 0) {
    return(ggplot() + 
           geom_text(aes(x = 0.5, y = 0.5, label = "No p-value data available"), 
                     size = 6, color = "gray50") +
           theme_void())
  }
  
  # Transform p-values to -log10 scale
  plot_data <- plot_data %>%
    mutate(
      gene_neg_log10p = -log10(pmax(gene_pval, 1e-10, na.rm = TRUE)),
      pathway_neg_log10p = -log10(pmax(pathway_pval, 1e-10, na.rm = TRUE)),
      significant_gene = gene_pval < 0.05,
      significant_pathway = pathway_pval < 0.05,
      significance_category = case_when(
        significant_gene & significant_pathway ~ "Both Significant",
        significant_gene & !significant_pathway ~ "Gene Only",
        !significant_gene & significant_pathway ~ "Pathway Only", 
        TRUE ~ "Neither Significant"
      )
    )
  
  # Create base plot
  p <- ggplot(plot_data, aes(x = gene_neg_log10p, y = pathway_neg_log10p)) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_vline(xintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_point(aes(color = significance_category, size = signature_strength), alpha = 0.7) +
    scale_color_manual(values = c(
      "Both Significant" = "#2E8B57",     # Sea green
      "Gene Only" = "#4682B4",           # Steel blue  
      "Pathway Only" = "#FF6347",        # Tomato
      "Neither Significant" = "#D3D3D3"  # Light gray
    )) +
    scale_size_continuous(range = c(2, 8), name = "Signature\nStrength") +
    labs(
      title = "Gene vs Pathway Overlap Significance",
      subtitle = "Each point represents a gene pair comparison (MAST vs CRISPRi)",
      x = "Gene Overlap Significance (-log10 p-value)",
      y = "Pathway Overlap Significance (-log10 p-value)",
      color = "Significance Category",
      caption = "Dashed lines show p=0.05 threshold. Pathway overlap often more meaningful for cross-method comparisons."
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40"),
      legend.position = "right"
    )
  
  if (interactive) {
    # Create interactive version with hover information
    plot_data$hover_text <- paste0(
      "<b>", plot_data$gene_pair, "</b><br>",
      "Cluster: ", ifelse(is.na(plot_data$cluster_info), "Pan-cluster", plot_data$cluster_info), "<br>",
      "Gene p-value: ", format(plot_data$gene_pval, digits = 3), "<br>",
      "Pathway p-value: ", format(plot_data$pathway_pval, digits = 3), "<br>",
      "Signature Strength: ", format(plot_data$signature_strength, digits = 2)
    )
    
    p <- p + aes(text = hover_text)
    return(ggplotly(p, tooltip = "text"))
  }
  
  return(p)
}

#' Create Interactive Signature Heatmap
#'
#' @param signature_data Data frame with signature results
#' @param metric Character, metric to display ("signature_strength", "gene_fisher_p", "pathway_fisher_p")
#' @param cluster_filter Character vector of clusters to include (NULL for all)
#' @return plotly heatmap object
#' @export
create_interactive_signature_heatmap <- function(signature_data, 
                                                metric = "signature_strength",
                                                cluster_filter = NULL) {
  
  # Validate input data
  validate_signature_data(signature_data)
  
  if (nrow(signature_data) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No signature data available"))
  }
  
  # Filter by clusters if specified
  if (!is.null(cluster_filter)) {
    cluster_col <- get_cluster_info(signature_data)
    signature_data <- signature_data[cluster_col %in% cluster_filter, ]
  }
  
  # Prepare data for heatmap
  plot_data <- signature_data %>%
    mutate(
      cluster_info = ifelse("cluster" %in% colnames(signature_data), cluster, 
                           ifelse("cluster_id" %in% colnames(signature_data), cluster_id, "Unknown")),
      metric_value = ifelse("signature_strength" %in% colnames(signature_data), signature_strength, 
                           ifelse("strength" %in% colnames(signature_data), strength, 1))
    )
  
  # Handle different metrics
  if (metric == "gene_fisher_p" && "gene_fisher_p" %in% colnames(signature_data)) {
    plot_data$metric_value <- -log10(pmax(signature_data$gene_fisher_p, 1e-10))
    metric_label <- "Gene Overlap -log10(p-value)"
  } else if (metric == "pathway_fisher_p" && "pathway_fisher_p" %in% colnames(signature_data)) {
    plot_data$metric_value <- -log10(pmax(signature_data$pathway_fisher_p, 1e-10))
    metric_label <- "Pathway Overlap -log10(p-value)"
  } else {
    metric_label <- "Signature Strength"
  }
  
  # Create matrix for heatmap
  heatmap_matrix <- plot_data %>%
    dplyr::select(gene_pair, cluster_info, metric_value) %>%
    tidyr::pivot_wider(names_from = cluster_info, values_from = metric_value, values_fill = 0) %>%
    tibble::column_to_rownames("gene_pair") %>%
    as.matrix()
  
  if (nrow(heatmap_matrix) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No data available for heatmap"))
  }
  
  # Create interactive heatmap
  plotly::plot_ly(
    z = heatmap_matrix,
    x = colnames(heatmap_matrix),
    y = rownames(heatmap_matrix),
    type = "heatmap",
    colorscale = "Viridis",
    hovertemplate = paste0(
      "<b>Gene Pair:</b> %{y}<br>",
      "<b>Cluster:</b> %{x}<br>",
      "<b>", metric_label, ":</b> %{z:.2f}<br>",
      "<extra></extra>"
    )
  ) %>%
  plotly::layout(
    title = paste("Signature", metric_label, "Across Gene Pairs and Clusters"),
    xaxis = list(title = "Cluster"),
    yaxis = list(title = "Gene Pair (MAST vs CRISPRi)")
  )
}

#' Enhanced Interactive Signature Heatmap with Full UI Controls
#'
#' @param signature_data Data frame with signature analysis results
#' @param metric Character, metric to display (signature_strength, gene_overlap_count, gene_fisher_p, gene_jaccard)
#' @param cluster_filter Character vector, clusters to include (NULL for all)
#' @param clustering Character, clustering option: "both", "row", "column", "none"
#' @param color_scale Character, color scale: "viridis", "RdBu", "Reds", "Blues"
#' @return plotly object
#' @export
create_interactive_signature_heatmap_enhanced <- function(signature_data, 
                                                         metric = "signature_strength",
                                                         cluster_filter = NULL,
                                                         clustering = "both",
                                                         color_scale = "viridis") {
  
  # Validate input data
  validate_signature_data(signature_data)
  
  if (nrow(signature_data) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No signature data available"))
  }
  
  # Filter by clusters if specified
  if (!is.null(cluster_filter)) {
    cluster_col <- get_cluster_info(signature_data)
    signature_data <- signature_data[cluster_col %in% cluster_filter, ]
  }
  
  # Prepare data for heatmap
  plot_data <- signature_data %>%
    mutate(
      cluster_info = ifelse("cluster" %in% colnames(signature_data), cluster, 
                           ifelse("cluster_id" %in% colnames(signature_data), cluster_id, "Unknown")),
      metric_value = ifelse("signature_strength" %in% colnames(signature_data), signature_strength, 
                           ifelse("strength" %in% colnames(signature_data), strength, 1))
    )
  
  # Handle different metrics
  if (metric == "gene_fisher_p" && "gene_fisher_p" %in% colnames(signature_data)) {
    plot_data$metric_value <- -log10(pmax(signature_data$gene_fisher_p, 1e-10))
    metric_label <- "Gene Overlap -log10(p-value)"
  } else if (metric == "pathway_fisher_p" && "pathway_fisher_p" %in% colnames(signature_data)) {
    plot_data$metric_value <- -log10(pmax(signature_data$pathway_fisher_p, 1e-10))
    metric_label <- "Pathway Overlap -log10(p-value)"
  } else if (metric == "gene_jaccard" && "gene_jaccard" %in% colnames(signature_data)) {
    plot_data$metric_value <- signature_data$gene_jaccard
    metric_label <- "Gene Jaccard Index"
  } else if (metric == "gene_overlap_count" && "gene_overlap_count" %in% colnames(signature_data)) {
    plot_data$metric_value <- signature_data$gene_overlap_count
    metric_label <- "Gene Overlap Count"
  } else {
    # Default to signature strength
    if ("signature_strength" %in% colnames(signature_data)) {
      plot_data$metric_value <- signature_data$signature_strength
    } else if ("strength" %in% colnames(signature_data)) {
      plot_data$metric_value <- signature_data$strength
    } else {
      plot_data$metric_value <- 1  # fallback value
    }
    metric_label <- "Signature Strength"
  }
  
  # Create matrix
  heatmap_matrix <- plot_data %>%
    dplyr::select(gene_pair, cluster_info, metric_value) %>%
    tidyr::pivot_wider(names_from = cluster_info, values_from = metric_value, values_fill = 0) %>%
    tibble::column_to_rownames("gene_pair") %>%
    as.matrix()
  
  if (nrow(heatmap_matrix) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No data available for heatmap"))
  }
  
  # Handle clustering
  row_order <- rownames(heatmap_matrix)
  col_order <- colnames(heatmap_matrix)
  
  if (clustering %in% c("both", "row")) {
    if (nrow(heatmap_matrix) > 1) {
      tryCatch({
        row_dist <- dist(heatmap_matrix)
        row_hclust <- hclust(row_dist)
        row_order <- rownames(heatmap_matrix)[row_hclust$order]
      }, error = function(e) {
        cat("[HEATMAP] Row clustering failed:", e$message, "\n")
      })
    }
  }
  
  if (clustering %in% c("both", "column")) {
    if (ncol(heatmap_matrix) > 1) {
      tryCatch({
        col_dist <- dist(t(heatmap_matrix))
        col_hclust <- hclust(col_dist)
        col_order <- colnames(heatmap_matrix)[col_hclust$order]
      }, error = function(e) {
        cat("[HEATMAP] Column clustering failed:", e$message, "\n")
      })
    }
  }
  
  # Reorder matrix
  heatmap_matrix <- heatmap_matrix[row_order, col_order, drop = FALSE]
  
  # Set color scale
  plotly_colorscale <- switch(color_scale,
    "viridis" = "Viridis",
    "RdBu" = list(c(0, "blue"), c(0.5, "white"), c(1, "red")),
    "Reds" = "Reds", 
    "Blues" = "Blues",
    "Viridis"  # fallback
  )
  
  # Create interactive heatmap
  plotly::plot_ly(
    z = heatmap_matrix,
    x = colnames(heatmap_matrix),
    y = rownames(heatmap_matrix),
    type = "heatmap",
    colorscale = plotly_colorscale,
    hovertemplate = paste0(
      "<b>Gene Pair:</b> %{y}<br>",
      "<b>Cluster:</b> %{x}<br>",
      "<b>", metric_label, ":</b> %{z:.2f}<br>",
      "<extra></extra>"
    )
  ) %>%
  plotly::layout(
    title = paste("Signature", metric_label, "Across Gene Pairs and Clusters"),
    xaxis = list(title = "Cluster"),
    yaxis = list(title = "Gene Pair (MAST vs CRISPRi)")
  )
}

#' Create Gene Pair Multi-Metric Dashboard
#'
#' @param signature_data Data frame with signature results
#' @param selected_gene_pair Character, specific gene pair to focus on (NULL for all)
#' @return plotly subplot object
#' @export
create_gene_pair_multi_metric_dashboard <- function(signature_data, selected_gene_pair = NULL) {
  
  # Validate input data
  validate_signature_data(signature_data)
  
  if (!is.null(selected_gene_pair)) {
    signature_data <- signature_data[signature_data$gene_pair == selected_gene_pair, ]
  }
  
  if (nrow(signature_data) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No data available for selected gene pair"))
  }
  
  # Prepare data
  plot_data <- signature_data %>%
    mutate(
      signature_strength = get_signature_strength(.),
      cluster_info = get_cluster_info(.),
      gene_pval = ifelse("gene_fisher_p" %in% colnames(.), gene_fisher_p, NA),
      pathway_pval = ifelse("pathway_fisher_p" %in% colnames(.), pathway_fisher_p, NA),
      gene_jaccard = ifelse("gene_jaccard" %in% colnames(.), gene_jaccard, NA),
      pathway_jaccard = ifelse("pathway_jaccard" %in% colnames(.), pathway_jaccard, NA)
    )
  
  # Plot 1: Signature Strength by Cluster
  p1 <- plotly::plot_ly(plot_data, x = ~cluster_info, y = ~signature_strength, type = "bar",
                       name = "Signature Strength", marker = list(color = "#2E8B57")) %>%
    plotly::layout(title = "Signature Strength by Cluster", xaxis = list(title = "Cluster"))
  
  # Plot 2: P-values comparison
  pval_data <- plot_data %>%
    select(cluster_info, gene_pval, pathway_pval) %>%
    tidyr::pivot_longer(cols = c(gene_pval, pathway_pval), names_to = "pval_type", values_to = "pvalue") %>%
    filter(!is.na(pvalue)) %>%
    mutate(neg_log10p = -log10(pmax(pvalue, 1e-10)))
  
  p2 <- plotly::plot_ly(pval_data, x = ~cluster_info, y = ~neg_log10p, color = ~pval_type,
                       type = "scatter", mode = "markers+lines") %>%
    plotly::layout(title = "P-value Comparison", xaxis = list(title = "Cluster"), 
                  yaxis = list(title = "-log10(p-value)"))
  
  # Plot 3: Jaccard Index comparison  
  jaccard_data <- plot_data %>%
    select(cluster_info, gene_jaccard, pathway_jaccard) %>%
    tidyr::pivot_longer(cols = c(gene_jaccard, pathway_jaccard), names_to = "jaccard_type", values_to = "jaccard") %>%
    filter(!is.na(jaccard))
  
  p3 <- plotly::plot_ly(jaccard_data, x = ~cluster_info, y = ~jaccard, color = ~jaccard_type,
                       type = "scatter", mode = "markers+lines") %>%
    plotly::layout(title = "Jaccard Similarity Index", xaxis = list(title = "Cluster"), 
                  yaxis = list(title = "Jaccard Index"))
  
  # Combine plots
  plotly::subplot(p1, p2, p3, nrows = 3, shareX = TRUE, titleY = TRUE) %>%
    plotly::layout(title = paste("Multi-Metric Dashboard:", selected_gene_pair %||% "All Gene Pairs"))
}

#' Create Pathway Category Bubble Chart
#'
#' @param signature_data Data frame with signature results  
#' @param enrichment_data Data frame with enrichment results (optional, for pathway categorization)
#' @return plotly bubble chart object
#' @export
create_pathway_category_bubble_chart <- function(signature_data, enrichment_data = NULL) {
  
  # Validate input data
  validate_signature_data(signature_data)
  
  if (nrow(signature_data) == 0) {
    return(plotly::plot_ly() %>% 
           plotly::add_text(x = 0.5, y = 0.5, text = "No signature data available"))
  }
  
  # Prepare data for bubble chart
  bubble_data <- signature_data %>%
    mutate(
      signature_strength = get_signature_strength(.),
      cluster_info = get_cluster_info(.),
      pathway_pval = ifelse("pathway_fisher_p" %in% colnames(.), pathway_fisher_p, 1),
      pathway_significance = ifelse(pathway_pval < 0.05, "Significant", "Not Significant")
    ) %>%
    group_by(gene_pair, pathway_significance) %>%
    summarise(
      avg_strength = mean(signature_strength, na.rm = TRUE),
      cluster_count = n(),
      min_pval = min(pathway_pval, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      bubble_size = sqrt(cluster_count) * 10,  # Scale for visibility
      neg_log10p = -log10(pmax(min_pval, 1e-10))
    )
  
  # Create bubble chart
  plotly::plot_ly(bubble_data, 
                 x = ~avg_strength, 
                 y = ~neg_log10p,
                 size = ~bubble_size,
                 color = ~pathway_significance,
                 colors = c("Significant" = "#FF6347", "Not Significant" = "#D3D3D3"),
                 text = ~gene_pair,
                 hovertemplate = paste0(
                   "<b>Gene Pair:</b> %{text}<br>",
                   "<b>Avg Signature Strength:</b> %{x:.2f}<br>",
                   "<b>Min Pathway p-value:</b> %{customdata:.2e}<br>",
                   "<b>Clusters:</b> %{marker.size}<br>",
                   "<extra></extra>"
                 ),
                 customdata = ~min_pval) %>%
  plotly::layout(
    title = "Pathway Significance vs Signature Strength",
    xaxis = list(title = "Average Signature Strength"),
    yaxis = list(title = "Pathway Overlap Significance (-log10 p-value)"),
    showlegend = TRUE
  )
}

#' Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a