#!/usr/bin/env Rscript

#' Unified Functional Enrichment Visualization Package
#' 
#' This file combines all heatmap visualization functions into a single, unified implementation:
#' - Basic flexible heatmaps with customizable metrics (V1)
#' - Enhanced intersection logic for method comparisons (V2) 
#' - GSEA-specific visualization with NES support (V3)
#' - Direction annotations for all versions
#' - Convenience wrapper functions for common use cases
#'
#' The unified approach provides a consistent interface while maintaining the specialized
#' functionality of each version. Direction annotations are included by default to help
#' users quickly identify regulation patterns.

# =============================================================================
# LOAD REQUIRED LIBRARIES
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(stringr)
  library(grid)
})

# =============================================================================
# CONFIGURATION AND CONSTANTS
# =============================================================================

#' Global configuration for enrichment visualization
CONFIG <- list(
  HEATMAP_COLUMN_GENES = c("ATP13A2", "DNAJC6", "FBXO7", "GBA", "LRRK2", 
                         "PARK7", "PINK1", "PRKN", "SNCA_A30P", "SNCA_A53T", 
                         "SYNJ1", "VPS13C_A444P", "VPS13C_W395C"),
  PLOT_ENRICHMENT_TYPES = c("GO_BP", "GO_CC", "GO_MF", "GO_GOALL", 
                          "KEGG", "Reactome", "WikiPathways", "GSEA"),
  ENRICHMENT_TYPES_IN_DATA = c("GO_BP", "GO_CC", "GO_MF", "GO_GOALL", "KEGG", 
                             "Reactome", "WikiPathways", "STRING", "GSEA"),
  METHOD_DISPLAY_NAMES = c(MAST = "iSCORE-PD (MAST)", MixScale = "CRISPRi (MixScale)", 
                         all="All Methods", intersection="Intersection", union="Union"),
  DIRECTION_COLORS = c(UP = "#E41A1C", DOWN = "#377EB8", ALL = "#4DAF4A"),  # Red, Blue, Green
  MAX_HEATMAP_ROWS = 500,
  DEFAULT_P_THRESHOLD = 0.05,
  DEFAULT_GSEA_P_THRESHOLD = 0.25
)

# =============================================================================
# MAIN UNIFIED HEATMAP FUNCTION
# =============================================================================

#' Create Unified Enrichment Heatmap
#' 
#' A comprehensive heatmap function that combines all visualization approaches:
#' - Standard enrichment metrics (p.adjust, FoldEnrichment, zScore)  
#' - GSEA-specific NES visualization with divergent colors
#' - Method intersection/union logic with detailed breakdown
#' - Direction annotations showing UP/DOWN/ALL regulation
#' - Flexible customization options for colors, dimensions, and styling
#' 
#' @param data_full Data frame containing enrichment results with required columns
#' @param heatmap_type Character. "standard" for regular metrics, "gsea" for NES visualization
#' @param metric_type Character. "p.adjust", "FoldEnrichment", "zScore", or "NES" (for GSEA)
#' @param cluster_val Character. Cluster to analyze (e.g., "cluster_1") or NULL for all clusters
#' @param method_val Character. "all", "MAST", "MixScale", "intersection", "union"
#' @param enrichment_types_to_include Character vector. Enrichment types to include
#' @param term_selection_strategy Character. "top_n_overall" or "all_significant"
#' @param max_terms_overall Integer. Maximum number of terms when using "top_n_overall"
#' @param p_threshold Numeric. P-value threshold (default: 0.05, 0.25 for GSEA)
#' @param min_nes_magnitude Numeric. Minimum absolute NES value for GSEA (default: 1.0)
#' @param show_direction_annotation Logical. Whether to show direction annotation (default: TRUE)
#' @param show_method_breakdown Logical. Whether to log detailed method breakdown (default: TRUE)
#' @param color_scaling_strategy Character. Color scaling method
#' @param custom_title Character. Custom title override (default: auto-generated)
#' @param log_messages Logical. Whether to print log messages (default: TRUE)
#' @param ... Additional styling parameters
#'
#' @return ComplexHeatmap object that can be drawn with draw()

create_unified_enrichment_heatmap <- function(data_full,
                                            heatmap_type = "standard",  # "standard" or "gsea"
                                            metric_type = "p.adjust",
                                            cluster_val = "cluster_1",
                                            method_val = "all",
                                            enrichment_types_to_include = "GO_BP",
                                            term_selection_strategy = "top_n_overall",
                                            max_terms_overall = 25,
                                            p_threshold = NULL,  # Auto-determined based on heatmap_type
                                            min_nes_magnitude = 1.0,
                                            show_direction_annotation = TRUE,
                                            show_method_breakdown = TRUE,
                                            color_scaling_strategy = "default_complexheatmap",
                                            custom_title = NULL,
                                            log_messages = TRUE,
                                            ...) {
  
  # =============================================================================
  # PARAMETER VALIDATION AND SETUP
  # =============================================================================
  
  # Logging function
  log_message_fn <- function(msg, level = "INFO") {
    if (log_messages) {
      cat(paste0("[UNIFIED-", level, "] ", msg, "\n"))
    }
  }
  
  # Validate heatmap type
  valid_heatmap_types <- c("standard", "gsea")
  if (!heatmap_type %in% valid_heatmap_types) {
    stop("heatmap_type must be one of: ", paste(valid_heatmap_types, collapse = ", "))
  }
  
  # Set default p_threshold based on heatmap type
  if (is.null(p_threshold)) {
    p_threshold <- if (heatmap_type == "gsea") CONFIG$DEFAULT_GSEA_P_THRESHOLD else CONFIG$DEFAULT_P_THRESHOLD
  }
  
  # Validate metric type based on heatmap type
  if (heatmap_type == "standard") {
    valid_metrics <- c("p.adjust", "FoldEnrichment", "zScore")
    if (!metric_type %in% valid_metrics) {
      stop("For standard heatmaps, metric_type must be one of: ", paste(valid_metrics, collapse = ", "))
    }
  } else if (heatmap_type == "gsea") {
    metric_type <- "NES"  # Force NES for GSEA heatmaps
    enrichment_types_to_include <- "GSEA"  # Force GSEA enrichment type
  }
  
  # Check for required direction column if direction annotation is requested
  if (show_direction_annotation && !"direction" %in% colnames(data_full)) {
    warning("Direction annotation requested but 'direction' column not found. Disabling direction annotation.")
    show_direction_annotation <- FALSE
  }
  
  log_message_fn(paste("Creating", heatmap_type, "heatmap using metric:", metric_type))
  log_message_fn(paste("Method:", method_val, "- Enrichment:", paste(enrichment_types_to_include, collapse = ", ")))
  log_message_fn(paste("P-value threshold:", p_threshold))
  
  # =============================================================================
  # DATA FILTERING
  # =============================================================================
  
  # Initial filtering by enrichment type
  heatmap_data_filtered <- data_full %>% 
    filter(enrichment_type %in% enrichment_types_to_include)
  
  # Cluster filtering
  if (!is.null(cluster_val)) { 
    heatmap_data_filtered <- heatmap_data_filtered %>% filter(cluster == cluster_val) 
  }
  
  # Method filtering with enhanced intersection logic
  if (method_val == "intersection") {
    log_message_fn("=== APPLYING TRUE INTERSECTION LOGIC ===")
    
    # Find terms existing in both methods
    term_methods <- heatmap_data_filtered %>% 
      group_by(ID) %>% 
      summarise(methods_count = n_distinct(method), 
               methods = paste(sort(unique(method)), collapse = "+"), .groups = "drop") %>% 
      filter(methods_count == 2)
    
    log_message_fn(paste("Found", nrow(term_methods), "terms in both methods"))
    heatmap_data_filtered <- heatmap_data_filtered %>% filter(ID %in% term_methods$ID)
    
    # Find genes with data in both methods
    gene_methods <- heatmap_data_filtered %>% 
      group_by(mutation_perturbation) %>%
      summarise(methods_count = n_distinct(method), .groups = "drop") %>%
      filter(methods_count == 2)
    
    log_message_fn(paste("Found", nrow(gene_methods), "genes with data in both methods"))
    if (show_method_breakdown) {
      log_message_fn(paste("Intersection genes:", paste(gene_methods$mutation_perturbation, collapse = ", ")))
    }
    
    heatmap_data_filtered <- heatmap_data_filtered %>% 
      filter(mutation_perturbation %in% gene_methods$mutation_perturbation)
    
  } else if (method_val == "union") {
    log_message_fn("=== APPLYING UNION LOGIC ===")
    # Union keeps all data - no additional filtering needed
    
  } else if (method_val != "all") {
    log_message_fn(paste("=== FILTERING TO SINGLE METHOD:", method_val, "==="))
    heatmap_data_filtered <- heatmap_data_filtered %>% filter(method == method_val)
  }
  
  # Significance and metric filtering
  if (heatmap_type == "gsea") {
    # GSEA-specific filtering
    heatmap_data_filtered <- heatmap_data_filtered %>% 
      filter(p.adjust < p_threshold, !is.na(p.adjust), !is.na(NES)) %>%
      filter(abs(NES) >= min_nes_magnitude)
  } else {
    # Standard enrichment filtering
    heatmap_data_filtered <- heatmap_data_filtered %>% 
      filter(p.adjust < p_threshold, !is.na(p.adjust), !is.na(.data[[metric_type]]))
  }
  
  if (nrow(heatmap_data_filtered) == 0) {
    log_message_fn("No data available after filtering", "WARNING")
    return(NULL)
  }
  
  log_message_fn(paste("Filtered data contains", nrow(heatmap_data_filtered), "rows"))
  
  # =============================================================================
  # TERM SELECTION
  # =============================================================================
  
  if (term_selection_strategy == "top_n_overall") {
    max_terms <- min(max_terms_overall, CONFIG$MAX_HEATMAP_ROWS)
    
    if (heatmap_type == "gsea") {
      # Select by NES magnitude for GSEA
      selected_terms <- heatmap_data_filtered %>% 
        group_by(ID) %>% 
        summarise(max_nes_magnitude = max(abs(NES), na.rm = TRUE), .groups = "drop") %>% 
        arrange(desc(max_nes_magnitude)) %>% 
        head(max_terms) %>% 
        pull(ID)
      
      significant_terms_for_heatmap <- heatmap_data_filtered %>% 
        filter(ID %in% selected_terms)
      
    } else {
      # Select by metric value for standard enrichment
      if (metric_type == "p.adjust") {
        selected_descriptions <- heatmap_data_filtered %>% 
          group_by(Description) %>% 
          summarise(min_metric = min(.data[[metric_type]], na.rm = TRUE), .groups = "drop") %>% 
          arrange(min_metric) %>% 
          head(max_terms) %>% 
          pull(Description)
      } else {
        selected_descriptions <- heatmap_data_filtered %>% 
          group_by(Description) %>% 
          summarise(max_metric = max(.data[[metric_type]], na.rm = TRUE), .groups = "drop") %>% 
          arrange(desc(max_metric)) %>% 
          head(max_terms) %>% 
          pull(Description)
      }
      
      significant_terms_for_heatmap <- heatmap_data_filtered %>% 
        filter(Description %in% selected_descriptions)
    }
    
  } else if (term_selection_strategy == "all_significant") {
    significant_terms_for_heatmap <- heatmap_data_filtered
    
    # Apply row limit if needed
    if (heatmap_type == "gsea") {
      unique_terms <- n_distinct(significant_terms_for_heatmap$ID)
    } else {
      unique_terms <- n_distinct(significant_terms_for_heatmap$Description)
    }
    
    if (unique_terms > CONFIG$MAX_HEATMAP_ROWS) {
      log_message_fn(paste("Too many terms (", unique_terms, "), limiting to", CONFIG$MAX_HEATMAP_ROWS))
      # Apply same selection logic as top_n_overall
      significant_terms_for_heatmap <- significant_terms_for_heatmap %>%
        head(CONFIG$MAX_HEATMAP_ROWS * 10)  # Rough estimate to get enough rows
    }
  } else {
    stop("Unknown term_selection_strategy: ", term_selection_strategy)
  }
  
  if (nrow(significant_terms_for_heatmap) < 1) {
    log_message_fn("Not enough data for heatmap after term selection", "WARNING")
    return(NULL)
  }
  
  # Count unique terms
  unique_term_count <- if (heatmap_type == "gsea") {
    n_distinct(significant_terms_for_heatmap$ID)
  } else {
    n_distinct(significant_terms_for_heatmap$Description)
  }
  log_message_fn(paste("Selected", unique_term_count, "terms for heatmap"))
  
  # =============================================================================
  # MATRIX CREATION WITH DIRECTION TRACKING
  # =============================================================================
  
  # Setup column order (all CONFIG genes)
  all_heatmap_cols_ordered <- if (is.null(cluster_val)) {
    interaction_df <- expand.grid(mutation_perturbation = CONFIG$HEATMAP_COLUMN_GENES, 
                                cluster = sort(unique(data_full$cluster)), 
                                KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
    paste(interaction_df$mutation_perturbation, interaction_df$cluster, sep="_")
  } else {
    CONFIG$HEATMAP_COLUMN_GENES
  }
  all_heatmap_cols_ordered <- unique(all_heatmap_cols_ordered)
  
  # Transform values for plotting
  transform_value <- function(metric_val, metric_name) {
    if (metric_name == "p.adjust") {
      return(-log10(metric_val + 1e-10))
    } else if (metric_name == "FoldEnrichment") {
      return(metric_val)
    } else if (metric_name == "zScore") {
      return(metric_val)
    } else if (metric_name == "NES") {
      return(metric_val)  # Use NES directly
    }
  }
  
  # Create aggregated values
  if (heatmap_type == "gsea") {
    # GSEA uses ID and NES
    aggregated_values_for_pivot <- significant_terms_for_heatmap %>%
      mutate(
        matrix_col_label = if (is.null(cluster_val)) paste(mutation_perturbation, cluster, sep = "_") else mutation_perturbation,
        value_col = NES,
        term_label = ID
      ) %>%
      filter(is.finite(value_col)) %>% 
      group_by(term_label, matrix_col_label) %>% 
      summarise(agg_value_col = ifelse(all(value_col < 0), min(value_col, na.rm = TRUE), max(value_col, na.rm = TRUE)), 
               .groups = "drop") %>%
      mutate(agg_value_col = ifelse(is.infinite(agg_value_col) | is.na(agg_value_col), 0, agg_value_col))
  } else {
    # Standard enrichment uses Description and transformed metric
    aggregated_values_for_pivot <- significant_terms_for_heatmap %>%
      mutate(
        matrix_col_label = if (is.null(cluster_val)) paste(mutation_perturbation, cluster, sep = "_") else mutation_perturbation,
        value_col = transform_value(.data[[metric_type]], metric_type),
        term_label = Description
      ) %>%
      filter(is.finite(value_col)) %>% 
      group_by(term_label, matrix_col_label) %>% 
      summarise(agg_value_col = max(value_col, na.rm = TRUE), .groups = "drop") %>%
      mutate(agg_value_col = ifelse(is.infinite(agg_value_col) | is.na(agg_value_col), 0, agg_value_col))
  }
  
  # Create direction mapping if direction annotation is enabled
  direction_mapping <- NULL
  row_directions <- NULL
  if (show_direction_annotation) {
    if (heatmap_type == "gsea") {
      direction_mapping <- significant_terms_for_heatmap %>%
        group_by(ID) %>%
        summarise(
          direction = if (n_distinct(direction) == 1) {
            first(direction)
          } else {
            dir_counts <- table(direction)
            if (length(dir_counts) > 0) names(which.max(dir_counts)) else "ALL"
          },
          .groups = "drop"
        )
      names(direction_mapping)[names(direction_mapping) == "ID"] <- "term_label"
    } else {
      direction_mapping <- significant_terms_for_heatmap %>%
        group_by(Description) %>%
        summarise(
          direction = if (n_distinct(direction) == 1) {
            first(direction)
          } else {
            dir_counts <- table(direction)
            if (length(dir_counts) > 0) names(which.max(dir_counts)) else "ALL"
          },
          .groups = "drop"
        )
      names(direction_mapping)[names(direction_mapping) == "Description"] <- "term_label"
    }
    
    log_message_fn(paste("Created direction mapping for", nrow(direction_mapping), "terms"))
  }
  
  # Create matrix
  all_terms_present <- unique(aggregated_values_for_pivot$term_label)
  matrix_values_df <- tidyr::crossing(term_label = all_terms_present, 
                                    matrix_col_label = all_heatmap_cols_ordered) %>%
    left_join(aggregated_values_for_pivot, by = c("term_label", "matrix_col_label")) %>%
    mutate(agg_value_col = ifelse(is.na(agg_value_col), 0, agg_value_col)) %>%
    pivot_wider(names_from = matrix_col_label, values_from = agg_value_col)
  
  # Join direction information if available
  if (show_direction_annotation && !is.null(direction_mapping)) {
    matrix_values_df <- matrix_values_df %>%
      left_join(direction_mapping, by = "term_label") %>%
      mutate(direction = ifelse(is.na(direction), "ALL", direction))
    row_directions <- matrix_values_df$direction
  }
  
  # Create display names and matrix
  if (heatmap_type == "gsea") {
    # For GSEA, try to get Description if available, otherwise use ID
    display_names <- matrix_values_df$term_label
    if ("Description" %in% colnames(significant_terms_for_heatmap)) {
      id_to_desc <- significant_terms_for_heatmap %>%
        select(ID, Description) %>%
        distinct() %>%
        group_by(ID) %>%
        slice(1) %>%
        ungroup()
      
      temp_mapping <- data.frame(term_label = matrix_values_df$term_label) %>%
        left_join(id_to_desc, by = c("term_label" = "ID"))
      
      display_names <- ifelse(!is.na(temp_mapping$Description), 
                            temp_mapping$Description, 
                            temp_mapping$term_label)
    }
  } else {
    display_names <- matrix_values_df$term_label
  }
  
  display_term_row_names <- make.unique(str_trunc(display_names, 80), sep = "_")
  cols_for_matrix <- setdiff(colnames(matrix_values_df), c("term_label", "direction"))
  heatmap_matrix <- as.matrix(matrix_values_df[, cols_for_matrix])
  rownames(heatmap_matrix) <- display_term_row_names
  heatmap_matrix[is.na(heatmap_matrix) | is.infinite(heatmap_matrix)] <- 0
  
  if (nrow(heatmap_matrix) < 2 || ncol(heatmap_matrix) < 1) {
    log_message_fn(sprintf("Matrix too small (R=%d,C=%d)", nrow(heatmap_matrix), ncol(heatmap_matrix)), "WARNING")
    return(NULL)
  }
  
  log_message_fn(paste("Final matrix dimensions:", nrow(heatmap_matrix), "x", ncol(heatmap_matrix)))
  
  # =============================================================================
  # COLOR FUNCTION SETUP
  # =============================================================================
  
  if (heatmap_type == "gsea") {
    # GSEA uses divergent colors for NES
    max_abs_nes <- max(abs(heatmap_matrix), na.rm = TRUE)
    if (max_abs_nes == 0) max_abs_nes <- 1
    
    col_fun <- colorRamp2(
      c(-max_abs_nes, -min_nes_magnitude/2, 0, min_nes_magnitude/2, max_abs_nes),
      c("blue", "lightblue", "white", "pink", "red")
    )
  } else {
    # Standard enrichment uses sequential colors
    if (color_scaling_strategy == "fixed_user") {
      default_breaks <- c(0, 5, 10, 20)
      default_colors <- c("white", "yellow", "orange", "red")
      col_fun <- colorRamp2(default_breaks, default_colors)
    } else if (color_scaling_strategy == "adaptive_quantile") {
      significant_values_for_scale <- as.vector(heatmap_matrix[heatmap_matrix > 0 & is.finite(heatmap_matrix)])
      if (length(significant_values_for_scale) > 10) {
        breaks_adaptive <- quantile(significant_values_for_scale, probs = c(0, 0.25, 0.5, 0.75, 0.95, 1.0), na.rm = TRUE)
        breaks_adaptive <- unique(round(sort(breaks_adaptive), 1))
        if (length(breaks_adaptive) == 0 || breaks_adaptive[1] > 0) {
          breaks_adaptive <- c(0, breaks_adaptive)
        }
        breaks_adaptive <- unique(sort(breaks_adaptive))
        if (length(breaks_adaptive) < 2) {
          breaks_adaptive <- c(0, max(1, max(significant_values_for_scale, 0, na.rm=TRUE)))
        }
        num_colors_needed_adaptive <- length(breaks_adaptive)
        palette_name <- "YlOrRd"
        base_palette_colors <- RColorBrewer::brewer.pal(max(3, min(9, num_colors_needed_adaptive)), palette_name)
        if (num_colors_needed_adaptive > length(base_palette_colors)) {
          colors_for_breaks_adaptive <- colorRampPalette(base_palette_colors)(num_colors_needed_adaptive)
        } else {
          colors_for_breaks_adaptive <- base_palette_colors[1:num_colors_needed_adaptive]
        }
        col_fun <- colorRamp2(breaks_adaptive, colors_for_breaks_adaptive)
      } else {
        col_fun <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
      }
    } else {
      col_fun <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(100)
    }
  }
  
  # =============================================================================
  # TITLE GENERATION
  # =============================================================================
  
  if (!is.null(custom_title)) {
    column_title <- custom_title
  } else {
    if (heatmap_type == "gsea") {
      metric_display <- "GSEA NES"
    } else {
      metric_display <- switch(metric_type,
        "p.adjust" = "-log10(p.adjust)",
        "FoldEnrichment" = "Fold Enrichment",
        "zScore" = "Z-Score"
      )
    }
    
    etype_display <- paste(enrichment_types_to_include, collapse = ", ")
    cluster_display <- ifelse(is.null(cluster_val), "All Clusters", cluster_val)
    
    method_suffix <- if (method_val == "intersection") {
      paste("- True Intersection (", ncol(heatmap_matrix), " genes)")
    } else if (method_val == "union") {
      paste("- Union (", ncol(heatmap_matrix), " genes)")
    } else if (method_val != "all") {
      paste("-", method_val, "(", ncol(heatmap_matrix), " genes)")
    } else {
      ""
    }
    
    column_title <- paste(metric_display, "-", etype_display, cluster_display, method_suffix)
  }
  
  # Create heatmap legend name
  heatmap_name <- if (heatmap_type == "gsea") {
    "NES"
  } else {
    switch(metric_type,
      "p.adjust" = "-log10(p.adjust)",
      "FoldEnrichment" = "Fold Enrichment",
      "zScore" = "Z-Score"
    )
  }
  
  # =============================================================================
  # CREATE DIRECTION ANNOTATION
  # =============================================================================
  
  row_ha <- NULL
  if (show_direction_annotation && !is.null(row_directions)) {
    direction_colors <- CONFIG$DIRECTION_COLORS
    row_directions_factor <- factor(row_directions, levels = c("UP", "DOWN", "ALL"))
    
    row_ha <- rowAnnotation(
      Direction = row_directions_factor,
      col = list(Direction = direction_colors),
      annotation_legend_param = list(
        Direction = list(
          title = "Direction",
          title_gp = gpar(fontsize = 10),
          labels_gp = gpar(fontsize = 8)
        )
      ),
      width = unit(0.5, "cm")
    )
    
    log_message_fn("Added direction annotation to heatmap")
  }
  
  # =============================================================================
  # CREATE HEATMAP
  # =============================================================================
  
  # Calculate dimensions
  dynamic_heatmap_height_cm <- max(10, min(60, nrow(heatmap_matrix) * 0.20 + 5))
  dynamic_heatmap_width_cm <- max(8, min(30, ncol(heatmap_matrix) * 0.8 + 3))
  dynamic_font_size_rows <- max(3, 10 - nrow(heatmap_matrix) %/% 30)
  
  # Create legend parameters
  legend_params <- if (heatmap_type == "gsea") {
    max_abs_nes <- max(abs(heatmap_matrix), na.rm = TRUE)
    list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8),
      title_position = "topleft",
      at = c(-max_abs_nes, -min_nes_magnitude, 0, min_nes_magnitude, max_abs_nes),
      labels = c(paste0("↓", round(-max_abs_nes, 1)), paste0("↓", min_nes_magnitude), "0",
               paste0("↑", min_nes_magnitude), paste0("↑", round(max_abs_nes, 1)))
    )
  } else {
    list(
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8),
      title_position = "topleft"
    )
  }
  
  # Create heatmap
  ht_plot <- Heatmap(
    heatmap_matrix,
    name = heatmap_name,
    col = col_fun,
    cluster_rows = (nrow(heatmap_matrix) >= 2),
    cluster_columns = (ncol(heatmap_matrix) >= 2),
    show_row_names = TRUE,
    show_column_names = TRUE,
    row_names_gp = gpar(fontsize = dynamic_font_size_rows),
    column_names_gp = gpar(fontsize = 12, fontface = "bold"),
    column_names_rot = 45,
    column_names_side = "bottom",
    column_dend_side = "top",
    column_title = column_title,
    column_title_gp = gpar(fontsize = 10, lineheight = 0.9),
    row_dend_side = "right",
    row_names_side = "left",
    height = unit(dynamic_heatmap_height_cm, "cm"),
    width = unit(dynamic_heatmap_width_cm, "cm"),
    use_raster = (nrow(heatmap_matrix) * ncol(heatmap_matrix) > 5000),
    heatmap_legend_param = legend_params,
    right_annotation = row_ha
  )
  
  log_message_fn(paste("Successfully created", heatmap_type, "heatmap"))
  return(ht_plot)
}

# =============================================================================
# CONVENIENCE WRAPPER FUNCTIONS
# =============================================================================

#' Create P-value Heatmap with Direction
#' @description Convenience wrapper for p-value heatmaps with direction annotation
#' @param ... Parameters passed to create_unified_enrichment_heatmap()
#' @return ComplexHeatmap object
create_pvalue_heatmap <- function(...) {
  create_unified_enrichment_heatmap(heatmap_type = "standard", metric_type = "p.adjust", ...)
}

#' Create Fold Enrichment Heatmap with Direction
#' @description Convenience wrapper for fold enrichment heatmaps with direction annotation
#' @param ... Parameters passed to create_unified_enrichment_heatmap()
#' @return ComplexHeatmap object
create_fold_enrichment_heatmap <- function(...) {
  create_unified_enrichment_heatmap(heatmap_type = "standard", metric_type = "FoldEnrichment", ...)
}

#' Create Z-Score Heatmap with Direction
#' @description Convenience wrapper for z-score heatmaps with direction annotation
#' @param ... Parameters passed to create_unified_enrichment_heatmap()
#' @return ComplexHeatmap object
create_zscore_heatmap <- function(...) {
  create_unified_enrichment_heatmap(heatmap_type = "standard", metric_type = "zScore", ...)
}

#' Create GSEA NES Heatmap with Direction
#' @description Convenience wrapper for GSEA NES heatmaps with direction annotation
#' @param ... Parameters passed to create_unified_enrichment_heatmap()
#' @return ComplexHeatmap object
create_gsea_nes_heatmap <- function(...) {
  create_unified_enrichment_heatmap(heatmap_type = "gsea", ...)
}

# =============================================================================
# BATCH CREATION FUNCTIONS
# =============================================================================

#' Create Multiple Heatmaps for Comparison
#' @description Create multiple heatmaps side-by-side for comparison
#' @param data_full Data frame containing enrichment results
#' @param heatmap_configs List of configurations for each heatmap
#' @param output_file Character. Output file path (optional)
#' @param ... Additional parameters for draw()
#' @return Combined ComplexHeatmap object
create_comparison_heatmaps <- function(data_full, heatmap_configs, output_file = NULL, ...) {
  
  # Create heatmaps
  heatmaps <- list()
  for (i in seq_along(heatmap_configs)) {
    config <- heatmap_configs[[i]]
    config$data_full <- data_full
    
    ht <- do.call(create_unified_enrichment_heatmap, config)
    if (!is.null(ht)) {
      heatmaps[[i]] <- ht
    }
  }
  
  if (length(heatmaps) == 0) {
    warning("No valid heatmaps could be created")
    return(NULL)
  }
  
  # Combine heatmaps
  combined_ht <- Reduce("+", heatmaps)
  
  # Save if output file specified
  if (!is.null(output_file)) {
    if (grepl("\\.pdf$", output_file)) {
      pdf(output_file, width = length(heatmaps) * 8, height = 12)
      draw(combined_ht, ht_gap = unit(1, "cm"), ...)
      dev.off()
    } else if (grepl("\\.png$", output_file)) {
      png(output_file, width = length(heatmaps) * 800, height = 1200, res = 150)
      draw(combined_ht, ht_gap = unit(1, "cm"), ...)
      dev.off()
    }
  }
  
  return(combined_ht)
}

# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if (FALSE) {
  # Example usage (not run)
  
  # Load data
  data <- readRDS("all_enrichment_padj005_complete_with_direction.rds")
  
  # Basic usage examples
  
  # P-value heatmap with direction annotation
  ht_pval <- create_pvalue_heatmap(data, 
                                 enrichment_types_to_include = "GO_BP",
                                 cluster_val = "cluster_1",
                                 method_val = "MAST")
  draw(ht_pval)
  
  # GSEA NES heatmap with direction annotation
  ht_gsea <- create_gsea_nes_heatmap(data,
                                   cluster_val = "cluster_1", 
                                   method_val = "intersection")
  draw(ht_gsea)
  
  # Advanced usage with custom parameters
  
  # Custom heatmap with specific parameters
  ht_custom <- create_unified_enrichment_heatmap(
    data,
    heatmap_type = "standard",
    metric_type = "zScore",
    cluster_val = "cluster_0",
    method_val = "union",
    enrichment_types_to_include = c("GO_BP", "KEGG"),
    max_terms_overall = 40,
    show_direction_annotation = TRUE,
    color_scaling_strategy = "adaptive_quantile",
    custom_title = "Custom Z-Score Analysis"
  )
  draw(ht_custom)
  
  # Create comparison heatmaps
  
  comparison_configs <- list(
    list(
      heatmap_type = "standard",
      metric_type = "p.adjust",
      enrichment_types_to_include = "GO_BP",
      cluster_val = "cluster_1",
      method_val = "MAST",
      custom_title = "Statistical Significance\n-log10(p.adjust)"
    ),
    list(
      heatmap_type = "standard", 
      metric_type = "zScore",
      enrichment_types_to_include = "GO_BP",
      cluster_val = "cluster_1",
      method_val = "MAST",
      custom_title = "Effect Size\nZ-Score"
    ),
    list(
      heatmap_type = "gsea",
      cluster_val = "cluster_1",
      method_val = "MAST",
      custom_title = "GSEA Enrichment\nNES"
    )
  )
  
  combined_ht <- create_comparison_heatmaps(data, comparison_configs, "comparison_heatmaps.pdf")
  draw(combined_ht, 
       ht_gap = unit(1, "cm"),
       row_title = "Enrichment Analysis Comparison with Direction Annotations",
       row_title_gp = gpar(fontsize = 14, fontface = "bold"))
}