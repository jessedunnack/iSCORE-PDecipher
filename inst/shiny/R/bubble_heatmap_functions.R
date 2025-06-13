# Bubble Heatmap Functions for Shiny App
# Adapted from main codebase implementation

# Load packages conditionally
if (requireNamespace("ComplexHeatmap", quietly = TRUE) && 
    requireNamespace("circlize", quietly = TRUE)) {
  suppressPackageStartupMessages({
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
  })
}

#' Create a clustered bubble heatmap for enrichment analysis
#'
#' @param data Data frame containing enrichment results
#' @param max_terms Maximum number of terms to display
#' @param cluster_rows Whether to cluster rows (terms)
#' @param cluster_cols Whether to cluster columns (mutations)
#' @param color_scale Color scheme: "red", "blue", "green", or "viridis"
#' @param size_encoding What to encode with bubble size: "count" or "pvalue"
#' @return A ComplexHeatmap object
create_bubble_heatmap <- function(data,
                                 max_terms = 30,
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 color_scale = "red",
                                 size_encoding = "count",
                                 title = "Pathway Enrichment") {
  
  # Ensure we have the right columns
  if ("mutation_perturbation" %in% names(data)) {
    data$gene <- data$mutation_perturbation
  }
  
  # Check required columns
  required_cols <- c("gene", "ID", "Description", "p.adjust")
  missing_cols <- setdiff(required_cols, names(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process data
  data <- data %>%
    mutate(
      term_id = paste0(substr(Description, 1, 40), " (", ID, ")"),
      neg_log10_padj = -log10(p.adjust + 1e-50)
    )
  
  # Extract gene count
  if ("Count" %in% names(data) && !all(is.na(data$Count))) {
    data$gene_count <- as.numeric(data$Count)
  } else if ("number_of_genes" %in% names(data)) {
    data$gene_count <- as.numeric(data$number_of_genes)
  } else if ("geneID" %in% names(data)) {
    data$gene_count <- lengths(strsplit(as.character(data$geneID), ","))
  } else {
    data$gene_count <- 1
  }
  
  # Select top terms if needed
  if (n_distinct(data$ID) > max_terms) {
    top_terms <- data %>%
      group_by(ID) %>%
      summarise(mean_sig = mean(neg_log10_padj)) %>%
      arrange(desc(mean_sig)) %>%
      head(max_terms) %>%
      pull(ID)
    
    data <- data %>%
      filter(ID %in% top_terms)
  }
  
  # Create matrices
  genes <- unique(data$gene)
  terms <- unique(data$term_id)
  
  p_mat <- matrix(0, nrow = length(terms), ncol = length(genes))
  rownames(p_mat) <- terms
  colnames(p_mat) <- genes
  
  count_mat <- matrix(0, nrow = length(terms), ncol = length(genes))
  rownames(count_mat) <- terms
  colnames(count_mat) <- genes
  
  # Fill matrices
  for (i in seq_len(nrow(data))) {
    row_idx <- which(terms == data$term_id[i])
    col_idx <- which(genes == data$gene[i])
    
    if (length(row_idx) > 0 && length(col_idx) > 0) {
      p_mat[row_idx, col_idx] <- data$neg_log10_padj[i]
      count_mat[row_idx, col_idx] <- data$gene_count[i]
    }
  }
  
  # Remove empty rows and columns
  keep_rows <- rowSums(p_mat) > 0
  keep_cols <- colSums(p_mat) > 0
  
  p_mat <- p_mat[keep_rows, keep_cols, drop = FALSE]
  count_mat <- count_mat[keep_rows, keep_cols, drop = FALSE]
  
  # Clustering
  row_clust <- NULL
  col_clust <- NULL
  
  if (cluster_rows && nrow(p_mat) > 2) {
    row_clust <- hclust(dist(p_mat), method = "ward.D2")
  }
  
  if (cluster_cols && ncol(p_mat) > 2) {
    col_clust <- hclust(dist(t(p_mat)), method = "ward.D2")
  }
  
  # Color schemes
  if (size_encoding == "count") {
    color_schemes <- list(
      red = colorRamp2(c(0, 2, 5, 10), c("white", "gold", "orange", "darkred")),
      blue = colorRamp2(c(0, 2, 5, 10), c("white", "lightblue", "blue", "darkblue")),
      green = colorRamp2(c(0, 2, 5, 10), c("white", "lightgreen", "green", "darkgreen")),
      viridis = colorRamp2(c(0, 2, 5, 10), viridis::viridis(4))
    )
    col_fun <- color_schemes[[color_scale]]
  } else {
    count_colors <- list(
      red = colorRamp2(c(0, 30, 60, 100), c("white", "pink", "red", "darkred")),
      blue = colorRamp2(c(0, 30, 60, 100), c("white", "lightcyan", "dodgerblue", "navy")),
      green = colorRamp2(c(0, 30, 60, 100), c("white", "palegreen", "forestgreen", "darkgreen")),
      viridis = colorRamp2(c(0, 30, 60, 100), viridis::plasma(4))
    )
    col_fun <- count_colors[[color_scale]]
  }
  
  # Cell function
  if (size_encoding == "count") {
    cell_fun <- function(j, i, x, y, width, height, fill) {
      count_val <- count_mat[i, j]
      
      if (count_val > 0) {
        min_count <- min(count_mat[count_mat > 0])
        max_count <- max(count_mat)
        
        size_scale <- 0.7 + 0.3 * (count_val - min_count) / (max_count - min_count)
        
        grid.circle(
          x = x, y = y,
          r = min(unit.c(width, height)) * size_scale * 0.85,
          gp = gpar(fill = fill, col = NA)
        )
      }
    }
  } else {
    cell_fun <- function(j, i, x, y, width, height, fill) {
      p_val <- p_mat[i, j]
      
      if (p_val > 0) {
        min_p <- min(p_mat[p_mat > 0])
        max_p <- max(p_mat)
        
        size_scale <- 0.7 + 0.3 * (p_val - min_p) / (max_p - min_p)
        
        grid.circle(
          x = x, y = y,
          r = min(unit.c(width, height)) * size_scale * 0.85,
          gp = gpar(fill = fill, col = NA)
        )
      }
    }
  }
  
  # Choose matrix for color
  color_matrix <- if (size_encoding == "count") p_mat else count_mat
  color_name <- if (size_encoding == "count") "-log10(p.adjust)" else "Gene Count"
  
  # Create heatmap
  ht <- Heatmap(
    color_matrix,
    name = color_name,
    col = col_fun,
    rect_gp = gpar(fill = "transparent", col = "lightgray", lwd = 0.5),
    cell_fun = cell_fun,
    cluster_rows = row_clust,
    cluster_columns = col_clust,
    show_row_dend = cluster_rows,
    show_column_dend = cluster_cols,
    width = unit(ncol(p_mat) * 1.2, "cm"),
    height = unit(nrow(p_mat) * 0.6, "cm"),
    row_names_side = "left",
    row_dend_side = "right",
    row_names_gp = gpar(fontsize = 10),
    column_names_gp = gpar(fontsize = 10),
    column_names_rot = 45,
    column_title = title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_gap = unit(0.5, "mm"),
    column_gap = unit(0.5, "mm"),
    heatmap_legend_param = list(
      title = color_name,
      legend_direction = "vertical",
      title_gp = gpar(fontsize = 10),
      labels_gp = gpar(fontsize = 8)
    )
  )
  
  return(ht)
}