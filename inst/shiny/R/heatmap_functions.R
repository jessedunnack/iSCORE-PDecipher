# Heatmap functions for PerturbSeq enrichment analysis
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(pheatmap)
library(RColorBrewer)

# Prepare data for heatmap visualization
prepare_enrichment_heatmap <- function(data, 
                                     genes = NULL,
                                     enrichment_types = c("GO_BP", "KEGG", "Reactome"),
                                     direction = "ALL",
                                     max_terms = 50,
                                     min_frequency = 0.2,
                                     p_cutoff = 0.05) {
  
  # Filter data
  filtered_data <- data
  
  # Filter by genes if specified
  if (!is.null(genes) && length(genes) > 0) {
    filtered_data <- filtered_data %>%
      filter(gene %in% genes)
  }
  
  # Filter by enrichment types
  filtered_data <- filtered_data %>%
    filter(enrichment_type %in% enrichment_types)
  
  # Filter by direction
  if (direction != "ALL") {
    filtered_data <- filtered_data %>%
      filter(direction == !!direction)
  }
  
  # Filter by p-value
  filtered_data <- filtered_data %>%
    filter(p.adjust <= p_cutoff)
  
  # Create condition labels
  filtered_data <- filtered_data %>%
    mutate(condition = paste(gene, cluster, direction, sep = "_"))
  
  # Calculate term frequency across conditions
  term_freq <- filtered_data %>%
    group_by(Description) %>%
    summarise(
      n_conditions = n_distinct(condition),
      total_conditions = n_distinct(filtered_data$condition),
      frequency = n_conditions / total_conditions,
      mean_p = mean(p.adjust)
    ) %>%
    filter(frequency >= min_frequency) %>%
    arrange(desc(frequency), mean_p)
  
  # Select top terms
  top_terms <- head(term_freq$Description, max_terms)
  
  # Create matrix for heatmap
  heatmap_data <- filtered_data %>%
    filter(Description %in% top_terms) %>%
    select(Description, condition, p.adjust) %>%
    mutate(neg_log10_p = -log10(p.adjust)) %>%
    select(-p.adjust) %>%
    pivot_wider(names_from = condition, values_from = neg_log10_p, values_fill = 0)
  
  # Convert to matrix
  mat <- as.matrix(heatmap_data[, -1])
  rownames(mat) <- heatmap_data$Description
  
  return(list(
    matrix = mat,
    data = filtered_data,
    term_freq = term_freq
  ))
}

# Create interactive heatmap using plotly
create_interactive_heatmap <- function(heatmap_data, title = "Enrichment Heatmap") {
  library(plotly)
  
  mat <- heatmap_data$matrix
  
  # Create hover text
  hover_text <- matrix(
    paste("Term:", rownames(mat)[row(mat)], 
          "<br>Condition:", colnames(mat)[col(mat)],
          "<br>-log10(p):", round(mat, 3)),
    nrow = nrow(mat),
    ncol = ncol(mat)
  )
  
  p <- plot_ly(
    z = mat,
    x = colnames(mat),
    y = rownames(mat),
    type = "heatmap",
    colorscale = "Viridis",
    hovertext = hover_text,
    hoverinfo = "text"
  ) %>%
    layout(
      title = title,
      xaxis = list(title = "Condition", tickangle = -45),
      yaxis = list(title = "Enrichment Term"),
      height = 600 + nrow(mat) * 15,
      margin = list(l = 300, b = 150)
    )
  
  return(p)
}

# Create static heatmap using pheatmap
create_static_heatmap <- function(heatmap_data, 
                                 title = "Enrichment Heatmap",
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 show_rownames = TRUE,
                                 show_colnames = TRUE) {
  
  mat <- heatmap_data$matrix
  
  # Handle empty matrix
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    return(NULL)
  }
  
  # Create color palette
  colors <- colorRampPalette(c("white", "lightblue", "darkblue"))(100)
  
  # Create annotation for columns (conditions)
  col_annotation <- data.frame(
    Gene = sapply(strsplit(colnames(mat), "_"), `[`, 1),
    Cluster = sapply(strsplit(colnames(mat), "_"), `[`, 2),
    Direction = sapply(strsplit(colnames(mat), "_"), `[`, 3),
    row.names = colnames(mat)
  )
  
  # Create annotation colors
  n_genes <- length(unique(col_annotation$Gene))
  n_clusters <- length(unique(col_annotation$Cluster))
  
  ann_colors <- list(
    Gene = setNames(rainbow(n_genes), unique(col_annotation$Gene)),
    Cluster = setNames(brewer.pal(min(n_clusters, 12), "Set3"), unique(col_annotation$Cluster)[1:min(n_clusters, 12)]),
    Direction = c("UP" = "red", "DOWN" = "blue", "ALL" = "gray")
  )
  
  # Create heatmap
  p <- pheatmap(
    mat,
    main = title,
    color = colors,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    show_rownames = show_rownames,
    show_colnames = show_colnames,
    annotation_col = col_annotation,
    annotation_colors = ann_colors,
    fontsize_row = 8,
    fontsize_col = 8,
    angle_col = 45,
    border_color = NA
  )
  
  return(p)
}

# Create comparison heatmap between modalities
create_modality_comparison_heatmap <- function(data,
                                             genes = NULL,
                                             enrichment_type = "GO_BP",
                                             max_terms = 30,
                                             p_cutoff = 0.05) {
  
  # Filter data
  filtered_data <- data %>%
    filter(enrichment_type == !!enrichment_type,
           p.adjust <= p_cutoff)
  
  if (!is.null(genes)) {
    filtered_data <- filtered_data %>%
      filter(gene %in% genes)
  }
  
  # Get terms that appear in multiple modalities
  term_modality_count <- filtered_data %>%
    group_by(Description, modality) %>%
    summarise(
      n_occurrences = n(),
      mean_p = mean(p.adjust),
      .groups = 'drop'
    )
  
  # Find terms present in at least 2 modalities
  multi_modal_terms <- term_modality_count %>%
    group_by(Description) %>%
    summarise(n_modalities = n_distinct(modality)) %>%
    filter(n_modalities >= 2) %>%
    pull(Description)
  
  # Select top terms
  top_terms <- filtered_data %>%
    filter(Description %in% multi_modal_terms) %>%
    group_by(Description) %>%
    summarise(
      mean_p = mean(p.adjust),
      n_total = n()
    ) %>%
    arrange(mean_p) %>%
    head(max_terms) %>%
    pull(Description)
  
  # Create matrix for each modality
  modality_summary <- filtered_data %>%
    filter(Description %in% top_terms) %>%
    group_by(Description, modality) %>%
    summarise(
      mean_neg_log10_p = mean(-log10(p.adjust)),
      n_significant = n(),
      .groups = 'drop'
    ) %>%
    pivot_wider(
      names_from = modality,
      values_from = mean_neg_log10_p,
      values_fill = 0
    )
  
  # Convert to matrix
  mat <- as.matrix(modality_summary[, -1])
  rownames(mat) <- modality_summary$Description
  
  return(list(
    matrix = mat,
    data = filtered_data,
    summary = modality_summary
  ))
}

# Create gene-centric heatmap
create_gene_enrichment_heatmap <- function(data,
                                         enrichment_type = "GO_BP",
                                         cluster = NULL,
                                         direction = "ALL",
                                         max_terms = 40,
                                         p_cutoff = 0.05) {
  
  # Filter data
  filtered_data <- data %>%
    filter(enrichment_type == !!enrichment_type,
           p.adjust <= p_cutoff)
  
  if (!is.null(cluster)) {
    filtered_data <- filtered_data %>%
      filter(cluster == !!cluster)
  }
  
  if (direction != "ALL") {
    filtered_data <- filtered_data %>%
      filter(direction == !!direction)
  }
  
  # Find terms enriched in multiple genes
  term_gene_summary <- filtered_data %>%
    group_by(Description, gene) %>%
    summarise(
      min_p = min(p.adjust),
      .groups = 'drop'
    )
  
  # Select terms present in multiple genes
  multi_gene_terms <- term_gene_summary %>%
    group_by(Description) %>%
    summarise(n_genes = n_distinct(gene)) %>%
    filter(n_genes >= 3) %>%
    arrange(desc(n_genes)) %>%
    head(max_terms) %>%
    pull(Description)
  
  # Create matrix
  heatmap_data <- term_gene_summary %>%
    filter(Description %in% multi_gene_terms) %>%
    mutate(neg_log10_p = -log10(min_p)) %>%
    select(Description, gene, neg_log10_p) %>%
    pivot_wider(names_from = gene, values_from = neg_log10_p, values_fill = 0)
  
  # Convert to matrix
  mat <- as.matrix(heatmap_data[, -1])
  rownames(mat) <- heatmap_data$Description
  
  return(list(
    matrix = mat,
    data = filtered_data,
    summary = term_gene_summary
  ))
}