# Source Labeling Functions
# Provides functions to distinguish between iSCORE-PD mutations and PerturbSeq perturbations

#' Create source-aware gene labels
#' 
#' @param gene Gene name
#' @param method Analysis method (MAST, MixScale, MixScale_CRISPRa)
#' @param modality Optional modality for MixScale (CRISPRi, CRISPRa)
#' @return Gene label with source annotation
create_source_label <- function(gene, method, modality = NULL) {
  # Handle vectors
  if (length(gene) > 1) {
    return(mapply(create_source_label, gene, method, modality, USE.NAMES = FALSE))
  }
  
  # Single value logic
  if (is.na(gene) || is.null(gene)) return(NA_character_)
  
  label <- dplyr::case_when(
    method == "MAST" ~ paste0(gene, " (iSCORE-PD)"),
    method == "MixScale" & !is.null(modality) & modality == "CRISPRi" ~ paste0(gene, " (CRISPRi)"),
    method == "MixScale" & !is.null(modality) & modality == "CRISPRa" ~ paste0(gene, " (CRISPRa)"),
    method == "MixScale_CRISPRa" ~ paste0(gene, " (CRISPRa)"),
    method == "MixScale" ~ paste0(gene, " (CRISPRi)"), # Default MixScale to CRISPRi
    TRUE ~ as.character(gene)
  )
  
  return(label)
}

#' Add source labels to data frame
#' 
#' @param data Data frame with gene/mutation_perturbation and method columns
#' @param gene_col Name of gene column (default: "mutation_perturbation")
#' @param method_col Name of method column (default: "method")
#' @param modality_col Name of modality column (default: "modality")
#' @return Data frame with added source_label column
add_source_labels <- function(data, 
                            gene_col = "mutation_perturbation", 
                            method_col = "method",
                            modality_col = "modality") {
  
  if (!gene_col %in% names(data)) {
    warning(paste("Column", gene_col, "not found in data"))
    return(data)
  }
  
  if (!method_col %in% names(data)) {
    warning(paste("Column", method_col, "not found in data"))
    return(data)
  }
  
  # Get modality if available
  modality <- if (modality_col %in% names(data)) data[[modality_col]] else NULL
  
  # Create source labels
  data$source_label <- create_source_label(
    gene = data[[gene_col]],
    method = data[[method_col]],
    modality = modality
  )
  
  return(data)
}

#' Extract source from labeled gene
#' 
#' @param labeled_gene Gene with source label e.g., "ATP13A2 (iSCORE-PD)"
#' @return List with gene and source components
parse_source_label <- function(labeled_gene) {
  if (is.na(labeled_gene) || is.null(labeled_gene)) {
    return(list(gene = NA_character_, source = NA_character_))
  }
  
  # Extract gene and source using regex
  match <- regexpr("^(.+)\\s+\\((.+)\\)$", labeled_gene, perl = TRUE)
  
  if (match == -1) {
    # No source label found
    return(list(gene = labeled_gene, source = NA_character_))
  }
  
  # Extract captured groups
  gene <- sub("^(.+)\\s+\\((.+)\\)$", "\\1", labeled_gene)
  source <- sub("^(.+)\\s+\\((.+)\\)$", "\\2", labeled_gene)
  
  return(list(gene = gene, source = source))
}

#' Get unique genes with sources from data
#' 
#' @param data Data frame with source labels
#' @param source_label_col Column containing source labels
#' @return Data frame with unique gene-source combinations
get_unique_gene_sources <- function(data, source_label_col = "source_label") {
  if (!source_label_col %in% names(data)) {
    warning(paste("Column", source_label_col, "not found"))
    return(data.frame(gene = character(), source = character()))
  }
  
  unique_labels <- unique(data[[source_label_col]])
  unique_labels <- unique_labels[!is.na(unique_labels)]
  
  # Parse each label
  parsed <- lapply(unique_labels, parse_source_label)
  
  result <- data.frame(
    gene = sapply(parsed, `[[`, "gene"),
    source = sapply(parsed, `[[`, "source"),
    source_label = unique_labels,
    stringsAsFactors = FALSE
  )
  
  return(result[order(result$gene, result$source), ])
}

#' Group genes by source
#' 
#' @param data Data frame with source labels
#' @return List of genes grouped by source
group_genes_by_source <- function(data) {
  unique_sources <- get_unique_gene_sources(data)
  
  sources <- unique(unique_sources$source[!is.na(unique_sources$source)])
  
  result <- list()
  for (src in sources) {
    result[[src]] <- unique_sources$gene[unique_sources$source == src]
  }
  
  return(result)
}