# Source Labeling Functions
# Provides functions to distinguish between iSCORE-PD mutations and PerturbSeq perturbations

#' Create source-aware gene labels
#' 
#' @param gene Gene name
#' @param method Analysis method (MAST, MixScale, MixScale_CRISPRa)
#' @param modality Optional modality for MixScale (CRISPRi, CRISPRa)
#' @return Gene label with source annotation
create_source_label <- function(gene, method, modality = NULL) {
  # Handle empty inputs
  if (length(gene) == 0 || length(method) == 0) {
    return(character(0))
  }
  
  # Handle vectors
  if (length(gene) > 1) {
    # Ensure all vectors are same length
    max_len <- max(length(gene), length(method), ifelse(is.null(modality), 1, length(modality)))
    
    # Recycle vectors to same length
    gene <- rep_len(gene, max_len)
    method <- rep_len(method, max_len)
    if (!is.null(modality)) {
      modality <- rep_len(modality, max_len)
    } else {
      modality <- rep(NA_character_, max_len)
    }
    
    return(mapply(create_source_label, gene, method, modality, USE.NAMES = FALSE))
  }
  
  # Single value logic
  if (is.na(gene) || is.null(gene) || gene == "") return(NA_character_)
  if (is.na(method) || is.null(method) || method == "") return(as.character(gene))
  
  # Safe case_when replacement that handles NAs
  if (method == "MAST") {
    label <- paste0(gene, " (iSCORE-PD)")
  } else if (method == "MixScale" && !is.null(modality) && !is.na(modality) && modality == "CRISPRi") {
    label <- paste0(gene, " (CRISPRi)")
  } else if (method == "MixScale" && !is.null(modality) && !is.na(modality) && modality == "CRISPRa") {
    label <- paste0(gene, " (CRISPRa)")
  } else if (method == "MixScale_CRISPRa") {
    label <- paste0(gene, " (CRISPRa)")
  } else if (method == "MixScale") {
    label <- paste0(gene, " (CRISPRi)")  # Default MixScale to CRISPRi
  } else {
    label <- as.character(gene)
  }
  
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
  
  # Input validation
  if (is.null(data) || nrow(data) == 0) {
    warning("Data is NULL or empty")
    return(data)
  }
  
  if (!gene_col %in% names(data)) {
    warning(paste("Column", gene_col, "not found in data"))
    return(data)
  }
  
  if (!method_col %in% names(data)) {
    warning(paste("Column", method_col, "not found in data"))
    return(data)
  }
  
  # Get modality if available, otherwise NULL
  modality <- if (!is.null(modality_col) && modality_col %in% names(data)) {
    data[[modality_col]]
  } else {
    NULL
  }
  
  # Validate gene and method columns have values
  if (all(is.na(data[[gene_col]])) || all(data[[gene_col]] == "")) {
    warning("Gene column contains only NA or empty values")
    return(data)
  }
  
  if (all(is.na(data[[method_col]])) || all(data[[method_col]] == "")) {
    warning("Method column contains only NA or empty values")
    return(data)
  }
  
  # Create source labels safely
  tryCatch({
    data$source_label <- create_source_label(
      gene = data[[gene_col]],
      method = data[[method_col]],
      modality = modality
    )
  }, error = function(e) {
    warning(paste("Error creating source labels:", e$message))
    # Return data without source labels
    return(data)
  })
  
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