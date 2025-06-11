# Term extraction and visualization preparation functions
# These functions extract terms from enrichment results and prepare them for Shiny visualization

#' Handle different enrichment result formats
#'
#' @param result_obj The enrichment result object (can be S4, list, or data.frame)
#' @param min_padj Adjusted p-value threshold for filtering
#' @return A data frame of enrichment results with standardized columns
handle_enrichment_result <- function(result_obj, min_padj = 0.05) {
  # Initialize empty data frame for results
  empty_df <- data.frame(
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    Count = integer(0),
    stringsAsFactors = FALSE
  )
  
  # Case 1: If result is NULL or NA
  if (is.null(result_obj) || length(result_obj) == 0 || all(is.na(result_obj))) {
    return(empty_df)
  }
  
  # Check for GSEA format
  is_gsea_format <- FALSE
  if (is.data.frame(result_obj) && 
      all(c("pathway", "pval", "padj", "NES") %in% colnames(result_obj))) {
    is_gsea_format <- TRUE
  } else if (methods::isS4(result_obj) && "result" %in% methods::slotNames(result_obj)) {
    result_slot <- tryCatch({
      methods::slot(result_obj, "result")
    }, error = function(e) NULL)
    
    if (is.data.frame(result_slot) && 
        all(c("pathway", "pval", "padj", "NES") %in% colnames(result_slot))) {
      is_gsea_format <- TRUE
    }
  }
  
  if (is_gsea_format) {
    return(extract_gsea_terms(result_obj, min_padj))
  }
  
  # Check for STRING results
  is_string <- FALSE
  if (is.list(result_obj) && "enrichment" %in% names(result_obj) && 
      is.data.frame(result_obj$enrichment)) {
    is_string <- TRUE
  } else if (is.data.frame(result_obj) && 
             any(c("term", "description", "p_value", "fdr", "category") %in% colnames(result_obj))) {
    is_string <- TRUE
  }
  
  if (is_string) {
    return(extract_string_terms(result_obj, min_padj))
  }
  
  # Handle standard enrichResult objects (GO, KEGG, Reactome, WikiPathways)
  if (methods::isS4(result_obj) && methods::is(result_obj, "enrichResult")) {
    result_df <- methods::slot(result_obj, "result")
    
    # Filter by adjusted p-value
    if ("p.adjust" %in% colnames(result_df)) {
      result_df <- result_df[result_df$p.adjust <= min_padj, ]
    }
    
    return(result_df)
  }
  
  # Handle data.frame directly
  if (is.data.frame(result_obj)) {
    # Filter by adjusted p-value if available
    if ("p.adjust" %in% colnames(result_obj)) {
      result_df <- result_obj[result_obj$p.adjust <= min_padj, ]
      return(result_df)
    }
  }
  
  # Return empty data frame if format not recognized
  return(empty_df)
}

#' Extract terms from STRING results
#'
#' @param result_obj STRING result object
#' @param min_padj Maximum adjusted p-value threshold
#' @return Data frame with standardized columns
extract_string_terms <- function(result_obj, min_padj = 0.05) {
  # Extract STRING data
  if (is.list(result_obj) && "enrichment" %in% names(result_obj) && 
      is.data.frame(result_obj$enrichment)) {
    string_df <- result_obj$enrichment
  } else if (is.data.frame(result_obj)) {
    string_df <- result_obj
  } else {
    return(data.frame())
  }
  
  # Standardize column names
  if ("term" %in% colnames(string_df) && !"ID" %in% colnames(string_df)) {
    string_df$ID <- string_df$term
  }
  if ("description" %in% colnames(string_df) && !"Description" %in% colnames(string_df)) {
    string_df$Description <- string_df$description
  }
  if ("p_value" %in% colnames(string_df) && !"pvalue" %in% colnames(string_df)) {
    string_df$pvalue <- string_df$p_value
  }
  if ("fdr" %in% colnames(string_df) && !"p.adjust" %in% colnames(string_df)) {
    string_df$p.adjust <- string_df$fdr
  }
  if ("number_of_genes" %in% colnames(string_df) && !"Count" %in% colnames(string_df)) {
    string_df$Count <- string_df$number_of_genes
  }
  
  # Apply p-value cutoff
  if (min_padj < 1 && "p.adjust" %in% colnames(string_df)) {
    string_df <- string_df[string_df$p.adjust <= min_padj, ]
  }
  
  return(string_df)
}

#' Extract terms from GSEA results
#'
#' @param result_obj GSEA result object
#' @param min_padj Maximum adjusted p-value threshold
#' @return Data frame with standardized columns
extract_gsea_terms <- function(result_obj, min_padj = 0.05) {
  # Handle S4 objects with result slot
  if (methods::isS4(result_obj) && "result" %in% methods::slotNames(result_obj)) {
    result_df <- methods::slot(result_obj, "result")
  } else if (is.data.frame(result_obj)) {
    result_df <- result_obj
  } else {
    return(data.frame())
  }
  
  # Check for GSEA columns
  if (!all(c("pathway", "pval", "padj") %in% colnames(result_df))) {
    return(data.frame())
  }
  
  # Convert to standard format
  standardized <- data.frame(
    ID = result_df$pathway,
    Description = result_df$pathway,
    pvalue = result_df$pval,
    p.adjust = result_df$padj,
    stringsAsFactors = FALSE
  )
  
  # Add additional GSEA-specific columns if they exist
  if ("size" %in% colnames(result_df)) {
    standardized$Count <- result_df$size
  }
  if ("ES" %in% colnames(result_df)) {
    standardized$ES <- result_df$ES
  }
  if ("NES" %in% colnames(result_df)) {
    standardized$NES <- result_df$NES
  }
  
  # Filter by adjusted p-value
  standardized <- standardized[standardized$p.adjust <= min_padj, ]
  
  return(standardized)
}

#' Load and extract terms from an enrichment result file
#'
#' @param file_path Path to the RDS file
#' @param min_padj Maximum adjusted p-value threshold
#' @return Data frame with enrichment terms
load_and_extract_terms <- function(file_path, min_padj = 0.05) {
  if (!file.exists(file_path)) {
    warning("File not found: ", file_path)
    return(NULL)
  }
  
  tryCatch({
    # Load the RDS file
    result_obj <- readRDS(file_path)
    
    # Extract terms using appropriate handler
    terms_df <- handle_enrichment_result(result_obj, min_padj)
    
    # Add source file information
    if (!is.null(terms_df) && nrow(terms_df) > 0) {
      terms_df$source_file <- file_path
    }
    
    return(terms_df)
    
  }, error = function(e) {
    warning("Error loading file ", file_path, ": ", e$message)
    return(NULL)
  })
}

#' Get enrichment file path
#'
#' @param base_dir Base directory for enrichment results
#' @param analysis_type "MAST" or "MixScale" 
#' @param gene Gene/mutation name
#' @param cluster Cluster ID
#' @param experiment Experiment ID (default for MAST, specific for MixScale)
#' @param enrichment_type Type of enrichment (GO_BP, KEGG, etc.)
#' @param direction UP, DOWN, or ALL
#' @return Full path to the enrichment result file
get_enrichment_file_path <- function(base_dir, analysis_type, gene, cluster, 
                                   experiment, enrichment_type, direction) {
  
  # Special handling for GSEA
  if (enrichment_type == "GSEA") {
    # GSEA has subdirectories for different databases
    gsea_path <- file.path(base_dir, analysis_type, gene, cluster, experiment, "GSEA")
    
    if (dir.exists(gsea_path)) {
      # Look for gsea_results.rds in subdirectories
      gsea_files <- list.files(gsea_path, pattern = "gsea_results\\.rds$", 
                              recursive = TRUE, full.names = TRUE)
      if (length(gsea_files) > 0) {
        return(gsea_files[1])  # Return first found
      }
    }
    return(NULL)
  }
  
  # Standard enrichment types
  file_name <- paste0(enrichment_type, "_", direction, ".rds")
  file_path <- file.path(base_dir, analysis_type, gene, cluster, experiment, 
                        enrichment_type, file_name)
  
  if (file.exists(file_path)) {
    return(file_path)
  } else {
    return(NULL)
  }
}

#' Compare terms between conditions
#'
#' @param terms_list List of term data frames to compare
#' @param comparison_type "intersection" or "union"
#' @return Data frame with comparison results
compare_terms <- function(terms_list, comparison_type = "intersection") {
  if (length(terms_list) < 2) {
    stop("Need at least 2 term sets to compare")
  }
  
  # Extract term IDs from each set
  term_ids_list <- lapply(terms_list, function(df) {
    if (!is.null(df) && "ID" %in% colnames(df)) {
      unique(df$ID)
    } else {
      character(0)
    }
  })
  
  # Remove empty sets
  term_ids_list <- term_ids_list[sapply(term_ids_list, length) > 0]
  
  if (length(term_ids_list) == 0) {
    return(data.frame())
  }
  
  # Calculate intersection or union
  if (comparison_type == "intersection") {
    common_terms <- Reduce(intersect, term_ids_list)
  } else {
    common_terms <- Reduce(union, term_ids_list)
  }
  
  # Create result data frame with terms from all sets
  result_list <- list()
  
  for (i in seq_along(terms_list)) {
    df <- terms_list[[i]]
    if (!is.null(df) && nrow(df) > 0) {
      # Filter to common terms
      if (comparison_type == "intersection") {
        df_filtered <- df[df$ID %in% common_terms, ]
      } else {
        df_filtered <- df
      }
      
      if (nrow(df_filtered) > 0) {
        # Add source information
        df_filtered$source_index <- i
        result_list[[i]] <- df_filtered
      }
    }
  }
  
  # Combine all results
  if (length(result_list) > 0) {
    combined <- do.call(rbind, result_list)
    
    # Sort by p-value
    if ("p.adjust" %in% colnames(combined)) {
      combined <- combined[order(combined$p.adjust), ]
    }
    
    return(combined)
  } else {
    return(data.frame())
  }
}

#' Find frequently occurring terms across multiple analyses
#'
#' @param terms_list List of term data frames
#' @param min_frequency Minimum frequency (proportion) to be considered frequent
#' @return Data frame with frequent terms and their statistics
find_frequent_terms <- function(terms_list, min_frequency = 0.5) {
  # Remove NULL entries
  terms_list <- terms_list[!sapply(terms_list, is.null)]
  
  if (length(terms_list) == 0) {
    return(data.frame())
  }
  
  # Count term occurrences
  all_terms <- list()
  
  for (i in seq_along(terms_list)) {
    df <- terms_list[[i]]
    if (!is.null(df) && "ID" %in% colnames(df) && nrow(df) > 0) {
      # Get unique terms from this analysis
      terms <- unique(df$ID)
      all_terms[[i]] <- terms
    }
  }
  
  # Count frequencies
  all_term_ids <- unlist(all_terms)
  term_counts <- table(all_term_ids)
  
  # Calculate frequency
  n_analyses <- length(all_terms)
  term_freq <- term_counts / n_analyses
  
  # Filter by minimum frequency
  frequent_terms <- names(term_freq)[term_freq >= min_frequency]
  
  if (length(frequent_terms) == 0) {
    return(data.frame())
  }
  
  # Create summary data frame
  result <- data.frame(
    ID = frequent_terms,
    frequency = as.numeric(term_freq[frequent_terms]),
    count = as.numeric(term_counts[frequent_terms]),
    total_analyses = n_analyses,
    stringsAsFactors = FALSE
  )
  
  # Add descriptions from the first occurrence
  descriptions <- character(length(frequent_terms))
  
  for (i in seq_along(frequent_terms)) {
    term_id <- frequent_terms[i]
    
    # Find first occurrence with description
    for (df in terms_list) {
      if (!is.null(df) && "ID" %in% colnames(df) && "Description" %in% colnames(df)) {
        idx <- which(df$ID == term_id)[1]
        if (!is.na(idx)) {
          descriptions[i] <- df$Description[idx]
          break
        }
      }
    }
  }
  
  result$Description <- descriptions
  
  # Sort by frequency
  result <- result[order(result$frequency, decreasing = TRUE), ]
  
  return(result)
}

#' Prepare terms for heatmap visualization
#'
#' @param consolidated_data Data frame with all enrichment results
#' @param enrichment_types Vector of enrichment types to include
#' @param max_terms Maximum number of terms to include
#' @return Matrix suitable for heatmap visualization
prepare_heatmap_data <- function(consolidated_data, 
                               enrichment_types = c("GO_BP", "KEGG", "Reactome"),
                               max_terms = 50) {
  
  # Filter by enrichment types
  data_filtered <- consolidated_data[consolidated_data$enrichment_type %in% enrichment_types, ]
  
  if (nrow(data_filtered) == 0) {
    return(NULL)
  }
  
  # Get top terms by frequency across conditions
  term_counts <- table(data_filtered$ID)
  top_terms <- names(sort(term_counts, decreasing = TRUE))[1:min(max_terms, length(term_counts))]
  
  # Filter to top terms
  data_top <- data_filtered[data_filtered$ID %in% top_terms, ]
  
  # Create condition labels
  data_top$condition <- paste(data_top$mutation_perturbation, 
                             data_top$cluster, 
                             data_top$direction, 
                             sep = "_")
  
  # Convert to wide format for heatmap
  # Using -log10(p.adjust) as the value
  data_top$neg_log_padj <- -log10(data_top$p.adjust + 1e-100)  # Add small value to avoid Inf
  
  # Create matrix
  heatmap_matrix <- reshape2::dcast(data_top, 
                                   ID ~ condition, 
                                   value.var = "neg_log_padj",
                                   fun.aggregate = max,
                                   fill = 0)
  
  # Convert to matrix
  rownames(heatmap_matrix) <- heatmap_matrix$ID
  heatmap_matrix$ID <- NULL
  heatmap_matrix <- as.matrix(heatmap_matrix)
  
  return(heatmap_matrix)
}

# Example usage:
# # Load consolidated data
# consolidated_data <- readRDS("all_enrichment_padj005_complete_with_direction.rds")
# 
# # Extract terms for specific analysis
# terms <- load_and_extract_terms("enrichment_results/MAST/LRRK2/cluster_0/default/GO_BP/GO_BP_UP.rds")
# 
# # Compare terms between conditions
# terms1 <- consolidated_data[consolidated_data$mutation_perturbation == "LRRK2" & 
#                            consolidated_data$direction == "UP", ]
# terms2 <- consolidated_data[consolidated_data$mutation_perturbation == "PINK1" & 
#                            consolidated_data$direction == "UP", ]
# comparison <- compare_terms(list(LRRK2 = terms1, PINK1 = terms2), "intersection")
# 
# # Find frequent terms
# frequent <- find_frequent_terms(list(terms1, terms2), min_frequency = 0.5)
# 
# # Prepare heatmap data
# heatmap_data <- prepare_heatmap_data(consolidated_data)