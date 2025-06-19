# DE Data Processing Functions
# Enhanced functions for processing differential expression data

#' Process DE data with metadata integration
#' 
#' @param de_results DE results from full_DE_results.rds
#' @param metadata Cell metadata from final_dataset_metadata.rds
#' @return Processed data frame with DE statistics and metadata
process_de_data_with_metadata <- function(de_results, metadata = NULL, max_genes_per_condition = 1000, cluster_filter = NULL) {
  
  processed_data <- data.frame()
  
  # Process MAST results with filtering
  if ("iSCORE_PD_MAST" %in% names(de_results)) {
    mast_data <- de_results$iSCORE_PD_MAST
    
    for (gene in names(mast_data)) {
      for (cluster in names(mast_data[[gene]])) {
        # Skip clusters not in filter (if specified)
        if (!is.null(cluster_filter) && !cluster %in% cluster_filter) {
          next
        }
        
        if (!is.null(mast_data[[gene]][[cluster]]$results)) {
          results <- mast_data[[gene]][[cluster]]$results
          
          if (nrow(results) > 0) {
            # Filter to significant genes only to reduce data size
            sig_results <- results[!is.na(results$p_val_adj) & results$p_val_adj < 0.05, ]
            
            # Further limit to top genes by significance if still too many
            if (nrow(sig_results) > max_genes_per_condition) {
              sig_results <- sig_results[order(sig_results$p_val_adj), ]
              sig_results <- sig_results[1:max_genes_per_condition, ]
            }
            
            if (nrow(sig_results) > 0) {
              cluster_data <- data.frame(
                gene = gene,
                cluster = cluster,
                gene_name = rownames(sig_results),
                log2FC = sig_results$avg_log2FC,
                pvalue = sig_results$p_val_adj,
                method = "MAST",
                source = "iSCORE-PD",
                experiment = "default",
                stringsAsFactors = FALSE
              )
              processed_data <- rbind(processed_data, cluster_data)
            }
          }
        }
      }
    }
  }
  
  # Process CRISPRi MixScale results with filtering
  if ("CRISPRi_Mixscale" %in% names(de_results)) {
    crispi_data <- de_results$CRISPRi_Mixscale
    
    for (gene in names(crispi_data)) {
      for (cluster in names(crispi_data[[gene]])) {
        # Skip clusters not in filter (if specified)
        if (!is.null(cluster_filter) && !cluster %in% cluster_filter) {
          next
        }
        
        if (!is.null(crispi_data[[gene]][[cluster]]$results)) {
          results <- crispi_data[[gene]][[cluster]]$results
          
          # Find experiment-specific columns
          log2fc_cols <- grep("^log2FC_", names(results), value = TRUE)
          
          for (log2fc_col in log2fc_cols) {
            exp <- gsub("^log2FC_", "", log2fc_col)
            pval_col <- paste0("p_cell_type", exp, ":weight")
            
            if (pval_col %in% names(results) && nrow(results) > 0) {
              # Filter to significant genes only
              valid_rows <- !is.na(results[[pval_col]]) & results[[pval_col]] < 0.05
              sig_results <- results[valid_rows, ]
              
              # Limit to top genes if still too many
              if (nrow(sig_results) > max_genes_per_condition) {
                sig_results <- sig_results[order(sig_results[[pval_col]]), ]
                sig_results <- sig_results[1:max_genes_per_condition, ]
              }
              
              if (nrow(sig_results) > 0) {
                cluster_data <- data.frame(
                  gene = gene,
                  cluster = cluster,
                  gene_name = rownames(sig_results),
                  log2FC = sig_results[[log2fc_col]],
                  pvalue = sig_results[[pval_col]],
                  method = "MixScale",
                  source = "CRISPRi", 
                  experiment = exp,
                  stringsAsFactors = FALSE
                )
                processed_data <- rbind(processed_data, cluster_data)
              }
            }
          }
        }
      }
    }
  }
  
  # Process CRISPRa MixScale results
  if ("CRISPRa_Mixscale" %in% names(de_results)) {
    crispa_data <- de_results$CRISPRa_Mixscale
    
    for (gene in names(crispa_data)) {
      for (cluster in names(crispa_data[[gene]])) {
        # Skip clusters not in filter (if specified)
        if (!is.null(cluster_filter) && !cluster %in% cluster_filter) {
          next
        }
        
        if (!is.null(crispa_data[[gene]][[cluster]]$results)) {
          results <- crispa_data[[gene]][[cluster]]$results
          
          log2fc_cols <- grep("^log2FC_", names(results), value = TRUE)
          
          for (log2fc_col in log2fc_cols) {
            exp <- gsub("^log2FC_", "", log2fc_col)
            pval_col <- paste0("p_cell_type", exp, ":weight")
            
            if (pval_col %in% names(results) && nrow(results) > 0) {
              cluster_data <- data.frame(
                gene = gene,
                cluster = cluster,
                gene_name = rownames(results),
                log2FC = results[[log2fc_col]],
                pvalue = results[[pval_col]],
                method = "MixScale",
                source = "CRISPRa",
                experiment = exp,
                stringsAsFactors = FALSE
              )
              processed_data <- rbind(processed_data, cluster_data)
            }
          }
        }
      }
    }
  }
  
  # Clean up data
  processed_data <- processed_data[!is.na(processed_data$log2FC) & !is.na(processed_data$pvalue), ]
  
  # Add source labels
  processed_data$source_label <- paste0(processed_data$gene, " (", processed_data$source, ")")
  
  return(processed_data)
}

#' Prepare multi-condition DE comparison
#' 
#' @param de_data Processed DE data
#' @param conditions List of conditions to compare
#' @param cutoffs List with log2FC and pvalue thresholds
#' @return Filtered data ready for comparison
prepare_multi_condition_de <- function(de_data, conditions, cutoffs = list(log2FC = 0.25, pvalue = 0.05)) {
  
  # Apply significance cutoffs
  filtered_data <- de_data %>%
    dplyr::filter(
      abs(log2FC) >= cutoffs$log2FC,
      pvalue <= cutoffs$pvalue
    )
  
  # Filter by conditions if specified
  if (!is.null(conditions) && length(conditions) > 0) {
    # Create condition identifier
    filtered_data$condition_id <- paste(filtered_data$gene, filtered_data$cluster, filtered_data$source, sep = "_")
    
    # Filter to selected conditions
    condition_ids <- paste(conditions$gene, conditions$cluster, conditions$source, sep = "_")
    filtered_data <- filtered_data[filtered_data$condition_id %in% condition_ids, ]
  }
  
  return(filtered_data)
}

#' Sample DE data for preview
#' 
#' @param de_data DE data frame
#' @param max_genes Maximum number of genes to include
#' @param method Sampling method ("top_significant", "random", "balanced")
#' @return Sampled data frame
sample_de_for_preview <- function(de_data, max_genes = 500, method = "top_significant") {
  
  if (nrow(de_data) <= max_genes) {
    return(de_data)
  }
  
  if (method == "top_significant") {
    # Sample top genes by significance
    sampled <- de_data %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice_head(n = max_genes)
      
  } else if (method == "balanced") {
    # Sample equally from each source
    sources <- unique(de_data$source)
    genes_per_source <- floor(max_genes / length(sources))
    
    sampled <- de_data %>%
      dplyr::group_by(source) %>%
      dplyr::arrange(pvalue) %>%
      dplyr::slice_head(n = genes_per_source) %>%
      dplyr::ungroup()
      
  } else {
    # Random sampling
    sampled <- de_data %>%
      dplyr::sample_n(min(max_genes, nrow(de_data)))
  }
  
  return(sampled)
}

#' Calculate DE gene overlaps
#' 
#' @param de_data DE data frame
#' @param comparison_groups List of conditions to compare
#' @param cutoffs Significance cutoffs
#' @return List with overlap information
calculate_de_overlaps <- function(de_data, comparison_groups, cutoffs = list(log2FC = 0.25, pvalue = 0.05)) {
  
  # Get significant genes for each group
  gene_lists <- list()
  
  for (group_name in names(comparison_groups)) {
    group_condition <- comparison_groups[[group_name]]
    
    # Filter data for this group
    group_data <- de_data %>%
      dplyr::filter(
        gene == group_condition$gene,
        cluster == group_condition$cluster,
        source == group_condition$source,
        abs(log2FC) >= cutoffs$log2FC,
        pvalue <= cutoffs$pvalue
      )
    
    gene_lists[[group_name]] <- unique(group_data$gene_name)
  }
  
  # Calculate intersections
  if (length(gene_lists) >= 2) {
    # Pairwise intersections
    pairs <- combn(names(gene_lists), 2, simplify = FALSE)
    intersections <- list()
    
    for (pair in pairs) {
      pair_name <- paste(pair, collapse = " ∩ ")
      intersections[[pair_name]] <- intersect(gene_lists[[pair[1]]], gene_lists[[pair[2]]])
    }
    
    # Full intersection (if more than 2 groups)
    if (length(gene_lists) > 2) {
      intersections[["All"]] <- Reduce(intersect, gene_lists)
    }
    
    # Union
    intersections[["Union"]] <- Reduce(union, gene_lists)
    
    return(list(
      gene_lists = gene_lists,
      intersections = intersections,
      summary = data.frame(
        group = names(gene_lists),
        gene_count = sapply(gene_lists, length),
        stringsAsFactors = FALSE
      )
    ))
  } else {
    return(list(
      gene_lists = gene_lists,
      intersections = list(),
      summary = data.frame(
        group = names(gene_lists),
        gene_count = sapply(gene_lists, length),
        stringsAsFactors = FALSE
      )
    ))
  }
}

#' Calculate DE correlations between conditions
#' 
#' @param de_data DE data frame
#' @param comparison_groups List of conditions
#' @return Correlation matrix
calculate_de_correlations <- function(de_data, comparison_groups) {
  
  # Pivot data to wide format for correlation
  correlation_data <- de_data %>%
    dplyr::select(gene_name, gene, cluster, source, log2FC) %>%
    dplyr::mutate(condition = paste(gene, cluster, source, sep = "_")) %>%
    dplyr::select(gene_name, condition, log2FC) %>%
    tidyr::pivot_wider(names_from = condition, values_from = log2FC, values_fill = 0)
  
  # Calculate correlation matrix
  correlation_matrix <- correlation_data %>%
    dplyr::select(-gene_name) %>%
    cor(use = "complete.obs")
  
  return(correlation_matrix)
}