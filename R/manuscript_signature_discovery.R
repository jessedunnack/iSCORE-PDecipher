#' Manuscript-Focused Signature Discovery Pipeline
#'
#' This module implements rapid signature discovery prioritized for manuscript development,
#' focusing on the strongest shared signatures between MAST and CRISPRi methods.

# Required packages (loaded via namespace imports)
# Dependencies: dplyr

#' Discover top signatures for manuscript prioritization
#'
#' @param enrichment_data Consolidated enrichment data
#' @param de_data Original differential expression data (optional, for proper background genes)
#' @param top_n Integer, number of top signatures to return
#' @param min_cluster_breadth Integer, minimum clusters showing signature (for pan-cluster analysis)
#' @param combine_variants Logical, whether to combine gene variants
#' @param progress_callback Function for progress updates
#' @param use_enhanced_analysis Logical, whether to use enhanced direction-aware analysis (default TRUE)
#' @param direction Character, which direction to analyze ("ALL", "UP", "DOWN"). Default "ALL" prevents directionality inflation.
#' @return List with ranked signature discoveries
#' @export
discover_top_signatures <- function(enrichment_data, de_data = NULL, top_n = 10, min_cluster_breadth = 8,
                                   combine_variants = TRUE, progress_callback = NULL, 
                                   use_enhanced_analysis = TRUE, direction = "ALL") {
  
  if (!is.null(progress_callback)) {
    progress_callback("Initializing signature discovery...", value = 0.1)
  }
  
  # Get comparable gene pairs (MAST vs CRISPRi only)
  gene_pairs <- get_comparable_gene_pairs(
    combine_snca_variants = combine_variants,
    combine_vps13c_variants = combine_variants,
    include_mast_only = FALSE  # Exclude GBA for cross-method comparison
  )
  
  if (!is.null(progress_callback)) {
    progress_callback("Analyzing gene pair signatures...", value = 0.2)
  }
  
  # Analyze each gene pair
  all_signatures <- list()
  
  # Calculate total work for progress tracking
  total_work <- nrow(gene_pairs)
  
  for (i in seq_len(nrow(gene_pairs))) {
    gene_pair <- gene_pairs[i, ]
    
    if (!is.null(progress_callback)) {
      progress_callback(
        paste("Analyzing", gene_pair$mast_gene, "vs", gene_pair$crispri_gene), 
        value = 0.2 + (i / nrow(gene_pairs)) * 0.5,
        detail = paste("Gene pair", i, "of", total_work)
      )
    }
    
    # Analyze this gene pair across all clusters
    pair_start_time <- Sys.time()
    
    pair_analysis <- analyze_gene_pair_signatures(
      list(mast_gene = gene_pair$mast_gene, crispri_gene = gene_pair$crispri_gene),
      enrichment_data,
      de_data = de_data,
      include_pathways = TRUE,
      use_enhanced_analysis = use_enhanced_analysis,
      direction = direction,
      progress_callback = function(msg) {
        if (!is.null(progress_callback)) {
          pair_elapsed <- as.numeric(difftime(Sys.time(), pair_start_time, units = "secs"))
          progress_callback(
            paste("Analyzing", gene_pair$mast_gene, "vs", gene_pair$crispri_gene, "-", msg),
            value = 0.2 + ((i - 1 + 0.5) / nrow(gene_pairs)) * 0.5,
            detail = paste("Gene pair", i, "of", total_work, sprintf("(%.0fs)", pair_elapsed))
          )
        }
      }
    )
    
    if (!"error" %in% names(pair_analysis)) {
      all_signatures[[paste(gene_pair$mast_gene, gene_pair$crispri_gene, sep = "_vs_")]] <- pair_analysis
    }
    
    if (!is.null(progress_callback)) {
      pair_elapsed <- as.numeric(difftime(Sys.time(), pair_start_time, units = "secs"))
      progress_callback(
        paste("Completed", gene_pair$mast_gene, "vs", gene_pair$crispri_gene), 
        value = 0.2 + (i / nrow(gene_pairs)) * 0.5,
        detail = paste("Gene pair", i, "of", total_work, sprintf("completed in %.0fs", pair_elapsed))
      )
    }
  }
  
  if (!is.null(progress_callback)) {
    progress_callback("Computing signature strength scores...", value = 0.7)
  }
  
  # Calculate signature strength for each gene pair and cluster
  signature_rankings <- compute_signature_rankings(all_signatures)
  
  if (!is.null(progress_callback)) {
    progress_callback("Identifying pan-cluster signatures...", value = 0.8)
  }
  
  # Identify pan-cluster signatures (appear across many clusters)
  pan_cluster_signatures <- identify_pan_cluster_signatures(
    signature_rankings, min_cluster_breadth = min_cluster_breadth
  )
  
  if (!is.null(progress_callback)) {
    progress_callback("Identifying cluster-specific signatures...", value = 0.9)
  }
  
  # Identify cluster-specific signatures
  cluster_specific_signatures <- identify_cluster_specific_signatures(signature_rankings)
  
  if (!is.null(progress_callback)) {
    progress_callback("Analysis complete!", value = 1.0)
  }
  
  # Compile final results
  results <- list(
    top_signatures = head(signature_rankings, top_n),
    pan_cluster_signatures = pan_cluster_signatures,
    cluster_specific_signatures = cluster_specific_signatures,
    all_signatures = signature_rankings,
    gene_pairs_analyzed = gene_pairs,
    analysis_summary = list(
      total_gene_pairs = nrow(gene_pairs),
      total_signatures = nrow(signature_rankings),
      pan_cluster_count = nrow(pan_cluster_signatures),
      cluster_specific_count = length(cluster_specific_signatures),
      strongest_gene_pair = if(nrow(signature_rankings) > 0) signature_rankings$gene_pair[1] else "None",
      top_signature_strength = if(nrow(signature_rankings) > 0) max(signature_rankings$signature_strength) else 0
    )
  )
  
  return(results)
}

#' Apply hierarchical FDR correction to Fisher's exact test p-values
#'
#' @param signature_rankings Data frame with signature rankings and p-values
#' @return Data frame with FDR-corrected p-values added
apply_hierarchical_fdr_correction <- function(signature_rankings) {
  
  if (nrow(signature_rankings) == 0) {
    return(signature_rankings)
  }
  
  cat("[FDR DEBUG] Starting hierarchical FDR correction for", nrow(signature_rankings), "signatures\n")
  
  # Step 1: Within-gene-pair FDR correction
  signature_rankings$gene_fisher_p_fdr <- NA
  signature_rankings$intersection_fisher_p_fdr <- NA  
  signature_rankings$union_fisher_p_fdr <- NA
  signature_rankings$pathway_fisher_p_fdr <- NA
  
  # Add database-specific FDR columns
  databases <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING")
  for (db in databases) {
    signature_rankings[[paste0(db, "_fisher_p_fdr")]] <- NA
  }
  
  # Correct within each gene pair
  gene_pairs <- unique(signature_rankings$gene_pair)
  cat("[FDR DEBUG] Applying within-gene-pair correction for", length(gene_pairs), "gene pairs\n")
  
  for (gene_pair in gene_pairs) {
    pair_indices <- which(signature_rankings$gene_pair == gene_pair)
    pair_data <- signature_rankings[pair_indices, ]
    
    # Collect all p-values for this gene pair
    all_p_values <- c()
    p_value_info <- list()
    
    # Gene-level Fisher's tests
    if ("gene_fisher_p" %in% names(pair_data)) {
      valid_gene_p <- !is.na(pair_data$gene_fisher_p)
      if (any(valid_gene_p)) {
        all_p_values <- c(all_p_values, pair_data$gene_fisher_p[valid_gene_p])
        p_value_info <- c(p_value_info, lapply(which(valid_gene_p), function(i) list(type = "gene_fisher_p", index = pair_indices[i])))
      }
    }
    
    # Intersection approach Fisher's tests
    if ("intersection_fisher_p" %in% names(pair_data)) {
      valid_int_p <- !is.na(pair_data$intersection_fisher_p)
      if (any(valid_int_p)) {
        all_p_values <- c(all_p_values, pair_data$intersection_fisher_p[valid_int_p])
        p_value_info <- c(p_value_info, lapply(which(valid_int_p), function(i) list(type = "intersection_fisher_p", index = pair_indices[i])))
      }
    }
    
    # Union approach Fisher's tests
    if ("union_fisher_p" %in% names(pair_data)) {
      valid_union_p <- !is.na(pair_data$union_fisher_p)
      if (any(valid_union_p)) {
        all_p_values <- c(all_p_values, pair_data$union_fisher_p[valid_union_p])
        p_value_info <- c(p_value_info, lapply(which(valid_union_p), function(i) list(type = "union_fisher_p", index = pair_indices[i])))
      }
    }
    
    # Legacy pathway Fisher's tests
    if ("pathway_fisher_p" %in% names(pair_data)) {
      valid_pathway_p <- !is.na(pair_data$pathway_fisher_p)
      if (any(valid_pathway_p)) {
        all_p_values <- c(all_p_values, pair_data$pathway_fisher_p[valid_pathway_p])
        p_value_info <- c(p_value_info, lapply(which(valid_pathway_p), function(i) list(type = "pathway_fisher_p", index = pair_indices[i])))
      }
    }
    
    # Database-specific pathway Fisher's tests
    for (db in databases) {
      db_col <- paste0(db, "_fisher_p")
      if (db_col %in% names(pair_data)) {
        valid_db_p <- !is.na(pair_data[[db_col]])
        if (any(valid_db_p)) {
          all_p_values <- c(all_p_values, pair_data[[db_col]][valid_db_p])
          p_value_info <- c(p_value_info, lapply(which(valid_db_p), function(i) list(type = db_col, index = pair_indices[i])))
        }
      }
    }
    
    # Apply FDR correction within this gene pair
    if (length(all_p_values) > 0) {
      fdr_corrected <- p.adjust(all_p_values, method = "BH")
      
      cat("[FDR DEBUG]", gene_pair, "- correcting", length(all_p_values), "p-values within gene pair\n")
      
      # Assign FDR-corrected values back to appropriate columns
      for (i in seq_along(fdr_corrected)) {
        info <- p_value_info[[i]]
        fdr_col <- gsub("_p$", "_p_fdr", info$type)
        signature_rankings[info$index, fdr_col] <- fdr_corrected[i]
      }
    }
  }
  
  # Step 2: Across-gene-pair FDR correction 
  cat("[FDR DEBUG] Applying across-gene-pair correction\n")
  
  # Collect all FDR-corrected p-values for second level correction
  all_fdr_p_values <- c()
  fdr_p_value_info <- list()
  
  fdr_columns <- c("gene_fisher_p_fdr", "intersection_fisher_p_fdr", "union_fisher_p_fdr", "pathway_fisher_p_fdr",
                   paste0(databases, "_fisher_p_fdr"))
  
  for (col in fdr_columns) {
    if (col %in% names(signature_rankings)) {
      valid_fdr_p <- !is.na(signature_rankings[[col]])
      if (any(valid_fdr_p)) {
        all_fdr_p_values <- c(all_fdr_p_values, signature_rankings[[col]][valid_fdr_p])
        fdr_p_value_info <- c(fdr_p_value_info, lapply(which(valid_fdr_p), function(i) list(column = col, index = i)))
      }
    }
  }
  
  # Apply second level FDR correction across all gene pairs
  if (length(all_fdr_p_values) > 0) {
    final_fdr_corrected <- p.adjust(all_fdr_p_values, method = "BH")
    
    cat("[FDR DEBUG] Final hierarchical correction on", length(all_fdr_p_values), "FDR-corrected p-values\n")
    
    # Create final FDR columns
    for (col in fdr_columns) {
      final_col <- gsub("_fdr$", "_fdr_hierarchical", col)
      signature_rankings[[final_col]] <- NA
    }
    
    # Assign final FDR-corrected values
    for (i in seq_along(final_fdr_corrected)) {
      info <- fdr_p_value_info[[i]]
      final_col <- gsub("_fdr$", "_fdr_hierarchical", info$column)
      signature_rankings[info$index, final_col] <- final_fdr_corrected[i]
    }
  }
  
  cat("[FDR DEBUG] Hierarchical FDR correction completed\n")
  return(signature_rankings)
}

#' Apply enhanced hierarchical FDR correction for factorial designs (v0.2.6)
#'
#' This function implements a research-based three-level FDR correction approach
#' optimized for factorial designs with experiments × directions × gene pairs.
#' 
#' Method: 
#' - Level 1: Benjamini-Yekutieli (BY) for dependent tests within gene pairs
#' - Level 2: Benjamini-Hochberg (BH) for independent tests across gene pairs
#' - Expected inflation factor: δ* ≈ 1.44 (documented in literature)
#'
#' @param signature_rankings Data frame with signature rankings and p-values
#' @param use_enhanced_method Logical, whether to use enhanced method (default TRUE)
#' @return Data frame with enhanced FDR-corrected p-values added
#' @export
apply_enhanced_fdr_correction_v026 <- function(signature_rankings, use_enhanced_method = TRUE) {
  
  if (!use_enhanced_method) {
    # Fallback to legacy method for backward compatibility
    return(apply_hierarchical_fdr_correction(signature_rankings))
  }
  
  if (nrow(signature_rankings) == 0) {
    return(signature_rankings)
  }
  
  cat("[ENHANCED FDR] Starting enhanced hierarchical FDR correction for", nrow(signature_rankings), "signatures\n")
  cat("[ENHANCED FDR] Method: Benjamini-Yekutieli (within gene pairs) + Benjamini-Hochberg (across gene pairs)\n")
  
  # Initialize enhanced FDR columns
  signature_rankings$gene_fisher_p_fdr_enhanced <- NA
  signature_rankings$intersection_fisher_p_fdr_enhanced <- NA  
  signature_rankings$union_fisher_p_fdr_enhanced <- NA
  signature_rankings$pathway_fisher_p_fdr_enhanced <- NA
  
  # Add database-specific enhanced FDR columns
  databases <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING")
  for (db in databases) {
    signature_rankings[[paste0(db, "_fisher_p_fdr_enhanced")]] <- NA
  }
  
  # Add experiment and direction-specific enhanced FDR columns
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  directions <- c("same_direction", "opposite_direction")
  
  for (exp in experiments) {
    for (direction in directions) {
      col_name <- paste0(exp, "_", direction, "_fisher_p_fdr_enhanced")
      signature_rankings[[col_name]] <- NA
    }
  }
  
  # Level 1: Within-gene-pair correction using Benjamini-Yekutieli (BY) method
  gene_pairs <- unique(signature_rankings$gene_pair)
  cat("[ENHANCED FDR] Level 1: Applying Benjamini-Yekutieli correction within", length(gene_pairs), "gene pairs\n")
  
  for (gene_pair in gene_pairs) {
    pair_indices <- which(signature_rankings$gene_pair == gene_pair)
    pair_data <- signature_rankings[pair_indices, ]
    
    # Collect all p-values for this gene pair (including experiment and direction specific)
    all_p_values <- c()
    p_value_info <- list()
    
    # Collect standard Fisher's test p-values
    standard_columns <- c("gene_fisher_p", "intersection_fisher_p", "union_fisher_p", "pathway_fisher_p",
                         paste0(databases, "_fisher_p"))
    
    for (col in standard_columns) {
      if (col %in% names(pair_data)) {
        valid_p <- !is.na(pair_data[[col]])
        if (any(valid_p)) {
          all_p_values <- c(all_p_values, pair_data[[col]][valid_p])
          p_value_info <- c(p_value_info, lapply(which(valid_p), function(i) list(type = col, index = pair_indices[i])))
        }
      }
    }
    
    # Collect experiment and direction-specific p-values
    for (exp in experiments) {
      for (direction in directions) {
        col_name <- paste0(exp, "_", direction, "_fisher_p")
        if (col_name %in% names(pair_data)) {
          valid_p <- !is.na(pair_data[[col_name]])
          if (any(valid_p)) {
            all_p_values <- c(all_p_values, pair_data[[col_name]][valid_p])
            p_value_info <- c(p_value_info, lapply(which(valid_p), function(i) list(type = col_name, index = pair_indices[i])))
          }
        }
      }
    }
    
    # Apply Benjamini-Yekutieli correction within gene pair (handles dependencies)
    if (length(all_p_values) > 0) {
      # BY method: more conservative than BH, appropriate for dependent tests
      fdr_corrected_by <- p.adjust(all_p_values, method = "BY")
      
      cat("[ENHANCED FDR]", gene_pair, "- BY correction on", length(all_p_values), "dependent tests within gene pair\n")
      
      # Assign BY-corrected values back to enhanced columns
      for (i in seq_along(fdr_corrected_by)) {
        info <- p_value_info[[i]]
        enhanced_col <- paste0(info$type, "_fdr_enhanced")
        if (grepl("_fisher_p$", info$type)) {
          enhanced_col <- gsub("_fisher_p$", "_fisher_p_fdr_enhanced", info$type)
        }
        signature_rankings[info$index, enhanced_col] <- fdr_corrected_by[i]
      }
    }
  }
  
  # Level 2: Across-gene-pair correction using Benjamini-Hochberg (BH) method
  cat("[ENHANCED FDR] Level 2: Applying Benjamini-Hochberg correction across gene pairs\n")
  
  # Collect all BY-corrected p-values for second level correction
  all_by_corrected_p_values <- c()
  by_p_value_info <- list()
  
  enhanced_fdr_columns <- c("gene_fisher_p_fdr_enhanced", "intersection_fisher_p_fdr_enhanced", 
                           "union_fisher_p_fdr_enhanced", "pathway_fisher_p_fdr_enhanced",
                           paste0(databases, "_fisher_p_fdr_enhanced"))
  
  # Add experiment and direction-specific enhanced columns
  for (exp in experiments) {
    for (direction in directions) {
      enhanced_fdr_columns <- c(enhanced_fdr_columns, paste0(exp, "_", direction, "_fisher_p_fdr_enhanced"))
    }
  }
  
  for (col in enhanced_fdr_columns) {
    if (col %in% names(signature_rankings)) {
      valid_by_p <- !is.na(signature_rankings[[col]])
      if (any(valid_by_p)) {
        all_by_corrected_p_values <- c(all_by_corrected_p_values, signature_rankings[[col]][valid_by_p])
        by_p_value_info <- c(by_p_value_info, lapply(which(valid_by_p), function(i) list(column = col, index = i)))
      }
    }
  }
  
  # Apply BH correction across gene pairs (gene pairs are independent)
  if (length(all_by_corrected_p_values) > 0) {
    final_fdr_corrected_bh <- p.adjust(all_by_corrected_p_values, method = "BH")
    
    cat("[ENHANCED FDR] Level 2: BH correction on", length(all_by_corrected_p_values), "independent tests across gene pairs\n")
    
    # Create final enhanced hierarchical FDR columns
    for (col in enhanced_fdr_columns) {
      final_col <- gsub("_fdr_enhanced$", "_fdr_enhanced_hierarchical", col)
      signature_rankings[[final_col]] <- NA
    }
    
    # Assign final BH-corrected values
    for (i in seq_along(final_fdr_corrected_bh)) {
      info <- by_p_value_info[[i]]
      final_col <- gsub("_fdr_enhanced$", "_fdr_enhanced_hierarchical", info$column)
      signature_rankings[info$index, final_col] <- final_fdr_corrected_bh[i]
    }
  }
  
  # Add method documentation as attributes
  attr(signature_rankings, "fdr_method_enhanced") <- "BY_within_pairs_BH_across_pairs"
  attr(signature_rankings, "correction_levels_enhanced") <- c("Benjamini_Yekutieli_within_gene_pairs", "Benjamini_Hochberg_across_gene_pairs")
  attr(signature_rankings, "expected_inflation_factor") <- 1.44
  attr(signature_rankings, "dependency_handling") <- "BY_method_for_dependent_tests_within_gene_pairs"
  attr(signature_rankings, "independence_assumption") <- "BH_method_for_independent_gene_pairs"
  
  cat("[ENHANCED FDR] Enhanced hierarchical FDR correction completed\n")
  cat("[ENHANCED FDR] Expected inflation factor: δ* ≈ 1.44 (literature-based)\n")
  cat("[ENHANCED FDR] Method documentation saved as data frame attributes\n")
  
  return(signature_rankings)
}

#' Compute signature strength rankings
#'
#' @param all_signatures List of signature analysis results
#' @return Data frame with ranked signatures
compute_signature_rankings <- function(all_signatures) {
  
  signature_rankings <- data.frame()
  
  for (gene_pair_name in names(all_signatures)) {
    analysis <- all_signatures[[gene_pair_name]]
    
    if ("cluster_results" %in% names(analysis)) {
      for (cluster in names(analysis$cluster_results)) {
        cluster_result <- analysis$cluster_results[[cluster]]
        
        if (!"error" %in% names(cluster_result) && !is.null(cluster_result$overlap_stats)) {
          overlap_stats <- cluster_result$overlap_stats
          pathway_stats <- cluster_result$pathway_overlap_stats
          
          # Handle both new (intersection/union) and legacy overlap stats
          if ("intersection_approach" %in% names(overlap_stats) && "union_approach" %in% names(overlap_stats)) {
            # New proper Fisher's exact test results
            intersection_stats <- overlap_stats$intersection_approach
            union_stats <- overlap_stats$union_approach
            
            cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- NEW PROPER FISHER'S TEST:\n")
            cat("  INTERSECTION: count=", intersection_stats$overlap_count %||% "NULL",
                ", fisher_p=", intersection_stats$fisher_p %||% "NULL",
                ", background=", intersection_stats$background_size %||% "NULL", "\n")
            cat("  UNION: count=", union_stats$overlap_count %||% "NULL",
                ", fisher_p=", union_stats$fisher_p %||% "NULL", 
                ", background=", union_stats$background_size %||% "NULL", "\n")
            
            # Use intersection approach for conservative results (prioritize for UI)
            overlap_stats_display <- intersection_stats
          } else {
            # Legacy overlap stats
            cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- LEGACY overlap_stats:",
                "count=", overlap_stats$overlap_count %||% "NULL",
                ", jaccard=", overlap_stats$jaccard_index %||% "NULL",
                ", fisher_p=", overlap_stats$fisher_p %||% "NULL", "\n")
            overlap_stats_display <- overlap_stats
          }
          
          # Skip if overlap_stats contains error
          if ("error" %in% names(overlap_stats_display)) {
            cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- SKIPPED: Error in overlap stats\n")
            next
          }
          
          # Extract direction statistics if available from enhanced analysis
          direction_stats <- NULL
          if ("enhanced_statistics" %in% names(cluster_result)) {
            # Enhanced analysis results available
            direction_stats <- cluster_result$enhanced_statistics
            cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- ENHANCED DIRECTION STATS available\n")
          } else if ("direction_analysis" %in% names(cluster_result)) {
            # Direction analysis results available
            direction_stats <- cluster_result$direction_analysis
            cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- DIRECTION STATS available\n")
          } else {
            # Check if we can derive direction stats from overlap stats
            if ("same_direction" %in% names(overlap_stats) || "opposite_direction" %in% names(overlap_stats)) {
              direction_stats <- list(
                same_direction_p = overlap_stats$same_direction$fisher_p %||% NA,
                opposite_direction_p = overlap_stats$opposite_direction$fisher_p %||% NA,
                biological_expectation = overlap_stats$biological_expectation %||% "unknown",
                primary_pattern = overlap_stats$primary_pattern %||% "unknown"
              )
              cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- DERIVED direction stats from overlap\n")
            } else {
              cat("[RANKING DEBUG]", gene_pair_name, "-", cluster, "- NO direction stats available\n")
            }
          }
          
          # Calculate composite score using display stats with enhanced direction information
          composite_score <- calculate_composite_signature_score(
            overlap_stats = overlap_stats_display,
            correlation_stats = NULL,  # Not available from enrichment data directly
            direction_stats = direction_stats,    # ← ENABLED: Enhanced direction statistics
            pathway_overlap_stats = pathway_stats
          )
          
          # Create signature entry with BOTH intersection and union info if available
          signature_entry <- data.frame(
            gene_pair = gene_pair_name,
            mast_gene = analysis$gene_pair$mast_gene,
            crispri_gene = analysis$gene_pair$crispri_gene,
            cluster = cluster,
            signature_strength = composite_score$composite_score,
            gene_overlap_count = overlap_stats_display$overlap_count %||% 0,
            gene_fisher_p = overlap_stats_display$fisher_p %||% NA,
            gene_jaccard = overlap_stats_display$jaccard_index %||% 0,
            # Add background size information (current: intersection approach for main display)
            background_size = overlap_stats_display$background_size %||% NA,
            background_type = overlap_stats_display$background_type %||% "legacy",
            # Add INTERSECTION approach details
            intersection_overlap_count = if("intersection_approach" %in% names(overlap_stats)) overlap_stats$intersection_approach$overlap_count %||% 0 else overlap_stats_display$overlap_count %||% 0,
            intersection_fisher_p = if("intersection_approach" %in% names(overlap_stats)) overlap_stats$intersection_approach$fisher_p %||% NA else overlap_stats_display$fisher_p %||% NA,
            intersection_background_size = if("intersection_approach" %in% names(overlap_stats)) overlap_stats$intersection_approach$background_size %||% NA else overlap_stats_display$background_size %||% NA,
            # Add UNION approach details
            union_overlap_count = if("union_approach" %in% names(overlap_stats)) overlap_stats$union_approach$overlap_count %||% 0 else 0,
            union_fisher_p = if("union_approach" %in% names(overlap_stats)) overlap_stats$union_approach$fisher_p %||% NA else NA,
            union_background_size = if("union_approach" %in% names(overlap_stats)) overlap_stats$union_approach$background_size %||% NA else NA,
            pathway_overlap_count = if(!is.null(pathway_stats)) pathway_stats$overlap_count %||% 0 else 0,
            pathway_fisher_p = if(!is.null(pathway_stats)) pathway_stats$fisher_p %||% NA else NA,
            mast_term_count = cluster_result$mast_term_count %||% 0,
            crispri_term_count = cluster_result$crispri_term_count %||% 0,
            stringsAsFactors = FALSE
          )
          
          # Add database-specific pathway Fisher's test results (NEW)
          database_pathway_stats <- cluster_result$database_specific_pathway_stats
          databases <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING")
          
          for (db in databases) {
            if (!is.null(database_pathway_stats) && db %in% names(database_pathway_stats)) {
              db_result <- database_pathway_stats[[db]]
              signature_entry[[paste0(db, "_overlap_count")]] <- db_result$overlap_count %||% 0
              signature_entry[[paste0(db, "_fisher_p")]] <- db_result$fisher_p %||% NA
              signature_entry[[paste0(db, "_jaccard")]] <- db_result$jaccard_index %||% 0
              signature_entry[[paste0(db, "_mast_pathway_count")]] <- db_result$mast_pathway_count %||% 0
              signature_entry[[paste0(db, "_crispri_pathway_count")]] <- db_result$crispri_pathway_count %||% 0
            } else {
              # No data for this database
              signature_entry[[paste0(db, "_overlap_count")]] <- 0
              signature_entry[[paste0(db, "_fisher_p")]] <- NA
              signature_entry[[paste0(db, "_jaccard")]] <- 0
              signature_entry[[paste0(db, "_mast_pathway_count")]] <- 0
              signature_entry[[paste0(db, "_crispri_pathway_count")]] <- 0
            }
          }
          
          signature_rankings <- rbind(signature_rankings, signature_entry)
        }
      }
    }
  }
  
  # Apply enhanced hierarchical FDR correction (v0.2.6)
  if (nrow(signature_rankings) > 0) {
    cat("[COMPUTE RANKINGS] Applying enhanced hierarchical FDR correction to", nrow(signature_rankings), "signatures\n")
    
    # Check if we should use enhanced method (default TRUE for v0.2.6)
    use_enhanced_fdr <- getOption("iscore.use_enhanced_analysis", TRUE)
    
    if (use_enhanced_fdr) {
      signature_rankings <- apply_enhanced_fdr_correction_v026(signature_rankings, use_enhanced_method = TRUE)
    } else {
      # Fallback to legacy method
      signature_rankings <- apply_hierarchical_fdr_correction(signature_rankings)
    }
  }
  
  # Sort by signature strength (descending)
  if (nrow(signature_rankings) > 0) {
    signature_rankings <- signature_rankings[order(signature_rankings$signature_strength, decreasing = TRUE), ]
    signature_rankings$rank <- seq_len(nrow(signature_rankings))
  }
  
  return(signature_rankings)
}

#' Identify pan-cluster signatures
#'
#' @param signature_rankings Data frame with signature rankings
#' @param min_cluster_breadth Minimum number of clusters
#' @return Data frame with pan-cluster signatures
identify_pan_cluster_signatures <- function(signature_rankings, min_cluster_breadth = 8) {
  
  if (nrow(signature_rankings) == 0) {
    return(data.frame())
  }
  
  # Count clusters per gene pair
  cluster_breadth <- signature_rankings %>%
    group_by(gene_pair) %>%
    summarise(
      cluster_count = n(),
      mean_signature_strength = mean(signature_strength, na.rm = TRUE),
      max_signature_strength = max(signature_strength, na.rm = TRUE),
      total_gene_overlaps = sum(gene_overlap_count, na.rm = TRUE),
      total_pathway_overlaps = sum(pathway_overlap_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Filter for pan-cluster signatures
  pan_cluster <- cluster_breadth[cluster_breadth$cluster_count >= min_cluster_breadth, ]
  
  if (nrow(pan_cluster) > 0) {
    # Sort by mean signature strength
    pan_cluster <- pan_cluster[order(pan_cluster$mean_signature_strength, decreasing = TRUE), ]
    pan_cluster$pan_cluster_rank <- seq_len(nrow(pan_cluster))
  }
  
  return(pan_cluster)
}

#' Identify cluster-specific signatures
#'
#' @param signature_rankings Data frame with signature rankings
#' @param min_strength_threshold Minimum signature strength for consideration
#' @return List of cluster-specific signatures by cluster
identify_cluster_specific_signatures <- function(signature_rankings, min_strength_threshold = 0.5) {
  
  if (nrow(signature_rankings) == 0) {
    return(list())
  }
  
  # Filter for strong signatures
  strong_signatures <- signature_rankings[signature_rankings$signature_strength >= min_strength_threshold, ]
  
  if (nrow(strong_signatures) == 0) {
    return(list())
  }
  
  # Group by cluster
  cluster_specific <- split(strong_signatures, strong_signatures$cluster)
  
  # Sort within each cluster by signature strength
  cluster_specific <- lapply(cluster_specific, function(cluster_data) {
    cluster_data[order(cluster_data$signature_strength, decreasing = TRUE), ]
  })
  
  return(cluster_specific)
}

#' Generate manuscript-ready signature summary
#'
#' @param signature_discovery_results Results from discover_top_signatures
#' @param output_dir Directory to save results
#' @return List with summary statistics and file paths
#' @export
generate_manuscript_signature_summary <- function(signature_discovery_results, output_dir = "manuscript_signatures") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Extract key results
  top_signatures <- signature_discovery_results$top_signatures
  pan_cluster <- signature_discovery_results$pan_cluster_signatures
  cluster_specific <- signature_discovery_results$cluster_specific_signatures
  
  # Generate summary statistics
  summary_stats <- list(
    total_signatures_analyzed = nrow(signature_discovery_results$all_signatures),
    top_signature_strength = if(nrow(top_signatures) > 0) max(top_signatures$signature_strength) else 0,
    pan_cluster_signatures = nrow(pan_cluster),
    cluster_specific_clusters = length(cluster_specific),
    strongest_gene_pair = if(nrow(top_signatures) > 0) top_signatures$gene_pair[1] else "None"
  )
  
  # Save top signatures table
  top_signatures_file <- file.path(output_dir, "top_signatures_for_manuscript.csv")
  if (nrow(top_signatures) > 0) {
    write.csv(top_signatures, top_signatures_file, row.names = FALSE)
  }
  
  # Save pan-cluster signatures
  pan_cluster_file <- file.path(output_dir, "pan_cluster_signatures.csv")
  if (nrow(pan_cluster) > 0) {
    write.csv(pan_cluster, pan_cluster_file, row.names = FALSE)
  }
  
  # Save cluster-specific signatures
  if (length(cluster_specific) > 0) {
    cluster_files <- list()
    for (cluster in names(cluster_specific)) {
      cluster_file <- file.path(output_dir, paste0("cluster_", cluster, "_specific_signatures.csv"))
      write.csv(cluster_specific[[cluster]], cluster_file, row.names = FALSE)
      cluster_files[[cluster]] <- cluster_file
    }
  }
  
  # Generate summary report
  summary_file <- file.path(output_dir, "manuscript_signature_summary.txt")
  summary_text <- paste(
    "=== Manuscript Signature Discovery Summary ===",
    "",
    paste("Analysis Date:", Sys.Date()),
    paste("Total Signatures Analyzed:", summary_stats$total_signatures_analyzed),
    paste("Top Signature Strength:", round(summary_stats$top_signature_strength, 2)),
    paste("Pan-Cluster Signatures Found:", summary_stats$pan_cluster_signatures),
    paste("Clusters with Specific Signatures:", summary_stats$cluster_specific_clusters),
    paste("Strongest Gene Pair:", summary_stats$strongest_gene_pair),
    "",
    "=== Files Generated ===",
    paste("- Top signatures:", top_signatures_file),
    paste("- Pan-cluster signatures:", pan_cluster_file),
    if(exists("cluster_files")) paste("- Cluster-specific files:", length(cluster_files), "files") else "",
    "",
    "=== Next Steps for Manuscript ===",
    "1. Review top signatures for biological relevance",
    "2. Focus on pan-cluster signatures for main findings", 
    "3. Investigate cluster-specific signatures for detailed analysis",
    "4. Generate heatmap visualizations for selected signatures",
    sep = "\n"
  )
  
  writeLines(summary_text, summary_file)
  
  return(list(
    summary_stats = summary_stats,
    files_generated = list(
      top_signatures = top_signatures_file,
      pan_cluster = pan_cluster_file,
      cluster_specific = if(exists("cluster_files")) cluster_files else list(),
      summary = summary_file
    ),
    manuscript_priorities = list(
      focus_on_pan_cluster = nrow(pan_cluster) > 0,
      strongest_pairs = if(nrow(top_signatures) > 0) head(top_signatures$gene_pair, 5) else character(),
      recommend_clusters = names(cluster_specific)[1:3]  # Top 3 clusters with specific signatures
    )
  ))
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a