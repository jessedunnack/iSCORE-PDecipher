#' Comprehensive Fisher's Exact Test Analysis
#'
#' @description
#' Performs systematic analysis of cross-method significance for all MAST vs CRISPRi 
#' gene-cluster combinations using both intersection and union approaches.
#' This function replicates the exact logic from the DE Results module to provide
#' a comprehensive ranking of significant cross-method overlaps.
#'
#' @param de_data_path Path to full_DE_results.rds file
#' @param output_prefix Prefix for output CSV files (default: "comprehensive_fishers")
#' @param verbose Logical, whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \item{results}{Data frame with Fisher's exact test results for all combinations}
#' \item{gene_summary}{Data frame with gene-level summary statistics}
#' @export
#'
#' @examples
#' \dontrun{
#' analysis <- run_comprehensive_fishers_analysis("path/to/full_DE_results.rds")
#' View(analysis$results)
#' View(analysis$gene_summary)
#' }
run_comprehensive_fishers_analysis <- function(de_data_path = "full_DE_results.rds",
                                             output_prefix = "comprehensive_fishers",
                                             verbose = TRUE) {
  
  # Helper function for null coalescing
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  if (verbose) {
    cat("=== Comprehensive Fisher's Exact Test Analysis ===\n")
    cat("Analyzing MAST vs CRISPRi cross-method significance\n\n")
  }
  
  # Load DE data
  if (verbose) cat("Loading DE data...\n")
  if (!file.exists(de_data_path)) {
    stop("DE data file not found: ", de_data_path)
  }
  
  de_results <- readRDS(de_data_path)
  
  # Process volcano data
  if (verbose) cat("Processing MAST data...\n")
  mast_data <- process_mast_for_volcano(de_results$iSCORE_PD_MAST)
  if (verbose) {
    cat("  MAST volcano data:", nrow(mast_data), "rows\n")
    cat("  Available MAST genes:", paste(unique(mast_data$gene), collapse = ", "), "\n")
  }
  
  if (verbose) cat("Processing CRISPRi data...\n")
  crispri_data <- process_mixscale_for_volcano(de_results$CRISPRi_Mixscale)
  if (verbose) {
    cat("  CRISPRi volcano data:", nrow(crispri_data), "rows\n")
    cat("  Available CRISPRi genes:", paste(unique(crispri_data$gene), collapse = ", "), "\n")
  }
  
  # Define gene pairs for analysis (exclude GBA - MAST only)
  mast_genes <- unique(mast_data$gene)
  mast_genes <- mast_genes[mast_genes != "GBA"]  # Exclude GBA (no CRISPRi counterpart)
  
  # Get all clusters that have data for both methods
  mast_clusters <- unique(mast_data$cluster)
  crispri_clusters <- unique(crispri_data$cluster)
  common_clusters <- intersect(mast_clusters, crispri_clusters)
  
  if (verbose) {
    cat("\nAnalysis scope:\n")
    cat("  MAST genes to analyze:", length(mast_genes), "genes\n")
    cat("  Common clusters:", length(common_clusters), "clusters\n")
    cat("  Total combinations to test:", length(mast_genes) * length(common_clusters), "\n\n")
  }
  
  # Initialize results storage
  all_results <- data.frame()
  
  # Progress tracking
  total_combinations <- length(mast_genes) * length(common_clusters)
  current_combination <- 0
  
  # Analyze each gene-cluster combination
  for (mast_gene in mast_genes) {
    # Apply gene harmonization
    crispri_gene <- apply_gene_harmonization(mast_gene)
    
    # Check if CRISPRi data exists for this gene
    if (!crispri_gene %in% unique(crispri_data$gene)) {
      if (verbose) cat("  Skipping", mast_gene, "- no CRISPRi data for", crispri_gene, "\n")
      next
    }
    
    for (cluster in common_clusters) {
      current_combination <- current_combination + 1
      
      if (verbose && current_combination %% 10 == 0) {
        cat("  Progress:", current_combination, "/", total_combinations, 
            sprintf("(%.1f%%)", 100 * current_combination / total_combinations), "\n")
      }
      
      # Filter data for this combination
      mast_filtered <- mast_data[mast_data$gene == mast_gene & mast_data$cluster == cluster, ]
      crispri_filtered <- crispri_data[crispri_data$gene == crispri_gene & crispri_data$cluster == cluster, ]
      
      # Check if we have data for both methods
      if (nrow(mast_filtered) == 0 || nrow(crispri_filtered) == 0) {
        next
      }
      
      # Calculate significant genes (p < 0.05, |log2FC| > 0.25)
      mast_sig_idx <- !is.na(mast_filtered$pvalue) & !is.na(mast_filtered$log2FC) &
                     mast_filtered$pvalue < 0.05 & abs(mast_filtered$log2FC) > 0.25
      crispri_sig_idx <- !is.na(crispri_filtered$pvalue) & !is.na(crispri_filtered$log2FC) &
                        crispri_filtered$pvalue < 0.05 & abs(crispri_filtered$log2FC) > 0.25
      
      mast_sig_genes <- mast_filtered$gene_name[mast_sig_idx]
      crispri_sig_genes <- crispri_filtered$gene_name[crispri_sig_idx]
      
      mast_sig_count <- length(mast_sig_genes)
      crispri_sig_count <- length(crispri_sig_genes)
      overlap_count <- length(intersect(mast_sig_genes, crispri_sig_genes))
      
      # Skip if no significant genes in either method
      if (mast_sig_count == 0 || crispri_sig_count == 0) {
        next
      }
      
      # Get background genes (all genes tested)
      mast_background_genes <- unique(mast_filtered$gene_name)
      crispri_background_genes <- unique(crispri_filtered$gene_name)
      
      # Calculate Fisher's exact test for both approaches
      fisher_results <- list(
        intersection_approach = list(fisher_p = NA, fisher_or = NA, background_size = 0),
        union_approach = list(fisher_p = NA, fisher_or = NA, background_size = 0)
      )
      
      # INTERSECTION APPROACH (Conservative)
      intersection_background <- intersect(mast_background_genes, crispri_background_genes)
      intersection_size <- length(intersection_background)
      
      if (intersection_size > max(mast_sig_count, crispri_sig_count) && intersection_size > 0) {
        # Filter significant genes to intersection background
        mast_sig_filtered <- intersect(mast_sig_genes, intersection_background)
        crispri_sig_filtered <- intersect(crispri_sig_genes, intersection_background)
        overlap_filtered <- length(intersect(mast_sig_filtered, crispri_sig_filtered))
        
        mast_only <- length(mast_sig_filtered) - overlap_filtered
        crispri_only <- length(crispri_sig_filtered) - overlap_filtered  
        neither <- intersection_size - length(mast_sig_filtered) - crispri_only
        
        if (mast_only >= 0 && crispri_only >= 0 && neither >= 0) {
          contingency_matrix <- matrix(c(overlap_filtered, mast_only, crispri_only, neither), nrow = 2, byrow = TRUE)
          tryCatch({
            fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
            fisher_results$intersection_approach$fisher_p <- fisher_result$p.value
            fisher_results$intersection_approach$fisher_or <- fisher_result$estimate
            fisher_results$intersection_approach$background_size <- intersection_size
          }, error = function(e) {
            # Fisher's test failed
          })
        }
      }
      
      # UNION APPROACH (Liberal)
      union_background <- unique(c(mast_background_genes, crispri_background_genes))
      union_size <- length(union_background)
      
      if (union_size > max(mast_sig_count, crispri_sig_count) && union_size > 0) {
        mast_only <- mast_sig_count - overlap_count
        crispri_only <- crispri_sig_count - overlap_count
        neither <- union_size - mast_sig_count - crispri_only
        
        if (mast_only >= 0 && crispri_only >= 0 && neither >= 0) {
          contingency_matrix <- matrix(c(overlap_count, mast_only, crispri_only, neither), nrow = 2, byrow = TRUE)
          tryCatch({
            fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
            fisher_results$union_approach$fisher_p <- fisher_result$p.value
            fisher_results$union_approach$fisher_or <- fisher_result$estimate  
            fisher_results$union_approach$background_size <- union_size
          }, error = function(e) {
            # Fisher's test failed
          })
        }
      }
      
      # Store results
      result_row <- data.frame(
        mast_gene = mast_gene,
        crispri_gene = crispri_gene,
        cluster = cluster,
        mast_sig_count = mast_sig_count,
        crispri_sig_count = crispri_sig_count,
        overlap_count = overlap_count,
        overlap_percentage = round(100 * overlap_count / min(mast_sig_count, crispri_sig_count), 1),
        
        # Intersection approach
        intersection_p = fisher_results$intersection_approach$fisher_p,
        intersection_or = fisher_results$intersection_approach$fisher_or,
        intersection_background = fisher_results$intersection_approach$background_size,
        
        # Union approach  
        union_p = fisher_results$union_approach$fisher_p,
        union_or = fisher_results$union_approach$fisher_or,
        union_background = fisher_results$union_approach$background_size,
        
        stringsAsFactors = FALSE
      )
      
      all_results <- rbind(all_results, result_row)
    }
  }
  
  if (verbose) {
    cat("\nAnalysis complete!\n")
    cat("Total valid combinations analyzed:", nrow(all_results), "\n")
  }
  
  if (nrow(all_results) > 0) {
    # Add significance categories
    all_results$intersection_sig_level <- ifelse(is.na(all_results$intersection_p), "Not tested",
                                            ifelse(all_results$intersection_p < 0.001, "p < 0.001",
                                            ifelse(all_results$intersection_p < 0.01, "p < 0.01", 
                                            ifelse(all_results$intersection_p < 0.05, "p < 0.05", "Not significant"))))
    
    all_results$union_sig_level <- ifelse(is.na(all_results$union_p), "Not tested",
                                     ifelse(all_results$union_p < 0.001, "p < 0.001",
                                     ifelse(all_results$union_p < 0.01, "p < 0.01", 
                                     ifelse(all_results$union_p < 0.05, "p < 0.05", "Not significant"))))
    
    # Sort by most significant intersection p-value
    results_sorted <- all_results[order(all_results$intersection_p, na.last = TRUE), ]
    
    # Generate summary statistics
    if (verbose) {
      cat("\n=== SUMMARY STATISTICS ===\n")
      cat("Total combinations analyzed:", nrow(all_results), "\n")
      cat("Combinations with intersection p < 0.05:", sum(all_results$intersection_p < 0.05, na.rm = TRUE), "\n")
      cat("Combinations with intersection p < 0.01:", sum(all_results$intersection_p < 0.01, na.rm = TRUE), "\n")
      cat("Combinations with intersection p < 0.001:", sum(all_results$intersection_p < 0.001, na.rm = TRUE), "\n")
      cat("Combinations with union p < 0.05:", sum(all_results$union_p < 0.05, na.rm = TRUE), "\n")
      cat("Combinations with union p < 0.01:", sum(all_results$union_p < 0.01, na.rm = TRUE), "\n")
      cat("Combinations with union p < 0.001:", sum(all_results$union_p < 0.001, na.rm = TRUE), "\n")
    }
    
    # Gene-level summary
    gene_summary <- results_sorted %>%
      dplyr::filter(!is.na(intersection_p) & intersection_p < 0.05) %>%
      dplyr::group_by(mast_gene, crispri_gene) %>%
      dplyr::summarise(
        significant_clusters = n(),
        best_intersection_p = min(intersection_p, na.rm = TRUE),
        best_union_p = min(union_p, na.rm = TRUE),
        total_overlap = sum(overlap_count),
        avg_overlap_pct = mean(overlap_percentage),
        .groups = "drop"
      ) %>%
      dplyr::arrange(best_intersection_p)
    
    # Save results if output prefix provided
    if (!is.null(output_prefix) && output_prefix != "") {
      write.csv(results_sorted, paste0(output_prefix, "_results.csv"), row.names = FALSE)
      write.csv(gene_summary, paste0(output_prefix, "_gene_summary.csv"), row.names = FALSE)
      
      if (verbose) {
        cat("\nResults saved to:\n")
        cat("  ", paste0(output_prefix, "_results.csv"), "\n")
        cat("  ", paste0(output_prefix, "_gene_summary.csv"), "\n")
      }
    }
    
    return(list(
      results = results_sorted,
      gene_summary = gene_summary
    ))
    
  } else {
    if (verbose) cat("No valid combinations found for analysis.\n")
    return(list(
      results = data.frame(),
      gene_summary = data.frame()
    ))
  }
}

#' Helper function for p-value formatting
#' @noRd
format_pvalue <- function(p_val) {
  if (is.na(p_val)) return("NA")
  if (p_val >= 0.001) {
    return(sprintf("%.3f", p_val))
  } else {
    return(format(p_val, digits = 2, scientific = TRUE))
  }
}

#' Apply gene harmonization for MAST to CRISPRi mapping
#' @noRd
apply_gene_harmonization <- function(mast_gene) {
  # Map MAST gene names to corresponding CRISPRi names
  mixscale_gene <- mast_gene
  
  if (mast_gene == "PRKN") {
    mixscale_gene <- "PARK2"
  } else if (mast_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
    mixscale_gene <- "SNCA" 
  } else if (mast_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
    mixscale_gene <- "VPS13C"
  }
  
  return(mixscale_gene)
}