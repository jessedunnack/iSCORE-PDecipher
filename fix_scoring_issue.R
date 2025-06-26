# Quick fix for scoring issue - patch the function in current session
# This directly updates the function that's causing overlap_score = 0

# Define the fixed function
calculate_composite_signature_score <- function(overlap_stats, correlation_stats, direction_stats,
                                               pathway_overlap_stats = NULL,
                                               weights = list(overlap = 0.3, correlation = 0.3, 
                                                            direction = 0.2, pathway = 0.2)) {
  
  # Initialize score components
  overlap_score <- 0
  correlation_score <- 0
  direction_score <- 0
  pathway_score <- 0
  
  # Helper function for null coalescing
  `%||%` <- function(a, b) if (is.null(a)) b else a
  
  # Overlap component (based on Fisher's p-value and Jaccard index)
  if (!is.null(overlap_stats)) {
    cat("[SCORE FIX DEBUG] overlap_stats available, jaccard:", overlap_stats$jaccard_index %||% "NULL", 
        ", overlap_count:", overlap_stats$overlap_count %||% "NULL", 
        ", fisher_p:", overlap_stats$fisher_p %||% "NULL", "\n")
    
    if (!is.null(overlap_stats$fisher_p) && !is.na(overlap_stats$fisher_p)) {
      # Convert p-value to score (higher score for lower p-value)
      overlap_score <- -log10(max(overlap_stats$fisher_p, 1e-10)) * overlap_stats$jaccard_index
      cat("[SCORE FIX DEBUG] Using Fisher p-value scoring: ", overlap_score, "\n")
    } else {
      # Use Jaccard index alone if no p-value (more robust fallback)
      jaccard_val <- overlap_stats$jaccard_index %||% 0
      overlap_count <- overlap_stats$overlap_count %||% 0
      # Score based on both Jaccard similarity and absolute overlap count
      overlap_score <- (jaccard_val * 10) + (log10(max(overlap_count, 1)) * 2)
      cat("[SCORE FIX DEBUG] Using fallback scoring: jaccard=", jaccard_val, ", count=", overlap_count, ", score=", overlap_score, "\n")
    }
  } else {
    cat("[SCORE FIX DEBUG] overlap_stats is NULL\n")
  }
  
  # Correlation component
  if (!is.null(correlation_stats) && !is.na(correlation_stats$correlation)) {
    # Use absolute correlation (strength regardless of direction)
    correlation_score <- abs(correlation_stats$correlation) * 10
  }
  
  # Direction consistency component
  if (!is.null(direction_stats)) {
    direction_score <- direction_stats$consistency_percent / 10  # Scale to 0-10
  }
  
  # Pathway overlap component (if provided)
  if (!is.null(pathway_overlap_stats) && !is.na(pathway_overlap_stats$fisher_p)) {
    pathway_score <- -log10(max(pathway_overlap_stats$fisher_p, 1e-10)) * 
                    pathway_overlap_stats$jaccard_index
  }
  
  # Calculate weighted composite score
  composite_score <- (overlap_score * weights$overlap +
                     correlation_score * weights$correlation +
                     direction_score * weights$direction +
                     pathway_score * weights$pathway)
  
  # Debug output for scoring
  cat("[SCORE FIX DEBUG] Final - Overlap score:", round(overlap_score, 3), 
      ", Correlation:", round(correlation_score, 3),
      ", Direction:", round(direction_score, 3), 
      ", Pathway:", round(pathway_score, 3),
      ", Composite:", round(composite_score, 3), "\n")
  
  return(list(
    composite_score = composite_score,
    components = list(
      overlap = overlap_score,
      correlation = correlation_score,
      direction = direction_score,
      pathway = pathway_score
    ),
    weights = weights
  ))
}

# Also fix the variant combining in the analysis function
analyze_gene_pair_signatures <- function(gene_pair, enrichment_data, clusters = NULL,
                                        include_pathways = TRUE, progress_callback = NULL) {
  
  if (!is.null(progress_callback)) {
    progress_callback(paste("Analyzing", gene_pair$mast_gene, "vs", gene_pair$crispri_gene))
  }
  
  # Filter data for this gene pair (handle variant combining)
  if (gene_pair$mast_gene == "SNCA_combined") {
    # Combine SNCA variants
    mast_data <- enrichment_data[enrichment_data$method == "MAST" & 
                                enrichment_data$mutation_perturbation %in% c("SNCA_A30P", "SNCA_A53T"), ]
  } else if (gene_pair$mast_gene == "VPS13C_combined") {
    # Combine VPS13C variants  
    mast_data <- enrichment_data[enrichment_data$method == "MAST" & 
                                enrichment_data$mutation_perturbation %in% c("VPS13C_A444P", "VPS13C_W395C"), ]
  } else {
    # Regular single gene lookup
    mast_data <- enrichment_data[enrichment_data$method == "MAST" & 
                                enrichment_data$mutation_perturbation == gene_pair$mast_gene, ]
  }
  
  crispri_data <- enrichment_data[enrichment_data$method == "MixScale" & 
                                 enrichment_data$mutation_perturbation == gene_pair$crispri_gene, ]
  
  cat("[GENE PAIR DEBUG] Gene pair:", gene_pair$mast_gene, "vs", gene_pair$crispri_gene, "\n")
  cat("[GENE PAIR DEBUG] MAST data found:", nrow(mast_data), "rows\n")
  cat("[GENE PAIR DEBUG] CRISPRi data found:", nrow(crispri_data), "rows\n")
  
  if (nrow(mast_data) > 0) {
    cat("[GENE PAIR DEBUG] MAST clusters:", paste(unique(mast_data$cluster), collapse = ", "), "\n")
  }
  if (nrow(crispri_data) > 0) {
    cat("[GENE PAIR DEBUG] CRISPRi clusters:", paste(unique(crispri_data$cluster), collapse = ", "), "\n")
  }
  
  # Filter by clusters if specified
  if (!is.null(clusters)) {
    mast_data <- mast_data[mast_data$cluster %in% clusters, ]
    crispri_data <- crispri_data[crispri_data$cluster %in% clusters, ]
  }
  
  if (nrow(mast_data) == 0 || nrow(crispri_data) == 0) {
    return(list(
      gene_pair = gene_pair,
      error = "No data available for this gene pair and cluster selection"
    ))
  }
  
  # Analyze by cluster
  cluster_results <- list()
  all_clusters <- unique(c(mast_data$cluster, crispri_data$cluster))
  
  for (j in seq_along(all_clusters)) {
    cluster <- all_clusters[j]
    
    if (!is.null(progress_callback)) {
      progress_callback(paste("cluster", cluster, "(", j, "of", length(all_clusters), ")"))
    }
    
    cluster_mast <- mast_data[mast_data$cluster == cluster, ]
    cluster_crispri <- crispri_data[crispri_data$cluster == cluster, ]
    
    cat("[CLUSTER DEBUG]", cluster, "- MAST terms:", nrow(cluster_mast), ", CRISPRi terms:", nrow(cluster_crispri), "\n")
    
    if (nrow(cluster_mast) == 0 || nrow(cluster_crispri) == 0) {
      cluster_results[[cluster]] <- list(error = "Missing data for one method")
      cat("[CLUSTER DEBUG]", cluster, "- SKIPPED: Missing data for one method\n")
      next
    }
    
    # Gene overlap analysis (using significant genes from enrichment terms)
    mast_genes <- unique(unlist(strsplit(cluster_mast$geneID, "/")))
    crispri_genes <- unique(unlist(strsplit(cluster_crispri$geneID, "/")))
    
    # Remove NA and empty strings
    mast_genes <- mast_genes[!is.na(mast_genes) & mast_genes != ""]
    crispri_genes <- crispri_genes[!is.na(crispri_genes) & crispri_genes != ""]
    
    cat("[CLUSTER DEBUG]", cluster, "- Extracted MAST genes:", length(mast_genes), "\n")
    cat("[CLUSTER DEBUG]", cluster, "- Extracted CRISPRi genes:", length(crispri_genes), "\n")
    if (length(mast_genes) > 0) {
      cat("[CLUSTER DEBUG]", cluster, "- Sample MAST genes:", paste(head(mast_genes, 5), collapse = ", "), "\n")
    }
    if (length(crispri_genes) > 0) {
      cat("[CLUSTER DEBUG]", cluster, "- Sample CRISPRi genes:", paste(head(crispri_genes, 5), collapse = ", "), "\n")
    }
    
    if (length(mast_genes) > 0 && length(crispri_genes) > 0) {
      background_genes <- unique(c(mast_genes, crispri_genes))
      overlap_stats <- calculate_gene_overlap_significance(mast_genes, crispri_genes, background_genes)
      
      cat("[OVERLAP DEBUG]", cluster, "- Overlap count:", overlap_stats$overlap_count, "\n")
      cat("[OVERLAP DEBUG]", cluster, "- Jaccard index:", round(overlap_stats$jaccard_index, 3), "\n")
      if (overlap_stats$overlap_count > 0) {
        cat("[OVERLAP DEBUG]", cluster, "- Overlapping genes:", paste(head(overlap_stats$overlap_genes, 5), collapse = ", "), "\n")
      }
    } else {
      overlap_stats <- list(error = "No valid gene lists extracted")
      cat("[OVERLAP DEBUG]", cluster, "- ERROR: No valid gene lists extracted\n")
    }
    
    # Pathway overlap analysis
    pathway_overlap_stats <- NULL
    if (include_pathways) {
      mast_pathways <- cluster_mast$Description
      crispri_pathways <- cluster_crispri$Description
      
      if (length(mast_pathways) > 0 && length(crispri_pathways) > 0) {
        pathway_overlap_stats <- calculate_gene_overlap_significance(
          mast_pathways, crispri_pathways, 
          background_genes = unique(c(mast_pathways, crispri_pathways))
        )
      }
    }
    
    cluster_results[[cluster]] <- list(
      cluster = cluster,
      overlap_stats = overlap_stats,
      pathway_overlap_stats = pathway_overlap_stats,
      mast_term_count = nrow(cluster_mast),
      crispri_term_count = nrow(cluster_crispri)
    )
  }
  
  return(list(
    gene_pair = gene_pair,
    cluster_results = cluster_results,
    total_clusters_analyzed = length(cluster_results)
  ))
}

cat("Fixed functions loaded. Now please run the signature analysis again and you should see [SCORE FIX DEBUG] messages.\n")