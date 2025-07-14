#' Signature Analysis Functions for Cross-Method Comparison
#'
#' This module implements statistical methods for identifying shared signatures
#' between MAST (mutations) and CRISPRi (knockdowns) experiments.

#' Calculate gene overlap significance with proper background gene handling
#'
#' @param mast_genes Character vector of significant genes from MAST
#' @param crispri_genes Character vector of significant genes from CRISPRi  
#' @param mast_background_genes Character vector of all genes tested in MAST
#' @param crispri_background_genes Character vector of all genes tested in CRISPRi
#' @param alternative Character, test direction ("greater", "two.sided", "less")
#' @return List with overlap statistics for both intersection and union approaches
#' @export
calculate_gene_overlap_significance_proper <- function(mast_genes, crispri_genes, 
                                                      mast_background_genes = NULL,
                                                      crispri_background_genes = NULL,
                                                      alternative = "greater") {
  
  # Input validation
  if (length(mast_genes) == 0 || length(crispri_genes) == 0) {
    return(list(
      intersection_approach = list(
        overlap_count = 0,
        mast_count = length(mast_genes),
        crispri_count = length(crispri_genes),
        fisher_p = NA,
        fisher_or = NA,
        jaccard_index = 0,
        background_size = 0,
        background_type = "intersection",
        error = "Empty gene lists"
      ),
      union_approach = list(
        overlap_count = 0,
        mast_count = length(mast_genes),
        crispri_count = length(crispri_genes),
        fisher_p = NA,
        fisher_or = NA,
        jaccard_index = 0,
        background_size = 0,
        background_type = "union",
        error = "Empty gene lists"
      )
    ))
  }
  
  # Calculate basic overlap statistics (same for both approaches)
  overlap_genes <- intersect(mast_genes, crispri_genes)
  overlap_count <- length(overlap_genes)
  union_count <- length(union(mast_genes, crispri_genes))
  jaccard_index <- ifelse(union_count > 0, overlap_count / union_count, 0)
  
  # Initialize results for both approaches
  results <- list(
    intersection_approach = NULL,
    union_approach = NULL
  )
  
  # Approach 1: Intersection background (genes tested in BOTH methods)
  if (!is.null(mast_background_genes) && !is.null(crispri_background_genes)) {
    intersection_background <- intersect(mast_background_genes, crispri_background_genes)
    
    if (length(intersection_background) > 0) {
      # Filter gene lists to intersection background
      mast_genes_filtered <- intersect(mast_genes, intersection_background)
      crispri_genes_filtered <- intersect(crispri_genes, intersection_background)
      overlap_genes_filtered <- intersect(mast_genes_filtered, crispri_genes_filtered)
      overlap_count_filtered <- length(overlap_genes_filtered)
      
      # Calculate Fisher's exact test for intersection approach
      fisher_stats_int <- calculate_fisher_test(
        mast_genes_filtered, crispri_genes_filtered, intersection_background, alternative
      )
      
      results$intersection_approach <- list(
        overlap_count = overlap_count_filtered,
        overlap_genes = overlap_genes_filtered,
        mast_count = length(mast_genes_filtered),
        crispri_count = length(crispri_genes_filtered),
        fisher_p = fisher_stats_int$fisher_p,
        fisher_or = fisher_stats_int$fisher_or,
        jaccard_index = ifelse(length(union(mast_genes_filtered, crispri_genes_filtered)) > 0,
                              overlap_count_filtered / length(union(mast_genes_filtered, crispri_genes_filtered)), 0),
        background_size = length(intersection_background),
        background_type = "intersection",
        contingency_table = fisher_stats_int$contingency_table
      )
    } else {
      results$intersection_approach <- list(
        overlap_count = 0,
        mast_count = 0,
        crispri_count = 0,
        fisher_p = NA,
        fisher_or = NA,
        jaccard_index = 0,
        background_size = 0,
        background_type = "intersection",
        error = "No genes in intersection background"
      )
    }
    
    # Approach 2: Union background (genes tested in EITHER method)
    union_background <- unique(c(mast_background_genes, crispri_background_genes))
    
    # Filter gene lists to union background
    mast_genes_union <- intersect(mast_genes, union_background)
    crispri_genes_union <- intersect(crispri_genes, union_background)
    overlap_genes_union <- intersect(mast_genes_union, crispri_genes_union)
    overlap_count_union <- length(overlap_genes_union)
    
    # Calculate Fisher's exact test for union approach
    fisher_stats_union <- calculate_fisher_test(
      mast_genes_union, crispri_genes_union, union_background, alternative
    )
    
    results$union_approach <- list(
      overlap_count = overlap_count_union,
      overlap_genes = overlap_genes_union,
      mast_count = length(mast_genes_union),
      crispri_count = length(crispri_genes_union),
      fisher_p = fisher_stats_union$fisher_p,
      fisher_or = fisher_stats_union$fisher_or,
      jaccard_index = ifelse(length(union(mast_genes_union, crispri_genes_union)) > 0,
                            overlap_count_union / length(union(mast_genes_union, crispri_genes_union)), 0),
      background_size = length(union_background),
      background_type = "union",
      contingency_table = fisher_stats_union$contingency_table
    )
  } else {
    # Fallback if no proper background provided
    results$intersection_approach <- list(
      overlap_count = overlap_count,
      overlap_genes = overlap_genes,
      mast_count = length(mast_genes),
      crispri_count = length(crispri_genes),
      fisher_p = NA,
      fisher_or = NA,
      jaccard_index = jaccard_index,
      background_size = 0,
      background_type = "intersection",
      error = "No background gene information provided"
    )
    
    results$union_approach <- list(
      overlap_count = overlap_count,
      overlap_genes = overlap_genes,
      mast_count = length(mast_genes),
      crispri_count = length(crispri_genes),
      fisher_p = NA,
      fisher_or = NA,
      jaccard_index = jaccard_index,
      background_size = 0,
      background_type = "union",
      error = "No background gene information provided"
    )
  }
  
  return(results)
}

#' Helper function to calculate Fisher's exact test
#'
#' @param genes1 First gene set
#' @param genes2 Second gene set
#' @param background_genes Background gene universe
#' @param alternative Test alternative
#' @return List with Fisher's test results
calculate_fisher_test <- function(genes1, genes2, background_genes, alternative = "greater") {
  
  # Calculate overlap
  overlap_genes <- intersect(genes1, genes2)
  overlap_count <- length(overlap_genes)
  
  # Contingency table for Fisher's test
  genes1_only <- length(genes1) - overlap_count
  genes2_only <- length(genes2) - overlap_count  
  neither <- length(background_genes) - length(genes1) - length(genes2) + overlap_count
  
  fisher_p <- NA
  fisher_or <- NA
  contingency_matrix <- NULL
  
  if (genes1_only >= 0 && genes2_only >= 0 && neither >= 0) {
    contingency_matrix <- matrix(
      c(overlap_count, genes1_only, genes2_only, neither),
      nrow = 2, byrow = TRUE
    )
    
    tryCatch({
      fisher_result <- fisher.test(contingency_matrix, alternative = alternative)
      fisher_p <- fisher_result$p.value
      fisher_or <- fisher_result$estimate
    }, error = function(e) {
      warning(paste("Fisher's exact test failed:", e$message))
    })
  }
  
  return(list(
    fisher_p = fisher_p,
    fisher_or = fisher_or,
    contingency_table = contingency_matrix
  ))
}

#' Calculate gene overlap significance between two methods (LEGACY - keep for compatibility)
#'
#' @param mast_genes Character vector of significant genes from MAST
#' @param crispri_genes Character vector of significant genes from CRISPRi  
#' @param background_genes Character vector of background genes tested in both methods
#' @param alternative Character, test direction ("greater", "two.sided", "less")
#' @return List with overlap statistics
#' @export
calculate_gene_overlap_significance <- function(mast_genes, crispri_genes, 
                                               background_genes = NULL,
                                               alternative = "greater") {
  
  # Input validation
  if (length(mast_genes) == 0 || length(crispri_genes) == 0) {
    return(list(
      overlap_count = 0,
      mast_count = length(mast_genes),
      crispri_count = length(crispri_genes),
      fisher_p = NA,
      fisher_or = NA,
      jaccard_index = 0,
      background_size = length(background_genes),
      error = "Empty gene lists"
    ))
  }
  
  # Calculate overlap
  overlap_genes <- intersect(mast_genes, crispri_genes)
  overlap_count <- length(overlap_genes)
  
  # Calculate Jaccard similarity index
  union_count <- length(union(mast_genes, crispri_genes))
  jaccard_index <- ifelse(union_count > 0, overlap_count / union_count, 0)
  
  # Fisher's exact test (if background provided)
  fisher_p <- NA
  fisher_or <- NA
  
  if (!is.null(background_genes) && length(background_genes) > 0) {
    # Ensure all genes are in background
    mast_genes <- intersect(mast_genes, background_genes)
    crispri_genes <- intersect(crispri_genes, background_genes)
    overlap_genes <- intersect(mast_genes, crispri_genes)
    overlap_count <- length(overlap_genes)
    
    # Contingency table for Fisher's test
    mast_only <- length(mast_genes) - overlap_count
    crispri_only <- length(crispri_genes) - overlap_count  
    neither <- length(background_genes) - length(mast_genes) - length(crispri_genes) + overlap_count
    
    if (mast_only >= 0 && crispri_only >= 0 && neither >= 0) {
      contingency_matrix <- matrix(
        c(overlap_count, mast_only, crispri_only, neither),
        nrow = 2, byrow = TRUE
      )
      
      tryCatch({
        fisher_result <- fisher.test(contingency_matrix, alternative = alternative)
        fisher_p <- fisher_result$p.value
        fisher_or <- fisher_result$estimate
      }, error = function(e) {
        warning(paste("Fisher's exact test failed:", e$message))
      })
    }
  }
  
  return(list(
    overlap_count = overlap_count,
    overlap_genes = overlap_genes,
    mast_count = length(mast_genes),
    crispri_count = length(crispri_genes),
    fisher_p = fisher_p,
    fisher_or = fisher_or,
    jaccard_index = jaccard_index,
    background_size = length(background_genes),
    contingency_table = if(exists("contingency_matrix")) contingency_matrix else NULL
  ))
}

#' Calculate effect size correlation between methods
#'
#' @param mast_data Data frame with MAST results (gene_name, log2FC, pvalue)
#' @param crispri_data Data frame with CRISPRi results (gene_name, log2FC, pvalue)
#' @param method Character, correlation method ("pearson", "spearman") 
#' @return List with correlation statistics
#' @export
calculate_effect_size_correlation <- function(mast_data, crispri_data, method = "pearson") {
  
  # Input validation
  if (is.null(mast_data) || is.null(crispri_data) || 
      nrow(mast_data) == 0 || nrow(crispri_data) == 0) {
    return(list(
      correlation = NA,
      p_value = NA,
      n_genes = 0,
      error = "Empty data frames"
    ))
  }
  
  # Ensure required columns exist
  required_cols <- c("gene_name", "log2FC")
  if (!all(required_cols %in% names(mast_data)) || 
      !all(required_cols %in% names(crispri_data))) {
    return(list(
      correlation = NA,
      p_value = NA,
      n_genes = 0,
      error = "Missing required columns"
    ))
  }
  
  # Merge data by gene name
  merged_data <- merge(mast_data[, c("gene_name", "log2FC")], 
                       crispri_data[, c("gene_name", "log2FC")],
                       by = "gene_name", suffixes = c("_mast", "_crispri"))
  
  if (nrow(merged_data) < 3) {
    return(list(
      correlation = NA,
      p_value = NA,
      n_genes = nrow(merged_data),
      error = "Insufficient overlapping genes for correlation"
    ))
  }
  
  # Remove rows with missing values
  merged_data <- merged_data[complete.cases(merged_data), ]
  
  if (nrow(merged_data) < 3) {
    return(list(
      correlation = NA,
      p_value = NA,
      n_genes = nrow(merged_data),
      error = "Insufficient complete cases for correlation"
    ))
  }
  
  # Calculate correlation
  tryCatch({
    cor_result <- cor.test(merged_data$log2FC_mast, merged_data$log2FC_crispri, method = method)
    
    return(list(
      correlation = cor_result$estimate,
      p_value = cor_result$p.value,
      confidence_interval = cor_result$conf.int,
      n_genes = nrow(merged_data),
      merged_data = merged_data,
      method = method
    ))
  }, error = function(e) {
    return(list(
      correlation = NA,
      p_value = NA,
      n_genes = nrow(merged_data),
      error = paste("Correlation calculation failed:", e$message)
    ))
  })
}

#' Calculate direction consistency between methods
#'
#' @param mast_data Data frame with MAST results
#' @param crispri_data Data frame with CRISPRi results
#' @param fc_threshold Numeric, log2FC threshold for direction determination
#' @return List with direction consistency statistics
#' @export
calculate_direction_consistency <- function(mast_data, crispri_data, fc_threshold = 0) {
  
  # Merge data by gene name
  merged_data <- merge(mast_data[, c("gene_name", "log2FC")], 
                       crispri_data[, c("gene_name", "log2FC")],
                       by = "gene_name", suffixes = c("_mast", "_crispri"))
  
  if (nrow(merged_data) == 0) {
    return(list(
      same_direction_up = 0,
      same_direction_down = 0,
      opposite_direction = 0,
      total_genes = 0,
      consistency_percent = 0,
      error = "No overlapping genes"
    ))
  }
  
  # Classify directions
  merged_data$mast_direction <- ifelse(merged_data$log2FC_mast > fc_threshold, "up",
                                      ifelse(merged_data$log2FC_mast < -fc_threshold, "down", "neutral"))
  merged_data$crispri_direction <- ifelse(merged_data$log2FC_crispri > fc_threshold, "up",
                                         ifelse(merged_data$log2FC_crispri < -fc_threshold, "down", "neutral"))
  
  # Count direction patterns
  same_direction_up <- sum(merged_data$mast_direction == "up" & merged_data$crispri_direction == "up")
  same_direction_down <- sum(merged_data$mast_direction == "down" & merged_data$crispri_direction == "down")
  opposite_direction <- sum((merged_data$mast_direction == "up" & merged_data$crispri_direction == "down") |
                           (merged_data$mast_direction == "down" & merged_data$crispri_direction == "up"))
  
  total_directional <- same_direction_up + same_direction_down + opposite_direction
  consistency_percent <- ifelse(total_directional > 0, 
                               (same_direction_up + same_direction_down) / total_directional * 100, 0)
  
  return(list(
    same_direction_up = same_direction_up,
    same_direction_down = same_direction_down,
    opposite_direction = opposite_direction,
    total_genes = nrow(merged_data),
    total_directional = total_directional,
    consistency_percent = consistency_percent,
    detailed_data = merged_data
  ))
}

#' Calculate composite signature strength score
#'
#' @param overlap_stats List from calculate_gene_overlap_significance
#' @param correlation_stats List from calculate_effect_size_correlation  
#' @param direction_stats List from calculate_direction_consistency
#' @param pathway_overlap_stats List with pathway overlap information (optional)
#' @param weights List with weights for different components
#' @return Numeric composite score
#' @export
calculate_composite_signature_score <- function(overlap_stats, correlation_stats, direction_stats,
                                               pathway_overlap_stats = NULL,
                                               weights = list(overlap = 0.3, correlation = 0.3, 
                                                            direction = 0.2, pathway = 0.2)) {
  
  # Initialize score components
  overlap_score <- 0
  correlation_score <- 0
  direction_score <- 0
  pathway_score <- 0
  
  # Overlap component (based on Fisher's p-value and Jaccard index)
  if (!is.null(overlap_stats)) {
    jaccard_val <- overlap_stats$jaccard_index %||% 0
    overlap_count <- overlap_stats$overlap_count %||% 0
    fisher_p <- overlap_stats$fisher_p %||% 1
    
    cat("[SCORE FIX DEBUG] overlap_stats available, jaccard:", jaccard_val, 
        ", overlap_count:", overlap_count, 
        ", fisher_p:", fisher_p, "\n")
    
    # Use Fisher's p-value only if it's significant (p < 0.1)
    if (!is.null(fisher_p) && !is.na(fisher_p) && fisher_p < 0.1) {
      # Convert p-value to score (higher score for lower p-value)
      overlap_score <- -log10(max(fisher_p, 1e-10)) * jaccard_val
      cat("[SCORE FIX DEBUG] Using Fisher p-value scoring: ", overlap_score, "\n")
    } else {
      # Use Jaccard index + count when Fisher's test is not significant
      # Scale based on both similarity and absolute evidence
      overlap_score <- (jaccard_val * 20) + (log10(max(overlap_count, 1)) * 3)
      cat("[SCORE FIX DEBUG] Using Jaccard+count scoring: jaccard=", jaccard_val, 
          ", count=", overlap_count, ", score=", overlap_score, "\n")
    }
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
  cat("[SCORE DEBUG] Overlap score:", round(overlap_score, 3), 
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

#' Identify PD-relevant pathway enrichments
#'
#' @param enrichment_data Enrichment results data frame
#' @param pd_pathway_terms Character vector of PD-relevant terms
#' @return Filtered enrichment data with PD-relevant pathways
#' @export
identify_pd_relevant_enrichments <- function(enrichment_data, pd_pathway_terms = NULL) {
  
  if (is.null(pd_pathway_terms)) {
    # Load from gene_harmonization.R
    pd_pathway_terms <- get_pd_relevant_pathways()
  }
  
  if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
    return(data.frame())
  }
  
  # Create pattern for matching PD-relevant terms
  pd_pattern <- paste(pd_pathway_terms, collapse = "|")
  
  # Filter enrichment data for PD-relevant terms
  pd_enrichments <- enrichment_data[
    grepl(pd_pattern, enrichment_data$Description, ignore.case = TRUE) |
    grepl(pd_pattern, enrichment_data$ID, ignore.case = TRUE), 
  ]
  
  # Add PD relevance score based on number of matching keywords
  if (nrow(pd_enrichments) > 0) {
    pd_enrichments$pd_relevance_score <- sapply(pd_enrichments$Description, function(desc) {
      sum(sapply(pd_pathway_terms, function(term) {
        grepl(term, desc, ignore.case = TRUE)
      }))
    })
  }
  
  return(pd_enrichments)
}

#' Calculate database-specific pathway overlap significance
#'
#' @param mast_data Enrichment data for MAST method
#' @param crispri_data Enrichment data for CRISPRi method
#' @param databases Character vector of databases to analyze (e.g., "GO_BP", "KEGG", "Reactome")
#' @return List with database-specific pathway overlap statistics
#' @export
calculate_pathway_overlap_by_database <- function(mast_data, crispri_data, 
                                                 databases = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING")) {
  
  results <- list()
  
  for (database in databases) {
    # Filter data for this specific database
    mast_db <- mast_data[mast_data$enrichment_type == database, ]
    crispri_db <- crispri_data[crispri_data$enrichment_type == database, ]
    
    if (nrow(mast_db) > 0 && nrow(crispri_db) > 0) {
      # Extract pathway terms for this database
      mast_pathways <- mast_db$Description
      crispri_pathways <- crispri_db$Description
      
      # Calculate overlap using existing function
      overlap_stats <- calculate_gene_overlap_significance(
        mast_pathways, crispri_pathways, 
        background_genes = unique(c(mast_pathways, crispri_pathways))
      )
      
      results[[database]] <- list(
        database = database,
        mast_pathway_count = length(mast_pathways),
        crispri_pathway_count = length(crispri_pathways),
        overlap_count = overlap_stats$overlap_count,
        fisher_p = overlap_stats$fisher_p,
        fisher_or = overlap_stats$fisher_or,
        jaccard_index = overlap_stats$jaccard_index,
        overlapping_pathways = overlap_stats$overlap_genes
      )
    } else {
      # No data for this database in this comparison
      results[[database]] <- list(
        database = database,
        mast_pathway_count = nrow(mast_db),
        crispri_pathway_count = nrow(crispri_db),
        overlap_count = 0,
        fisher_p = NA,
        fisher_or = NA,
        jaccard_index = 0,
        overlapping_pathways = character(0),
        error = paste("Insufficient data: MAST =", nrow(mast_db), "CRISPRi =", nrow(crispri_db))
      )
    }
  }
  
  return(results)
}

#' Perform comprehensive signature analysis for a gene pair
#'
#' @param gene_pair List with mast_gene and crispri_gene names
#' @param enrichment_data Consolidated enrichment data
#' @param de_data Original differential expression data structure (optional, for proper background genes)
#' @param clusters Character vector of clusters to analyze (NULL for all)
#' @param include_pathways Logical, whether to include pathway analysis
#' @param progress_callback Function to call for progress updates (optional)
#' @return List with comprehensive signature analysis results
#' @export
analyze_gene_pair_signatures <- function(gene_pair, enrichment_data, de_data = NULL, clusters = NULL,
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
      
      # Extract proper background genes if DE data is available
      mast_background_genes <- NULL
      crispri_background_genes <- NULL
      
      if (!is.null(de_data)) {
        # For MAST data - handle variant combining
        if (gene_pair$mast_gene == "SNCA_combined") {
          # Try both SNCA variants and take the largest background
          snca_backgrounds <- list()
          for (variant in c("SNCA_A30P", "SNCA_A53T")) {
            if (!is.null(de_data$iSCORE_PD_MAST[[variant]][[cluster]])) {
              snca_backgrounds[[variant]] <- de_data$iSCORE_PD_MAST[[variant]][[cluster]]$background_genes
            }
          }
          if (length(snca_backgrounds) > 0) {
            mast_background_genes <- unique(unlist(snca_backgrounds))
          }
        } else if (gene_pair$mast_gene == "VPS13C_combined") {
          # Try both VPS13C variants and take the largest background
          vps13c_backgrounds <- list()
          for (variant in c("VPS13C_A444P", "VPS13C_W395C")) {
            if (!is.null(de_data$iSCORE_PD_MAST[[variant]][[cluster]])) {
              vps13c_backgrounds[[variant]] <- de_data$iSCORE_PD_MAST[[variant]][[cluster]]$background_genes
            }
          }
          if (length(vps13c_backgrounds) > 0) {
            mast_background_genes <- unique(unlist(vps13c_backgrounds))
          }
        } else {
          # Regular single gene lookup
          if (!is.null(de_data$iSCORE_PD_MAST[[gene_pair$mast_gene]][[cluster]])) {
            mast_background_genes <- de_data$iSCORE_PD_MAST[[gene_pair$mast_gene]][[cluster]]$background_genes
          }
        }
        
        # For CRISPRi data 
        if (!is.null(de_data$CRISPRi_Mixscale[[gene_pair$crispri_gene]][[cluster]])) {
          crispri_background_genes <- de_data$CRISPRi_Mixscale[[gene_pair$crispri_gene]][[cluster]]$background_genes
        }
        
        cat("[BACKGROUND DEBUG]", cluster, "- MAST background genes:", length(mast_background_genes %||% character(0)), "\n")
        cat("[BACKGROUND DEBUG]", cluster, "- CRISPRi background genes:", length(crispri_background_genes %||% character(0)), "\n")
      }
      
      # Use proper background calculation if available, otherwise fall back to old method
      if (!is.null(mast_background_genes) && !is.null(crispri_background_genes)) {
        overlap_stats <- calculate_gene_overlap_significance_proper(
          mast_genes, crispri_genes, 
          mast_background_genes, crispri_background_genes
        )
        
        cat("[OVERLAP DEBUG]", cluster, "- INTERSECTION approach:\n")
        cat("  Background size:", overlap_stats$intersection_approach$background_size, "\n")
        cat("  Overlap count:", overlap_stats$intersection_approach$overlap_count, "\n")
        cat("  Fisher p-value:", overlap_stats$intersection_approach$fisher_p, "\n")
        cat("  Jaccard index:", round(overlap_stats$intersection_approach$jaccard_index, 3), "\n")
        
        cat("[OVERLAP DEBUG]", cluster, "- UNION approach:\n")
        cat("  Background size:", overlap_stats$union_approach$background_size, "\n")
        cat("  Overlap count:", overlap_stats$union_approach$overlap_count, "\n")
        cat("  Fisher p-value:", overlap_stats$union_approach$fisher_p, "\n")
        cat("  Jaccard index:", round(overlap_stats$union_approach$jaccard_index, 3), "\n")
        
      } else {
        # Fallback to old method with warning
        cat("[OVERLAP DEBUG]", cluster, "- WARNING: Using legacy background calculation (union of significant genes)\n")
        background_genes <- unique(c(mast_genes, crispri_genes))
        overlap_stats <- calculate_gene_overlap_significance(mast_genes, crispri_genes, background_genes)
        
        cat("[OVERLAP DEBUG]", cluster, "- Overlap count:", overlap_stats$overlap_count, "\n")
        cat("[OVERLAP DEBUG]", cluster, "- Jaccard index:", round(overlap_stats$jaccard_index, 3), "\n")
        if (overlap_stats$overlap_count > 0) {
          cat("[OVERLAP DEBUG]", cluster, "- Overlapping genes:", paste(head(overlap_stats$overlap_genes, 5), collapse = ", "), "\n")
        }
      }
    } else {
      overlap_stats <- list(error = "No valid gene lists extracted")
      cat("[OVERLAP DEBUG]", cluster, "- ERROR: No valid gene lists extracted\n")
    }
    
    # Pathway overlap analysis (database-specific)
    pathway_overlap_stats <- NULL
    database_specific_pathway_stats <- NULL
    if (include_pathways) {
      # Legacy pathway overlap (all pathways combined)
      mast_pathways <- cluster_mast$Description
      crispri_pathways <- cluster_crispri$Description
      
      if (length(mast_pathways) > 0 && length(crispri_pathways) > 0) {
        pathway_overlap_stats <- calculate_gene_overlap_significance(
          mast_pathways, crispri_pathways, 
          background_genes = unique(c(mast_pathways, crispri_pathways))
        )
      }
      
      # Database-specific pathway overlap analysis (NEW)
      if (nrow(cluster_mast) > 0 && nrow(cluster_crispri) > 0) {
        database_specific_pathway_stats <- calculate_pathway_overlap_by_database(
          cluster_mast, cluster_crispri
        )
        
        cat("[PATHWAY DEBUG]", cluster, "- Database-specific pathway analysis completed for",
            length(database_specific_pathway_stats), "databases\n")
        
        # Log results for each database
        for (db in names(database_specific_pathway_stats)) {
          db_result <- database_specific_pathway_stats[[db]]
          cat("  ", db, ": MAST=", db_result$mast_pathway_count, 
              ", CRISPRi=", db_result$crispri_pathway_count,
              ", overlap=", db_result$overlap_count,
              ", p=", db_result$fisher_p %||% "NA", "\n")
        }
      }
    }
    
    cluster_results[[cluster]] <- list(
      cluster = cluster,
      overlap_stats = overlap_stats,
      pathway_overlap_stats = pathway_overlap_stats,
      database_specific_pathway_stats = database_specific_pathway_stats,
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

#' Enhanced gene overlap significance with direction-aware analysis and experiment weighting (v0.2.6)
#'
#' This function implements the enhanced statistical framework combining:
#' - Direction-aware Fisher's exact tests (same vs opposite direction)
#' - Experiment weighting based on cell counts
#' - Biological expectation weighting (LRRK2 vs SNCA)
#'
#' @param mast_data Data frame with MAST DE results (must have avg_log2FC, p_val_adj columns)
#' @param crispri_experiments_data List of CRISPRi DE results by experiment
#' @param gene_name Gene name for biological context
#' @param cluster_id Cluster identifier
#' @param experiment_weights Experiment weights from calculate_experiment_weights()
#' @param background_genes Character vector of background genes (optional)
#' @param lfc_threshold Log2 fold change threshold (default 0.25)
#' @param p_threshold P-value threshold (default 0.05)
#' @param use_enhanced_analysis Logical, whether to use enhanced method (default TRUE)
#' @return List with enhanced overlap significance results
#' @export
calculate_enhanced_overlap_significance <- function(mast_data, crispri_experiments_data, 
                                                   gene_name, cluster_id, experiment_weights,
                                                   background_genes = NULL,
                                                   lfc_threshold = 0.25, p_threshold = 0.05,
                                                   use_enhanced_analysis = TRUE) {
  
  if (!use_enhanced_analysis) {
    # Fallback to legacy method for backward compatibility
    if (is.list(crispri_experiments_data) && length(crispri_experiments_data) > 0) {
      # Use first experiment for legacy compatibility
      first_exp_data <- crispri_experiments_data[[1]]
      mast_genes <- rownames(mast_data)[mast_data$p_val_adj < p_threshold & 
                                        abs(mast_data$avg_log2FC) > lfc_threshold & 
                                        !is.na(mast_data$p_val_adj)]
      
      # Handle CRISPRi gene extraction (experiment-specific columns)
      crispri_lfc_col <- names(first_exp_data)[grepl("^log2FC", names(first_exp_data))][1]
      if (is.na(crispri_lfc_col)) crispri_lfc_col <- "log2FC"
      
      crispri_genes <- rownames(first_exp_data)[first_exp_data$p_val_adj < p_threshold & 
                                               abs(first_exp_data[[crispri_lfc_col]]) > lfc_threshold & 
                                               !is.na(first_exp_data$p_val_adj)]
      
      return(calculate_gene_overlap_significance_proper(
        mast_genes, crispri_genes, 
        rownames(mast_data), rownames(first_exp_data)
      ))
    }
  }
  
  cat("[ENHANCED OVERLAP] Analyzing", gene_name, "in", cluster_id, "with direction-aware and weighted approach\n")
  
  # Determine biological direction expectation
  if (!exists("get_biological_direction_expectation")) {
    source("R/enhanced_direction_analysis.R")
  }
  direction_expectation <- get_biological_direction_expectation(gene_name)
  cat("[ENHANCED OVERLAP] Biological expectation for", gene_name, ":", direction_expectation, "\n")
  
  # Initialize results structure
  enhanced_results <- list(
    gene_name = gene_name,
    cluster_id = cluster_id,
    biological_expectation = direction_expectation,
    experiment_results = list(),
    weighted_meta_analysis = NULL,
    enhanced_statistics = NULL
  )
  
  # Analyze each experiment separately with direction awareness
  experiments <- names(crispri_experiments_data)
  experiment_direction_results <- list()
  
  for (exp in experiments) {
    exp_data <- crispri_experiments_data[[exp]]
    if (is.null(exp_data) || nrow(exp_data) == 0) {
      cat("[ENHANCED OVERLAP] Skipping", exp, "- no data available\n")
      next
    }
    
    cat("[ENHANCED OVERLAP] Analyzing experiment:", exp, "\n")
    
    # Perform direction-aware analysis for this experiment
    if (!exists("enhanced_direction_analysis")) {
      source("R/enhanced_direction_analysis.R")
    }
    exp_direction_result <- enhanced_direction_analysis(
      mast_data, exp_data, gene_name, background_genes, lfc_threshold, p_threshold
    )
    
    experiment_direction_results[[exp]] <- exp_direction_result
    enhanced_results$experiment_results[[exp]] <- exp_direction_result
  }
  
  # Weighted meta-analysis across experiments
  if (length(experiment_direction_results) > 0) {
    cat("[ENHANCED OVERLAP] Performing weighted meta-analysis across", length(experiment_direction_results), "experiments\n")
    
    # Extract weights for this cluster
    meta_analysis_results <- list()
    
    # Combine same-direction results across experiments
    same_direction_effects <- list()
    same_direction_pvalues <- list()
    same_direction_weights <- list()
    
    # Combine opposite-direction results across experiments  
    opposite_direction_effects <- list()
    opposite_direction_pvalues <- list()
    opposite_direction_weights <- list()
    
    for (exp in names(experiment_direction_results)) {
      weight_key <- paste0(exp, "_", cluster_id)
      if (weight_key %in% names(experiment_weights$weights)) {
        exp_weight <- experiment_weights$weights[[weight_key]]
        
        if (exp_weight > 0) {
          exp_result <- experiment_direction_results[[exp]]
          
          # Same direction
          if (!is.null(exp_result$same_direction) && !is.na(exp_result$same_direction$fisher_p)) {
            same_direction_effects[[exp]] <- exp_result$same_direction$overlap_count
            same_direction_pvalues[[exp]] <- exp_result$same_direction$fisher_p
            same_direction_weights[[exp]] <- exp_weight
          }
          
          # Opposite direction
          if (!is.null(exp_result$opposite_direction) && !is.na(exp_result$opposite_direction$fisher_p)) {
            opposite_direction_effects[[exp]] <- exp_result$opposite_direction$overlap_count
            opposite_direction_pvalues[[exp]] <- exp_result$opposite_direction$fisher_p
            opposite_direction_weights[[exp]] <- exp_weight
          }
        }
      }
    }
    
    # Weighted combination of same-direction results
    same_direction_meta <- combine_weighted_results(
      same_direction_effects, same_direction_pvalues, same_direction_weights, "same_direction"
    )
    
    # Weighted combination of opposite-direction results
    opposite_direction_meta <- combine_weighted_results(
      opposite_direction_effects, opposite_direction_pvalues, opposite_direction_weights, "opposite_direction"
    )
    
    # Apply biological weighting to combine same vs opposite direction results
    if (direction_expectation == "opposing") {
      primary_meta <- opposite_direction_meta
      secondary_meta <- same_direction_meta
      primary_weight <- 0.8
      secondary_weight <- 0.2
    } else if (direction_expectation == "same") {
      primary_meta <- same_direction_meta
      secondary_meta <- opposite_direction_meta
      primary_weight <- 0.8
      secondary_weight <- 0.2
    } else {
      # Mixed/unknown: equal weighting
      primary_meta <- same_direction_meta
      secondary_meta <- opposite_direction_meta
      primary_weight <- 0.5
      secondary_weight <- 0.5
    }
    
    # Final combined p-value using weighted Fisher's method
    final_combined_p <- combine_meta_pvalues(
      primary_meta$combined_p, secondary_meta$combined_p,
      primary_weight, secondary_weight
    )
    
    enhanced_results$weighted_meta_analysis <- list(
      same_direction_meta = same_direction_meta,
      opposite_direction_meta = opposite_direction_meta,
      primary_direction = ifelse(direction_expectation == "opposing", "opposite", "same"),
      primary_meta = primary_meta,
      secondary_meta = secondary_meta,
      final_combined_p = final_combined_p,
      biological_weighting = c(primary_weight, secondary_weight),
      experiments_included = length(experiment_direction_results)
    )
  }
  
  # Enhanced statistics summary
  enhanced_results$enhanced_statistics <- list(
    direction_expectation = direction_expectation,
    experiments_analyzed = length(experiment_direction_results),
    primary_experiment = experiment_weights$primary_experiment,
    weighting_method = experiment_weights$weighting_method,
    analysis_timestamp = Sys.time(),
    enhanced_method_version = "v0.2.6"
  )
  
  cat("[ENHANCED OVERLAP] Enhanced analysis completed for", gene_name, "\n")
  cat("[ENHANCED OVERLAP] Final combined p-value:", enhanced_results$weighted_meta_analysis$final_combined_p %||% "NA", "\n")
  
  return(enhanced_results)
}

#' Helper function to combine weighted results across experiments
#'
#' @param effects List of effect sizes by experiment
#' @param pvalues List of p-values by experiment  
#' @param weights List of weights by experiment
#' @param direction_type Character, "same_direction" or "opposite_direction"
#' @return List with combined weighted results
combine_weighted_results <- function(effects, pvalues, weights, direction_type) {
  
  if (length(pvalues) == 0) {
    return(list(
      combined_p = NA,
      weighted_effect = NA,
      experiments_included = 0,
      direction_type = direction_type,
      error = "No valid experiments"
    ))
  }
  
  # Convert to vectors
  effect_vec <- unlist(effects)
  pvalue_vec <- unlist(pvalues)
  weight_vec <- unlist(weights)
  
  # Weighted average effect size
  weighted_effect <- sum(effect_vec * weight_vec) / sum(weight_vec)
  
  # Weighted combination of p-values using Fisher's method
  chi_square_stats <- -2 * log(pvalue_vec)
  weighted_chi_square <- sum(chi_square_stats * weight_vec) / sum(weight_vec)
  
  # Convert back to p-value
  effective_df <- 2 * length(pvalue_vec)
  combined_p <- pchisq(weighted_chi_square * 2 * length(pvalue_vec), 
                      df = effective_df, lower.tail = FALSE)
  
  return(list(
    combined_p = combined_p,
    weighted_effect = weighted_effect,
    experiments_included = length(pvalue_vec),
    direction_type = direction_type,
    experiment_pvalues = pvalue_vec,
    experiment_weights = weight_vec,
    weighted_chi_square = weighted_chi_square
  ))
}

#' Helper function to combine meta-analysis p-values with biological weighting
#'
#' @param primary_p P-value from primary (expected) direction
#' @param secondary_p P-value from secondary (alternative) direction
#' @param primary_weight Weight for primary direction
#' @param secondary_weight Weight for secondary direction
#' @return Combined p-value
combine_meta_pvalues <- function(primary_p, secondary_p, primary_weight, secondary_weight) {
  
  if (is.na(primary_p) && is.na(secondary_p)) {
    return(NA)
  }
  
  if (is.na(primary_p)) {
    return(secondary_p)
  }
  
  if (is.na(secondary_p)) {
    return(primary_p)
  }
  
  # Weighted Fisher's method
  chi_square_primary <- -2 * log(primary_p)
  chi_square_secondary <- -2 * log(secondary_p)
  
  weighted_chi_square <- (chi_square_primary * primary_weight + 
                         chi_square_secondary * secondary_weight)
  
  # Convert back to p-value
  effective_df <- 2 * (primary_weight + secondary_weight)
  combined_p <- pchisq(weighted_chi_square, df = effective_df, lower.tail = FALSE)
  
  return(combined_p)
}

# Helper function for null coalescing
if (!exists("%||%")) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
}