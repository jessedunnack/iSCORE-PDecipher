#' Signature Analysis Functions for Cross-Method Comparison
#'
#' This module implements statistical methods for identifying shared signatures
#' between MAST (mutations) and CRISPRi (knockdowns) experiments.

#' Calculate gene overlap significance between two methods
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
    neither <- length(background_genes) - length(mast_genes) - crispri_only
    
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
  if (!is.null(overlap_stats) && !is.na(overlap_stats$fisher_p)) {
    # Convert p-value to score (higher score for lower p-value)
    overlap_score <- -log10(max(overlap_stats$fisher_p, 1e-10)) * overlap_stats$jaccard_index
  } else if (!is.null(overlap_stats)) {
    # Use Jaccard index alone if no p-value
    overlap_score <- overlap_stats$jaccard_index * 5  # Scale to be comparable
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

#' Perform comprehensive signature analysis for a gene pair
#'
#' @param gene_pair List with mast_gene and crispri_gene names
#' @param enrichment_data Consolidated enrichment data
#' @param clusters Character vector of clusters to analyze (NULL for all)
#' @param include_pathways Logical, whether to include pathway analysis
#' @param progress_callback Function to call for progress updates (optional)
#' @return List with comprehensive signature analysis results
#' @export
analyze_gene_pair_signatures <- function(gene_pair, enrichment_data, clusters = NULL,
                                        include_pathways = TRUE, progress_callback = NULL) {
  
  if (!is.null(progress_callback)) {
    progress_callback(paste("Analyzing", gene_pair$mast_gene, "vs", gene_pair$crispri_gene))
  }
  
  # Filter data for this gene pair
  mast_data <- enrichment_data[enrichment_data$method == "MAST" & 
                              enrichment_data$mutation_perturbation == gene_pair$mast_gene, ]
  crispri_data <- enrichment_data[enrichment_data$method == "MixScale" & 
                                 enrichment_data$mutation_perturbation == gene_pair$crispri_gene, ]
  
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
  
  for (cluster in all_clusters) {
    cluster_mast <- mast_data[mast_data$cluster == cluster, ]
    cluster_crispri <- crispri_data[crispri_data$cluster == cluster, ]
    
    if (nrow(cluster_mast) == 0 || nrow(cluster_crispri) == 0) {
      cluster_results[[cluster]] <- list(error = "Missing data for one method")
      next
    }
    
    # Gene overlap analysis (using significant genes from enrichment terms)
    mast_genes <- unique(unlist(strsplit(cluster_mast$geneID, "/")))
    crispri_genes <- unique(unlist(strsplit(cluster_crispri$geneID, "/")))
    
    # Remove NA and empty strings
    mast_genes <- mast_genes[!is.na(mast_genes) & mast_genes != ""]
    crispri_genes <- crispri_genes[!is.na(crispri_genes) & crispri_genes != ""]
    
    if (length(mast_genes) > 0 && length(crispri_genes) > 0) {
      background_genes <- unique(c(mast_genes, crispri_genes))
      overlap_stats <- calculate_gene_overlap_significance(mast_genes, crispri_genes, background_genes)
    } else {
      overlap_stats <- list(error = "No valid gene lists extracted")
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