#' Manuscript-Focused Signature Discovery Pipeline
#'
#' This module implements rapid signature discovery prioritized for manuscript development,
#' focusing on the strongest shared signatures between MAST and CRISPRi methods.

# Required packages
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr package is required for signature analysis")
}
library(dplyr)

#' Discover top signatures for manuscript prioritization
#'
#' @param enrichment_data Consolidated enrichment data
#' @param top_n Integer, number of top signatures to return
#' @param min_cluster_breadth Integer, minimum clusters showing signature (for pan-cluster analysis)
#' @param combine_variants Logical, whether to combine gene variants
#' @param progress_callback Function for progress updates
#' @return List with ranked signature discoveries
#' @export
discover_top_signatures <- function(enrichment_data, top_n = 10, min_cluster_breadth = 8,
                                   combine_variants = TRUE, progress_callback = NULL) {
  
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
      include_pathways = TRUE,
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
      cluster_specific_count = length(cluster_specific_signatures)
    )
  )
  
  return(results)
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
          
          # Calculate composite score
          composite_score <- calculate_composite_signature_score(
            overlap_stats = overlap_stats,
            correlation_stats = NULL,  # Not available from enrichment data directly
            direction_stats = NULL,    # Would need DE data
            pathway_overlap_stats = pathway_stats
          )
          
          # Create signature entry
          signature_entry <- data.frame(
            gene_pair = gene_pair_name,
            mast_gene = analysis$gene_pair$mast_gene,
            crispri_gene = analysis$gene_pair$crispri_gene,
            cluster = cluster,
            signature_strength = composite_score$composite_score,
            gene_overlap_count = overlap_stats$overlap_count %||% 0,
            gene_fisher_p = overlap_stats$fisher_p %||% NA,
            gene_jaccard = overlap_stats$jaccard_index %||% 0,
            pathway_overlap_count = if(!is.null(pathway_stats)) pathway_stats$overlap_count %||% 0 else 0,
            pathway_fisher_p = if(!is.null(pathway_stats)) pathway_stats$fisher_p %||% NA else NA,
            mast_term_count = cluster_result$mast_term_count %||% 0,
            crispri_term_count = cluster_result$crispri_term_count %||% 0,
            stringsAsFactors = FALSE
          )
          
          signature_rankings <- rbind(signature_rankings, signature_entry)
        }
      }
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
identify_cluster_specific_signatures <- function(signature_rankings, min_strength_threshold = 5) {
  
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