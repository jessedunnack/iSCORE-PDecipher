#' Gene Association Lookup Functions
#'
#' Functions to load and lookup gene associations from the consolidated dataset
#'
#' @author Claude Code Assistant
#' @date 2025-01-13

#' Global variable to store loaded gene associations
.gene_associations <- NULL

#' Load gene associations data
#'
#' @param force_reload Logical, whether to force reload even if already loaded
#' @return Data frame with gene associations or NULL if not available
#' @export
load_gene_associations <- function(force_reload = FALSE) {
  if (is.null(.gene_associations) || force_reload) {
    # Try to load from package extdata
    extdata_path <- system.file("extdata", "gene_term_associations.rds", package = "iSCORE.PDecipher")
    
    if (file.exists(extdata_path)) {
      message("Loading gene associations from package data...")
      .gene_associations <<- readRDS(extdata_path)
      message("Loaded ", nrow(.gene_associations), " gene-term associations")
    } else {
      # Fallback: try local file for development
      local_path <- "inst/extdata/gene_term_associations.rds"
      if (file.exists(local_path)) {
        message("Loading gene associations from local development file...")
        .gene_associations <<- readRDS(local_path)
        message("Loaded ", nrow(.gene_associations), " gene-term associations")
      } else {
        warning("Gene associations file not found. Gene display functionality will be limited.")
        .gene_associations <<- NULL
      }
    }
  }
  
  return(.gene_associations)
}

#' Check if gene associations are available
#'
#' @return Logical indicating if gene associations are loaded
#' @export
gene_associations_available <- function() {
  return(!is.null(.gene_associations) && nrow(.gene_associations) > 0)
}

#' Get genes for a specific enrichment term
#'
#' @param term_id Character, the term ID to lookup
#' @param analysis_type Character, analysis type (MAST, MixScale)
#' @param gene Character, gene name
#' @param cluster Character, cluster name
#' @param enrichment_type Character, enrichment type (GO_BP, KEGG, etc.)
#' @param direction Character, direction (UP, DOWN, ALL, etc.)
#' @param experiment Character, experiment name (default for most)
#' @return List with genes vector and metadata, or NULL if not found
#' @export
get_genes_for_term <- function(term_id, analysis_type, gene, cluster, 
                              enrichment_type, direction, experiment = "default") {
  
  # Ensure data is loaded
  if (!gene_associations_available()) {
    load_gene_associations()
  }
  
  if (!gene_associations_available()) {
    return(list(genes = NULL, error = "Gene associations not available"))
  }
  
  # Create composite key for lookup
  composite_key <- paste(analysis_type, gene, cluster, enrichment_type, 
                        direction, experiment, term_id, sep = "|")
  
  # Find matching association
  match_idx <- which(.gene_associations$composite_key == composite_key)
  
  if (length(match_idx) == 0) {
    # Try without experiment specification (fallback)
    composite_key_fallback <- paste(analysis_type, gene, cluster, enrichment_type, 
                                   direction, "default", term_id, sep = "|")
    match_idx <- which(.gene_associations$composite_key == composite_key_fallback)
  }
  
  if (length(match_idx) == 0) {
    return(list(genes = NULL, error = "Term not found in associations"))
  }
  
  # Get the association
  association <- .gene_associations[match_idx[1], ]
  
  # Split gene string into vector
  genes <- unlist(strsplit(association$associated_genes, "/"))
  
  return(list(
    genes = genes,
    gene_count = association$gene_count,
    term_id = association$term_id,
    description = association$description,
    error = NULL
  ))
}

#' Get all genes for multiple terms (bulk lookup)
#'
#' @param term_ids Character vector of term IDs
#' @param analysis_type Character, analysis type 
#' @param gene Character, gene name
#' @param cluster Character, cluster name
#' @param enrichment_type Character, enrichment type
#' @param direction Character, direction
#' @param experiment Character, experiment name
#' @return Named list where names are term_ids and values are gene vectors
#' @export
get_genes_for_terms <- function(term_ids, analysis_type, gene, cluster,
                               enrichment_type, direction, experiment = "default") {
  
  results <- list()
  
  for (term_id in term_ids) {
    result <- get_genes_for_term(term_id, analysis_type, gene, cluster,
                                enrichment_type, direction, experiment)
    if (!is.null(result$genes)) {
      results[[term_id]] <- result$genes
    }
  }
  
  return(results)
}

#' Search gene associations by term description
#'
#' @param pattern Character, pattern to search in descriptions
#' @param analysis_type Character, optional filter by analysis type
#' @param gene Character, optional filter by gene
#' @param enrichment_type Character, optional filter by enrichment type
#' @return Data frame with matching associations
#' @export
search_gene_associations <- function(pattern, analysis_type = NULL, gene = NULL,
                                   enrichment_type = NULL) {
  
  if (!gene_associations_available()) {
    load_gene_associations()
  }
  
  if (!gene_associations_available()) {
    return(data.frame())
  }
  
  # Start with all associations
  results <- .gene_associations
  
  # Apply filters
  if (!is.null(analysis_type)) {
    results <- results[results$analysis_type == analysis_type, ]
  }
  
  if (!is.null(gene)) {
    results <- results[results$gene == gene, ]
  }
  
  if (!is.null(enrichment_type)) {
    results <- results[results$enrichment_type == enrichment_type, ]
  }
  
  # Search description
  if (!is.null(pattern) && pattern != "") {
    matches <- grepl(pattern, results$description, ignore.case = TRUE)
    results <- results[matches, ]
  }
  
  return(results)
}

#' Get association statistics
#'
#' @return List with summary statistics
#' @export
get_association_stats <- function() {
  if (!gene_associations_available()) {
    load_gene_associations()
  }
  
  if (!gene_associations_available()) {
    return(list(error = "Gene associations not available"))
  }
  
  stats <- list(
    total_associations = nrow(.gene_associations),
    unique_terms = length(unique(.gene_associations$term_id)),
    unique_genes = length(unique(.gene_associations$gene)),
    analysis_types = sort(unique(.gene_associations$analysis_type)),
    enrichment_types = sort(unique(.gene_associations$enrichment_type)),
    directions = sort(unique(.gene_associations$direction)),
    clusters = sort(unique(.gene_associations$cluster))
  )
  
  return(stats)
}