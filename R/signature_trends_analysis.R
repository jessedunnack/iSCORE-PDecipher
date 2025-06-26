#' Signature Trends Analysis Functions
#'
#' Data-driven discovery of signature patterns, frequency analysis, and impact assessment
#' for cross-method signature discovery without manual curation bias.

#' Analyze signature trends across all discovered signatures
#'
#' @param signature_results Results from discover_top_signatures()
#' @param enrichment_data Original consolidated enrichment data
#' @param min_frequency Minimum frequency for inclusion (default: 2)
#' @param top_n Number of top results to return for each analysis (default: 50)
#' @return List with frequency analysis, impact analysis, and term patterns
#' @export
analyze_signature_trends <- function(signature_results, enrichment_data, 
                                   min_frequency = 2, top_n = 50) {
  
  cat("[TRENDS] Starting comprehensive signature trends analysis...\n")
  
  # Validate inputs with comprehensive error handling
  validation_result <- validate_trends_inputs(signature_results, enrichment_data)
  if (!validation_result$valid) {
    warning("Input validation failed: ", validation_result$message)
    return(create_empty_trends_result(validation_result$message))
  }
  
  tryCatch({
    # Extract all signatures for analysis
    all_signatures <- signature_results$all_signatures
    pan_cluster_signatures <- signature_results$pan_cluster_signatures
    
    cat("[TRENDS] Analyzing", nrow(all_signatures), "individual signatures and", 
        nrow(pan_cluster_signatures), "pan-cluster signatures\n")
    
    # 1. Signature Frequency Analysis
    cat("[TRENDS] Computing signature frequency patterns...\n")
    frequency_analysis <- compute_signature_frequency_analysis(
      all_signatures, pan_cluster_signatures, min_frequency, top_n
    )
    
    # 2. Impact Analysis
    cat("[TRENDS] Computing impact-based rankings...\n") 
    impact_analysis <- compute_signature_impact_analysis(
      all_signatures, pan_cluster_signatures, top_n
    )
    
    # 3. Term Pattern Discovery
    cat("[TRENDS] Discovering enrichment term patterns...\n")
    term_patterns <- discover_enrichment_term_patterns(
      all_signatures, enrichment_data, top_n
    )
    
    # 4. Create comprehensive summary
    cat("[TRENDS] Creating trends summary...\n")
    trends_summary <- create_trends_summary(
      frequency_analysis, impact_analysis, term_patterns
    )
    
    cat("[TRENDS] Analysis complete! Found", nrow(frequency_analysis$top_frequent_signatures), 
        "frequent signatures and", nrow(impact_analysis$top_impact_signatures), "high-impact signatures\n")
    
    return(list(
      frequency_analysis = frequency_analysis,
      impact_analysis = impact_analysis,
      term_patterns = term_patterns,
      trends_summary = trends_summary,
      analysis_metadata = list(
        total_signatures = nrow(all_signatures),
        pan_cluster_count = nrow(pan_cluster_signatures),
        min_frequency_threshold = min_frequency,
        top_n_limit = top_n,
        analysis_timestamp = Sys.time()
      )
    ))
    
  }, error = function(e) {
    error_msg <- paste("Signature trends analysis failed:", e$message)
    cat("[ERROR]", error_msg, "\n")
    return(create_empty_trends_result(error_msg))
  })
}

#' Validate inputs for trends analysis with comprehensive checking
#'
#' @param signature_results Signature discovery results
#' @param enrichment_data Enrichment data
#' @return List with validation status and message
validate_trends_inputs <- function(signature_results, enrichment_data) {
  
  # Check signature_results structure
  if (is.null(signature_results)) {
    return(list(valid = FALSE, message = "signature_results is NULL"))
  }
  
  required_components <- c("all_signatures", "pan_cluster_signatures")
  missing_components <- setdiff(required_components, names(signature_results))
  if (length(missing_components) > 0) {
    return(list(valid = FALSE, 
               message = paste("Missing components:", paste(missing_components, collapse = ", "))))
  }
  
  # Check data availability
  if (is.null(signature_results$all_signatures) || nrow(signature_results$all_signatures) == 0) {
    return(list(valid = FALSE, message = "No individual signatures available"))
  }
  
  # Check enrichment data
  if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
    return(list(valid = FALSE, message = "No enrichment data available"))
  }
  
  # Check required columns in signatures
  required_sig_cols <- c("gene_pair", "signature_strength")
  missing_sig_cols <- setdiff(required_sig_cols, colnames(signature_results$all_signatures))
  if (length(missing_sig_cols) > 0) {
    return(list(valid = FALSE, 
               message = paste("Missing signature columns:", paste(missing_sig_cols, collapse = ", "))))
  }
  
  return(list(valid = TRUE, message = "Validation passed"))
}

#' Compute signature frequency analysis
#'
#' @param all_signatures Individual signatures data
#' @param pan_cluster_signatures Pan-cluster signatures data  
#' @param min_frequency Minimum frequency threshold
#' @param top_n Number of top results to return
#' @return List with frequency analysis results
compute_signature_frequency_analysis <- function(all_signatures, pan_cluster_signatures, 
                                               min_frequency, top_n) {
  
  tryCatch({
    # Calculate gene pair frequencies
    gene_pair_frequency <- table(all_signatures$gene_pair)
    gene_pair_frequency <- gene_pair_frequency[gene_pair_frequency >= min_frequency]
    gene_pair_frequency <- sort(gene_pair_frequency, decreasing = TRUE)
    
    # Create frequency summary
    frequency_summary <- data.frame(
      gene_pair = names(gene_pair_frequency),
      frequency_count = as.numeric(gene_pair_frequency),
      stringsAsFactors = FALSE
    )
    
    # Add signature strength statistics for frequent pairs
    frequency_summary$mean_signature_strength <- sapply(frequency_summary$gene_pair, function(pair) {
      pair_sigs <- all_signatures[all_signatures$gene_pair == pair, ]
      mean(pair_sigs$signature_strength, na.rm = TRUE)
    })
    
    frequency_summary$max_signature_strength <- sapply(frequency_summary$gene_pair, function(pair) {
      pair_sigs <- all_signatures[all_signatures$gene_pair == pair, ]
      max(pair_sigs$signature_strength, na.rm = TRUE)
    })
    
    # Calculate frequency score (normalized)
    frequency_summary$frequency_score <- frequency_summary$frequency_count / max(frequency_summary$frequency_count)
    
    # Add breadth analysis (how many clusters each pair appears in)
    frequency_summary$cluster_breadth <- sapply(frequency_summary$gene_pair, function(pair) {
      pair_sigs <- all_signatures[all_signatures$gene_pair == pair, ]
      if ("cluster" %in% colnames(pair_sigs)) {
        length(unique(pair_sigs$cluster))
      } else {
        NA
      }
    })
    
    # Create top frequent signatures list
    top_frequent_signatures <- head(frequency_summary, top_n)
    
    # Calculate frequency distribution
    frequency_distribution <- data.frame(
      frequency_bin = names(table(frequency_summary$frequency_count)),
      signature_count = as.numeric(table(frequency_summary$frequency_count)),
      stringsAsFactors = FALSE
    )
    
    return(list(
      top_frequent_signatures = top_frequent_signatures,
      frequency_distribution = frequency_distribution,
      total_unique_pairs = length(gene_pair_frequency),
      frequent_pairs_count = nrow(frequency_summary)
    ))
    
  }, error = function(e) {
    cat("[ERROR] Frequency analysis failed:", e$message, "\n")
    return(list(
      top_frequent_signatures = data.frame(),
      frequency_distribution = data.frame(),
      total_unique_pairs = 0,
      frequent_pairs_count = 0
    ))
  })
}

#' Compute signature impact analysis based on statistical strength
#'
#' @param all_signatures Individual signatures data
#' @param pan_cluster_signatures Pan-cluster signatures data
#' @param top_n Number of top results to return
#' @return List with impact analysis results
compute_signature_impact_analysis <- function(all_signatures, pan_cluster_signatures, top_n) {
  
  tryCatch({
    # Rank individual signatures by strength
    individual_ranked <- all_signatures[order(all_signatures$signature_strength, decreasing = TRUE), ]
    top_individual <- head(individual_ranked, top_n)
    
    # Add impact score (normalized signature strength)
    if (nrow(individual_ranked) > 0) {
      max_strength <- max(individual_ranked$signature_strength, na.rm = TRUE)
      top_individual$impact_score <- top_individual$signature_strength / max_strength
    } else {
      top_individual$impact_score <- numeric(0)
    }
    
    # Pan-cluster impact analysis
    top_pan_cluster <- data.frame()
    if (!is.null(pan_cluster_signatures) && nrow(pan_cluster_signatures) > 0) {
      if ("mean_signature_strength" %in% colnames(pan_cluster_signatures)) {
        pan_cluster_ranked <- pan_cluster_signatures[order(pan_cluster_signatures$mean_signature_strength, 
                                                          decreasing = TRUE), ]
        top_pan_cluster <- head(pan_cluster_ranked, top_n)
        
        # Add impact score for pan-cluster
        max_pan_strength <- max(pan_cluster_ranked$mean_signature_strength, na.rm = TRUE)
        top_pan_cluster$impact_score <- top_pan_cluster$mean_signature_strength / max_pan_strength
      }
    }
    
    # Calculate impact distribution
    strength_breaks <- seq(0, max(all_signatures$signature_strength, na.rm = TRUE), length.out = 10)
    impact_distribution <- data.frame(
      strength_bin = head(strength_breaks, -1),
      signature_count = as.numeric(table(cut(all_signatures$signature_strength, breaks = strength_breaks))),
      stringsAsFactors = FALSE
    )
    
    return(list(
      top_impact_signatures = top_individual,
      top_pan_cluster_impact = top_pan_cluster,
      impact_distribution = impact_distribution,
      max_individual_strength = max(all_signatures$signature_strength, na.rm = TRUE),
      mean_individual_strength = mean(all_signatures$signature_strength, na.rm = TRUE)
    ))
    
  }, error = function(e) {
    cat("[ERROR] Impact analysis failed:", e$message, "\n")
    return(list(
      top_impact_signatures = data.frame(),
      top_pan_cluster_impact = data.frame(), 
      impact_distribution = data.frame(),
      max_individual_strength = 0,
      mean_individual_strength = 0
    ))
  })
}

#' Discover enrichment term patterns across signatures
#'
#' @param all_signatures Individual signatures data
#' @param enrichment_data Original enrichment data
#' @param top_n Number of top terms to return
#' @return List with term pattern analysis
discover_enrichment_term_patterns <- function(all_signatures, enrichment_data, top_n) {
  
  tryCatch({
    # This is a simplified pattern discovery - in full implementation,
    # this would extract actual enrichment terms for each signature
    # and find the most commonly occurring terms across signatures
    
    # For now, create a framework for term frequency analysis
    if (!"Description" %in% colnames(enrichment_data)) {
      cat("[WARNING] No Description column in enrichment data for term analysis\n")
      return(create_empty_term_patterns())
    }
    
    # Sample implementation - would need to be expanded based on
    # actual signature-term mapping
    term_frequency <- table(enrichment_data$Description)
    term_frequency <- sort(term_frequency, decreasing = TRUE)
    
    # Create top terms summary
    top_terms <- data.frame(
      term = names(head(term_frequency, top_n)),
      frequency = as.numeric(head(term_frequency, top_n)),
      stringsAsFactors = FALSE
    )
    
    # Add term pattern categories (simplified)
    top_terms$pattern_category <- sapply(top_terms$term, classify_term_pattern)
    
    return(list(
      top_frequent_terms = top_terms,
      total_unique_terms = length(term_frequency),
      pattern_categories = table(top_terms$pattern_category)
    ))
    
  }, error = function(e) {
    cat("[ERROR] Term pattern discovery failed:", e$message, "\n")
    return(create_empty_term_patterns())
  })
}

#' Classify term into pattern categories
#'
#' @param term Term description
#' @return Pattern category
classify_term_pattern <- function(term) {
  if (is.na(term) || !is.character(term)) return("Unknown")
  
  term_lower <- tolower(term)
  
  if (grepl("mitochondri|respiratory|oxidative phosphorylation", term_lower)) {
    return("Mitochondrial")
  } else if (grepl("autophagy|lysosome|phagosome", term_lower)) {
    return("Autophagy") 
  } else if (grepl("protein fold|unfolded protein|proteasome", term_lower)) {
    return("Protein Quality")
  } else if (grepl("dopamine|neurotransmitter", term_lower)) {
    return("Neurotransmission")
  } else if (grepl("synapse|synaptic", term_lower)) {
    return("Synaptic")
  } else {
    return("Other")
  }
}

#' Create empty term patterns result
#'
#' @return Empty term patterns structure
create_empty_term_patterns <- function() {
  list(
    top_frequent_terms = data.frame(),
    total_unique_terms = 0,
    pattern_categories = table(character(0))
  )
}

#' Create comprehensive trends summary
#'
#' @param frequency_analysis Frequency analysis results
#' @param impact_analysis Impact analysis results  
#' @param term_patterns Term pattern results
#' @return Summary statistics
create_trends_summary <- function(frequency_analysis, impact_analysis, term_patterns) {
  
  list(
    most_frequent_signature = if(nrow(frequency_analysis$top_frequent_signatures) > 0) {
      frequency_analysis$top_frequent_signatures$gene_pair[1]
    } else "None",
    
    highest_impact_signature = if(nrow(impact_analysis$top_impact_signatures) > 0) {
      impact_analysis$top_impact_signatures$gene_pair[1]
    } else "None",
    
    most_common_term = if(nrow(term_patterns$top_frequent_terms) > 0) {
      term_patterns$top_frequent_terms$term[1]
    } else "None",
    
    total_patterns_analyzed = frequency_analysis$total_unique_pairs,
    average_signature_strength = impact_analysis$mean_individual_strength,
    max_signature_strength = impact_analysis$max_individual_strength
  )
}

#' Create empty trends result for error cases
#'
#' @param error_message Error message to include
#' @return Empty trends result structure
create_empty_trends_result <- function(error_message) {
  list(
    frequency_analysis = list(
      top_frequent_signatures = data.frame(),
      frequency_distribution = data.frame(),
      total_unique_pairs = 0,
      frequent_pairs_count = 0
    ),
    impact_analysis = list(
      top_impact_signatures = data.frame(),
      top_pan_cluster_impact = data.frame(),
      impact_distribution = data.frame(),
      max_individual_strength = 0,
      mean_individual_strength = 0
    ),
    term_patterns = create_empty_term_patterns(),
    trends_summary = list(
      most_frequent_signature = "Error",
      highest_impact_signature = "Error", 
      most_common_term = "Error",
      total_patterns_analyzed = 0,
      average_signature_strength = 0,
      max_signature_strength = 0
    ),
    analysis_metadata = list(
      error_message = error_message,
      analysis_timestamp = Sys.time()
    )
  )
}