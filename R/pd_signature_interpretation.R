#' Parkinson's Disease-Focused Signature Interpretation
#'
#' This module provides biological interpretation of signature analysis results,
#' focusing on PD-relevant pathways and creating manuscript-ready insights.

#' Analyze signatures for PD-relevant biological processes
#'
#' @param signature_results Results from discover_top_signatures
#' @param enrichment_data Original consolidated enrichment data
#' @param focus_on_pan_cluster Logical, prioritize pan-cluster signatures
#' @return Enhanced signature analysis with PD biological context
#' @export
analyze_pd_signatures <- function(signature_results, enrichment_data, focus_on_pan_cluster = TRUE) {
  
  cat("[PD ANALYSIS] Starting PD-focused signature interpretation...\n")
  
  # Get PD-relevant pathway terms
  pd_pathways <- get_pd_relevant_pathways()
  
  # Choose which signatures to analyze
  if (focus_on_pan_cluster && nrow(signature_results$pan_cluster_signatures) > 0) {
    target_signatures <- signature_results$pan_cluster_signatures
    analysis_type <- "pan-cluster"
    cat("[PD ANALYSIS] Focusing on", nrow(target_signatures), "pan-cluster signatures\n")
  } else {
    target_signatures <- head(signature_results$all_signatures, 20)  # Top 20 overall
    analysis_type <- "top-ranked"
    cat("[PD ANALYSIS] Focusing on top", nrow(target_signatures), "ranked signatures\n")
  }
  
  # Enhance each signature with PD biological context
  enhanced_signatures <- list()
  
  for (i in seq_len(nrow(target_signatures))) {
    signature <- target_signatures[i, ]
    
    cat("[PD ANALYSIS] Analyzing signature", i, ":", signature$gene_pair, "\n")
    
    # Get enrichment terms for this gene pair
    signature_context <- extract_signature_biological_context(
      signature, enrichment_data, pd_pathways
    )
    
    enhanced_signatures[[i]] <- signature_context
  }
  
  # Create comprehensive summary
  pd_summary <- create_pd_signature_summary(enhanced_signatures, analysis_type)
  
  return(list(
    enhanced_signatures = enhanced_signatures,
    pd_summary = pd_summary,
    analysis_type = analysis_type,
    signature_count = length(enhanced_signatures)
  ))
}

#' Extract biological context for a single signature
#'
#' @param signature Single row from signature rankings
#' @param enrichment_data Original enrichment data
#' @param pd_pathways PD-relevant pathway terms
#' @return List with biological context information
extract_signature_biological_context <- function(signature, enrichment_data, pd_pathways) {
  
  # Get enrichment terms for both methods in this signature
  mast_terms <- get_signature_enrichment_terms(
    signature, enrichment_data, method = "MAST"
  )
  
  crispri_terms <- get_signature_enrichment_terms(
    signature, enrichment_data, method = "MixScale"
  )
  
  # Identify PD-relevant terms
  mast_pd_terms <- filter_pd_relevant_terms(mast_terms, pd_pathways)
  crispri_pd_terms <- filter_pd_relevant_terms(crispri_terms, pd_pathways)
  
  # Find shared PD-relevant pathways
  shared_pd_pathways <- find_shared_pathways(mast_pd_terms, crispri_pd_terms)
  
  # Categorize biological processes
  biological_categories <- categorize_biological_processes(shared_pd_pathways, pd_pathways)
  
  # Generate interpretation
  interpretation <- generate_signature_interpretation(
    signature, biological_categories, shared_pd_pathways
  )
  
  return(list(
    signature = signature,
    mast_pd_terms = mast_pd_terms,
    crispri_pd_terms = crispri_pd_terms,
    shared_pd_pathways = shared_pd_pathways,
    biological_categories = biological_categories,
    interpretation = interpretation,
    pd_relevance_score = calculate_pd_relevance_score(biological_categories)
  ))
}

#' Get enrichment terms for a signature and method
#'
#' @param signature Signature row
#' @param enrichment_data Full enrichment data
#' @param method "MAST" or "MixScale"
#' @return Filtered enrichment terms
get_signature_enrichment_terms <- function(signature, enrichment_data, method) {
  
  # Handle gene name mapping
  if (method == "MAST") {
    if (signature$mast_gene == "SNCA_combined") {
      gene_filter <- c("SNCA_A30P", "SNCA_A53T")
    } else if (signature$mast_gene == "VPS13C_combined") {
      gene_filter <- c("VPS13C_A444P", "VPS13C_W395C")
    } else {
      gene_filter <- signature$mast_gene
    }
  } else {
    gene_filter <- signature$crispri_gene
  }
  
  # Filter enrichment data
  terms <- enrichment_data[
    enrichment_data$method == method &
    enrichment_data$mutation_perturbation %in% gene_filter &
    enrichment_data$cluster == signature$cluster,
  ]
  
  return(terms)
}

#' Filter terms for PD relevance
#'
#' @param terms Enrichment terms data frame
#' @param pd_pathways PD-relevant pathway keywords
#' @return PD-relevant terms with relevance scores
filter_pd_relevant_terms <- function(terms, pd_pathways) {
  
  if (nrow(terms) == 0) {
    return(data.frame())
  }
  
  # Create pattern for PD pathway matching
  pd_pattern <- paste(pd_pathways, collapse = "|")
  
  # Find PD-relevant terms
  pd_matches <- grepl(pd_pattern, terms$Description, ignore.case = TRUE)
  pd_terms <- terms[pd_matches, ]
  
  # Add PD relevance scoring
  if (nrow(pd_terms) > 0) {
    pd_terms$pd_relevance_score <- sapply(pd_terms$Description, function(desc) {
      sum(sapply(pd_pathways, function(pathway) {
        grepl(pathway, desc, ignore.case = TRUE)
      }))
    })
    
    # Sort by relevance and significance
    pd_terms <- pd_terms[order(-pd_terms$pd_relevance_score, pd_terms$p.adjust), ]
  }
  
  return(pd_terms)
}

#' Find shared pathways between methods
#'
#' @param mast_terms MAST PD-relevant terms
#' @param crispri_terms CRISPRi PD-relevant terms
#' @return Data frame of shared pathway information
find_shared_pathways <- function(mast_terms, crispri_terms) {
  
  if (nrow(mast_terms) == 0 || nrow(crispri_terms) == 0) {
    return(data.frame())
  }
  
  # Find overlapping pathway descriptions (exact matches)
  shared_exact <- intersect(mast_terms$Description, crispri_terms$Description)
  
  # Find conceptually similar pathways (partial matches)
  shared_partial <- find_conceptually_similar_pathways(mast_terms, crispri_terms)
  
  # Combine results
  shared_pathways <- data.frame(
    pathway = c(shared_exact, shared_partial$pathway),
    match_type = c(rep("exact", length(shared_exact)), 
                   rep("conceptual", nrow(shared_partial))),
    mast_pval = NA,
    crispri_pval = NA,
    stringsAsFactors = FALSE
  )
  
  # Add p-values for shared pathways
  for (i in seq_len(nrow(shared_pathways))) {
    pathway <- shared_pathways$pathway[i]
    
    mast_match <- mast_terms[mast_terms$Description == pathway, ]
    crispri_match <- crispri_terms[crispri_terms$Description == pathway, ]
    
    if (nrow(mast_match) > 0) {
      shared_pathways$mast_pval[i] <- min(mast_match$p.adjust)
    }
    if (nrow(crispri_match) > 0) {
      shared_pathways$crispri_pval[i] <- min(crispri_match$p.adjust)
    }
  }
  
  return(shared_pathways)
}

#' Find conceptually similar pathways between methods
#'
#' @param mast_terms MAST terms
#' @param crispri_terms CRISPRi terms
#' @return Data frame of conceptually similar pathways
find_conceptually_similar_pathways <- function(mast_terms, crispri_terms) {
  
  # Define pathway similarity keywords
  similarity_groups <- list(
    mitochondrial = c("mitochondri", "respiratory", "oxidative phosphorylation", "electron transport"),
    autophagy = c("autophagy", "lysosome", "phagosome", "vacuole"),
    protein_quality = c("protein folding", "unfolded protein", "proteasome", "ubiquitin", "chaperone"),
    dopamine = c("dopamine", "dopaminergic", "catecholamine", "neurotransmitter"),
    synaptic = c("synaptic", "synapse", "neurotransmitter", "vesicle"),
    oxidative_stress = c("oxidative stress", "reactive oxygen", "antioxidant")
  )
  
  similar_pathways <- data.frame()
  
  for (group_name in names(similarity_groups)) {
    keywords <- similarity_groups[[group_name]]
    pattern <- paste(keywords, collapse = "|")
    
    mast_group <- mast_terms[grepl(pattern, mast_terms$Description, ignore.case = TRUE), ]
    crispri_group <- crispri_terms[grepl(pattern, crispri_terms$Description, ignore.case = TRUE), ]
    
    if (nrow(mast_group) > 0 && nrow(crispri_group) > 0) {
      similar_pathways <- rbind(similar_pathways, data.frame(
        pathway = paste(group_name, "pathways"),
        mast_count = nrow(mast_group),
        crispri_count = nrow(crispri_group),
        stringsAsFactors = FALSE
      ))
    }
  }
  
  return(similar_pathways)
}

#' Categorize biological processes
#'
#' @param shared_pathways Shared pathway data
#' @param pd_pathways PD pathway keywords
#' @return List of biological process categories
categorize_biological_processes <- function(shared_pathways, pd_pathways) {
  
  if (nrow(shared_pathways) == 0) {
    return(list(
      mitochondrial = 0,
      protein_quality = 0,
      autophagy = 0,
      dopamine = 0,
      synaptic = 0,
      oxidative_stress = 0,
      neuronal = 0
    ))
  }
  
  # Count pathways in each category
  categories <- list(
    mitochondrial = sum(grepl("mitochondri|respiratory|oxidative phosphorylation|electron transport", 
                             shared_pathways$pathway, ignore.case = TRUE)),
    protein_quality = sum(grepl("protein folding|unfolded protein|proteasome|ubiquitin|chaperone", 
                                shared_pathways$pathway, ignore.case = TRUE)),
    autophagy = sum(grepl("autophagy|lysosome|phagosome|vacuole", 
                         shared_pathways$pathway, ignore.case = TRUE)),
    dopamine = sum(grepl("dopamine|dopaminergic|catecholamine", 
                        shared_pathways$pathway, ignore.case = TRUE)),
    synaptic = sum(grepl("synaptic|synapse|neurotransmitter|vesicle", 
                        shared_pathways$pathway, ignore.case = TRUE)),
    oxidative_stress = sum(grepl("oxidative stress|reactive oxygen|antioxidant", 
                                shared_pathways$pathway, ignore.case = TRUE)),
    neuronal = sum(grepl("axon|dendrite|neuron|neurogenesis", 
                        shared_pathways$pathway, ignore.case = TRUE))
  )
  
  return(categories)
}

#' Generate interpretation text for a signature
#'
#' @param signature Signature information
#' @param biological_categories Categorized biological processes
#' @param shared_pathways Shared pathway data
#' @return Character string with biological interpretation
generate_signature_interpretation <- function(signature, biological_categories, shared_pathways) {
  
  # Start with gene pair information
  interpretation <- paste0(
    "Gene Pair: ", signature$mast_gene, " (mutation) vs ", signature$crispri_gene, " (knockdown)\n",
    "Cluster: ", signature$cluster, " | Signature Strength: ", round(signature$signature_strength, 2), "\n\n"
  )
  
  # Identify dominant biological themes
  top_categories <- names(biological_categories)[biological_categories > 0]
  top_categories <- top_categories[order(biological_categories[top_categories], decreasing = TRUE)]
  
  if (length(top_categories) == 0) {
    interpretation <- paste0(interpretation, 
      "Limited overlap in PD-relevant pathways between mutation and knockdown approaches.")
    return(interpretation)
  }
  
  # Generate biological insights
  interpretation <- paste0(interpretation, "Key Biological Themes:\n")
  
  for (category in head(top_categories, 3)) {  # Top 3 categories
    count <- biological_categories[[category]]
    category_text <- generate_category_interpretation(category, count)
    interpretation <- paste0(interpretation, "• ", category_text, "\n")
  }
  
  # Add pathway examples if available
  if (nrow(shared_pathways) > 0) {
    interpretation <- paste0(interpretation, "\nShared Pathway Examples:\n")
    top_pathways <- head(shared_pathways[order(shared_pathways$mast_pval, shared_pathways$crispri_pval), ], 3)
    
    for (i in seq_len(nrow(top_pathways))) {
      pathway <- top_pathways$pathway[i]
      interpretation <- paste0(interpretation, "• ", pathway, "\n")
    }
  }
  
  # Add biological significance
  interpretation <- paste0(interpretation, "\nBiological Significance:\n")
  interpretation <- paste0(interpretation, generate_biological_significance(top_categories, signature))
  
  return(interpretation)
}

#' Generate category-specific interpretation
#'
#' @param category Biological category name
#' @param count Number of pathways in category
#' @return Interpretation text for the category
generate_category_interpretation <- function(category, count) {
  
  category_meanings <- list(
    mitochondrial = paste0("Mitochondrial dysfunction (", count, " pathways) - suggests impaired cellular energy production"),
    protein_quality = paste0("Protein quality control (", count, " pathways) - indicates protein misfolding/aggregation issues"),
    autophagy = paste0("Autophagy/lysosomal pathways (", count, " pathways) - suggests impaired cellular cleanup mechanisms"),
    dopamine = paste0("Dopamine metabolism (", count, " pathways) - direct relevance to PD pathophysiology"),
    synaptic = paste0("Synaptic function (", count, " pathways) - indicates neurotransmission defects"),
    oxidative_stress = paste0("Oxidative stress response (", count, " pathways) - suggests cellular damage mechanisms"),
    neuronal = paste0("Neuronal development/function (", count, " pathways) - indicates neurodegenerative processes")
  )
  
  return(category_meanings[[category]] %||% paste0(category, " (", count, " pathways)"))
}

#' Generate biological significance summary
#'
#' @param top_categories Top biological categories
#' @param signature Signature information
#' @return Biological significance text
generate_biological_significance <- function(top_categories, signature) {
  
  if (length(top_categories) == 0) {
    return("Limited biological convergence detected between mutation and knockdown approaches.")
  }
  
  # Generate significance based on dominant categories
  if ("mitochondrial" %in% head(top_categories, 2)) {
    significance <- "Strong evidence for mitochondrial dysfunction, a core PD mechanism. "
  } else if ("dopamine" %in% head(top_categories, 2)) {
    significance <- "Direct impact on dopaminergic pathways, central to PD pathophysiology. "
  } else if ("protein_quality" %in% head(top_categories, 2)) {
    significance <- "Protein quality control defects, consistent with PD protein aggregation pathology. "
  } else {
    significance <- "Multiple PD-relevant pathways affected, suggesting complex disease mechanisms. "
  }
  
  # Add context about signature strength
  strength <- signature$signature_strength
  if (strength > 3) {
    significance <- paste0(significance, "High signature strength suggests robust cross-method agreement.")
  } else if (strength > 2) {
    significance <- paste0(significance, "Moderate signature strength indicates meaningful biological convergence.")
  } else {
    significance <- paste0(significance, "Modest signature strength suggests subtle but consistent effects.")
  }
  
  return(significance)
}

#' Calculate PD relevance score
#'
#' @param biological_categories Categorized biological processes
#' @return Numeric PD relevance score
calculate_pd_relevance_score <- function(biological_categories) {
  
  # Weight categories by PD relevance
  pd_weights <- list(
    mitochondrial = 3,      # Core PD mechanism
    dopamine = 3,           # Central to PD
    protein_quality = 2.5,  # Major PD feature
    autophagy = 2,          # Important PD mechanism
    synaptic = 2,           # PD-relevant
    oxidative_stress = 1.5, # Contributing factor
    neuronal = 1.5          # General neurodegeneration
  )
  
  score <- 0
  for (category in names(biological_categories)) {
    if (category %in% names(pd_weights)) {
      score <- score + (biological_categories[[category]] * pd_weights[[category]])
    }
  }
  
  return(score)
}

#' Create comprehensive PD signature summary
#'
#' @param enhanced_signatures List of enhanced signatures
#' @param analysis_type Type of analysis performed
#' @return Comprehensive summary for manuscript/interpretation
create_pd_signature_summary <- function(enhanced_signatures, analysis_type) {
  
  cat("[PD ANALYSIS] Creating comprehensive summary...\n")
  
  # Extract key metrics
  signature_count <- length(enhanced_signatures)
  
  # Aggregate biological categories across all signatures
  all_categories <- list(
    mitochondrial = 0, protein_quality = 0, autophagy = 0,
    dopamine = 0, synaptic = 0, oxidative_stress = 0, neuronal = 0
  )
  
  pd_relevance_scores <- numeric(signature_count)
  gene_pairs <- character(signature_count)
  
  for (i in seq_along(enhanced_signatures)) {
    sig <- enhanced_signatures[[i]]
    
    # Aggregate categories
    for (cat in names(all_categories)) {
      all_categories[[cat]] <- all_categories[[cat]] + sig$biological_categories[[cat]]
    }
    
    pd_relevance_scores[i] <- sig$pd_relevance_score
    gene_pairs[i] <- sig$signature$gene_pair
  }
  
  # Create summary statistics
  summary_stats <- list(
    analysis_type = analysis_type,
    total_signatures = signature_count,
    mean_pd_relevance = mean(pd_relevance_scores),
    top_biological_category = names(all_categories)[which.max(unlist(all_categories))],
    most_relevant_signature = gene_pairs[which.max(pd_relevance_scores)]
  )
  
  # Create biological insights summary
  biological_insights <- generate_overall_biological_insights(all_categories, summary_stats)
  
  return(list(
    summary_stats = summary_stats,
    biological_categories = all_categories,
    biological_insights = biological_insights,
    enhanced_signatures = enhanced_signatures
  ))
}

#' Generate overall biological insights
#'
#' @param all_categories Aggregated biological categories
#' @param summary_stats Summary statistics
#' @return Biological insights text
generate_overall_biological_insights <- function(all_categories, summary_stats) {
  
  # Rank categories by frequency
  category_ranking <- all_categories[order(unlist(all_categories), decreasing = TRUE)]
  top_3_categories <- head(names(category_ranking)[unlist(category_ranking) > 0], 3)
  
  insights <- paste0(
    "Analysis of ", summary_stats$total_signatures, " ", summary_stats$analysis_type, 
    " signatures reveals:\n\n"
  )
  
  if (length(top_3_categories) > 0) {
    insights <- paste0(insights, "Dominant Biological Themes:\n")
    for (i in seq_along(top_3_categories)) {
      cat_name <- top_3_categories[i]
      count <- all_categories[[cat_name]]
      insights <- paste0(insights, i, ". ", str_to_title(gsub("_", " ", cat_name)), 
                        " (", count, " pathway occurrences)\n")
    }
    
    insights <- paste0(insights, "\nBiological Interpretation:\n")
    insights <- paste0(insights, generate_cross_signature_interpretation(top_3_categories))
  } else {
    insights <- paste0(insights, "Limited PD-specific pathway convergence detected.")
  }
  
  return(insights)
}

#' Generate cross-signature biological interpretation
#'
#' @param top_categories Top biological categories
#' @return Cross-signature interpretation text
generate_cross_signature_interpretation <- function(top_categories) {
  
  interpretation <- ""
  
  if ("mitochondrial" %in% top_categories && "protein_quality" %in% top_categories) {
    interpretation <- paste0(interpretation, 
      "• Convergent evidence for mitochondrial dysfunction coupled with protein quality control defects, ",
      "suggesting a unified mechanism where energy deficits and protein misfolding reinforce each other.\n")
  }
  
  if ("dopamine" %in% top_categories) {
    interpretation <- paste0(interpretation,
      "• Direct disruption of dopaminergic pathways provides mechanistic insight into PD-specific vulnerability.\n")
  }
  
  if ("autophagy" %in% top_categories) {
    interpretation <- paste0(interpretation,
      "• Impaired autophagy/lysosomal function suggests defective cellular clearance mechanisms, ",
      "consistent with accumulation of toxic protein aggregates.\n")
  }
  
  if ("synaptic" %in% top_categories && "neuronal" %in% top_categories) {
    interpretation <- paste0(interpretation,
      "• Combined synaptic and neuronal dysfunction indicates progression from cellular defects ",
      "to network-level impairment.\n")
  }
  
  if (interpretation == "") {
    interpretation <- "• Multiple PD-relevant pathways affected, suggesting complex multi-system disease mechanisms.\n"
  }
  
  return(interpretation)
}


# Helper function for string manipulation
str_to_title <- function(x) {
  gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x, perl = TRUE)
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a