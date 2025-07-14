#' Enhanced Direction-Aware Analysis Functions (v0.2.6)
#'
#' This module implements direction-aware statistical methods for cross-method
#' comparison between MAST (mutations) and CRISPRi (knockdowns), with biological
#' context weighting for different mutation types (LRRK2 vs SNCA variants).

#' Get biological direction expectation for each gene
#'
#' Based on user guidance:
#' - LRRK2: Opposing effects (gain-of-function mutation vs loss-of-function CRISPRi)
#' - SNCA variants: Same-direction effects (aggregation-related in both cases)  
#' - Other mutations: Generally same-direction (loss-of-function mutations)
#'
#' @param gene_name Gene name (MAST or CRISPRi format)
#' @return Character indicating expected direction pattern
#' @export
get_biological_direction_expectation <- function(gene_name) {
  
  # Normalize gene name variations (using base R instead of dplyr)
  if (gene_name %in% c("LRRK2")) {
    normalized_gene <- "LRRK2"
  } else if (gene_name %in% c("SNCA_A30P", "SNCA_A53T", "SNCA")) {
    normalized_gene <- "SNCA"
  } else if (gene_name %in% c("PRKN", "PARK2")) {
    normalized_gene <- "PARK2_PRKN"
  } else if (gene_name %in% c("VPS13C_A444P", "VPS13C_W395C", "VPS13C")) {
    normalized_gene <- "VPS13C"
  } else {
    normalized_gene <- gene_name
  }
  
  # Determine biological expectation (using base R)
  if (normalized_gene == "LRRK2") {
    # LRRK2: Hyperactive kinase (gain-of-function) vs CRISPRi knockdown (loss-of-function)
    # Expect OPPOSING effects on downstream gene expression
    direction_expectation <- "opposing"
  } else if (normalized_gene == "SNCA") {
    # SNCA variants: Aggregation-prone variants (Lewy body formation) vs CRISPRi (less aggregation)
    # Both may reduce aggregation pathways, expect SAME-DIRECTION effects
    direction_expectation <- "same"
  } else if (normalized_gene %in% c("PINK1", "PARK2_PRKN", "PARK7", "ATP13A2", "DNAJC6", "FBXO7", "GBA")) {
    # Loss-of-function mutations: Should show similar effects to CRISPRi knockdown
    # Expect SAME-DIRECTION effects
    direction_expectation <- "same"
  } else if (normalized_gene == "VPS13C") {
    # VPS13C: Point mutations, behavior may be complex
    direction_expectation <- "mixed"
  } else if (normalized_gene == "SYNJ1") {
    # SYNJ1: Point mutation, may have mixed effects
    direction_expectation <- "mixed"
  } else {
    # Unknown genes: Test both patterns equally
    direction_expectation <- "unknown"
  }
  
  return(direction_expectation)
}

#' Calculate same-direction gene overlap with Fisher's exact test
#'
#' Tests for genes that are differentially expressed in the SAME direction
#' in both MAST and CRISPRi experiments (both up or both down)
#'
#' @param mast_data Data frame with MAST DE results (must have avg_log2FC and p_val_adj columns)
#' @param crispri_data Data frame with CRISPRi DE results (must have log2FC and p_val_adj columns)
#' @param background_genes Character vector of genes tested in both methods
#' @param lfc_threshold Log2 fold change threshold for significance (default 0.25)
#' @param p_threshold P-value threshold for significance (default 0.05)
#' @return List with same-direction overlap statistics
#' @export
calculate_same_direction_overlap <- function(mast_data, crispri_data, background_genes = NULL,
                                           lfc_threshold = 0.25, p_threshold = 0.05) {
  
  # Input validation
  if (is.null(mast_data) || is.null(crispri_data)) {
    return(list(
      overlap_count = 0,
      same_up_count = 0,
      same_down_count = 0,
      fisher_p = NA,
      fisher_or = NA,
      contingency_table = NULL,
      error = "NULL input data"
    ))
  }
  
  # Ensure required columns exist
  required_mast_cols <- c("avg_log2FC", "p_val_adj")
  required_crispri_cols <- c("log2FC", "p_val_adj")  # May be experiment-specific like log2FC_C12_FPD-24
  
  # Handle experiment-specific column names for CRISPRi
  crispri_lfc_col <- names(crispri_data)[grepl("^log2FC", names(crispri_data))][1]
  if (is.na(crispri_lfc_col)) {
    crispri_lfc_col <- "log2FC"  # Fallback
  }
  
  # Extract significantly upregulated genes
  mast_up <- rownames(mast_data)[mast_data$avg_log2FC > lfc_threshold & 
                                 mast_data$p_val_adj < p_threshold & 
                                 !is.na(mast_data$avg_log2FC) & 
                                 !is.na(mast_data$p_val_adj)]
  
  mast_down <- rownames(mast_data)[mast_data$avg_log2FC < -lfc_threshold & 
                                   mast_data$p_val_adj < p_threshold & 
                                   !is.na(mast_data$avg_log2FC) & 
                                   !is.na(mast_data$p_val_adj)]
  
  # Handle CRISPRi data (may have different row structure)
  if (is.null(rownames(crispri_data)) || all(rownames(crispri_data) == 1:nrow(crispri_data))) {
    # Use gene column if rownames not set
    if ("gene" %in% names(crispri_data)) {
      crispri_genes <- crispri_data$gene
    } else {
      crispri_genes <- crispri_data[[1]]  # Assume first column is genes
    }
  } else {
    crispri_genes <- rownames(crispri_data)
  }
  
  crispri_up <- crispri_genes[crispri_data[[crispri_lfc_col]] > lfc_threshold & 
                              crispri_data$p_val_adj < p_threshold & 
                              !is.na(crispri_data[[crispri_lfc_col]]) & 
                              !is.na(crispri_data$p_val_adj)]
  
  crispri_down <- crispri_genes[crispri_data[[crispri_lfc_col]] < -lfc_threshold & 
                                crispri_data$p_val_adj < p_threshold & 
                                !is.na(crispri_data[[crispri_lfc_col]]) & 
                                !is.na(crispri_data$p_val_adj)]
  
  # Calculate same-direction overlaps
  same_up <- intersect(mast_up, crispri_up)  # Both methods upregulate
  same_down <- intersect(mast_down, crispri_down)  # Both methods downregulate
  same_direction_genes <- c(same_up, same_down)
  
  # Use background genes or infer from data
  if (is.null(background_genes)) {
    background_genes <- intersect(rownames(mast_data), crispri_genes)
  }
  
  # Calculate Fisher's exact test for same-direction overlap
  mast_any_direction <- c(mast_up, mast_down)
  crispri_any_direction <- c(crispri_up, crispri_down)
  
  # Create contingency table for same-direction genes
  same_and_both <- length(same_direction_genes)
  same_mast_not_crispri <- length(setdiff(mast_any_direction, crispri_any_direction))
  same_crispri_not_mast <- length(setdiff(crispri_any_direction, mast_any_direction))
  neither_same <- length(background_genes) - same_and_both - same_mast_not_crispri - same_crispri_not_mast
  
  # Ensure non-negative values
  neither_same <- max(0, neither_same)
  
  contingency_table <- matrix(c(same_and_both, same_mast_not_crispri,
                               same_crispri_not_mast, neither_same),
                             nrow = 2, ncol = 2,
                             dimnames = list(c("MAST_sig", "MAST_not_sig"),
                                           c("CRISPRi_sig", "CRISPRi_not_sig")))
  
  # Perform Fisher's exact test
  if (all(contingency_table >= 0) && sum(contingency_table) > 0) {
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    fisher_p <- fisher_result$p.value
    fisher_or <- fisher_result$estimate
  } else {
    fisher_p <- NA
    fisher_or <- NA
  }
  
  return(list(
    overlap_count = length(same_direction_genes),
    same_up_count = length(same_up),
    same_down_count = length(same_down),
    same_up_genes = same_up,
    same_down_genes = same_down,
    all_same_direction_genes = same_direction_genes,
    fisher_p = fisher_p,
    fisher_or = fisher_or,
    contingency_table = contingency_table,
    background_size = length(background_genes),
    mast_sig_count = length(mast_any_direction),
    crispri_sig_count = length(crispri_any_direction)
  ))
}

#' Calculate opposite-direction gene overlap with Fisher's exact test
#'
#' Tests for genes that are differentially expressed in OPPOSITE directions
#' in MAST vs CRISPRi experiments (one up, other down)
#'
#' @param mast_data Data frame with MAST DE results
#' @param crispri_data Data frame with CRISPRi DE results  
#' @param background_genes Character vector of genes tested in both methods
#' @param lfc_threshold Log2 fold change threshold for significance (default 0.25)
#' @param p_threshold P-value threshold for significance (default 0.05)
#' @return List with opposite-direction overlap statistics
#' @export
calculate_opposite_direction_overlap <- function(mast_data, crispri_data, background_genes = NULL,
                                               lfc_threshold = 0.25, p_threshold = 0.05) {
  
  # Input validation
  if (is.null(mast_data) || is.null(crispri_data)) {
    return(list(
      overlap_count = 0,
      mast_up_crispri_down_count = 0,
      mast_down_crispri_up_count = 0,
      fisher_p = NA,
      fisher_or = NA,
      contingency_table = NULL,
      error = "NULL input data"
    ))
  }
  
  # Handle experiment-specific column names for CRISPRi
  crispri_lfc_col <- names(crispri_data)[grepl("^log2FC", names(crispri_data))][1]
  if (is.na(crispri_lfc_col)) {
    crispri_lfc_col <- "log2FC"  # Fallback
  }
  
  # Extract significantly regulated genes by direction
  mast_up <- rownames(mast_data)[mast_data$avg_log2FC > lfc_threshold & 
                                 mast_data$p_val_adj < p_threshold & 
                                 !is.na(mast_data$avg_log2FC) & 
                                 !is.na(mast_data$p_val_adj)]
  
  mast_down <- rownames(mast_data)[mast_data$avg_log2FC < -lfc_threshold & 
                                   mast_data$p_val_adj < p_threshold & 
                                   !is.na(mast_data$avg_log2FC) & 
                                   !is.na(mast_data$p_val_adj)]
  
  # Handle CRISPRi data structure
  if (is.null(rownames(crispri_data)) || all(rownames(crispri_data) == 1:nrow(crispri_data))) {
    if ("gene" %in% names(crispri_data)) {
      crispri_genes <- crispri_data$gene
    } else {
      crispri_genes <- crispri_data[[1]]
    }
  } else {
    crispri_genes <- rownames(crispri_data)
  }
  
  crispri_up <- crispri_genes[crispri_data[[crispri_lfc_col]] > lfc_threshold & 
                              crispri_data$p_val_adj < p_threshold & 
                              !is.na(crispri_data[[crispri_lfc_col]]) & 
                              !is.na(crispri_data$p_val_adj)]
  
  crispri_down <- crispri_genes[crispri_data[[crispri_lfc_col]] < -lfc_threshold & 
                                crispri_data$p_val_adj < p_threshold & 
                                !is.na(crispri_data[[crispri_lfc_col]]) & 
                                !is.na(crispri_data$p_val_adj)]
  
  # Calculate opposite-direction overlaps
  mast_up_crispri_down <- intersect(mast_up, crispri_down)  # MAST up, CRISPRi down
  mast_down_crispri_up <- intersect(mast_down, crispri_up)  # MAST down, CRISPRi up
  opposite_direction_genes <- c(mast_up_crispri_down, mast_down_crispri_up)
  
  # Use background genes or infer from data
  if (is.null(background_genes)) {
    background_genes <- intersect(rownames(mast_data), crispri_genes)
  }
  
  # Calculate Fisher's exact test for opposite-direction overlap
  mast_any_direction <- c(mast_up, mast_down)
  crispri_any_direction <- c(crispri_up, crispri_down)
  
  # Create contingency table for opposite-direction genes
  opposite_and_both <- length(opposite_direction_genes)
  mast_sig_not_opposite <- length(mast_any_direction) - opposite_and_both
  crispri_sig_not_opposite <- length(crispri_any_direction) - opposite_and_both
  neither_nor_opposite <- length(background_genes) - opposite_and_both - mast_sig_not_opposite - crispri_sig_not_opposite
  
  # Ensure non-negative values
  neither_nor_opposite <- max(0, neither_nor_opposite)
  
  contingency_table <- matrix(c(opposite_and_both, mast_sig_not_opposite,
                               crispri_sig_not_opposite, neither_nor_opposite),
                             nrow = 2, ncol = 2,
                             dimnames = list(c("MAST_sig", "MAST_not_sig"),
                                           c("CRISPRi_opposite", "CRISPRi_not_opposite")))
  
  # Perform Fisher's exact test
  if (all(contingency_table >= 0) && sum(contingency_table) > 0) {
    fisher_result <- fisher.test(contingency_table, alternative = "greater")
    fisher_p <- fisher_result$p.value
    fisher_or <- fisher_result$estimate
  } else {
    fisher_p <- NA
    fisher_or <- NA
  }
  
  return(list(
    overlap_count = length(opposite_direction_genes),
    mast_up_crispri_down_count = length(mast_up_crispri_down),
    mast_down_crispri_up_count = length(mast_down_crispri_up),
    mast_up_crispri_down_genes = mast_up_crispri_down,
    mast_down_crispri_up_genes = mast_down_crispri_up,
    all_opposite_direction_genes = opposite_direction_genes,
    fisher_p = fisher_p,
    fisher_or = fisher_or,
    contingency_table = contingency_table,
    background_size = length(background_genes),
    mast_sig_count = length(mast_any_direction),
    crispri_sig_count = length(crispri_any_direction)
  ))
}

#' Combine direction-aware p-values using biological weighting
#'
#' Combines same-direction and opposite-direction Fisher's test p-values
#' using biological expectations for each gene type
#'
#' @param same_direction_p P-value from same-direction test
#' @param opposite_direction_p P-value from opposite-direction test
#' @param gene_name Gene name for biological context
#' @param primary_weight Weight for the primary (expected) direction test
#' @param secondary_weight Weight for the secondary (alternative) direction test
#' @return Combined p-value using weighted Fisher's method
#' @export
combine_direction_pvalues <- function(same_direction_p, opposite_direction_p, gene_name = NULL,
                                    primary_weight = 0.8, secondary_weight = 0.2) {
  
  # Determine biological expectation if gene name provided
  if (!is.null(gene_name)) {
    direction_expectation <- get_biological_direction_expectation(gene_name)
    
    if (direction_expectation == "opposing") {
      # LRRK2 case: prioritize opposite-direction test
      primary_p <- opposite_direction_p
      secondary_p <- same_direction_p
      primary_test <- "opposite"
    } else if (direction_expectation == "same") {
      # SNCA variants, loss-of-function mutations: prioritize same-direction test
      primary_p <- same_direction_p
      secondary_p <- opposite_direction_p
      primary_test <- "same"
    } else {
      # Unknown/mixed: equal weighting
      primary_p <- same_direction_p
      secondary_p <- opposite_direction_p
      primary_weight <- 0.5
      secondary_weight <- 0.5
      primary_test <- "equal"
    }
  } else {
    # No biological context: equal weighting
    primary_p <- same_direction_p
    secondary_p <- opposite_direction_p
    primary_weight <- 0.5
    secondary_weight <- 0.5
    primary_test <- "equal"
  }
  
  # Handle missing p-values
  if (is.na(primary_p) && is.na(secondary_p)) {
    return(list(
      combined_p = NA,
      primary_test = primary_test,
      primary_p = primary_p,
      secondary_p = secondary_p,
      weights = c(primary_weight, secondary_weight),
      error = "Both p-values are NA"
    ))
  }
  
  if (is.na(primary_p)) {
    # Use secondary test only
    combined_p <- secondary_p
    effective_test <- "secondary_only"
  } else if (is.na(secondary_p)) {
    # Use primary test only
    combined_p <- primary_p
    effective_test <- "primary_only"
  } else {
    # Combine using weighted Fisher's method
    # Convert p-values to chi-square statistics
    chi_square_primary <- -2 * log(primary_p)
    chi_square_secondary <- -2 * log(secondary_p)
    
    # Apply weights
    weighted_chi_square <- (chi_square_primary * primary_weight + 
                           chi_square_secondary * secondary_weight)
    
    # Convert back to p-value (approximation for weighted combination)
    # Using effective degrees of freedom based on weights
    effective_df <- 2 * (primary_weight + secondary_weight)
    combined_p <- pchisq(weighted_chi_square, df = effective_df, lower.tail = FALSE)
    effective_test <- "weighted_combination"
  }
  
  return(list(
    combined_p = combined_p,
    primary_test = primary_test,
    primary_p = primary_p,
    secondary_p = secondary_p,
    weights = c(primary_weight, secondary_weight),
    effective_test = effective_test,
    biological_expectation = if (!is.null(gene_name)) get_biological_direction_expectation(gene_name) else "unknown"
  ))
}

#' Enhanced direction-aware Fisher's test with biological expectations
#'
#' Main function that performs comprehensive direction-aware analysis
#' combining same and opposite direction tests with biological context
#'
#' @param mast_data Data frame with MAST DE results
#' @param crispri_data Data frame with CRISPRi DE results
#' @param gene_name Gene name for biological context
#' @param background_genes Character vector of background genes
#' @param lfc_threshold Log2 fold change threshold (default 0.25)
#' @param p_threshold P-value threshold (default 0.05)
#' @return List with comprehensive direction-aware analysis results
#' @export
enhanced_direction_analysis <- function(mast_data, crispri_data, gene_name, 
                                      background_genes = NULL,
                                      lfc_threshold = 0.25, p_threshold = 0.05) {
  
  cat("[DIRECTION ANALYSIS] Analyzing", gene_name, "with biological expectation\n")
  
  # Determine biological direction expectation
  direction_expectation <- get_biological_direction_expectation(gene_name)
  cat("[DIRECTION ANALYSIS] Biological expectation for", gene_name, ":", direction_expectation, "\n")
  
  # Calculate direction-specific overlaps
  same_direction_result <- calculate_same_direction_overlap(
    mast_data, crispri_data, background_genes, lfc_threshold, p_threshold
  )
  
  opposite_direction_result <- calculate_opposite_direction_overlap(
    mast_data, crispri_data, background_genes, lfc_threshold, p_threshold
  )
  
  # Combine p-values using biological weighting
  combined_result <- combine_direction_pvalues(
    same_direction_result$fisher_p,
    opposite_direction_result$fisher_p,
    gene_name
  )
  
  # Determine primary pattern based on biological expectation
  if (direction_expectation == "opposing") {
    primary_result <- opposite_direction_result
    secondary_result <- same_direction_result
    primary_pattern <- "opposite"
  } else if (direction_expectation == "same") {
    primary_result <- same_direction_result
    secondary_result <- opposite_direction_result
    primary_pattern <- "same"
  } else {
    # For mixed/unknown, choose the more significant result
    if (!is.na(same_direction_result$fisher_p) && !is.na(opposite_direction_result$fisher_p)) {
      if (same_direction_result$fisher_p <= opposite_direction_result$fisher_p) {
        primary_result <- same_direction_result
        secondary_result <- opposite_direction_result
        primary_pattern <- "same"
      } else {
        primary_result <- opposite_direction_result
        secondary_result <- same_direction_result
        primary_pattern <- "opposite"
      }
    } else if (!is.na(same_direction_result$fisher_p)) {
      primary_result <- same_direction_result
      secondary_result <- opposite_direction_result
      primary_pattern <- "same"
    } else {
      primary_result <- opposite_direction_result
      secondary_result <- same_direction_result
      primary_pattern <- "opposite"
    }
  }
  
  cat("[DIRECTION ANALYSIS] Results for", gene_name, ":\n")
  cat("  - Same direction p-value:", same_direction_result$fisher_p %||% "NA", "\n")
  cat("  - Opposite direction p-value:", opposite_direction_result$fisher_p %||% "NA", "\n") 
  cat("  - Combined p-value:", combined_result$combined_p %||% "NA", "\n")
  cat("  - Primary pattern:", primary_pattern, "\n")
  
  return(list(
    gene_name = gene_name,
    biological_expectation = direction_expectation,
    same_direction = same_direction_result,
    opposite_direction = opposite_direction_result,
    combined_analysis = combined_result,
    primary_pattern = primary_pattern,
    primary_result = primary_result,
    secondary_result = secondary_result,
    analysis_parameters = list(
      lfc_threshold = lfc_threshold,
      p_threshold = p_threshold,
      background_size = length(background_genes %||% c())
    )
  ))
}

# Helper function for null coalescing
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a