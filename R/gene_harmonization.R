#' Gene Harmonization Functions for MAST vs CRISPRi Comparisons
#'
#' This module handles gene name mapping between iSCORE-PD (MAST) and PerturbSeq (CRISPRi)
#' datasets, accounting for naming differences and variant handling.

#' Create gene mapping table for MAST vs CRISPRi comparisons
#'
#' @return Data frame with gene mappings
#' @export
create_gene_mapping_table <- function() {
  
  # Core gene mappings between MAST and CRISPRi
  gene_mappings <- data.frame(
    mast_gene = c(
      "ATP13A2", "DNAJC6", "FBXO7", "LRRK2", "PARK7", "PINK1", "SYNJ1",
      "PRKN",           # MAST name
      "SNCA_A30P",      # MAST variant 1
      "SNCA_A53T",      # MAST variant 2  
      "VPS13C_A444P",   # MAST variant 1
      "VPS13C_W395C"    # MAST variant 2
    ),
    crispri_gene = c(
      "ATP13A2", "DNAJC6", "FBXO7", "LRRK2", "PARK7", "PINK1", "SYNJ1",
      "PARK2",          # CRISPRi name for PRKN
      "SNCA",           # CRISPRi name (maps to both MAST variants)
      "SNCA",           # CRISPRi name (maps to both MAST variants)
      "VPS13C",         # CRISPRi name (maps to both MAST variants)
      "VPS13C"          # CRISPRi name (maps to both MAST variants)
    ),
    variant_group = c(
      "single", "single", "single", "single", "single", "single", "single",
      "single",         # PRKN/PARK2 
      "SNCA_variants",  # Group for SNCA variants
      "SNCA_variants",  # Group for SNCA variants
      "VPS13C_variants", # Group for VPS13C variants
      "VPS13C_variants"  # Group for VPS13C variants
    ),
    mast_available = TRUE,
    crispri_available = c(rep(TRUE, 7), TRUE, TRUE, TRUE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
  
  # Add GBA (MAST only - no CRISPRi counterpart)
  gba_row <- data.frame(
    mast_gene = "GBA",
    crispri_gene = NA_character_,
    variant_group = "mast_only",
    mast_available = TRUE,
    crispri_available = FALSE,
    stringsAsFactors = FALSE
  )
  
  gene_mappings <- rbind(gene_mappings, gba_row)
  
  return(gene_mappings)
}

#' Get comparable gene pairs for analysis
#'
#' @param combine_snca_variants Logical, whether to combine SNCA variants (default TRUE)
#' @param combine_vps13c_variants Logical, whether to combine VPS13C variants (default TRUE)
#' @param include_mast_only Logical, whether to include MAST-only genes like GBA (default FALSE)
#' @return Data frame with comparable gene pairs
#' @export
get_comparable_gene_pairs <- function(combine_snca_variants = TRUE, 
                                     combine_vps13c_variants = TRUE,
                                     include_mast_only = FALSE) {
  
  mapping_table <- create_gene_mapping_table()
  
  # Filter to genes available in both methods (unless including MAST-only)
  if (!include_mast_only) {
    mapping_table <- mapping_table[mapping_table$crispri_available, ]
  }
  
  # Initialize comparison_type column
  mapping_table$comparison_type <- NA_character_
  
  # Handle variant grouping
  if (combine_snca_variants) {
    # Combine SNCA variants - keep one representative
    snca_rows <- mapping_table[mapping_table$variant_group == "SNCA_variants", ]
    if (nrow(snca_rows) > 0) {
      # Create combined entry
      snca_combined <- snca_rows[1, ]
      snca_combined$mast_gene <- "SNCA_combined"
      snca_combined$comparison_type <- "variants_combined"
      
      # Remove individual variants and add combined
      mapping_table <- mapping_table[mapping_table$variant_group != "SNCA_variants", ]
      mapping_table <- rbind(mapping_table, snca_combined)
    }
  } else {
    # Keep variants separate - add comparison type
    mapping_table$comparison_type[mapping_table$variant_group == "SNCA_variants"] <- "variants_separate"
  }
  
  if (combine_vps13c_variants) {
    # Combine VPS13C variants - keep one representative
    vps13c_rows <- mapping_table[mapping_table$variant_group == "VPS13C_variants", ]
    if (nrow(vps13c_rows) > 0) {
      # Create combined entry
      vps13c_combined <- vps13c_rows[1, ]
      vps13c_combined$mast_gene <- "VPS13C_combined"
      vps13c_combined$comparison_type <- "variants_combined"
      
      # Remove individual variants and add combined
      mapping_table <- mapping_table[mapping_table$variant_group != "VPS13C_variants", ]
      mapping_table <- rbind(mapping_table, vps13c_combined)
    }
  } else {
    # Keep variants separate - add comparison type
    mapping_table$comparison_type[mapping_table$variant_group == "VPS13C_variants"] <- "variants_separate"
  }
  
  # Add comparison type for remaining genes
  mapping_table$comparison_type[is.na(mapping_table$comparison_type)] <- "direct"
  
  # Add metadata
  mapping_table$has_both_methods <- mapping_table$mast_available & mapping_table$crispri_available
  mapping_table$analysis_priority <- ifelse(mapping_table$has_both_methods, "high", "low")
  
  return(mapping_table)
}

#' Get mutation categories for genes
#'
#' @return Data frame with gene mutation categories
#' @export
get_mutation_categories <- function() {
  
  mutation_categories <- data.frame(
    gene = c(
      # Point Mutations (missense variants)
      "SNCA_A30P", "SNCA_A53T", "LRRK2", "VPS13C_W395C", "VPS13C_A444P", "SYNJ1",
      
      # Nonsense/Truncating Mutations
      "PINK1", "FBXO7",
      
      # Large Deletions
      "PRKN", "PARK7", 
      
      # Frameshift Mutations
      "ATP13A2", "DNAJC6",
      
      # Splice Site Mutations
      "GBA"
    ),
    mutation_category = c(
      # Point Mutations
      rep("Point_Mutation", 6),
      
      # Nonsense/Truncating
      rep("Nonsense_Truncating", 2),
      
      # Large Deletions
      rep("Large_Deletion", 2),
      
      # Frameshift
      rep("Frameshift", 2),
      
      # Splice Site
      "Splice_Site"
    ),
    expected_expression_effect = c(
      # Point mutations - minimal target gene effect expected
      rep("Minimal_Change", 6),
      
      # Nonsense/Truncating - moderate to severe reduction
      rep("Moderate_Severe_Reduction", 2),
      
      # Large deletions - strong reduction expected
      rep("Strong_Reduction", 2),
      
      # Frameshift - moderate to severe reduction
      rep("Moderate_Severe_Reduction", 2),
      
      # Splice site - reduced expression
      "Reduced_Expression"
    ),
    pathway_focus = c(
      # Point mutations - focus on downstream effects
      rep("Downstream_Pathways", 6),
      
      # Loss-of-function mutations - both target and downstream
      rep("Target_And_Downstream", 7)
    ),
    stringsAsFactors = FALSE
  )
  
  return(mutation_categories)
}

#' Filter enrichment data for specific gene comparisons
#'
#' @param enrichment_data Consolidated enrichment data
#' @param mast_genes Character vector of MAST gene names to include
#' @param crispri_genes Character vector of CRISPRi gene names to include
#' @param combine_variants Logical, whether to combine variant data
#' @return Filtered and harmonized enrichment data
#' @export
filter_for_gene_comparison <- function(enrichment_data, 
                                      mast_genes = NULL, 
                                      crispri_genes = NULL,
                                      combine_variants = TRUE) {
  
  if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
    warning("Empty or NULL enrichment data provided")
    return(data.frame())
  }
  
  # Filter by method
  mast_data <- enrichment_data[enrichment_data$method == "MAST", ]
  crispri_data <- enrichment_data[enrichment_data$method == "MixScale", ]
  
  # Filter by specific genes if provided
  if (!is.null(mast_genes)) {
    if (combine_variants) {
      # Handle variant combining for MAST data
      if ("SNCA_combined" %in% mast_genes) {
        snca_genes <- c("SNCA_A30P", "SNCA_A53T")
        mast_genes <- c(mast_genes[mast_genes != "SNCA_combined"], snca_genes)
      }
      if ("VPS13C_combined" %in% mast_genes) {
        vps13c_genes <- c("VPS13C_A444P", "VPS13C_W395C")
        mast_genes <- c(mast_genes[mast_genes != "VPS13C_combined"], vps13c_genes)
      }
    }
    
    mast_data <- mast_data[mast_data$mutation_perturbation %in% mast_genes, ]
  }
  
  if (!is.null(crispri_genes)) {
    crispri_data <- crispri_data[crispri_data$mutation_perturbation %in% crispri_genes, ]
  }
  
  # Combine filtered data
  filtered_data <- rbind(mast_data, crispri_data)
  
  # Add harmonized gene names for easier comparison
  if (combine_variants && nrow(filtered_data) > 0) {
    filtered_data$harmonized_gene <- filtered_data$mutation_perturbation
    
    # Map MAST variants to common names
    filtered_data$harmonized_gene[filtered_data$mutation_perturbation %in% c("SNCA_A30P", "SNCA_A53T")] <- "SNCA"
    filtered_data$harmonized_gene[filtered_data$mutation_perturbation %in% c("VPS13C_A444P", "VPS13C_W395C")] <- "VPS13C"
    filtered_data$harmonized_gene[filtered_data$mutation_perturbation == "PRKN"] <- "PARK2"
  } else {
    filtered_data$harmonized_gene <- filtered_data$mutation_perturbation
  }
  
  return(filtered_data)
}

#' Get PD-relevant pathway terms for prioritization
#'
#' @return Character vector of PD-relevant pathway terms
#' @export
get_pd_relevant_pathways <- function() {
  
  # PD-relevant pathway keywords for term matching
  pd_pathways <- c(
    # Mitochondrial dysfunction
    "mitochondrion", "mitochondrial", "respiratory chain", "oxidative phosphorylation",
    "electron transport", "ATP synthesis", 
    
    # Protein aggregation and quality control
    "protein folding", "unfolded protein response", "endoplasmic reticulum stress",
    "proteasome", "ubiquitin", "protein degradation", "chaperone",
    
    # Dopamine metabolism and signaling
    "dopamine", "dopaminergic", "catecholamine", "neurotransmitter",
    "synaptic transmission", "synaptic vesicle",
    
    # Autophagy and lysosomal function
    "autophagy", "lysosome", "phagosome", "endosome", "vacuole",
    
    # Oxidative stress and cellular defense
    "oxidative stress", "reactive oxygen species", "antioxidant",
    "cellular response to oxidative stress",
    
    # Neuronal function and development
    "axon", "dendrite", "neuron development", "neurogenesis",
    "synaptic plasticity", "neuron projection"
  )
  
  return(pd_pathways)
}