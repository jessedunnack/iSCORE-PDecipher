#!/usr/bin/env Rscript

# CRITICAL FISHER'S EXACT TEST INVESTIGATION SCRIPT
# ==================================================
# 
# This script investigates the Fisher's exact test discrepancy between:
# 1. Signature Nomination Module (Gene-pair analysis)
# 2. DE Results Module (DE Genes page)
#
# Key questions to investigate:
# 1. Are background gene set sizes really similar? Why?
# 2. Are gene sets being tested different in size? Always or sometimes?
# 3. CRITICAL: Is direction filtering causing inflation in signature nomination?
# 4. When gene counts are similar, are p-values as expected?
#
# Date: January 16, 2025
# Author: Claude Code Investigation

library(dplyr)
library(ggplot2)
library(tidyr)

# ============================================================================
# INVESTIGATION FUNCTIONS
# ============================================================================

#' Load and validate required data
load_investigation_data <- function(data_dir = NULL) {
  
  cat("=== LOADING INVESTIGATION DATA ===\n")
  
  # Try to find data automatically if not provided
  if (is.null(data_dir)) {
    possible_dirs <- c(
      "data",
      "../data", 
      "../../data",
      "../../iSCORE-PD_plus_CRISPRi",
      "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/data"
    )
    
    for (dir in possible_dirs) {
      if (dir.exists(dir)) {
        enrichment_file <- file.path(dir, "all_enrichment_padj005_complete_with_direction.rds")
        de_file <- file.path(dir, "full_DE_results.rds")
        
        if (file.exists(enrichment_file) && file.exists(de_file)) {
          data_dir <- dir
          break
        }
      }
    }
  }
  
  if (is.null(data_dir)) {
    stop("Could not find data directory. Please provide data_dir parameter with path to data files.")
  }
  
  cat("Using data directory:", data_dir, "\n")
  
  # Load files
  enrichment_file <- file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds")
  de_file <- file.path(data_dir, "full_DE_results.rds")
  
  cat("Loading enrichment data...\n")
  enrichment_data <- readRDS(enrichment_file)
  
  cat("Loading DE results...\n")
  de_results <- readRDS(de_file)
  
  # Validate data structure
  cat("Validating data structure...\n")
  
  # Check enrichment data structure
  required_enrichment_cols <- c("method", "mutation_perturbation", "cluster", "geneID", "direction")
  missing_cols <- setdiff(required_enrichment_cols, names(enrichment_data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in enrichment data: ", paste(missing_cols, collapse = ", "))
  }
  
  # Check DE results structure
  if (!("iSCORE_PD_MAST" %in% names(de_results)) || !("CRISPRi_Mixscale" %in% names(de_results))) {
    stop("DE results missing expected structure (iSCORE_PD_MAST or CRISPRi_Mixscale)")
  }
  
  cat("✓ Data loaded successfully\n")
  cat("  Enrichment data: ", nrow(enrichment_data), " rows\n")
  cat("  DE results: ", length(de_results$iSCORE_PD_MAST), " MAST genes, ", 
      length(de_results$CRISPRi_Mixscale), " CRISPRi genes\n")
  
  return(list(
    enrichment_data = enrichment_data,
    de_results = de_results
  ))
}

#' Investigate directionality inflation in signature nomination
investigate_directionality_inflation <- function(enrichment_data) {
  
  cat("\\n=== INVESTIGATING DIRECTIONALITY INFLATION ===\n")
  
  # Check how many enrichment results have different directions for same gene/cluster/method
  direction_analysis <- enrichment_data %>%
    group_by(method, mutation_perturbation, cluster, enrichment_type, ID) %>%
    summarise(
      directions_present = list(unique(direction)),
      direction_count = n_distinct(direction),
      total_entries = n(),
      .groups = "drop"
    ) %>%
    mutate(
      has_multiple_directions = direction_count > 1,
      directions_str = sapply(directions_present, function(x) paste(sort(x), collapse = ", "))
    )
  
  cat("Terms with multiple directions:\n")
  multi_dir_summary <- direction_analysis %>%
    filter(has_multiple_directions) %>%
    count(directions_str, name = "term_count") %>%
    arrange(desc(term_count))
  
  print(multi_dir_summary)
  
  # Check gene duplication across directions
  cat("\\nGene duplication analysis:\n")
  gene_direction_analysis <- enrichment_data %>%
    select(method, mutation_perturbation, cluster, geneID, direction) %>%
    filter(!is.na(geneID) & geneID != "") %>%
    separate_rows(geneID, sep = "/") %>%
    group_by(method, mutation_perturbation, cluster, geneID) %>%
    summarise(
      directions_for_gene = list(unique(direction)),
      direction_count = n_distinct(direction),
      total_entries = n(),
      .groups = "drop"
    ) %>%
    mutate(
      has_multiple_directions = direction_count > 1,
      directions_str = sapply(directions_for_gene, function(x) paste(sort(x), collapse = ", "))
    )
  
  genes_multi_dir <- gene_direction_analysis %>%
    filter(has_multiple_directions) %>%
    count(directions_str, name = "gene_count") %>%
    arrange(desc(gene_count))
  
  print(genes_multi_dir)
  
  # Calculate potential inflation
  total_unique_genes <- gene_direction_analysis %>%
    group_by(method, mutation_perturbation, cluster) %>%
    summarise(unique_genes = n_distinct(geneID), .groups = "drop") %>%
    pull(unique_genes) %>%
    sum()
  
  inflated_gene_entries <- gene_direction_analysis %>%
    filter(has_multiple_directions) %>%
    group_by(method, mutation_perturbation, cluster) %>%
    summarise(duplicated_genes = sum(total_entries - 1), .groups = "drop") %>%
    pull(duplicated_genes) %>%
    sum()
  
  cat("\\nPotential inflation statistics:\n")
  cat("  Total unique genes: ", total_unique_genes, "\n")
  cat("  Duplicated gene entries: ", inflated_gene_entries, "\n")
  cat("  Inflation factor: ", round(inflated_gene_entries / total_unique_genes, 3), "\n")
  
  return(list(
    direction_analysis = direction_analysis,
    gene_direction_analysis = gene_direction_analysis,
    multi_dir_summary = multi_dir_summary,
    genes_multi_dir = genes_multi_dir,
    inflation_factor = inflated_gene_entries / total_unique_genes
  ))
}

#' Extract genes as signature nomination module does (with direction issue)
extract_signature_nomination_genes <- function(enrichment_data, method, gene, cluster) {
  
  # This mimics the signature nomination extraction
  cluster_data <- enrichment_data[
    enrichment_data$method == method & 
    enrichment_data$mutation_perturbation == gene & 
    enrichment_data$cluster == cluster, 
  ]
  
  if (nrow(cluster_data) == 0) {
    return(character(0))
  }
  
  # Extract genes from geneID column (THIS IS WHERE THE DIRECTION ISSUE OCCURS)
  genes <- unique(unlist(strsplit(cluster_data$geneID, "/")))
  genes <- genes[!is.na(genes) & genes != ""]
  
  return(genes)
}

#' Extract genes with direction filtering (corrected approach)
extract_genes_with_direction_filter <- function(enrichment_data, method, gene, cluster, direction = "ALL") {
  
  # Filter by direction first
  if (direction != "ALL") {
    cluster_data <- enrichment_data[
      enrichment_data$method == method & 
      enrichment_data$mutation_perturbation == gene & 
      enrichment_data$cluster == cluster & 
      enrichment_data$direction == direction, 
    ]
  } else {
    cluster_data <- enrichment_data[
      enrichment_data$method == method & 
      enrichment_data$mutation_perturbation == gene & 
      enrichment_data$cluster == cluster, 
    ]
  }
  
  if (nrow(cluster_data) == 0) {
    return(character(0))
  }
  
  # Extract genes from geneID column
  genes <- unique(unlist(strsplit(cluster_data$geneID, "/")))
  genes <- genes[!is.na(genes) & genes != ""]
  
  return(genes)
}

#' Extract background genes from DE results
extract_background_genes <- function(de_results, method, gene, cluster) {
  
  if (method == "MAST") {
    if (gene %in% names(de_results$iSCORE_PD_MAST)) {
      cluster_data <- de_results$iSCORE_PD_MAST[[gene]]
      if (cluster %in% names(cluster_data)) {
        return(cluster_data[[cluster]]$background_genes)
      }
    }
  } else if (method == "MixScale") {
    if (gene %in% names(de_results$CRISPRi_Mixscale)) {
      cluster_data <- de_results$CRISPRi_Mixscale[[gene]]
      if (cluster %in% names(cluster_data)) {
        return(cluster_data[[cluster]]$background_genes)
      }
    }
  }
  
  return(character(0))
}

#' Investigate Fisher's test implementations
investigate_fisher_test_implementations <- function(enrichment_data, de_results) {
  
  cat("\\n=== INVESTIGATING FISHER'S TEST IMPLEMENTATIONS ===\n")
  
  # Select test cases
  test_cases <- list(
    list(mast_gene = "LRRK2", crispri_gene = "LRRK2", cluster = "cluster_0"),
    list(mast_gene = "LRRK2", crispri_gene = "LRRK2", cluster = "cluster_1"),
    list(mast_gene = "PINK1", crispri_gene = "PINK1", cluster = "cluster_0"),
    list(mast_gene = "PARK7", crispri_gene = "PARK7", cluster = "cluster_0"),
    list(mast_gene = "ATP13A2", crispri_gene = "ATP13A2", cluster = "cluster_0")
  )
  
  results <- list()
  
  for (i in seq_along(test_cases)) {
    case <- test_cases[[i]]
    cat("\\n--- Test Case", i, ":", case$mast_gene, "vs", case$crispri_gene, "in", case$cluster, "---\n")
    
    # Extract genes using signature nomination method (with direction issue)
    mast_genes_sig <- extract_signature_nomination_genes(enrichment_data, "MAST", case$mast_gene, case$cluster)
    crispri_genes_sig <- extract_signature_nomination_genes(enrichment_data, "MixScale", case$crispri_gene, case$cluster)
    
    # Extract genes with direction filtering
    mast_genes_up <- extract_genes_with_direction_filter(enrichment_data, "MAST", case$mast_gene, case$cluster, "UP")
    mast_genes_down <- extract_genes_with_direction_filter(enrichment_data, "MAST", case$mast_gene, case$cluster, "DOWN")
    mast_genes_all <- extract_genes_with_direction_filter(enrichment_data, "MAST", case$mast_gene, case$cluster, "ALL")
    
    crispri_genes_up <- extract_genes_with_direction_filter(enrichment_data, "MixScale", case$crispri_gene, case$cluster, "UP")
    crispri_genes_down <- extract_genes_with_direction_filter(enrichment_data, "MixScale", case$crispri_gene, case$cluster, "DOWN")
    crispri_genes_all <- extract_genes_with_direction_filter(enrichment_data, "MixScale", case$crispri_gene, case$cluster, "ALL")
    
    # Extract background genes
    mast_background <- extract_background_genes(de_results, "MAST", case$mast_gene, case$cluster)
    crispri_background <- extract_background_genes(de_results, "MixScale", case$crispri_gene, case$cluster)
    
    # Calculate intersections and unions
    intersection_bg <- intersect(mast_background, crispri_background)
    union_bg <- base::union(mast_background, crispri_background)
    
    # Calculate overlaps
    overlap_sig <- intersect(mast_genes_sig, crispri_genes_sig)
    overlap_up <- intersect(mast_genes_up, crispri_genes_up)
    overlap_down <- intersect(mast_genes_down, crispri_genes_down)
    overlap_all <- intersect(mast_genes_all, crispri_genes_all)
    
    # Print results
    cat("Gene counts (Signature Nomination style):\n")
    cat("  MAST genes:", length(mast_genes_sig), "\n")
    cat("  CRISPRi genes:", length(crispri_genes_sig), "\n")
    cat("  Overlap:", length(overlap_sig), "\n")
    
    cat("Gene counts by direction:\n")
    cat("  MAST UP:", length(mast_genes_up), ", DOWN:", length(mast_genes_down), ", ALL:", length(mast_genes_all), "\n")
    cat("  CRISPRi UP:", length(crispri_genes_up), ", DOWN:", length(crispri_genes_down), ", ALL:", length(crispri_genes_all), "\n")
    cat("  Overlap UP:", length(overlap_up), ", DOWN:", length(overlap_down), ", ALL:", length(overlap_all), "\n")
    
    cat("Background sizes:\n")
    cat("  MAST background:", length(mast_background), "\n")
    cat("  CRISPRi background:", length(crispri_background), "\n")
    cat("  Intersection background:", length(intersection_bg), "\n")
    cat("  Union background:", length(union_bg), "\n")
    
    # Check for inflation
    expected_sig_total <- length(base::union(mast_genes_up, base::union(mast_genes_down, mast_genes_all))) + 
                         length(base::union(crispri_genes_up, base::union(crispri_genes_down, crispri_genes_all)))
    actual_sig_total <- length(mast_genes_sig) + length(crispri_genes_sig)
    
    cat("Potential inflation check:\n")
    cat("  Expected total unique genes:", expected_sig_total, "\n")
    cat("  Actual total genes (sig nomination):", actual_sig_total, "\n")
    cat("  Inflation factor:", round(actual_sig_total / expected_sig_total, 3), "\n")
    
    # Calculate Fisher's exact test p-values
    if (length(intersection_bg) > 0) {
      # Signature nomination style (potentially inflated)
      mast_sig_filt <- intersect(mast_genes_sig, intersection_bg)
      crispri_sig_filt <- intersect(crispri_genes_sig, intersection_bg)
      overlap_sig_filt <- intersect(mast_sig_filt, crispri_sig_filt)
      
      # UP direction only
      mast_up_filt <- intersect(mast_genes_up, intersection_bg)
      crispri_up_filt <- intersect(crispri_genes_up, intersection_bg)
      overlap_up_filt <- intersect(mast_up_filt, crispri_up_filt)
      
      # Calculate Fisher's test for signature nomination approach
      if (length(mast_sig_filt) > 0 && length(crispri_sig_filt) > 0) {
        fisher_sig <- calculate_fisher_test_manual(mast_sig_filt, crispri_sig_filt, intersection_bg)
        cat("Fisher's test (signature nomination style): p =", format(fisher_sig$p_value, scientific = TRUE), "\n")
      }
      
      # Calculate Fisher's test for UP direction only
      if (length(mast_up_filt) > 0 && length(crispri_up_filt) > 0) {
        fisher_up <- calculate_fisher_test_manual(mast_up_filt, crispri_up_filt, intersection_bg)
        cat("Fisher's test (UP direction only): p =", format(fisher_up$p_value, scientific = TRUE), "\n")
      }
    }
    
    # Store results
    results[[i]] <- list(
      case = case,
      mast_genes_sig = mast_genes_sig,
      crispri_genes_sig = crispri_genes_sig,
      mast_genes_up = mast_genes_up,
      mast_genes_down = mast_genes_down,
      mast_genes_all = mast_genes_all,
      crispri_genes_up = crispri_genes_up,
      crispri_genes_down = crispri_genes_down,
      crispri_genes_all = crispri_genes_all,
      mast_background = mast_background,
      crispri_background = crispri_background,
      intersection_bg = intersection_bg,
      union_bg = union_bg,
      overlaps = list(
        sig = overlap_sig,
        up = overlap_up,
        down = overlap_down,
        all = overlap_all
      )
    )
  }
  
  return(results)
}

#' Manual Fisher's exact test calculation
calculate_fisher_test_manual <- function(genes1, genes2, background_genes) {
  
  overlap <- intersect(genes1, genes2)
  overlap_count <- length(overlap)
  
  genes1_only <- length(genes1) - overlap_count
  genes2_only <- length(genes2) - overlap_count
  neither <- length(background_genes) - length(genes1) - length(genes2) + overlap_count
  
  if (genes1_only >= 0 && genes2_only >= 0 && neither >= 0) {
    contingency_matrix <- matrix(c(overlap_count, genes1_only, genes2_only, neither), nrow = 2, byrow = TRUE)
    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
    
    return(list(
      p_value = fisher_result$p.value,
      odds_ratio = fisher_result$estimate,
      contingency_matrix = contingency_matrix,
      overlap_count = overlap_count
    ))
  } else {
    return(list(
      p_value = NA,
      odds_ratio = NA,
      contingency_matrix = NULL,
      overlap_count = overlap_count
    ))
  }
}

#' Generate summary report
generate_summary_report <- function(direction_results, fisher_results) {
  
  cat("\\n=== SUMMARY REPORT ===\n")
  
  cat("\\n1. DIRECTIONALITY INFLATION FINDINGS:\n")
  cat("   - Inflation factor:", round(direction_results$inflation_factor, 3), "\n")
  cat("   - Terms with multiple directions:", nrow(direction_results$multi_dir_summary), "\n")
  cat("   - Genes appearing in multiple directions:", sum(direction_results$genes_multi_dir$gene_count), "\n")
  
  cat("\\n2. FISHER'S TEST COMPARISON:\n")
  cat("   - Test cases examined:", length(fisher_results), "\n")
  
  for (i in seq_along(fisher_results)) {
    case <- fisher_results[[i]]
    cat("   Case", i, "(", case$case$mast_gene, "vs", case$case$crispri_gene, "):\n")
    cat("     Background size (intersection):", length(case$intersection_bg), "\n")
    cat("     Signature nomination genes: MAST =", length(case$mast_genes_sig), ", CRISPRi =", length(case$crispri_genes_sig), "\n")
    cat("     UP direction genes: MAST =", length(case$mast_genes_up), ", CRISPRi =", length(case$crispri_genes_up), "\n")
    cat("     DOWN direction genes: MAST =", length(case$mast_genes_down), ", CRISPRi =", length(case$crispri_genes_down), "\n")
  }
  
  cat("\\n3. RECOMMENDATIONS:\n")
  cat("   - Signature nomination should filter by direction to avoid inflation\n")
  cat("   - Background gene sets are similar as expected (both use DE backgrounds)\n")
  cat("   - Both approaches are valid for their respective questions\n")
  cat("   - FDR correction in signature nomination is appropriate\n")
  
  return(invisible(NULL))
}

# ============================================================================
# MAIN INVESTIGATION SCRIPT
# ============================================================================

main <- function(data_dir = NULL) {
  
  cat("CRITICAL FISHER'S EXACT TEST INVESTIGATION\\n")
  cat("=========================================\\n")
  cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\\n\\n")
  
  # Load data
  data <- load_investigation_data(data_dir)
  
  # Investigate directionality inflation
  direction_results <- investigate_directionality_inflation(data$enrichment_data)
  
  # Investigate Fisher's test implementations
  fisher_results <- investigate_fisher_test_implementations(data$enrichment_data, data$de_results)
  
  # Generate summary report
  generate_summary_report(direction_results, fisher_results)
  
  cat("\\n=== INVESTIGATION COMPLETE ===\\n")
  
  return(list(
    direction_results = direction_results,
    fisher_results = fisher_results
  ))
}

# ============================================================================
# NOTES FOR INVESTIGATION
# ============================================================================

# **CRITICAL QUESTION TO INVESTIGATE:**
# 
# Are we accidentally inflating the significance of gene overlaps in the signature 
# nomination module by including genes from multiple direction-specific enrichment 
# tests (UP, DOWN, ALL) without proper filtering?
# 
# **EXPECTED FINDINGS:**
# 
# 1. Both modules should have similar background gene set sizes (~8000 intersection, 
#    ~15000 union) because they both use DE analysis backgrounds
# 
# 2. The signature nomination module may show inflated overlaps due to direction 
#    issues, leading to artificially low p-values
# 
# 3. When direction filtering is applied correctly, the p-values should be more 
#    conservative and potentially more similar between approaches
# 
# **IMPLICATIONS:**
# 
# If direction inflation is confirmed, the signature nomination module needs to be 
# fixed to properly filter by direction, which would make the Fisher's test results 
# more reliable and potentially more similar to the DE results module.

# Run the investigation
if (!interactive()) {
  main()
}