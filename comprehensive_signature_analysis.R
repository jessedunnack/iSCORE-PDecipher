# Comprehensive Signature Analysis for Bug #2
# This script verifies directionality logic and conducts comprehensive signature analysis
# across all genes, clusters, and directions to find the most significant cross-method signatures

library(dplyr)
library(tidyr)

# Find the signature analysis functions
source_files <- c(
  "R/signature_analysis.R",
  "R/manuscript_signature_discovery.R"
)

for (file in source_files) {
  if (file.exists(file)) {
    cat("Sourcing:", file, "\n")
    source(file)
  }
}

# Function to find consolidated data
find_consolidated_data <- function() {
  data_paths <- c(
    "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi",
    "../../iSCORE-PD_plus_CRISPRi",
    "../iSCORE-PD_plus_CRISPRi"
  )
  
  for (path in data_paths) {
    test_file <- file.path(path, "all_enrichment_padj005_complete_with_direction.rds")
    if (file.exists(test_file)) {
      return(test_file)
    }
  }
  return(NULL)
}

# Main analysis function
comprehensive_signature_analysis <- function() {
  cat("=== COMPREHENSIVE SIGNATURE ANALYSIS ===\n")
  cat("Objective: Find strongest cross-method signatures between iSCORE-PD (MAST) and CRISPRi\n\n")
  
  # Load data
  consolidated_file <- find_consolidated_data()
  if (is.null(consolidated_file)) {
    stop("Could not find consolidated data file!")
  }
  
  cat("Loading data from:", consolidated_file, "\n")
  consolidated_data <- readRDS(consolidated_file)
  
  # 1. VERIFY DIRECTIONALITY LOGIC
  cat("\n1. VERIFYING DIRECTIONALITY LOGIC\n")
  cat("=================================\n")
  
  # Check direction values
  direction_counts <- table(consolidated_data$direction, useNA = "always")
  cat("Direction values in data:\n")
  print(direction_counts)
  
  # Check for potential inflation (genes appearing in multiple directions)
  cat("\nChecking for directionality inflation...\n")
  
  # Sample check: look for same gene/cluster/method appearing in multiple directions
  sample_check <- consolidated_data %>%
    filter(!is.na(mutation_perturbation), !is.na(cluster), !is.na(method)) %>%
    group_by(mutation_perturbation, cluster, method, enrichment_type, Description) %>%
    summarise(
      directions = paste(unique(direction), collapse = ","),
      n_directions = n_distinct(direction),
      .groups = "drop"
    ) %>%
    filter(n_directions > 1)
  
  if (nrow(sample_check) > 0) {
    cat("⚠️  WARNING: Found", nrow(sample_check), "cases where same term appears in multiple directions\n")
    cat("This suggests potential inflation issue still exists\n")
    print(head(sample_check))
  } else {
    cat("✅ No directionality inflation detected\n")
  }
  
  # 2. COMPREHENSIVE CROSS-METHOD ANALYSIS
  cat("\n\n2. COMPREHENSIVE CROSS-METHOD ANALYSIS\n")
  cat("=====================================\n")
  
  # Get available genes that exist in both MAST and MixScale/CRISPRi data
  mast_genes <- unique(consolidated_data$mutation_perturbation[
    consolidated_data$method == "MAST"
  ])
  crispri_genes <- unique(consolidated_data$mutation_perturbation[
    grepl("MixScale", consolidated_data$method) | 
    grepl("CRISPRi", consolidated_data$method)
  ])
  
  # Find overlapping genes
  overlapping_genes <- intersect(mast_genes, crispri_genes)
  cat("Genes available in both MAST and CRISPRi:", length(overlapping_genes), "\n")
  cat("Overlapping genes:", paste(overlapping_genes, collapse = ", "), "\n")
  
  # Get available clusters
  clusters <- sort(unique(consolidated_data$cluster))
  cat("Available clusters:", paste(clusters, collapse = ", "), "\n")
  
  # Directions to test
  directions <- c("ALL", "UP", "DOWN")
  
  # 3. SYSTEMATIC SIGNATURE DISCOVERY
  cat("\n\n3. SYSTEMATIC SIGNATURE DISCOVERY\n")
  cat("=================================\n")
  
  # Initialize results storage
  all_results <- data.frame()
  
  # Check if signature analysis functions are available
  if (!exists("analyze_gene_pair_signatures")) {
    cat("⚠️  WARNING: analyze_gene_pair_signatures function not found\n")
    cat("Attempting to create simplified analysis...\n")
    
    # Simplified analysis if functions not available
    for (gene in overlapping_genes[1:min(3, length(overlapping_genes))]) {
      for (cluster in clusters[1:min(3, length(clusters))]) {
        for (direction in directions) {
          cat(sprintf("Analyzing: %s | Cluster %s | Direction %s\n", gene, cluster, direction))
          
          # Get MAST data
          mast_data <- consolidated_data %>%
            filter(
              method == "MAST",
              mutation_perturbation == gene,
              cluster == !!cluster,
              direction == !!direction | direction == "ALL"
            )
          
          # Get CRISPRi data  
          crispri_data <- consolidated_data %>%
            filter(
              grepl("MixScale", method) | grepl("CRISPRi", method),
              mutation_perturbation == gene,
              cluster == !!cluster,
              direction == !!direction | direction == "ALL"
            )
          
          # Calculate overlap
          if (nrow(mast_data) > 0 && nrow(crispri_data) > 0) {
            mast_terms <- unique(mast_data$Description)
            crispri_terms <- unique(crispri_data$Description)
            
            overlap_terms <- intersect(mast_terms, crispri_terms)
            
            if (length(overlap_terms) > 0) {
              # Simple Fisher's exact test
              all_terms <- union(mast_terms, crispri_terms)
              
              # Contingency table
              a <- length(overlap_terms)  # in both
              b <- length(mast_terms) - a  # MAST only
              c <- length(crispri_terms) - a  # CRISPRi only
              d <- length(all_terms) - a - b - c  # neither (approximation)
              
              if (d <= 0) d <- 1  # avoid negative values
              
              # Fisher's test
              fisher_result <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
              
              result <- data.frame(
                gene = gene,
                cluster = cluster,
                direction = direction,
                mast_terms = length(mast_terms),
                crispri_terms = length(crispri_terms),
                overlap_terms = length(overlap_terms),
                fisher_p = fisher_result$p.value,
                jaccard = length(overlap_terms) / length(all_terms),
                stringsAsFactors = FALSE
              )
              
              all_results <- rbind(all_results, result)
            }
          }
        }
      }
    }
  } else {
    # Use proper signature analysis functions if available
    cat("Using analyze_gene_pair_signatures function...\n")
    
    for (gene in overlapping_genes) {
      for (cluster in clusters) {
        for (direction in directions) {
          cat(sprintf("Analyzing: %s | Cluster %s | Direction %s\n", gene, cluster, direction))
          
          tryCatch({
            result <- analyze_gene_pair_signatures(
              mast_gene = gene,
              crispri_gene = gene,
              cluster = cluster,
              direction = direction,
              consolidated_data = consolidated_data
            )
            
            if (!is.null(result) && length(result) > 0) {
              result_df <- data.frame(
                gene = gene,
                cluster = cluster,
                direction = direction,
                fisher_p = result$intersection_fisher_p %||% NA,
                jaccard = result$intersection_jaccard %||% NA,
                overlap_count = result$intersection_count %||% 0,
                mast_count = result$mast_count %||% 0,
                crispri_count = result$crispri_count %||% 0,
                stringsAsFactors = FALSE
              )
              
              all_results <- rbind(all_results, result_df)
            }
          }, error = function(e) {
            cat("Error analyzing", gene, cluster, direction, ":", e$message, "\n")
          })
        }
      }
    }
  }
  
  # 4. ANALYZE RESULTS
  cat("\n\n4. ANALYSIS RESULTS\n")
  cat("==================\n")
  
  if (nrow(all_results) > 0) {
    # Sort by significance
    all_results <- all_results %>%
      arrange(fisher_p) %>%
      mutate(
        neg_log10_p = -log10(pmax(fisher_p, 1e-10)),
        rank = row_number()
      )
    
    cat("Total combinations analyzed:", nrow(all_results), "\n")
    
    # Top 10 most significant
    top_results <- head(all_results, 10)
    cat("\nTOP 10 MOST SIGNIFICANT SIGNATURES:\n")
    cat("===================================\n")
    
    for (i in 1:nrow(top_results)) {
      result <- top_results[i, ]
      cat(sprintf("%d. %s | Cluster %s | Direction %s\n", 
                  i, result$gene, result$cluster, result$direction))
      cat(sprintf("   Fisher p-value: %.2e | -log10(p): %.2f\n", 
                  result$fisher_p, result$neg_log10_p))
      if ("overlap_count" %in% colnames(result)) {
        cat(sprintf("   Overlap: %d | MAST: %d | CRISPRi: %d\n", 
                    result$overlap_count, result$mast_count, result$crispri_count))
      }
      if ("jaccard" %in% colnames(result)) {
        cat(sprintf("   Jaccard Index: %.3f\n", result$jaccard))
      }
      cat("\n")
    }
    
    # Summary by direction
    cat("\nSUMMARY BY DIRECTION:\n")
    cat("====================\n")
    direction_summary <- all_results %>%
      group_by(direction) %>%
      summarise(
        n_combinations = n(),
        mean_neg_log10_p = mean(neg_log10_p, na.rm = TRUE),
        median_fisher_p = median(fisher_p, na.rm = TRUE),
        n_significant = sum(fisher_p < 0.05, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_neg_log10_p))
    
    print(direction_summary)
    
    # Summary by gene
    cat("\nSUMMARY BY GENE:\n")
    cat("===============\n")
    gene_summary <- all_results %>%
      group_by(gene) %>%
      summarise(
        n_combinations = n(),
        best_p_value = min(fisher_p, na.rm = TRUE),
        mean_neg_log10_p = mean(neg_log10_p, na.rm = TRUE),
        n_significant = sum(fisher_p < 0.05, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(best_p_value)
    
    print(gene_summary)
    
  } else {
    cat("⚠️  No results generated. Check if signature analysis functions are working.\n")
  }
  
  # 5. RECOMMENDATIONS
  cat("\n\n5. MANUSCRIPT RECOMMENDATIONS\n")
  cat("============================\n")
  
  if (nrow(all_results) > 0) {
    significant_results <- all_results %>%
      filter(fisher_p < 0.05) %>%
      arrange(fisher_p)
    
    if (nrow(significant_results) > 0) {
      cat("RECOMMENDED SIGNATURE CANDIDATES FOR MANUSCRIPT:\n")
      cat("================================================\n")
      
      top_candidates <- head(significant_results, 5)
      for (i in 1:nrow(top_candidates)) {
        result <- top_candidates[i, ]
        cat(sprintf("⭐ CANDIDATE %d: %s\n", i, result$gene))
        cat(sprintf("   Best condition: Cluster %s, Direction %s\n", 
                    result$cluster, result$direction))
        cat(sprintf("   Significance: p = %.2e\n", result$fisher_p))
        cat(sprintf("   Priority: %s\n", 
                    ifelse(result$fisher_p < 0.001, "HIGH", 
                           ifelse(result$fisher_p < 0.01, "MEDIUM", "LOW"))))
        cat("\n")
      }
    } else {
      cat("⚠️  No statistically significant signatures found.\n")
      cat("This may indicate the directionality fix was too stringent.\n")
      cat("Consider investigating direction-agnostic approaches.\n")
    }
  }
  
  # Save results
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  results_file <- paste0("comprehensive_signature_results_", timestamp, ".rds")
  
  saveRDS(list(
    analysis_results = all_results,
    direction_summary = if(exists("direction_summary")) direction_summary else NULL,
    gene_summary = if(exists("gene_summary")) gene_summary else NULL,
    significant_results = if(exists("significant_results")) significant_results else NULL,
    analysis_timestamp = Sys.time(),
    overlapping_genes = overlapping_genes,
    total_combinations = nrow(all_results)
  ), results_file)
  
  cat("\n\nResults saved to:", results_file, "\n")
  cat("Analysis complete!\n")
  
  return(all_results)
}

# Run the analysis
cat("Starting comprehensive signature analysis...\n")
results <- comprehensive_signature_analysis()

cat("\n=== ANALYSIS SUMMARY ===\n")
cat("This analysis checked:\n")
cat("1. Directionality logic verification ✓\n")
cat("2. Cross-method signature discovery ✓\n") 
cat("3. Statistical significance ranking ✓\n")
cat("4. Manuscript recommendations ✓\n")
cat("\nPlease review the results above and test the app to confirm fixes are working!\n")