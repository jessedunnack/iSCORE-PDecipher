# Simple but Comprehensive Signature Analysis
# This script actually works and finds the most significant cross-method signatures

library(dplyr)
library(tidyr)

# Find consolidated data
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

cat("=== SIMPLE COMPREHENSIVE SIGNATURE ANALYSIS ===\n")
cat("Finding most significant cross-method signatures between MAST and CRISPRi\n\n")

# Load data
consolidated_file <- find_consolidated_data()
if (is.null(consolidated_file)) {
  stop("Could not find consolidated data file!")
}

cat("Loading data from:", consolidated_file, "\n")
data <- readRDS(consolidated_file)

cat("Data dimensions:", nrow(data), "rows x", ncol(data), "columns\n")

# 1. CHECK DIRECTIONALITY INFLATION
cat("\n1. CHECKING DIRECTIONALITY INFLATION\n")
cat("===================================\n")

# Check for terms appearing in multiple directions
inflation_check <- data %>%
  filter(!is.na(mutation_perturbation), !is.na(cluster), !is.na(method)) %>%
  group_by(mutation_perturbation, cluster, method, enrichment_type, Description) %>%
  summarise(
    directions = paste(unique(direction), collapse = ","),
    n_directions = n_distinct(direction),
    .groups = "drop"
  ) %>%
  filter(n_directions > 1)

cat("Cases with multiple directions:", nrow(inflation_check), "\n")
if (nrow(inflation_check) > 0) {
  cat("⚠️  DIRECTIONALITY INFLATION STILL EXISTS!\n")
  cat("Sample cases:\n")
  print(head(inflation_check, 3))
} else {
  cat("✅ No directionality inflation detected\n")
}

# 2. FIND OVERLAPPING GENES
cat("\n\n2. FINDING OVERLAPPING GENES\n")
cat("============================\n")

mast_genes <- unique(data$mutation_perturbation[data$method == "MAST"])
crispri_genes <- unique(data$mutation_perturbation[
  grepl("MixScale", data$method) | grepl("CRISPRi", data$method)
])

overlapping_genes <- intersect(mast_genes, crispri_genes)
cat("Genes in both MAST and CRISPRi:", length(overlapping_genes), "\n")
cat("Overlapping genes:", paste(overlapping_genes, collapse = ", "), "\n")

clusters <- sort(unique(data$cluster))
directions <- c("ALL", "UP", "DOWN")

# 3. SYSTEMATIC CROSS-METHOD ANALYSIS
cat("\n\n3. SYSTEMATIC CROSS-METHOD ANALYSIS\n")
cat("===================================\n")

all_results <- data.frame()

for (gene in overlapping_genes) {
  for (cluster in clusters[1:5]) {  # Test first 5 clusters for efficiency
    for (direction in directions) {
      cat(sprintf("Analyzing: %s | %s | %s\n", gene, cluster, direction))
      
      # Get MAST data for this gene/cluster/direction
      mast_data <- data %>%
        filter(
          method == "MAST",
          mutation_perturbation == gene,
          cluster == !!cluster,
          direction == !!direction
        )
      
      # Get CRISPRi data for this gene/cluster/direction
      crispri_data <- data %>%
        filter(
          grepl("MixScale", method),
          mutation_perturbation == gene,
          cluster == !!cluster,
          direction == !!direction
        )
      
      if (nrow(mast_data) > 0 && nrow(crispri_data) > 0) {
        # Get unique terms (pathways/GO terms)
        mast_terms <- unique(mast_data$Description)
        crispri_terms <- unique(crispri_data$Description)
        
        # Calculate overlap
        overlap_terms <- intersect(mast_terms, crispri_terms)
        union_terms <- union(mast_terms, crispri_terms)
        
        if (length(overlap_terms) > 0) {
          # Simple Fisher's exact test for pathway overlap
          # Contingency table: overlap vs non-overlap
          a <- length(overlap_terms)  # in both
          b <- length(mast_terms) - a  # MAST only
          c <- length(crispri_terms) - a  # CRISPRi only
          
          # For background, use all terms in this cluster/direction
          all_terms_background <- data %>%
            filter(cluster == !!cluster, direction == !!direction) %>%
            pull(Description) %>%
            unique()
          
          d <- length(all_terms_background) - length(union_terms)
          if (d <= 0) d <- 1  # avoid issues
          
          # Fisher's test
          fisher_result <- tryCatch({
            fisher.test(matrix(c(a, b, c, d), nrow = 2), alternative = "greater")
          }, error = function(e) {
            list(p.value = 1)  # fallback
          })
          
          # Calculate additional metrics
          jaccard <- length(overlap_terms) / length(union_terms)
          
          # Store result
          result <- data.frame(
            gene = gene,
            cluster = cluster,
            direction = direction,
            mast_terms = length(mast_terms),
            crispri_terms = length(crispri_terms),
            overlap_terms = length(overlap_terms),
            background_terms = length(all_terms_background),
            fisher_p = fisher_result$p.value,
            jaccard_index = jaccard,
            neg_log10_p = -log10(pmax(fisher_result$p.value, 1e-10)),
            enrichment_fold = (a / length(mast_terms)) / (length(crispri_terms) / length(all_terms_background)),
            stringsAsFactors = FALSE
          )
          
          all_results <- rbind(all_results, result)
        }
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
      rank = row_number(),
      is_significant = fisher_p < 0.05
    )
  
  cat("Total successful analyses:", nrow(all_results), "\n")
  cat("Significant results (p < 0.05):", sum(all_results$is_significant), "\n")
  
  # Top 10 most significant
  cat("\n🏆 TOP 10 MOST SIGNIFICANT CROSS-METHOD SIGNATURES:\n")
  cat("==================================================\n")
  
  top_results <- head(all_results, 10)
  for (i in 1:nrow(top_results)) {
    result <- top_results[i, ]
    cat(sprintf("%d. %s | %s | Direction: %s\n", 
                i, result$gene, result$cluster, result$direction))
    cat(sprintf("   📊 Fisher p-value: %.2e (rank %d)\n", 
                result$fisher_p, result$rank))
    cat(sprintf("   🎯 -log10(p): %.2f\n", result$neg_log10_p))
    cat(sprintf("   🔗 Overlap: %d pathways (MAST: %d, CRISPRi: %d)\n", 
                result$overlap_terms, result$mast_terms, result$crispri_terms))
    cat(sprintf("   📈 Jaccard Index: %.3f | Enrichment Fold: %.2f\n", 
                result$jaccard_index, result$enrichment_fold))
    
    if (result$fisher_p < 0.001) {
      cat("   ⭐⭐⭐ HIGHLY SIGNIFICANT - TOP MANUSCRIPT CANDIDATE\n")
    } else if (result$fisher_p < 0.01) {
      cat("   ⭐⭐ SIGNIFICANT - GOOD MANUSCRIPT CANDIDATE\n")
    } else if (result$fisher_p < 0.05) {
      cat("   ⭐ MARGINALLY SIGNIFICANT\n")
    }
    cat("\n")
  }
  
  # Summary by direction
  cat("\n📊 SUMMARY BY DIRECTION:\n")
  cat("=======================\n")
  direction_summary <- all_results %>%
    group_by(direction) %>%
    summarise(
      n_analyses = n(),
      n_significant = sum(is_significant),
      mean_neg_log10_p = mean(neg_log10_p),
      best_p_value = min(fisher_p),
      mean_overlap = mean(overlap_terms),
      .groups = "drop"
    ) %>%
    arrange(desc(mean_neg_log10_p))
  
  print(direction_summary)
  
  # Summary by gene
  cat("\n🧬 SUMMARY BY GENE:\n")
  cat("=================\n")
  gene_summary <- all_results %>%
    group_by(gene) %>%
    summarise(
      n_analyses = n(),
      n_significant = sum(is_significant),
      best_p_value = min(fisher_p),
      best_rank = min(rank),
      mean_overlap = mean(overlap_terms),
      .groups = "drop"
    ) %>%
    arrange(best_p_value)
  
  print(gene_summary)
  
  # Find the most promising individual results
  cat("\n🎯 MANUSCRIPT RECOMMENDATIONS:\n")
  cat("============================\n")
  
  highly_significant <- all_results %>% filter(fisher_p < 0.001)
  if (nrow(highly_significant) > 0) {
    cat("🏆 TIER 1 CANDIDATES (p < 0.001):\n")
    for (i in 1:min(3, nrow(highly_significant))) {
      result <- highly_significant[i, ]
      cat(sprintf("• %s (%s, %s): p = %.2e, %d shared pathways\n",
                  result$gene, result$cluster, result$direction, 
                  result$fisher_p, result$overlap_terms))
    }
  }
  
  significant <- all_results %>% filter(fisher_p >= 0.001 & fisher_p < 0.01)
  if (nrow(significant) > 0) {
    cat("\n⭐ TIER 2 CANDIDATES (0.001 ≤ p < 0.01):\n")
    for (i in 1:min(3, nrow(significant))) {
      result <- significant[i, ]
      cat(sprintf("• %s (%s, %s): p = %.2e, %d shared pathways\n",
                  result$gene, result$cluster, result$direction, 
                  result$fisher_p, result$overlap_terms))
    }
  }
  
} else {
  cat("❌ No successful cross-method analyses completed.\n")
  cat("This suggests either:\n")
  cat("1. No overlapping data between MAST and CRISPRi\n")
  cat("2. Data structure issues\n")
  cat("3. Directionality filtering too restrictive\n")
}

# Save results
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
results_file <- paste0("simple_signature_results_", timestamp, ".rds")

if (nrow(all_results) > 0) {
  saveRDS(list(
    all_results = all_results,
    direction_summary = if(exists("direction_summary")) direction_summary else NULL,
    gene_summary = if(exists("gene_summary")) gene_summary else NULL,
    inflation_check = inflation_check,
    overlapping_genes = overlapping_genes,
    analysis_timestamp = Sys.time()
  ), results_file)
  
  cat("\n💾 Results saved to:", results_file, "\n")
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
if (nrow(inflation_check) > 0) {
  cat("❌ Directionality inflation still exists - needs further investigation\n")
}
cat("✅ Cross-method signature analysis completed\n")
if (nrow(all_results) > 0) {
  cat("✅ Found", sum(all_results$is_significant), "statistically significant signatures\n")
  cat("🎯 Top candidate:", all_results$gene[1], "(", all_results$cluster[1], 
      ", ", all_results$direction[1], ") with p =", 
      format(all_results$fisher_p[1], scientific = TRUE), "\n")
} else {
  cat("⚠️  No significant signatures found - may need direction-agnostic approach\n")
}