# Test Script for Signature Analysis Functions
# Quick validation that core functions are working properly

# Load required functions
source("R/gene_harmonization.R")
source("R/signature_analysis.R") 
source("R/manuscript_signature_discovery.R")

cat("=== Testing Gene Harmonization Functions ===\n")

# Test gene mapping
gene_mapping <- create_gene_mapping_table()
cat("Gene mapping table created with", nrow(gene_mapping), "entries\n")
print(head(gene_mapping))

# Test comparable gene pairs
comparable_pairs <- get_comparable_gene_pairs(combine_snca_variants = TRUE, 
                                             combine_vps13c_variants = TRUE)
cat("\nComparable gene pairs (variants combined):", nrow(comparable_pairs), "pairs\n")
print(comparable_pairs[, c("mast_gene", "crispri_gene", "comparison_type")])

# Test mutation categories
mutation_cats <- get_mutation_categories()
cat("\nMutation categories:", nrow(mutation_cats), "genes categorized\n")
print(table(mutation_cats$mutation_category))

cat("\n=== Testing Signature Analysis Functions ===\n")

# Test overlap analysis with dummy data
mast_genes <- c("GENE1", "GENE2", "GENE3", "GENE4")
crispri_genes <- c("GENE2", "GENE3", "GENE5", "GENE6")
background_genes <- c("GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7", "GENE8")

overlap_result <- calculate_gene_overlap_significance(mast_genes, crispri_genes, background_genes)
cat("Overlap analysis test:\n")
cat("- Overlap count:", overlap_result$overlap_count, "\n")
cat("- Jaccard index:", round(overlap_result$jaccard_index, 3), "\n")
cat("- Fisher p-value:", if(is.na(overlap_result$fisher_p)) "NA" else format(overlap_result$fisher_p, digits = 3), "\n")

# Test effect size correlation with dummy data
mast_data <- data.frame(
  gene_name = c("GENE1", "GENE2", "GENE3"),
  log2FC = c(1.5, -0.8, 2.1),
  pvalue = c(0.01, 0.05, 0.001)
)

crispri_data <- data.frame(
  gene_name = c("GENE1", "GENE2", "GENE4"),
  log2FC = c(1.2, -0.9, 1.8),
  pvalue = c(0.02, 0.03, 0.01)
)

correlation_result <- calculate_effect_size_correlation(mast_data, crispri_data)
cat("\nCorrelation analysis test:\n")
cat("- Correlation:", if(is.na(correlation_result$correlation)) "NA" else round(correlation_result$correlation, 3), "\n")
cat("- N genes:", correlation_result$n_genes, "\n")

# Test direction consistency
direction_result <- calculate_direction_consistency(mast_data, crispri_data)
cat("\nDirection consistency test:\n")
cat("- Same direction up:", direction_result$same_direction_up, "\n")
cat("- Same direction down:", direction_result$same_direction_down, "\n")
cat("- Opposite direction:", direction_result$opposite_direction, "\n")
cat("- Consistency %:", round(direction_result$consistency_percent, 1), "\n")

# Test composite score
composite_score <- calculate_composite_signature_score(overlap_result, correlation_result, direction_result)
cat("\nComposite signature score:", round(composite_score$composite_score, 2), "\n")

cat("\n=== Testing PD Pathway Functions ===\n")

# Test PD pathway identification
pd_pathways <- get_pd_relevant_pathways()
cat("PD-relevant pathway terms loaded:", length(pd_pathways), "terms\n")
cat("Sample terms:", paste(head(pd_pathways, 5), collapse = ", "), "\n")

cat("\n=== All Function Tests Completed Successfully! ===\n")
cat("Ready to proceed with implementation.\n")