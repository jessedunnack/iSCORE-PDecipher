# Test Module Integration
# Quick test to ensure signature nomination module can be loaded and basic functions work

# Load required functions
source("R/gene_harmonization.R")
source("R/signature_analysis.R") 
source("R/manuscript_signature_discovery.R")

cat("=== Testing Module Integration ===\n")

# Test a minimal enrichment analysis workflow
cat("Testing signature discovery workflow...\n")

# Create dummy enrichment data structure
dummy_data <- data.frame(
  method = c("MAST", "MixScale"),
  mutation_perturbation = c("LRRK2", "LRRK2"),
  cluster = c("cluster_0", "cluster_0"),
  Description = c("mitochondrial function", "protein folding"),
  geneID = c("GENE1/GENE2", "GENE2/GENE3"),
  p.adjust = c(0.01, 0.02),
  stringsAsFactors = FALSE
)

cat("Created dummy data with", nrow(dummy_data), "rows\n")

# Test gene pair analysis
gene_pair <- list(mast_gene = "LRRK2", crispri_gene = "LRRK2")
result <- analyze_gene_pair_signatures(gene_pair, dummy_data, include_pathways = TRUE)

if ("gene_pair" %in% names(result)) {
  cat("Gene pair analysis successful\n")
  cat("Result structure:", paste(names(result), collapse = ", "), "\n")
} else {
  cat("Gene pair analysis returned unexpected format\n")
  cat("Result keys:", paste(names(result), collapse = ", "), "\n")
}

# Test gene mapping functions
gene_pairs <- get_comparable_gene_pairs(combine_snca_variants = TRUE)
cat("Found", nrow(gene_pairs), "comparable gene pairs\n")

cat("=== Module Integration Test Complete ===\n")