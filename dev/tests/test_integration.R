# Test Gene Association Integration
# Quick test to verify everything works together

library(iSCORE.PDecipher)

cat("=== Integration Test ===\n")

# Test package functions
cat("1. Testing package functions...\n")
result <- load_gene_associations()
if (!is.null(result)) {
  cat("   Package functions: WORKING\n")
} else {
  cat("   Package functions: FAILED\n")
}

# Test data availability
available <- gene_associations_available()
cat("2. Gene associations available:", available, "\n")

if (available) {
  # Test lookup
  stats <- get_association_stats()
  cat("3. Dataset stats:\n")
  cat("   Total associations:", stats$total_associations, "\n")
  cat("   Unique terms:", stats$unique_terms, "\n")
  cat("   Unique genes:", stats$unique_genes, "\n")
  
  # Test search
  results <- search_gene_associations("mitochondrial")
  cat("4. Search test: Found", nrow(results), "mitochondrial terms\n")
  
  cat("\n✅ INTEGRATION SUCCESSFUL\n")
  cat("Ready for Shiny app testing!\n")
} else {
  cat("\n❌ INTEGRATION FAILED\n")
  cat("Gene associations not available\n")
}

