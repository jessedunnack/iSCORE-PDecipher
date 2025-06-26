# Test the PD analysis fix with real pan-cluster signature data
cat("🧪 TESTING PD ANALYSIS FIX\n")
cat("==========================\n\n")

# Load required functions
source("R/pd_signature_interpretation.R")
source("R/gene_harmonization.R")

# Create test pan-cluster signature that matches real data structure
test_signature <- data.frame(
  gene_pair = "FBXO7_vs_FBXO7",
  cluster_count = 11,
  mean_signature_strength = 1.5,
  max_signature_strength = 2.7,
  total_gene_overlaps = 983,
  total_pathway_overlaps = 5015,
  pan_cluster_rank = 1,
  stringsAsFactors = FALSE
)

cat("Test signature structure:\n")
print(str(test_signature))
cat("\n")

# Test the helper function directly
cat("Testing extract_signature_genes()...\n")
tryCatch({
  gene_info <- extract_signature_genes(test_signature[1, ])
  cat("✅ Gene extraction successful:\n")
  cat("  MAST gene:", gene_info$mast_gene, "\n")
  cat("  CRISPRi gene:", gene_info$crispri_gene, "\n")
}, error = function(e) {
  cat("❌ Gene extraction failed:", e$message, "\n")
})

# Test the enrichment term extraction with minimal data
cat("\nTesting get_signature_enrichment_terms()...\n")

# Load minimal enrichment data
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
cat("Loaded", nrow(enrichment_data), "enrichment terms\n")

# Test with MAST method
cat("\nTesting MAST method:\n")
tryCatch({
  mast_terms <- get_signature_enrichment_terms(
    test_signature[1, ], 
    enrichment_data, 
    method = "MAST"
  )
  cat("✅ MAST enrichment extraction successful:", nrow(mast_terms), "terms\n")
}, error = function(e) {
  cat("❌ MAST enrichment extraction failed:", e$message, "\n")
})

# Test with MixScale method  
cat("\nTesting MixScale method:\n")
tryCatch({
  crispri_terms <- get_signature_enrichment_terms(
    test_signature[1, ], 
    enrichment_data, 
    method = "MixScale"
  )
  cat("✅ MixScale enrichment extraction successful:", nrow(crispri_terms), "terms\n")
}, error = function(e) {
  cat("❌ MixScale enrichment extraction failed:", e$message, "\n")
})

cat("\n🎯 TEST COMPLETE\n")