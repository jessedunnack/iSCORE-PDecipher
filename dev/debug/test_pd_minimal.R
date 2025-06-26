# Minimal test to isolate the exact error location
cat("🔍 MINIMAL PD ANALYSIS TEST\n")
cat("============================\n\n")

# Load functions
source("R/pd_signature_interpretation.R")
source("R/gene_harmonization.R")

# Test just the PD pathways function
cat("Testing get_pd_relevant_pathways()...\n")
pd_pathways <- get_pd_relevant_pathways()
cat("✅ PD pathways loaded:", length(pd_pathways), "terms\n")
cat("First few pathways:", paste(head(pd_pathways, 3), collapse = ", "), "\n\n")

# Load minimal enrichment data
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
fbxo7_mast <- enrichment_data[
  enrichment_data$method == "MAST" & 
  enrichment_data$mutation_perturbation == "FBXO7",
][1:100, ]  # Just 100 terms for testing

cat("Testing filter_pd_relevant_terms()...\n")
tryCatch({
  pd_terms <- filter_pd_relevant_terms(fbxo7_mast, pd_pathways)
  cat("✅ PD term filtering works:", nrow(pd_terms), "terms found\n")
}, error = function(e) {
  cat("❌ PD term filtering failed:", e$message, "\n")
})

# Test with empty data
cat("\nTesting with empty input...\n")
tryCatch({
  empty_terms <- filter_pd_relevant_terms(data.frame(), pd_pathways)
  cat("✅ Empty input handling works\n")
}, error = function(e) {
  cat("❌ Empty input handling failed:", e$message, "\n")
})

cat("\n🎯 MINIMAL TEST COMPLETE\n")