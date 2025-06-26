# Test the complete PD analysis with the fixes
cat("🔬 TESTING COMPLETE PD ANALYSIS WITH FIXES\n")
cat("=========================================\n\n")

# Load required functions
source("R/pd_signature_interpretation.R")
source("R/gene_harmonization.R")
source("R/manuscript_signature_discovery.R")

# Create test signature data that matches the app's real structure
test_signatures <- data.frame(
  gene_pair = c("FBXO7_vs_FBXO7", "SNCA_combined_vs_SNCA"),
  cluster_count = c(11, 9),
  mean_signature_strength = c(1.5, 2.2),
  max_signature_strength = c(2.7, 2.6),
  total_gene_overlaps = c(983, 745),
  total_pathway_overlaps = c(5015, 4327),
  pan_cluster_rank = c(1, 2),
  stringsAsFactors = FALSE
)

# Load enrichment data
cat("Loading enrichment data...\n")
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
cat("Loaded", nrow(enrichment_data), "enrichment terms\n\n")

# Create test signature results structure
signature_results <- list(
  pan_cluster_signatures = test_signatures,
  all_signatures = data.frame()  # Empty for this test
)

# Test the full PD analysis
cat("Testing analyze_pd_signatures()...\n")
tryCatch({
  pd_analysis <- analyze_pd_signatures(
    signature_results = signature_results,
    enrichment_data = enrichment_data,
    focus_on_pan_cluster = TRUE
  )
  
  cat("✅ PD analysis completed successfully!\n")
  cat("Analysis type:", pd_analysis$analysis_type, "\n")
  cat("Signatures analyzed:", pd_analysis$signature_count, "\n")
  cat("Enhanced signatures created:", length(pd_analysis$enhanced_signatures), "\n")
  
  if (length(pd_analysis$enhanced_signatures) > 0) {
    cat("\nFirst signature summary:\n")
    first_sig <- pd_analysis$enhanced_signatures[[1]]
    cat("  Gene pair:", first_sig$signature$gene_pair, "\n")
    cat("  MAST PD terms:", nrow(first_sig$mast_pd_terms), "\n")
    cat("  CRISPRi PD terms:", nrow(first_sig$crispri_pd_terms), "\n")
    cat("  Shared pathways:", nrow(first_sig$shared_pd_pathways), "\n")
    cat("  PD relevance score:", first_sig$pd_relevance_score, "\n")
  }
  
}, error = function(e) {
  cat("❌ PD analysis failed:", e$message, "\n")
  cat("Traceback:\n")
  print(traceback())
})

cat("\n🎯 COMPLETE TEST FINISHED\n")