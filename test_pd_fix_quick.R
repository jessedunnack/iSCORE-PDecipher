# Quick test of the PD analysis fix
cat("🧪 QUICK PD ANALYSIS FIX TEST\n")
cat("=============================\n\n")

# Load functions
source("R/signature_analysis.R")
source("R/manuscript_signature_discovery.R") 
source("R/gene_harmonization.R")
source("R/pd_signature_interpretation.R")

# Load small test data
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
test_data <- enrichment_data[enrichment_data$method %in% c("MAST", "MixScale"), ]

cat("✅ Data loaded:", nrow(test_data), "enrichment terms\n\n")

# Create minimal test signature results
test_signature_results <- list(
  all_signatures = data.frame(
    gene_pair = "FBXO7_vs_FBXO7",
    mast_gene = "FBXO7", 
    crispri_gene = "FBXO7",
    cluster = "cluster_0",
    signature_strength = 2.5,
    stringsAsFactors = FALSE
  ),
  pan_cluster_signatures = data.frame(
    gene_pair = "FBXO7_vs_FBXO7",
    cluster_count = 10,
    mean_signature_strength = 2.0,
    stringsAsFactors = FALSE
  )
)

cat("Test data structure created\n")
cat("Individual signature columns:", paste(colnames(test_signature_results$all_signatures), collapse = ", "), "\n")
cat("Pan-cluster signature columns:", paste(colnames(test_signature_results$pan_cluster_signatures), collapse = ", "), "\n\n")

# Test helper function
cat("🔧 TESTING HELPER FUNCTION\n")
cat("===========================\n")

# Test with individual signature
individual_sig <- test_signature_results$all_signatures[1, ]
gene_info1 <- extract_signature_genes(individual_sig)
cat("Individual signature extraction:\n")
cat("  mast_gene:", gene_info1$mast_gene, "\n")
cat("  crispri_gene:", gene_info1$crispri_gene, "\n\n")

# Test with pan-cluster signature  
pan_sig <- test_signature_results$pan_cluster_signatures[1, ]
gene_info2 <- extract_signature_genes(pan_sig)
cat("Pan-cluster signature extraction:\n")
cat("  mast_gene:", gene_info2$mast_gene, "\n")
cat("  crispri_gene:", gene_info2$crispri_gene, "\n\n")

# Test PD analysis with pan-cluster focus
cat("🩺 TESTING PD ANALYSIS\n")
cat("======================\n")

tryCatch({
  pd_analysis <- analyze_pd_signatures(
    signature_results = test_signature_results,
    enrichment_data = test_data,
    focus_on_pan_cluster = TRUE
  )
  
  cat("✅ PD analysis completed successfully!\n")
  cat("Enhanced signatures:", length(pd_analysis$enhanced_signatures), "\n")
  cat("Analysis type:", pd_analysis$analysis_type, "\n")
  
  if (length(pd_analysis$enhanced_signatures) > 0) {
    first_enhanced <- pd_analysis$enhanced_signatures[[1]]
    cat("PD relevance score:", first_enhanced$pd_relevance_score, "\n")
    cat("Shared PD pathways:", nrow(first_enhanced$shared_pd_pathways), "\n")
  }
  
}, error = function(e) {
  cat("❌ PD analysis failed with error:", e$message, "\n")
})

cat("\n🎯 QUICK FIX TEST COMPLETE\n")
cat("===========================\n")