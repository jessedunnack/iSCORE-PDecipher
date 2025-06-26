# Test the PD analysis fix for column name mismatch
cat("🧪 TESTING PD ANALYSIS FIX\n")
cat("==========================\n\n")

# Load functions
source("R/signature_analysis.R")
source("R/manuscript_signature_discovery.R") 
source("R/gene_harmonization.R")
source("R/pd_signature_interpretation.R")

# Load data
data_paths <- c(
  "all_enrichment_padj005_complete_with_direction.rds",
  "../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds",
  "../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
)

data_file <- NULL
for (path in data_paths) {
  if (file.exists(path)) {
    data_file <- path
    break
  }
}

if (is.null(data_file)) {
  stop("Data file not found")
}

enrichment_data <- readRDS(data_file)
test_data <- enrichment_data[enrichment_data$method %in% c("MAST", "MixScale"), ]

cat("✅ Data loaded:", nrow(test_data), "enrichment terms\n\n")

# Run signature discovery
cat("Running quick signature discovery...\n")
signature_results <- discover_top_signatures(
  enrichment_data = test_data,
  top_n = 5,
  min_cluster_breadth = 6,
  combine_variants = TRUE
)

cat("✅ Signature discovery complete\n")
cat("All signatures found:", nrow(signature_results$all_signatures), "\n")
cat("Pan-cluster signatures found:", nrow(signature_results$pan_cluster_signatures), "\n\n")

# Test the helper function
cat("🔧 TESTING HELPER FUNCTION\n")
cat("===========================\n")

if (nrow(signature_results$all_signatures) > 0) {
  individual_sig <- signature_results$all_signatures[1, ]
  cat("Individual signature columns:", paste(colnames(individual_sig), collapse = ", "), "\n")
  
  gene_info1 <- extract_signature_genes(individual_sig)
  cat("Extracted from individual signature:\n")
  cat("  mast_gene:", gene_info1$mast_gene, "\n")
  cat("  crispri_gene:", gene_info1$crispri_gene, "\n\n")
}

if (nrow(signature_results$pan_cluster_signatures) > 0) {
  pan_sig <- signature_results$pan_cluster_signatures[1, ]
  cat("Pan-cluster signature columns:", paste(colnames(pan_sig), collapse = ", "), "\n")
  
  gene_info2 <- extract_signature_genes(pan_sig)
  cat("Extracted from pan-cluster signature:\n")
  cat("  mast_gene:", gene_info2$mast_gene, "\n")
  cat("  crispri_gene:", gene_info2$crispri_gene, "\n\n")
}

# Test the actual PD analysis
cat("🩺 TESTING PD ANALYSIS\n")
cat("======================\n")

tryCatch({
  pd_analysis <- analyze_pd_signatures(
    signature_results = signature_results,
    enrichment_data = test_data,
    focus_on_pan_cluster = TRUE
  )
  
  cat("✅ PD analysis completed successfully!\n")
  cat("Enhanced signatures:", length(pd_analysis$enhanced_signatures), "\n")
  cat("Analysis type:", pd_analysis$analysis_type, "\n")
  
  if (length(pd_analysis$enhanced_signatures) > 0) {
    first_enhanced <- pd_analysis$enhanced_signatures[[1]]
    cat("\nFirst enhanced signature keys:", paste(names(first_enhanced), collapse = ", "), "\n")
    cat("PD relevance score:", first_enhanced$pd_relevance_score, "\n")
  }
  
}, error = function(e) {
  cat("❌ PD analysis failed with error:", e$message, "\n")
  cat("Full error:", toString(e), "\n")
})

cat("\n🎯 FIX TEST COMPLETE\n")
cat("=====================\n")