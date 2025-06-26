# Debug enrichment data structure
cat("🔍 ENRICHMENT DATA STRUCTURE DEBUG\n")
cat("===================================\n\n")

# Load data
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
cat("Total enrichment terms:", nrow(enrichment_data), "\n")

cat("\nColumn names and types:\n")
for (col in colnames(enrichment_data)) {
  sample_val <- enrichment_data[[col]][1]
  cat("  ", col, ":", class(sample_val), "\n")
  
  # Check if it's a list column
  if (is.list(sample_val)) {
    cat("    -> LIST COLUMN! First element:", sample_val[[1]], "\n")
  }
}

# Check FBXO7 MAST data specifically
fbxo7_mast <- enrichment_data[
  enrichment_data$method == "MAST" & 
  enrichment_data$mutation_perturbation == "FBXO7",
]
cat("\nFBXO7 MAST data:", nrow(fbxo7_mast), "terms\n")

if (nrow(fbxo7_mast) > 0) {
  cat("Sample p.adjust values:\n")
  padj_vals <- fbxo7_mast$p.adjust[1:5]
  for (i in 1:length(padj_vals)) {
    cat("  [", i, "]:", padj_vals[i], "(", class(padj_vals[i]), ")\n")
  }
  
  cat("\nTesting order function on p.adjust:\n")
  tryCatch({
    test_order <- order(fbxo7_mast$p.adjust[1:10])
    cat("✅ Order function works\n")
  }, error = function(e) {
    cat("❌ Order function failed:", e$message, "\n")
  })
}