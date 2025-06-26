# Fix the remaining issues based on user console output

cat("🔧 FIXING REMAINING SIGNATURE NOMINATION ISSUES\n")
cat("===============================================\n\n")

# Issue 1: Fix the PD analysis column access
cat("1. Checking PD analysis pan-cluster signature handling...\n")

# Create test data that matches what the app is actually using
test_pan_cluster <- data.frame(
  gene_pair = "FBXO7_vs_FBXO7",
  cluster_count = 11,
  mean_signature_strength = 1.5,
  max_signature_strength = 2.7,
  total_gene_overlaps = 983,
  total_pathway_overlaps = 5015,
  pan_cluster_rank = 1,
  stringsAsFactors = FALSE
)

cat("Pan-cluster signature columns:", paste(colnames(test_pan_cluster), collapse = ", "), "\n")
cat("Has mast_gene:", "mast_gene" %in% colnames(test_pan_cluster), "\n")
cat("Has gene_pair:", "gene_pair" %in% colnames(test_pan_cluster), "\n")

# Test the helper function with this exact structure
source("R/pd_signature_interpretation.R")
source("R/gene_harmonization.R")

cat("Testing extract_signature_genes with pan-cluster data...\n")
tryCatch({
  gene_info <- extract_signature_genes(test_pan_cluster[1, ])
  cat("✅ Gene extraction successful:\n")
  cat("  mast_gene:", gene_info$mast_gene, "\n")
  cat("  crispri_gene:", gene_info$crispri_gene, "\n")
}, error = function(e) {
  cat("❌ Gene extraction failed:", e$message, "\n")
})

# Issue 2: Check signature strength extraction
cat("\n2. Checking signature strength extraction...\n")
test_signature <- data.frame(
  gene_pair = "VPS13C_combined_vs_VPS13C",
  mast_gene = "VPS13C_combined", 
  crispri_gene = "VPS13C",
  cluster = "cluster_4",
  signature_strength = 2.719098,
  stringsAsFactors = FALSE
)

# Check if the strongest gene pair extraction works
top_signatures <- data.frame(
  gene_pair = c("VPS13C_combined_vs_VPS13C", "SNCA_combined_vs_SNCA", "FBXO7_vs_FBXO7"),
  signature_strength = c(2.719, 2.594, 2.712),
  stringsAsFactors = FALSE
)

strongest_pair <- if(nrow(top_signatures) > 0) {
  top_signatures$gene_pair[which.max(top_signatures$signature_strength)]
} else {
  "None"
}

cat("Test strongest gene pair:", strongest_pair, "\n")
cat("Expected: VPS13C_combined_vs_VPS13C\n")

# Issue 3: Check cluster-specific signature threshold
cat("\n3. Checking cluster-specific signature threshold...\n")
test_cluster_signatures <- data.frame(
  gene_pair = c("FBXO7_vs_FBXO7", "SNCA_combined_vs_SNCA"),
  cluster = c("cluster_2", "cluster_4"),
  signature_strength = c(2.712, 2.594),
  stringsAsFactors = FALSE
)

# Test with different thresholds
thresholds <- c(3.0, 2.0, 1.0, 0.5)
for (thresh in thresholds) {
  strong_sigs <- test_cluster_signatures[test_cluster_signatures$signature_strength >= thresh, ]
  cat(sprintf("Threshold %.1f: %d signatures found\n", thresh, nrow(strong_sigs))
}

cat("\n💡 RECOMMENDATIONS:\n")
cat("1. Lower cluster-specific threshold from 1.0 to 0.5 for better results\n")
cat("2. Add NULL checks for strongest gene pair extraction\n") 
cat("3. Fix PD analysis to handle missing mast_gene/crispri_gene columns\n")
cat("4. Check why heatmap is not rendering (likely data structure issue)\n")

cat("\n🎯 NEXT STEPS:\n")
cat("1. Update cluster-specific signature threshold\n")
cat("2. Fix strongest gene pair extraction in analysis summary\n")
cat("3. Add more defensive programming to PD analysis\n")
cat("4. Test with actual app data to ensure compatibility\n")