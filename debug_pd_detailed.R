# Debug the specific PD analysis issue
cat("🔍 DETAILED PD ANALYSIS DEBUG\n")
cat("==============================\n\n")

# Load functions
source("R/signature_analysis.R")
source("R/manuscript_signature_discovery.R") 
source("R/gene_harmonization.R")
source("R/pd_signature_interpretation.R")

# Create test signature with all necessary columns
test_signature_results <- list(
  pan_cluster_signatures = data.frame(
    gene_pair = "FBXO7_vs_FBXO7",
    cluster_count = 10,
    mean_signature_strength = 2.0,
    max_signature_strength = 2.5,
    total_gene_overlaps = 100,
    total_pathway_overlaps = 200,
    pan_cluster_rank = 1,
    stringsAsFactors = FALSE
  )
)

cat("Pan-cluster signature structure:\n")
for (col in colnames(test_signature_results$pan_cluster_signatures)) {
  val <- test_signature_results$pan_cluster_signatures[[col]][1]
  cat("  ", col, ":", val, "(", class(val), ")\n")
}

# Test helper function with detailed debug
cat("\n🔧 TESTING HELPER FUNCTION\n")
pan_sig <- test_signature_results$pan_cluster_signatures[1, ]
cat("Available columns:", paste(colnames(pan_sig), collapse = ", "), "\n")

gene_info <- extract_signature_genes(pan_sig)
cat("Extracted genes:\n")
cat("  mast_gene:", gene_info$mast_gene, "\n")
cat("  crispri_gene:", gene_info$crispri_gene, "\n")

# Test what happens when we try to get enrichment terms
cat("\n🧪 TESTING ENRICHMENT TERM EXTRACTION\n")
enrichment_data <- readRDS("../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds")
test_data <- enrichment_data[enrichment_data$method %in% c("MAST", "MixScale"), ]

cat("Enrichment data loaded:", nrow(test_data), "terms\n")

# The issue: pan-cluster signatures don't have a cluster column!
cat("Pan-cluster signature cluster column exists:", "cluster" %in% colnames(pan_sig), "\n")

# This is the problem - we need to handle the missing cluster
cat("\n💡 SOLUTION NEEDED:\n")
cat("Pan-cluster signatures don't have individual cluster info\n")
cat("We need to either:\n")
cat("1. Skip cluster-specific filtering for pan-cluster analysis\n")
cat("2. Use a representative cluster or aggregate across all clusters\n")
cat("3. Modify the PD analysis to handle pan-cluster differently\n")