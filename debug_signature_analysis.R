# Debug script to understand what's happening with your signature analysis
# This will help identify why the "strongest signature" is showing as "None"

# Load required functions
source("R/signature_analysis.R")
source("R/manuscript_signature_discovery.R") 
source("R/gene_harmonization.R")

# Load your data - check multiple possible locations
data_paths <- c(
  "all_enrichment_padj005_complete_with_direction.rds",  # Current directory
  "../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds",  # Dataset 2
  "../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"  # Alternative
)

data_file <- NULL
for (path in data_paths) {
  if (file.exists(path)) {
    data_file <- path
    break
  }
}

if (is.null(data_file)) {
  cat("❌ Data file not found in expected locations:\n")
  for (path in data_paths) {
    cat("  ", path, "\n")
  }
  cat("\n💡 Please navigate to the correct directory or specify the path:\n")
  cat("setwd('/path/to/iSCORE-PD_plus_CRISPRi')\n")
  cat("# OR\n")
  cat("enrichment_data <- readRDS('/full/path/to/all_enrichment_padj005_complete_with_direction.rds')\n")
  stop("Data file not found")
}

enrichment_data <- readRDS(data_file)
cat("✅ Data loaded from:", data_file, "\n")
cat("✅ Data size:", nrow(enrichment_data), "enrichment terms\n\n")

# Test with a small subset to understand the issue
cat("🔍 DEBUGGING SIGNATURE ANALYSIS\n")
cat(paste(rep("=", 50), collapse=""), "\n\n")

# Filter data similar to what the app does
test_data <- enrichment_data[enrichment_data$method %in% c("MAST", "MixScale"), ]
cat("Filtered data (MAST + MixScale):", nrow(test_data), "terms\n")

# Check available gene pairs
gene_pairs <- get_comparable_gene_pairs(
  combine_snca_variants = TRUE,
  combine_vps13c_variants = TRUE,
  include_mast_only = FALSE
)
cat("Available gene pairs:", nrow(gene_pairs), "\n")
print(gene_pairs[, c("mast_gene", "crispri_gene")])

# Run discovery with lower thresholds
cat("\n🧪 RUNNING TEST ANALYSIS\n")
cat(paste(rep("=", 30), collapse=""), "\n")

signature_results <- discover_top_signatures(
  enrichment_data = test_data,
  top_n = 15,
  min_cluster_breadth = 6,  # Lower threshold
  combine_variants = TRUE,
  progress_callback = function(msg, value = NULL, detail = NULL) {
    cat("  ", msg, "\n")
  }
)

cat("\n📊 ANALYSIS RESULTS STRUCTURE\n")
cat(paste(rep("=", 35), collapse=""), "\n")

cat("Results components:", paste(names(signature_results), collapse = ", "), "\n\n")

# Check top signatures
if ("all_signatures" %in% names(signature_results)) {
  all_sigs <- signature_results$all_signatures
  cat("All signatures found:", nrow(all_sigs), "\n")
  
  if (nrow(all_sigs) > 0) {
    cat("Signature strength range:", 
        round(min(all_sigs$signature_strength, na.rm = TRUE), 3), "to",
        round(max(all_sigs$signature_strength, na.rm = TRUE), 3), "\n")
    
    cat("\nTop 5 signatures:\n")
    top_5 <- head(all_sigs[order(-all_sigs$signature_strength), ], 5)
    print(top_5[, c("gene_pair", "cluster", "signature_strength", "gene_overlap_count")])
  }
} else {
  cat("❌ No 'all_signatures' component found\n")
}

# Check analysis summary
if ("analysis_summary" %in% names(signature_results)) {
  summary_stats <- signature_results$analysis_summary
  cat("\n📈 ANALYSIS SUMMARY\n")
  cat(paste(rep("=", 20), collapse=""), "\n")
  
  for (stat_name in names(summary_stats)) {
    cat(paste0(stat_name, ": ", summary_stats[[stat_name]]), "\n")
  }
} else {
  cat("❌ No 'analysis_summary' component found\n")
}

# Check pan-cluster signatures
if ("pan_cluster_signatures" %in% names(signature_results)) {
  pan_cluster <- signature_results$pan_cluster_signatures
  cat("\n🌐 PAN-CLUSTER SIGNATURES\n")
  cat(paste(rep("=", 25), collapse=""), "\n")
  cat("Pan-cluster signatures found:", nrow(pan_cluster), "\n")
  
  if (nrow(pan_cluster) > 0) {
    print(pan_cluster[, c("gene_pair", "cluster_count", "mean_signature_strength")])
  }
} else {
  cat("❌ No 'pan_cluster_signatures' component found\n")
}

# Check cluster-specific signatures  
if ("cluster_specific_signatures" %in% names(signature_results)) {
  cluster_specific <- signature_results$cluster_specific_signatures
  cat("\n🎯 CLUSTER-SPECIFIC SIGNATURES\n")
  cat(paste(rep("=", 30), collapse=""), "\n")
  cat("Cluster-specific signature groups:", length(cluster_specific), "\n")
  
  if (length(cluster_specific) > 0) {
    for (cluster in names(cluster_specific)) {
      cat(paste0("  ", cluster, ": ", nrow(cluster_specific[[cluster]]), " signatures\n"))
    }
  }
} else {
  cat("❌ No 'cluster_specific_signatures' component found\n")
}

cat("\n🔧 RECOMMENDATIONS\n")
cat(paste(rep("=", 15), collapse=""), "\n")

if (exists("all_sigs") && nrow(all_sigs) > 0) {
  max_strength <- max(all_sigs$signature_strength, na.rm = TRUE)
  if (max_strength < 2.0) {
    cat("✅ Your signatures have modest but meaningful strength (", round(max_strength, 2), ")\n")
    cat("✅ This is normal for cross-method comparisons\n")
    cat("✅ Focus on relative ranking rather than absolute thresholds\n")
  } else {
    cat("✅ You have strong signatures (max:", round(max_strength, 2), ")\n")
  }
} else {
  cat("❌ No signatures found - check data filtering\n")
}

cat("\n💡 Next steps:\n")
cat("1. Install updated package with: remotes::install_github('jessedunnack/iSCORE-PDecipher', force = TRUE)\n")
cat("2. Lower cluster-specific threshold in app settings\n")
cat("3. Focus on pan-cluster results which show meaningful patterns\n")

cat("\n", paste(rep("=", 50), collapse=""), "\n")
cat("DEBUG ANALYSIS COMPLETE\n")
cat(paste(rep("=", 50), collapse=""), "\n")