# Comprehensive debug script for PD analysis column mismatch issue
# This writes detailed output to a file to capture everything

# Create output file
debug_file <- "pd_analysis_debug_output.txt"
cat("🔍 COMPREHENSIVE PD ANALYSIS DEBUG\n", file = debug_file)
cat("=====================================\n\n", file = debug_file, append = TRUE)

# Capture all console output to file
sink(debug_file, append = TRUE, split = TRUE)

cat("Debug started at:", as.character(Sys.time()), "\n\n")

# Load required functions
tryCatch({
  source("R/signature_analysis.R")
  source("R/manuscript_signature_discovery.R") 
  source("R/gene_harmonization.R")
  source("R/pd_signature_interpretation.R")
  cat("✅ All R functions loaded successfully\n\n")
}, error = function(e) {
  cat("❌ Error loading functions:", e$message, "\n\n")
})

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
  cat("❌ Data file not found in expected locations\n")
  stop("Data file not found")
}

enrichment_data <- readRDS(data_file)
cat("✅ Data loaded from:", data_file, "\n")
cat("✅ Data size:", nrow(enrichment_data), "enrichment terms\n\n")

# Filter data for MAST + MixScale only
test_data <- enrichment_data[enrichment_data$method %in% c("MAST", "MixScale"), ]
cat("Filtered data (MAST + MixScale):", nrow(test_data), "terms\n\n")

# Run signature discovery with minimal parameters
cat("🧪 RUNNING SIGNATURE DISCOVERY\n")
cat("==============================\n")

signature_results <- discover_top_signatures(
  enrichment_data = test_data,
  top_n = 10,
  min_cluster_breadth = 6,
  combine_variants = TRUE,
  progress_callback = function(msg, value = NULL, detail = NULL) {
    cat("  ", msg, "\n")
  }
)

cat("\n📊 SIGNATURE RESULTS STRUCTURE\n")
cat("===============================\n")

# Examine all_signatures structure
if ("all_signatures" %in% names(signature_results)) {
  all_sigs <- signature_results$all_signatures
  cat("All signatures found:", nrow(all_sigs), "\n")
  cat("All signatures columns:", paste(colnames(all_sigs), collapse = ", "), "\n")
  
  if (nrow(all_sigs) > 0) {
    cat("\nFirst signature row structure:\n")
    first_sig <- all_sigs[1, ]
    for (col in colnames(first_sig)) {
      cat("  ", col, ":", first_sig[[col]], "\n")
    }
  }
}

# Examine pan_cluster_signatures structure  
if ("pan_cluster_signatures" %in% names(signature_results)) {
  pan_cluster <- signature_results$pan_cluster_signatures
  cat("\nPan-cluster signatures found:", nrow(pan_cluster), "\n")
  cat("Pan-cluster columns:", paste(colnames(pan_cluster), collapse = ", "), "\n")
  
  if (nrow(pan_cluster) > 0) {
    cat("\nFirst pan-cluster signature structure:\n")
    first_pan <- pan_cluster[1, ]
    for (col in colnames(first_pan)) {
      cat("  ", col, ":", first_pan[[col]], "\n")
    }
  }
}

cat("\n🔬 TESTING PD ANALYSIS WITH DIFFERENT DATA TYPES\n")
cat("=================================================\n")

# Test 1: Try PD analysis with all_signatures (should work)
if ("all_signatures" %in% names(signature_results) && nrow(signature_results$all_signatures) > 0) {
  cat("\nTest 1: PD analysis with individual signatures\n")
  cat("----------------------------------------------\n")
  
  test_signature <- signature_results$all_signatures[1, ]
  cat("Test signature columns:", paste(colnames(test_signature), collapse = ", "), "\n")
  cat("Has mast_gene column:", "mast_gene" %in% colnames(test_signature), "\n")
  cat("Has crispri_gene column:", "crispri_gene" %in% colnames(test_signature), "\n")
  
  if ("mast_gene" %in% colnames(test_signature)) {
    cat("mast_gene value:", test_signature$mast_gene, "\n")
    cat("crispri_gene value:", test_signature$crispri_gene, "\n")
  }
}

# Test 2: Try PD analysis with pan_cluster_signatures (will likely fail)
if ("pan_cluster_signatures" %in% names(signature_results) && nrow(signature_results$pan_cluster_signatures) > 0) {
  cat("\nTest 2: PD analysis with pan-cluster signatures\n")
  cat("-----------------------------------------------\n")
  
  test_pan_signature <- signature_results$pan_cluster_signatures[1, ]
  cat("Pan-cluster signature columns:", paste(colnames(test_pan_signature), collapse = ", "), "\n")
  cat("Has mast_gene column:", "mast_gene" %in% colnames(test_pan_signature), "\n")
  cat("Has crispri_gene column:", "crispri_gene" %in% colnames(test_pan_signature), "\n")
  cat("Has gene_pair column:", "gene_pair" %in% colnames(test_pan_signature), "\n")
  
  if ("gene_pair" %in% colnames(test_pan_signature)) {
    cat("gene_pair value:", test_pan_signature$gene_pair, "\n")
    
    # Test parsing gene_pair
    gene_pair_parts <- strsplit(test_pan_signature$gene_pair, "_vs_")[[1]]
    if (length(gene_pair_parts) == 2) {
      cat("Parsed mast_gene:", gene_pair_parts[1], "\n")
      cat("Parsed crispri_gene:", gene_pair_parts[2], "\n")
    }
  }
}

cat("\n🧪 ATTEMPTING ACTUAL PD ANALYSIS\n")
cat("================================\n")

# Try running the actual PD analysis that's failing
if (exists("analyze_pd_signatures") && "pan_cluster_signatures" %in% names(signature_results)) {
  cat("Attempting analyze_pd_signatures with pan-cluster data...\n")
  
  tryCatch({
    pd_analysis <- analyze_pd_signatures(
      signature_results = signature_results,
      enrichment_data = test_data,
      focus_on_pan_cluster = TRUE
    )
    cat("✅ PD analysis completed successfully!\n")
    cat("PD analysis results:", length(pd_analysis$enhanced_signatures), "enhanced signatures\n")
  }, error = function(e) {
    cat("❌ PD analysis failed with error:", e$message, "\n")
    cat("Error details:", toString(e), "\n")
  })
}

cat("\n💡 DIAGNOSIS AND SOLUTION\n")
cat("=========================\n")

cat("The issue is clear:\n")
cat("1. Individual signatures have mast_gene/crispri_gene columns\n")
cat("2. Pan-cluster signatures only have gene_pair column\n") 
cat("3. PD analysis expects mast_gene/crispri_gene columns\n")
cat("4. Solution: Parse gene_pair when individual columns missing\n\n")

cat("Recommended fix:\n")
cat("- Add helper function to extract gene names from either format\n")
cat("- Update extract_signature_biological_context() to use helper\n")
cat("- Update get_signature_enrichment_terms() to use helper\n\n")

# Stop capturing output
sink()

cat("Debug output written to:", debug_file, "\n")
cat("Check the file for complete details about the column mismatch issue.\n")