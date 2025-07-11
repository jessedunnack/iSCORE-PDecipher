# Simple test of file detection and parsing
source("R/extract_gene_associations.R")

cat("=== Testing File Detection ===\n")

# Test file detection with correct path
test_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD"

# Use find command instead of Sys.glob
files_cmd <- paste0("find ", test_dir, "/enrichment_results -name '*.rds' | grep -E '(GO_|KEGG|Reactome|STRING|GSEA)' | head -5")
sample_files <- system(files_cmd, intern = TRUE)

cat("Found", length(sample_files), "sample files:\n")
for (file in sample_files) {
  cat(" ", file, "\n")
}

# Test metadata parsing on first file
if (length(sample_files) > 0) {
  test_file <- sample_files[1]
  cat("\nTesting metadata parsing on:", test_file, "\n")
  
  tryCatch({
    metadata <- parse_file_metadata(test_file)
    cat("Analysis type:", metadata$analysis_type, "\n")
    cat("Gene:", metadata$gene, "\n") 
    cat("Cluster:", metadata$cluster, "\n")
    cat("Enrichment type:", metadata$enrichment_type, "\n")
    cat("Direction:", metadata$direction, "\n")
    
    # Test extraction
    result <- extract_from_single_file(test_file, metadata)
    if (!is.null(result)) {
      cat("Successfully extracted", nrow(result), "associations\n")
      cat("Sample data:\n")
      print(head(result[, c("term_id", "description", "gene_count")], 3))
    }
    
  }, error = function(e) {
    cat("Error:", e$message, "\n")
  })
}