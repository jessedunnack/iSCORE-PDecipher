# Test the gene association extraction with a small subset
source("R/extract_gene_associations.R")

# Test with just a few files first
cat("=== Testing Gene Association Extraction ===\n")

# Get a few sample files
sample_files <- c(
  "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/enrichment_results/MAST/ATP13A2/cluster_0/default/GO_BP/GO_BP_UP.rds",
  "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/enrichment_results/MAST/ATP13A2/cluster_0/default/KEGG/KEGG_ALL.rds"
)

# Test metadata parsing
for (file in sample_files) {
  if (file.exists(file)) {
    cat("\nTesting file:", file, "\n")
    metadata <- parse_file_metadata(file)
    cat("Metadata:", str(metadata), "\n")
    
    # Test extraction
    result <- extract_from_single_file(file, metadata)
    if (!is.null(result)) {
      cat("Extracted", nrow(result), "associations\n")
      cat("Sample term:", result$term_id[1], "\n")
      cat("Sample genes:", substr(result$associated_genes[1], 1, 50), "...\n")
    } else {
      cat("No results extracted\n")
    }
  } else {
    cat("File does not exist:", file, "\n")
  }
}

cat("\n=== Test Complete ===\n")