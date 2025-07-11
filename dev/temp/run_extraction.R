# Run the full gene association extraction
source("R/extract_gene_associations.R")

cat("=== Starting Full Gene Association Extraction ===\n")
cat("This will process ~16,000 enrichment files...\n")

# Run the extraction
start_time <- Sys.time()
associations <- main()
end_time <- Sys.time()

cat("Total extraction time:", round(difftime(end_time, start_time, units = "mins"), 1), "minutes\n")

# Show summary statistics
if (!is.null(associations) && nrow(associations) > 0) {
  cat("\n=== Final Summary ===\n")
  cat("Total associations:", nrow(associations), "\n")
  cat("Unique terms:", n_distinct(associations$term_id), "\n")
  cat("Unique genes:", n_distinct(associations$gene), "\n")
  cat("Analysis types:", paste(sort(unique(associations$analysis_type)), collapse = ", "), "\n")
  cat("Enrichment types:", paste(sort(unique(associations$enrichment_type)), collapse = ", "), "\n")
  cat("Directions:", paste(sort(unique(associations$direction)), collapse = ", "), "\n")
  
  # Check file size
  output_file <- "inst/extdata/gene_term_associations.rds"
  if (file.exists(output_file)) {
    file_size_mb <- round(file.size(output_file) / (1024^2), 1)
    cat("Output file size:", file_size_mb, "MB\n")
    
    if (file_size_mb > 100) {
      cat("WARNING: File size may be too large for GitHub deployment!\n")
    } else {
      cat("File size acceptable for GitHub deployment.\n")
    }
  }
} else {
  cat("ERROR: No associations extracted!\n")
}

cat("=== Extraction Complete ===\n")