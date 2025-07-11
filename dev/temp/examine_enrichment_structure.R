# Script to examine enrichment file structure
library(dplyr)

# Examine a sample MAST file
sample_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD/enrichment_results/MAST/ATP13A2/cluster_0/default/GO_BP/GO_BP_UP.rds"

cat("=== Examining enrichment file structure ===\n")
cat("File:", sample_file, "\n")

if (file.exists(sample_file)) {
  result <- readRDS(sample_file)
  
  cat("Class:", class(result), "\n")
  cat("Type:", typeof(result), "\n")
  
  if (is.data.frame(result)) {
    cat("Data frame with", nrow(result), "rows and", ncol(result), "columns\n")
    cat("Column names:", paste(colnames(result), collapse = ", "), "\n")
    
    # Check for gene ID column
    if ("geneID" %in% colnames(result)) {
      cat("Sample geneID entries:\n")
      print(head(result$geneID, 3))
    }
    
    if ("ID" %in% colnames(result)) {
      cat("Sample term IDs:\n") 
      print(head(result$ID, 3))
    }
    
    if ("Description" %in% colnames(result)) {
      cat("Sample descriptions:\n")
      print(head(result$Description, 3))
    }
    
  } else if (inherits(result, "enrichResult")) {
    cat("clusterProfiler enrichResult object\n")
    cat("Result data frame has", nrow(result@result), "rows\n")
    cat("Columns:", paste(colnames(result@result), collapse = ", "), "\n")
    
    # Check for gene ID column in result slot
    if ("geneID" %in% colnames(result@result)) {
      cat("Sample geneID entries:\n")
      print(head(result@result$geneID, 3))
    }
    
  } else {
    cat("Unknown object type, structure:\n")
    str(result, max.level = 2)
  }
  
} else {
  cat("File does not exist\n")
}

cat("\n=== Done ===\n")