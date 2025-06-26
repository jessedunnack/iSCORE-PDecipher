# Debug Data Methods Detection
# Check what methods and modalities are actually in the loaded data

data_file <- "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"

cat("=== Debugging Data Methods Detection ===\n")
cat("Loading data from:", data_file, "\n")

if (file.exists(data_file)) {
  data <- readRDS(data_file)
  
  cat("Data loaded successfully with", nrow(data), "rows\n")
  cat("Column names:", paste(names(data), collapse = ", "), "\n\n")
  
  # Check method column
  if ("method" %in% names(data)) {
    methods <- unique(data$method)
    cat("Unique values in 'method' column:", paste(methods, collapse = ", "), "\n")
    
    # Count by method
    method_counts <- table(data$method)
    cat("Method counts:\n")
    print(method_counts)
  } else {
    cat("❌ 'method' column not found!\n")
  }
  
  # Check modality column
  if ("modality" %in% names(data)) {
    modalities <- unique(data$modality)
    cat("\nUnique values in 'modality' column:", paste(modalities, collapse = ", "), "\n")
    
    # Count by modality
    modality_counts <- table(data$modality)
    cat("Modality counts:\n")
    print(modality_counts)
    
    # Cross-tabulation
    if ("method" %in% names(data)) {
      cat("\nMethod × Modality cross-tabulation:\n")
      print(table(data$method, data$modality))
    }
  } else {
    cat("\n❌ 'modality' column not found!\n")
    cat("Available columns:", paste(names(data), collapse = ", "), "\n")
  }
  
  # Check what the detect function would return
  cat("\n=== Detection Logic Results ===\n")
  methods <- unique(data$method)
  modalities <- if ("modality" %in% names(data)) unique(data$modality) else character(0)
  
  cat("Methods found:", paste(methods, collapse = ", "), "\n")
  cat("Modalities found:", paste(modalities, collapse = ", "), "\n")
  
  detection_results <- list(
    MAST = "MAST" %in% methods,
    CRISPRi = "MixScale" %in% methods && "CRISPRi" %in% modalities,
    CRISPRa = "MixScale" %in% methods && "CRISPRa" %in% modalities
  )
  
  cat("Detection results:\n")
  for (method in names(detection_results)) {
    cat(" ", method, ":", detection_results[[method]], "\n")
  }
  
} else {
  cat("❌ Data file not found:", data_file, "\n")
}

cat("\n=== Debug Complete ===\n")