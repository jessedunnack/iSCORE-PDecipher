# Simple debug script for CRISPRi dotplot issue (Bug #5)
# Find and load the consolidated data

# Find data file
data_paths <- c(
  "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi",
  "../../iSCORE-PD_plus_CRISPRi"
)

consolidated_file <- NULL
for (path in data_paths) {
  test_file <- file.path(path, "all_enrichment_padj005_complete_with_direction.rds")
  if (file.exists(test_file)) {
    consolidated_file <- test_file
    cat("Found consolidated data at:", consolidated_file, "\n")
    break
  }
}

if (is.null(consolidated_file)) {
  stop("Could not find consolidated data file!")
}

# Load the data
cat("\nLoading consolidated data...\n")
consolidated_data <- readRDS(consolidated_file)

# Basic data inspection
cat("\nData dimensions:", nrow(consolidated_data), "rows x", ncol(consolidated_data), "columns\n")
cat("\nColumn names:\n")
print(colnames(consolidated_data))

# Check for modality column
cat("\n\n=== MODALITY COLUMN CHECK ===\n")
if ("modality" %in% colnames(consolidated_data)) {
  cat("✓ 'modality' column exists\n")
  cat("Unique modality values:\n")
  print(table(consolidated_data$modality, useNA = "always"))
} else {
  cat("✗ 'modality' column NOT FOUND\n")
}

# Check method column (alternative to modality)
cat("\n\n=== METHOD COLUMN CHECK ===\n")
if ("method" %in% colnames(consolidated_data)) {
  cat("✓ 'method' column exists\n")
  cat("Unique method values:\n")
  print(table(consolidated_data$method, useNA = "always"))
} else {
  cat("✗ 'method' column NOT FOUND\n")
}

# Check how CRISPRi data is identified
cat("\n\n=== CRISPRI DATA IDENTIFICATION ===\n")

# Check analysis_type patterns
if ("analysis_type" %in% colnames(consolidated_data)) {
  cat("\nUnique analysis_type values:\n")
  analysis_types <- unique(consolidated_data$analysis_type)
  print(analysis_types)
  
  # Check for CRISPRi patterns
  crispri_patterns <- grep("CRISPRi|MixScale|Perturb", analysis_types, value = TRUE, ignore.case = TRUE)
  if (length(crispri_patterns) > 0) {
    cat("\nFound CRISPRi-related analysis types:\n")
    print(crispri_patterns)
  }
}

# Look for LRRK2 CRISPRi data specifically
cat("\n\n=== SEARCHING FOR LRRK2 CRISPRI DATA ===\n")
lrrk2_data <- consolidated_data[grepl("LRRK2", consolidated_data$mutation_perturbation, ignore.case = TRUE), ]

cat("Total LRRK2 rows:", nrow(lrrk2_data), "\n")

# Check methods for LRRK2
if (nrow(lrrk2_data) > 0) {
  cat("\nMethods available for LRRK2:\n")
  print(table(lrrk2_data$method))
  
  # Check if there's CRISPRi data
  crispri_lrrk2 <- lrrk2_data[grepl("MixScale|CRISPRi", lrrk2_data$method, ignore.case = TRUE), ]
  
  cat("\nLRRK2 CRISPRi/MixScale rows:", nrow(crispri_lrrk2), "\n")
  
  if (nrow(crispri_lrrk2) > 0) {
    cat("\nSample of LRRK2 CRISPRi data:\n")
    print(crispri_lrrk2[1:min(5, nrow(crispri_lrrk2)), c("mutation_perturbation", "method", "cluster", "enrichment_type", "direction", "Description")])
  }
}

# Final diagnosis
cat("\n\n=== DIAGNOSIS ===\n")
if ("modality" %in% colnames(consolidated_data)) {
  if (!"CRISPRi" %in% consolidated_data$modality) {
    cat("❌ ISSUE: 'modality' column exists but does not contain 'CRISPRi' value\n")
    cat("   The app filters for modality='CRISPRi' but this value doesn't exist in the data\n")
  }
} else {
  cat("❌ ISSUE: 'modality' column is missing from consolidated data\n")
  cat("   The app tries to filter by modality='CRISPRi' but the column doesn't exist\n")
}

cat("\nSaving results to debug_results.txt\n")