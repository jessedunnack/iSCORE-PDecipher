# Debug script for CRISPRi dotplot issue (Bug #5)
# This script investigates why CRISPRi dotplots fail to appear

library(tidyverse)
library(iSCORE.PDecipher)

# Find and load the consolidated data
data_paths <- c(
  "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi",
  "../../iSCORE-PD_plus_CRISPRi",
  "../iSCORE-PD_plus_CRISPRi",
  "iSCORE-PD_plus_CRISPRi"
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

# Test the filtering function used in the app
cat("\n\n=== TESTING APP FILTERING LOGIC ===\n")

# Source the filtering function if not available
if (!exists("get_significant_terms_from_consolidated")) {
  # Define it inline for testing
  get_significant_terms_from_consolidated <- function(data, gene = NULL, cluster = NULL, 
                                                     enrichment_type = NULL, direction = "ALL",
                                                     modality = NULL, analysis_type = NULL,
                                                     experiment = NULL, pval_threshold = 0.05) {
    
    if (is.null(data) || nrow(data) == 0) {
      return(data.frame())
    }
    
    # Start with all data
    filtered_data <- data
    
    # Filter by gene/mutation
    if (!is.null(gene) && gene != "All" && gene != "") {
      # Check which column to use
      if ("gene" %in% names(filtered_data)) {
        filtered_data <- filtered_data[filtered_data$gene == gene, ]
      } else if ("mutation_perturbation" %in% names(filtered_data)) {
        filtered_data <- filtered_data[filtered_data$mutation_perturbation == gene, ]
      }
    }
    
    # Filter by cluster
    if (!is.null(cluster) && cluster != "All") {
      filtered_data <- filtered_data[filtered_data$cluster == cluster, ]
    }
    
    # Filter by enrichment type
    if (!is.null(enrichment_type) && enrichment_type != "All") {
      filtered_data <- filtered_data[filtered_data$enrichment_type == enrichment_type, ]
    }
    
    # Filter by direction (FIXED: only filter when not "ALL")
    if (!is.null(direction) && direction != "ALL") {
      filtered_data <- filtered_data[filtered_data$direction == direction, ]
    }
    
    # Filter by modality if specified
    if (!is.null(modality) && modality != "All" && "modality" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$modality == modality, ]
    }
    
    # Filter by analysis type (method) if specified
    if (!is.null(analysis_type) && analysis_type != "All" && "method" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$method == analysis_type, ]
    }
    
    # Filter by experiment if specified
    if (!is.null(experiment) && experiment != "All" && experiment != "default" && "experiment" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$experiment == experiment, ]
    }
    
    # Filter by p-value threshold
    if (!is.null(pval_threshold) && "p.adjust" %in% names(filtered_data)) {
      filtered_data <- filtered_data[filtered_data$p.adjust <= pval_threshold, ]
    }
    
    # Sort by p-value and limit for performance
    if (nrow(filtered_data) > 0) {
      filtered_data <- filtered_data[order(filtered_data$p.adjust), ]
      
      # Limit to top 1000 for performance
      if (nrow(filtered_data) > 1000) {
        filtered_data <- filtered_data[1:1000, ]
      }
    }
    
    return(filtered_data)
  }
}

# Test different filtering scenarios
cat("\nTest 1: Filter for MAST data (baseline)\n")
mast_data <- get_significant_terms_from_consolidated(
  consolidated_data,
  analysis_type = "MAST",
  gene = "LRRK2",
  cluster = "1",
  enrichment_type = "GO_BP",
  direction = "ALL"
)
cat("MAST filtered rows:", nrow(mast_data), "\n")

cat("\nTest 2: Filter for CRISPRi with modality='CRISPRi'\n")
crispri_data_mod <- get_significant_terms_from_consolidated(
  consolidated_data,
  modality = "CRISPRi",
  gene = "LRRK2",
  cluster = "1",
  enrichment_type = "GO_BP",
  direction = "ALL"
)
cat("CRISPRi filtered rows (modality):", nrow(crispri_data_mod), "\n")

cat("\nTest 3: Filter for MixScale analysis type\n")
mixscale_data <- get_significant_terms_from_consolidated(
  consolidated_data,
  analysis_type = "MixScale_CRISPRi",
  gene = "LRRK2",
  cluster = "1",
  enrichment_type = "GO_BP",
  direction = "ALL"
)
cat("MixScale filtered rows:", nrow(mixscale_data), "\n")

# Look for LRRK2 CRISPRi data specifically
cat("\n\n=== SEARCHING FOR LRRK2 CRISPRI DATA ===\n")
lrrk2_data <- consolidated_data %>%
  filter(grepl("LRRK2", mutation_perturbation, ignore.case = TRUE) | 
         grepl("LRRK2", gene, ignore.case = TRUE))

cat("Total LRRK2 rows:", nrow(lrrk2_data), "\n")

# Check methods for LRRK2
if (nrow(lrrk2_data) > 0) {
  cat("\nMethods available for LRRK2:\n")
  print(table(lrrk2_data$method))
  
  # Check if there's CRISPRi data
  crispri_lrrk2 <- lrrk2_data %>%
    filter(grepl("MixScale|CRISPRi", method, ignore.case = TRUE))
  
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

if (nrow(crispri_data_mod) == 0 && nrow(mixscale_data) == 0) {
  cat("❌ CONFIRMED: No CRISPRi data is returned by the filtering function\n")
  cat("   This explains why dotplots disappear when switching to CRISPRi\n")
}

# Save diagnostic info
saveRDS(list(
  consolidated_columns = colnames(consolidated_data),
  modality_values = if("modality" %in% colnames(consolidated_data)) unique(consolidated_data$modality) else NULL,
  method_values = if("method" %in% colnames(consolidated_data)) unique(consolidated_data$method) else NULL,
  mast_test_rows = nrow(mast_data),
  crispri_test_rows = nrow(crispri_data_mod),
  mixscale_test_rows = nrow(mixscale_data)
), "debug_crispri_diagnostics.rds")

cat("\nDiagnostic information saved to: debug_crispri_diagnostics.rds\n")