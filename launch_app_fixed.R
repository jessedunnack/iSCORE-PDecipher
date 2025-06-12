#!/usr/bin/env Rscript

# Fixed launcher for iSCORE-PDecipher Shiny app
# This bypasses the package loading issues and launches directly

# Get the dataset directory
args <- commandArgs(trailingOnly = TRUE)

if (length(args) > 0) {
  dataset_choice <- args[1]
} else {
  cat("\n=== Select Dataset ===\n")
  cat("Available datasets:\n")
  cat("[1] iSCORE-PD only\n")
  cat("    E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD\n")
  cat("[2] iSCORE-PD + CRISPRi\n") 
  cat("    E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi\n")
  cat("[3] iSCORE-PD + CRISPRi + CRISPRa\n")
  cat("    E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa\n")
  cat("[Q] Quit\n")
  
  dataset_choice <- readline("Enter your choice: ")
}

# Map choice to directory
data_dir <- switch(dataset_choice,
  "1" = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD",
  "2" = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi",
  "3" = "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi_and_CRISPRa",
  NULL
)

if (is.null(data_dir) || tolower(dataset_choice) == "q") {
  cat("Exiting...\n")
  quit("no")
}

# Validate the directory
if (!dir.exists(data_dir)) {
  stop("Dataset directory not found: ", data_dir)
}

enrichment_file <- file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds")
if (!file.exists(enrichment_file)) {
  stop("Enrichment file not found: ", enrichment_file)
}

cat("\n✓ Dataset validated. Launching app...\n\n")

# Set environment variables
Sys.setenv(ISCORE_DATA_DIR = normalizePath(data_dir))
Sys.setenv(ISCORE_HAS_DATA = "TRUE")
Sys.setenv(ISCORE_DE_FILE = file.path(data_dir, "full_DE_results.rds"))
Sys.setenv(ISCORE_ENRICHMENT_FILE = enrichment_file)
Sys.setenv(ISCORE_ENRICHMENT_DIR = file.path(data_dir, "enrichment_results"))

# Launch the app directly
app_dir <- "inst/shiny"
if (!dir.exists(app_dir)) {
  # Try from different location
  if (dir.exists("iSCORE-PDecipher/inst/shiny")) {
    app_dir <- "iSCORE-PDecipher/inst/shiny"
  } else {
    stop("Could not find Shiny app directory")
  }
}

# Load shiny and run the app
library(shiny)
runApp(app_dir, launch.browser = TRUE)