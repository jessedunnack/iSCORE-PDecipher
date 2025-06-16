#!/usr/bin/env Rscript

# Test script to diagnose launch issues

# Set working directory to package root
if (!file.exists("DESCRIPTION")) {
  if (file.exists("iSCORE-PDecipher/DESCRIPTION")) {
    setwd("iSCORE-PDecipher")
  }
}

# Load the package
library(iSCORE.PDecipher)

# Check what functions are available
cat("\nAvailable launch functions:\n")
if (exists("launch_app")) cat("- launch_app exists\n")
if (exists("launch_iscore_app")) cat("- launch_iscore_app exists\n")

# Check if they're the same
if (exists("launch_app") && exists("launch_iscore_app")) {
  cat("\nAre they identical?", identical(launch_app, launch_iscore_app), "\n")
}

# Try launching with dataset 1
cat("\nAttempting to launch with dataset 1...\n")
try({
  # Set to non-interactive mode for testing
  options(iscore.dataset.choice = 1)
  
  # Try a dry run first
  cat("\nChecking dataset directory...\n")
  data_dir <- "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD"
  
  if (dir.exists(data_dir)) {
    cat("✓ Dataset directory exists\n")
    
    # Check for required files
    enrichment_file <- file.path(data_dir, "all_enrichment_padj005_complete_with_direction.rds")
    if (file.exists(enrichment_file)) {
      cat("✓ Enrichment file exists\n")
      cat("  Size:", round(file.size(enrichment_file) / 1024^2, 1), "MB\n")
    } else {
      cat("✗ Enrichment file not found\n")
    }
  } else {
    cat("✗ Dataset directory not found\n")
  }
})