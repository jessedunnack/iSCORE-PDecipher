#!/usr/bin/env Rscript

# Test script to verify dataset selection behavior
# This script tests the select_dataset_directory function in isolation

# Source the required functions
source("R/launch_app.R")
source("R/dataset_validator.R")

cat("Testing dataset selection function...\n")
cat("This should show ONE prompt and wait for input.\n")
cat("Press Enter without typing anything to test empty input handling.\n")
cat("Then type 'q' to quit.\n\n")

# Test the function
result <- select_dataset_directory()

if (is.null(result)) {
  cat("Selection cancelled by user.\n")
} else {
  cat("Selected dataset:", result, "\n")
}