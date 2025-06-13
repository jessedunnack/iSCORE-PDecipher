#!/usr/bin/env Rscript

# Test script to verify the double prompt fix
# This tests the full launch_iscore_app() function

cat("=== Testing iSCORE-PDecipher Launch Fix ===\n")
cat("This should show ONE dataset selection prompt.\n")
cat("Enter 'q' to quit without launching the app.\n\n")

# Load the package functions
library(iSCORE.PDecipher)

# Test with no data_dir to trigger dataset selection
tryCatch({
  launch_iscore_app(data_dir = NULL)
}, error = function(e) {
  cat("Launch cancelled or failed:", e$message, "\n")
})

cat("Test completed.\n")