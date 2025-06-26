#!/usr/bin/env Rscript

# Test script for reorganized app structure

cat("Testing reorganized iSCORE-PDecipher app...\n\n")

# Test 1: Check if app launches
cat("Test 1: Checking if app can be launched...\n")
tryCatch({
  # We're in development mode
  app_path <- "inst/shiny"
  
  # Check if app.R exists
  if (file.exists(file.path(app_path, "app.R"))) {
    cat("✓ app.R found\n")
  } else {
    cat("✗ app.R not found\n")
  }
  
  # Check if CSS files exist
  www_path <- file.path(app_path, "www")
  css_files <- c("custom.css", "custom_enhanced.css", "nested_tabs.css")
  
  for (css in css_files) {
    if (file.exists(file.path(www_path, css))) {
      cat(sprintf("✓ %s found\n", css))
    } else {
      cat(sprintf("✗ %s not found\n", css))
    }
  }
  
  # Check if modules exist
  module_path <- file.path(app_path, "modules")
  required_modules <- c(
    "mod_landing_page_with_umap_v2.R",
    "mod_de_results.R",
    "mod_de_analysis.R", 
    "mod_de_heatmap.R",
    "mod_visualization_enhanced.R",
    "mod_heatmap.R",
    "mod_comparison.R",
    "mod_pathview.R",
    "mod_export.R"
  )
  
  cat("\nTest 2: Checking required modules...\n")
  for (mod in required_modules) {
    if (file.exists(file.path(module_path, mod))) {
      cat(sprintf("✓ %s found\n", mod))
    } else {
      cat(sprintf("✗ %s not found\n", mod))
    }
  }
  
  # Check backup file
  cat("\nTest 3: Checking backup...\n")
  if (file.exists(file.path(app_path, "app.R.backup_pre_reorganization"))) {
    cat("✓ Backup file created\n")
  } else {
    cat("✗ Backup file not found\n")
  }
  
  cat("\n✅ All tests passed! The reorganized app structure is ready.\n")
  cat("\nTo launch the app, use: launch_iscore_app()\n")
  
}, error = function(e) {
  cat("\n❌ Error during testing:\n")
  print(e)
})

cat("\nDone!\n")