#!/usr/bin/env Rscript
# Test script to validate app functions without launching Shiny

cat("=== Testing iSCORE-PDecipher App Functions ===\n\n")

# Set working directory to shiny app
setwd("inst/shiny")

# Test 1: Check if global.R loads properly
cat("Test 1: Loading global.R...\n")
tryCatch({
  source("global.R")
  cat("✓ global.R loaded successfully\n")
  cat("  - APP_CONFIG exists:", exists("APP_CONFIG"), "\n")
  cat("  - get_significant_terms_from_consolidated exists:", 
      exists("get_significant_terms_from_consolidated"), "\n\n")
}, error = function(e) {
  cat("✗ Error loading global.R:", e$message, "\n\n")
})

# Test 2: Check data manager
cat("Test 2: Testing data manager...\n")
tryCatch({
  source("R/data_manager.R")
  cat("✓ data_manager.R loaded\n")
  
  # Check if we can get enrichment data
  if (Sys.getenv("ISCORE_ENRICHMENT_FILE") != "") {
    data <- get_enrichment_data()
    if (!is.null(data)) {
      cat("✓ Enrichment data loaded:", nrow(data), "rows\n")
      cat("  - Columns:", paste(head(names(data), 10), collapse=", "), "...\n")
    } else {
      cat("✗ Failed to load enrichment data\n")
    }
  } else {
    cat("⚠ ISCORE_ENRICHMENT_FILE not set\n")
  }
  cat("\n")
}, error = function(e) {
  cat("✗ Error in data manager:", e$message, "\n\n")
})

# Test 3: Check startup manager
cat("Test 3: Testing startup manager...\n")
tryCatch({
  source("R/startup_manager.R")
  cat("✓ startup_manager.R loaded\n")
  cat("  - initialize_app_with_data exists:", exists("initialize_app_with_data"), "\n\n")
}, error = function(e) {
  cat("✗ Error loading startup manager:", e$message, "\n\n")
})

# Test 4: Test filtering function
cat("Test 4: Testing data filtering...\n")
if (exists("get_significant_terms_from_consolidated") && exists("data") && !is.null(data)) {
  tryCatch({
    # Test basic filtering
    filtered <- get_significant_terms_from_consolidated(
      data, 
      gene = "LRRK2",
      cluster = "cluster_0",
      enrichment_type = "GO_BP",
      direction = "UP",
      pval_threshold = 0.05
    )
    cat("✓ Filtering function works\n")
    cat("  - Input rows:", nrow(data), "\n")
    cat("  - Filtered rows:", nrow(filtered), "\n\n")
  }, error = function(e) {
    cat("✗ Error in filtering:", e$message, "\n\n")
  })
} else {
  cat("⚠ Cannot test filtering - missing function or data\n\n")
}

# Test 5: Load and test heatmap module functions
cat("Test 5: Testing heatmap module...\n")
tryCatch({
  source("modules/mod_heatmap.R")
  cat("✓ mod_heatmap.R loaded\n")
  
  # Check if helper functions exist
  env <- environment(mod_heatmap_server)
  cat("  - mod_heatmap_ui exists:", exists("mod_heatmap_ui"), "\n")
  cat("  - mod_heatmap_server exists:", exists("mod_heatmap_server"), "\n")
  
  # Test color validation if the module loaded
  if (exists("mod_heatmap_server")) {
    # The helper functions are defined inside the server function
    # We can't easily test them without running the server
    cat("  - Heatmap module structure appears valid\n")
  }
  cat("\n")
}, error = function(e) {
  cat("✗ Error loading heatmap module:", e$message, "\n\n")
})

# Test 6: Check all required modules
cat("Test 6: Checking all required modules...\n")
modules <- c(
  "modules/mod_landing_page_with_umap_v2.R",
  "modules/mod_precomputed_reactive.R",
  "modules/mod_visualization_enhanced.R",
  "modules/mod_comparison.R",
  "modules/mod_pathview.R",
  "modules/mod_export.R"
)

for (mod in modules) {
  if (file.exists(mod)) {
    cat("✓", mod, "exists\n")
  } else {
    cat("✗", mod, "NOT FOUND\n")
  }
}
cat("\n")

# Test 7: Package dependencies
cat("Test 7: Checking package dependencies...\n")
required_packages <- c("shiny", "shinyjs", "shinycssloaders", "shinyWidgets", 
                      "DT", "plotly", "ggplot2", "dplyr", "tidyr", "heatmaply")
for (pkg in required_packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat("✓", pkg, "installed\n")
  } else {
    cat("✗", pkg, "NOT INSTALLED\n")
  }
}

cat("\n=== Test Summary ===\n")
cat("Run this script with: Rscript test_app_functions.R\n")
cat("Fix any ✗ errors before running the app\n")