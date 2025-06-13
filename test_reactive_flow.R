#!/usr/bin/env Rscript
# Test reactive data flow without UI

cat("=== Testing Reactive Data Flow ===\n\n")

suppressPackageStartupMessages({
  library(shiny)
  library(dplyr)
})

# Change to shiny directory
setwd("inst/shiny")

# Source required files
source("global.R")
source("R/data_manager.R")

# Create mock data similar to real enrichment data
create_mock_data <- function(n = 1000) {
  genes <- c("LRRK2", "PINK1", "PARK7", "SNCA", "GBA", "FBXO7", "ATP13A2")
  clusters <- paste0("cluster_", 0:9)
  types <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "STRING")
  directions <- c("UP", "DOWN", "ALL")
  
  data.frame(
    method = sample(c("MAST", "MixScale"), n, replace = TRUE),
    mutation_perturbation = sample(genes, n, replace = TRUE),
    gene = sample(genes, n, replace = TRUE),
    cluster = sample(clusters, n, replace = TRUE),
    enrichment_type = sample(types, n, replace = TRUE),
    direction = sample(directions, n, replace = TRUE),
    Description = paste("Term", 1:n),
    p.adjust = runif(n, 0.0001, 0.1),
    FoldEnrichment = runif(n, 1.2, 5),
    Count = sample(5:100, n, replace = TRUE),
    experiment = "default",
    stringsAsFactors = FALSE
  )
}

# Test the filtering function with various inputs
test_filtering <- function() {
  cat("Creating mock data...\n")
  mock_data <- create_mock_data(5000)
  cat("✓ Created", nrow(mock_data), "mock enrichment terms\n\n")
  
  # Test 1: Basic filtering
  cat("Test 1: Basic filtering (LRRK2, cluster_0, GO_BP, UP)\n")
  result1 <- tryCatch({
    get_significant_terms_from_consolidated(
      mock_data,
      gene = "LRRK2",
      cluster = "cluster_0",
      enrichment_type = "GO_BP",
      direction = "UP",
      pval_threshold = 0.05
    )
  }, error = function(e) {
    cat("✗ Error:", e$message, "\n")
    NULL
  })
  
  if (!is.null(result1)) {
    cat("✓ Filtered to", nrow(result1), "terms\n")
    cat("  P-value range:", range(result1$p.adjust), "\n\n")
  }
  
  # Test 2: Filter with "ALL" direction
  cat("Test 2: Testing 'ALL' direction filter\n")
  result2 <- tryCatch({
    get_significant_terms_from_consolidated(
      mock_data,
      gene = "PINK1",
      cluster = "cluster_1",
      enrichment_type = "KEGG",
      direction = "ALL",
      pval_threshold = 0.05
    )
  }, error = function(e) {
    cat("✗ Error:", e$message, "\n")
    NULL
  })
  
  if (!is.null(result2)) {
    cat("✓ Filtered to", nrow(result2), "terms\n")
    cat("  Directions present:", unique(result2$direction), "\n\n")
  }
  
  # Test 3: Method-specific filtering
  cat("Test 3: Method-specific filtering\n")
  for (method in c("MAST", "MixScale")) {
    result <- get_significant_terms_from_consolidated(
      mock_data,
      analysis_type = method,
      pval_threshold = 0.05
    )
    cat("  ", method, ":", nrow(result), "terms\n")
  }
  cat("\n")
  
  # Test 4: Edge cases
  cat("Test 4: Edge cases\n")
  
  # Empty data
  empty_result <- get_significant_terms_from_consolidated(
    data.frame(),
    gene = "LRRK2"
  )
  cat("  Empty data returns:", nrow(empty_result), "rows ✓\n")
  
  # NULL inputs
  null_result <- get_significant_terms_from_consolidated(
    mock_data,
    gene = NULL,
    cluster = NULL
  )
  cat("  NULL filters return:", nrow(null_result), "rows\n")
  
  # Non-existent gene
  no_gene <- get_significant_terms_from_consolidated(
    mock_data,
    gene = "NONEXISTENT"
  )
  cat("  Non-existent gene returns:", nrow(no_gene), "rows ✓\n")
}

# Test reactive chain simulation
test_reactive_chain <- function() {
  cat("\nTest 5: Simulating reactive chain\n")
  
  # Simulate the reactive flow in the app
  mock_data <- create_mock_data(10000)
  
  # Step 1: User selects global settings
  selection <- list(
    analysis_type = "MAST",
    gene = "LRRK2", 
    cluster = "cluster_0",
    enrichment_type = "GO_BP",
    direction = "UP",
    pval_threshold = 0.05
  )
  
  cat("  User selection:", paste(names(selection), selection, sep="=", collapse=", "), "\n")
  
  # Step 2: Filter data (simulates filtered_data reactive)
  filtered <- get_significant_terms_from_consolidated(
    mock_data,
    gene = selection$gene,
    cluster = selection$cluster,
    enrichment_type = selection$enrichment_type,
    direction = selection$direction,
    analysis_type = selection$analysis_type,
    pval_threshold = selection$pval_threshold
  )
  
  cat("  Filtered data:", nrow(filtered), "terms\n")
  
  # Step 3: Prepare for visualization (top N terms)
  if (nrow(filtered) > 0) {
    top_terms <- filtered %>%
      arrange(p.adjust) %>%
      head(20)
    
    cat("  Top 20 terms for visualization\n")
    cat("  P-value range:", range(top_terms$p.adjust), "\n")
    cat("  Enrichment range:", range(top_terms$FoldEnrichment), "\n")
  }
}

# Run all tests
cat("Testing filtering function robustness...\n")
cat("=====================================\n\n")

test_filtering()
test_reactive_chain()

cat("\n=== Summary ===\n")
cat("This test validates the data flow without launching the Shiny app.\n")
cat("Any errors here would cause issues in the live app.\n")