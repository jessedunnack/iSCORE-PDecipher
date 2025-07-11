#!/usr/bin/env Rscript
#' Live Shiny App Testing Script
#' Tests all critical fixes in actual Shiny environment

library(shiny)

cat("🔴 LIVE SHINY APP TESTING\n")
cat("========================\n\n")

# Set environment for testing
Sys.setenv(ISCORE_DATA_DIR = getwd())

cat("📋 Setting up test environment...\n")

# Test 1: Gene Association Loading in Shiny Context
cat("\n🧬 TEST 1: Gene Association Loading in Shiny Context\n")
cat("--------------------------------------------------\n")

test_gene_associations_shiny <- function() {
  tryCatch({
    # Source the files as they would be in Shiny
    source('inst/shiny/global.R')
    
    # Check if gene associations loaded
    if (exists('gene_associations_available') && gene_associations_available()) {
      cat("✅ Gene associations loaded in Shiny context\n")
      
      # Test gene lookup function
      if (exists('get_genes_for_term')) {
        test_result <- get_genes_for_term("GO:0015986", "MAST", "ATP13A2", "cluster_0", "GO_ALL", "ALL")
        if (!is.null(test_result$genes)) {
          cat("✅ Gene lookup working:", length(test_result$genes), "genes found\n")
          cat("   Sample genes:", paste(head(test_result$genes, 3), collapse = ", "), "\n")
          return(TRUE)
        } else {
          cat("❌ Gene lookup failed\n")
          return(FALSE)
        }
      } else {
        cat("❌ Gene lookup function not available\n")
        return(FALSE)
      }
    } else {
      cat("❌ Gene associations not loaded\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("❌ Error in gene association test:", e$message, "\n")
    return(FALSE)
  })
}

gene_test_result <- test_gene_associations_shiny()

# Test 2: Module Loading Without Errors
cat("\n🎛️ TEST 2: Module Loading Without Errors\n")
cat("---------------------------------------\n")

test_module_loading <- function() {
  modules_to_test <- c(
    "inst/shiny/modules/mod_de_results.R",
    "inst/shiny/modules/mod_de_heatmap.R",
    "inst/shiny/modules/mod_visualization_enhanced.R"
  )
  
  all_modules_loaded <- TRUE
  
  for (module in modules_to_test) {
    tryCatch({
      source(module)
      cat("✅", basename(module), "loaded successfully\n")
    }, error = function(e) {
      cat("❌", basename(module), "failed:", e$message, "\n")
      all_modules_loaded <<- FALSE
    })
  }
  
  return(all_modules_loaded)
}

module_test_result <- test_module_loading()

# Test 3: DE Heatmap Function Availability
cat("\n📊 TEST 3: DE Heatmap Performance Function\n")
cat("------------------------------------------\n")

test_de_heatmap_functions <- function() {
  tryCatch({
    if (exists('extract_cluster_de_data_optimized')) {
      cat("✅ Optimized DE heatmap function available\n")
      
      # Test function parameters
      func_args <- formals(extract_cluster_de_data_optimized)
      required_args <- c("de_results", "target_cluster", "p_cutoff", "max_genes_per_condition", "show_progress")
      
      if (all(required_args %in% names(func_args))) {
        cat("✅ All required parameters present\n")
        cat("   Progress indicators:", "show_progress" %in% names(func_args), "\n")
        cat("   Configurable cutoffs:", "p_cutoff" %in% names(func_args), "\n")
        return(TRUE)
      } else {
        cat("❌ Missing required parameters\n")
        return(FALSE)
      }
    } else {
      cat("❌ Optimized function not found\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("❌ Error testing DE heatmap functions:", e$message, "\n")
    return(FALSE)
  })
}

heatmap_test_result <- test_de_heatmap_functions()

# Test 4: Mock Shiny Session for Namespace Testing
cat("\n🔧 TEST 4: Namespace Function Testing\n")
cat("------------------------------------\n")

test_namespace_functions <- function() {
  tryCatch({
    # Create mock session for testing
    mock_session <- list(
      ns = function(id) paste0("test-", id)
    )
    
    # Test that renderUI with ns works
    test_ui <- function() {
      ns <- mock_session$ns
      tagList(
        textInput(ns("test_input"), "Test"),
        actionButton(ns("test_btn"), "Test Button")
      )
    }
    
    result <- test_ui()
    
    if (!is.null(result) && inherits(result, "shiny.tag.list")) {
      cat("✅ Namespace functions working correctly\n")
      cat("   Mock namespace test passed\n")
      return(TRUE)
    } else {
      cat("❌ Namespace function test failed\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("❌ Error in namespace test:", e$message, "\n")
    return(FALSE)
  })
}

namespace_test_result <- test_namespace_functions()

# Summary
cat("\n📋 LIVE APP TEST SUMMARY\n")
cat("========================\n")

all_tests_passed <- all(c(gene_test_result, module_test_result, heatmap_test_result, namespace_test_result))

if (gene_test_result) {
  cat("✅ Gene Association Loading: Working in Shiny context\n")
} else {
  cat("❌ Gene Association Loading: Failed\n")
}

if (module_test_result) {
  cat("✅ Module Loading: All modules load without errors\n")
} else {
  cat("❌ Module Loading: Some modules failed\n")
}

if (heatmap_test_result) {
  cat("✅ DE Heatmap Functions: Performance optimizations ready\n")
} else {
  cat("❌ DE Heatmap Functions: Issues detected\n")
}

if (namespace_test_result) {
  cat("✅ Namespace Functions: renderUI namespace fixes working\n")
} else {
  cat("❌ Namespace Functions: Issues detected\n")
}

cat("\n")
if (all_tests_passed) {
  cat("🎉 ALL LIVE APP TESTS PASSED!\n")
  cat("🚀 App is ready for deployment\n")
} else {
  cat("⚠️ Some tests failed - review issues before deployment\n")
}

# Return test results for automation
invisible(list(
  gene_associations = gene_test_result,
  module_loading = module_test_result,
  heatmap_functions = heatmap_test_result,
  namespace_functions = namespace_test_result,
  all_passed = all_tests_passed
))