# Test package integration to see what's happening with the Shiny app
# This will help diagnose why the updates aren't showing up

cat("🔍 TESTING PACKAGE INTEGRATION\n")
cat(paste(rep("=", 40), collapse=""), "\n\n")

# Test 1: Check if package is loaded
cat("1. PACKAGE STATUS\n")
cat("Package loaded:", "iSCORE.PDecipher" %in% loadedNamespaces(), "\n")

# Test 2: Check if signature analysis functions are available
cat("\n2. FUNCTION AVAILABILITY\n")
functions_to_check <- c(
  "discover_top_signatures",
  "analyze_pd_signatures", 
  "get_comparable_gene_pairs",
  "calculate_composite_signature_score"
)

for (func in functions_to_check) {
  available <- exists(func, mode = "function")
  cat(paste0("  ", func, ": ", if(available) "✅ Available" else "❌ Missing"), "\n")
}

# Test 3: Check if we can load the data
cat("\n3. DATA LOADING TEST\n")
tryCatch({
  # Try to find the data file
  data_paths <- c(
    "all_enrichment_padj005_complete_with_direction.rds",
    "../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
  )
  
  data_file <- NULL
  for (path in data_paths) {
    if (file.exists(path)) {
      data_file <- path
      break
    }
  }
  
  if (!is.null(data_file)) {
    cat("  Data file found:", data_file, "\n")
    data <- readRDS(data_file)
    cat("  Data loaded:", nrow(data), "rows\n")
    
    # Test if we can run the function
    if (exists("get_comparable_gene_pairs", mode = "function")) {
      gene_pairs <- get_comparable_gene_pairs()
      cat("  Gene pairs available:", nrow(gene_pairs), "\n")
    } else {
      cat("  ❌ get_comparable_gene_pairs function not available\n")
    }
  } else {
    cat("  ❌ Data file not found\n")
  }
}, error = function(e) {
  cat("  ❌ Error loading data:", e$message, "\n")
})

# Test 4: Check Shiny module files
cat("\n4. SHINY MODULE FILE CHECK\n")
module_path <- "inst/shiny/modules/mod_signature_nomination.R"
if (file.exists(module_path)) {
  cat("  Module file exists:", module_path, "\n")
  
  # Check if PD Biology Focus tab is in the file
  module_content <- readLines(module_path)
  pd_tab_line <- grep("PD Biology Focus", module_content)
  if (length(pd_tab_line) > 0) {
    cat("  ✅ PD Biology Focus tab found at line", pd_tab_line[1], "\n")
  } else {
    cat("  ❌ PD Biology Focus tab not found in module file\n")
  }
  
  # Check if PD analysis functions are sourced
  pd_source_line <- grep("pd_signature_interpretation", module_content)
  if (length(pd_source_line) > 0) {
    cat("  ✅ PD analysis sourcing found at line", pd_source_line[1], "\n")
  } else {
    cat("  ❌ PD analysis sourcing not found\n")
  }
} else {
  cat("  ❌ Module file not found:", module_path, "\n")
}

# Test 5: Check if package installation includes the right files
cat("\n5. PACKAGE INSTALLATION CHECK\n")
if ("iSCORE.PDecipher" %in% loadedNamespaces()) {
  pkg_path <- system.file(package = "iSCORE.PDecipher")
  cat("  Package installed at:", pkg_path, "\n")
  
  # Check if the signature nomination module exists in the installed package
  installed_module <- file.path(pkg_path, "shiny", "modules", "mod_signature_nomination.R")
  if (file.exists(installed_module)) {
    cat("  ✅ Signature nomination module found in installed package\n")
    
    # Check if it has the PD Biology Focus tab
    installed_content <- readLines(installed_module)
    if (any(grepl("PD Biology Focus", installed_content))) {
      cat("  ✅ PD Biology Focus tab found in installed module\n")
    } else {
      cat("  ❌ PD Biology Focus tab NOT found in installed module\n")
      cat("  This suggests the package installation is using an old version!\n")
    }
  } else {
    cat("  ❌ Signature nomination module not found in installed package\n")
  }
  
  # Check if PD analysis functions are in the installed package
  pd_file <- file.path(pkg_path, "R", "pd_signature_interpretation.R") 
  if (file.exists(pd_file)) {
    cat("  ✅ PD analysis file found in installed package\n")
  } else {
    cat("  ❌ PD analysis file NOT found in installed package\n")
  }
} else {
  cat("  ❌ Package not loaded\n")
}

cat("\n6. DIAGNOSIS\n")
cat("If you see ❌ symbols above, that explains why the Shiny app updates aren't visible.\n")
cat("The most likely issues:\n")
cat("- Package installation didn't include the updated files\n") 
cat("- GitHub cache is serving old version\n")
cat("- Module sourcing is failing silently\n")
cat("- Functions exist but PD analysis is not running\n")

cat("\n", paste(rep("=", 40), collapse=""), "\n")
cat("INTEGRATION TEST COMPLETE\n")
cat(paste(rep("=", 40), collapse=""), "\n")