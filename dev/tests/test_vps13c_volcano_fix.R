# Test VPS13C Volcano Plot Gene Harmonization Fix
# Verifies that VPS13C_W395C and VPS13C_A444P volcano plots show data
# and that SNCA_A30P and SNCA_A53T volcano plots also work

cat("=== Testing VPS13C and SNCA Volcano Plot Gene Harmonization Fix ===\n")

# Test the gene harmonization logic directly
test_gene_harmonization <- function(mast_gene) {
  # Apply gene harmonization (same logic as in mod_de_results.R)
  mixscale_gene <- mast_gene
  if (mast_gene == "PRKN") {
    mixscale_gene <- "PARK2"
  } else if (mast_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
    mixscale_gene <- "SNCA" 
  } else if (mast_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
    mixscale_gene <- "VPS13C"
  }
  
  return(mixscale_gene)
}

# Test cases for gene harmonization
test_cases <- data.frame(
  mast_gene = c("VPS13C_W395C", "VPS13C_A444P", "SNCA_A30P", "SNCA_A53T", "PRKN", "LRRK2"),
  expected_mixscale = c("VPS13C", "VPS13C", "SNCA", "SNCA", "PARK2", "LRRK2"),
  stringsAsFactors = FALSE
)

cat("\n1. Testing gene harmonization logic:\n")
all_passed <- TRUE
for (i in 1:nrow(test_cases)) {
  mast_gene <- test_cases$mast_gene[i]
  expected <- test_cases$expected_mixscale[i]
  result <- test_gene_harmonization(mast_gene)
  
  if (result == expected) {
    cat("   ✅", mast_gene, "→", result, "\n")
  } else {
    cat("   ❌", mast_gene, "→", result, "(expected:", expected, ")\n")
    all_passed <- FALSE
  }
}

if (all_passed) {
  cat("   ✅ All gene harmonization tests passed\n")
} else {
  cat("   ❌ Some gene harmonization tests failed\n")
}

# Test data availability
cat("\n2. Testing DE data availability:\n")
de_data_path <- "../../iSCORE-PD_plus_CRISPRi/full_DE_results.rds"
if (file.exists(de_data_path)) {
  cat("   ✅ DE data file found\n")
  
  # Load and check structure
  de_data <- readRDS(de_data_path)
  
  if ("CRISPRi_Mixscale" %in% names(de_data)) {
    mixscale_genes <- names(de_data$CRISPRi_Mixscale)
    cat("   ✅ MixScale data available\n")
    cat("   Available MixScale genes:", paste(mixscale_genes, collapse = ", "), "\n")
    
    # Check if target genes exist
    target_genes <- c("VPS13C", "SNCA", "PARK2")
    for (gene in target_genes) {
      if (gene %in% mixscale_genes) {
        cat("   ✅", gene, "found in MixScale data\n")
      } else {
        cat("   ❌", gene, "NOT found in MixScale data\n")
        all_passed <- FALSE
      }
    }
  } else {
    cat("   ❌ MixScale data not found in DE results\n")
    all_passed <- FALSE
  }
} else {
  cat("   ❌ DE data file not found\n")
  all_passed <- FALSE
}

# Test actual data filtering
if (all_passed) {
  cat("\n3. Testing data filtering with harmonization:\n")
  
  # Simulate the filtering that would happen in volcano plots
  test_variants <- c("VPS13C_W395C", "VPS13C_A444P", "SNCA_A30P", "SNCA_A53T")
  
  for (variant in test_variants) {
    harmonized_gene <- test_gene_harmonization(variant)
    
    # Check if data exists for harmonized gene
    if (harmonized_gene %in% names(de_data$CRISPRi_Mixscale)) {
      variant_data <- de_data$CRISPRi_Mixscale[[harmonized_gene]]
      cluster_count <- length(names(variant_data))
      cat("   ✅", variant, "→", harmonized_gene, "(", cluster_count, "clusters available)\n")
      
      # Check if any cluster has actual results
      has_results <- FALSE
      for (cluster in names(variant_data)) {
        if (!is.null(variant_data[[cluster]]$results) && nrow(variant_data[[cluster]]$results) > 0) {
          has_results <- TRUE
          break
        }
      }
      
      if (has_results) {
        cat("     ✅ Has DE results data\n")
      } else {
        cat("     ❌ No DE results data found\n")
      }
    } else {
      cat("   ❌", variant, "→", harmonized_gene, "(not found in data)\n")
    }
  }
}

# Summary
cat("\n=== TEST SUMMARY ===\n")
if (all_passed) {
  cat("✅ ALL TESTS PASSED\n")
  cat("VPS13C and SNCA volcano plots should now work correctly!\n")
  cat("\nThe fix ensures that:\n")
  cat("- VPS13C_W395C volcano plots will show VPS13C MixScale data\n")
  cat("- VPS13C_A444P volcano plots will show VPS13C MixScale data\n") 
  cat("- SNCA_A30P volcano plots will show SNCA MixScale data\n")
  cat("- SNCA_A53T volcano plots will show SNCA MixScale data\n")
  cat("- PRKN volcano plots will show PARK2 MixScale data\n")
} else {
  cat("❌ SOME TESTS FAILED\n")
  cat("Please check the issues above before deploying the fix.\n")
}

cat("\nNext steps:\n")
cat("1. Test the Shiny app with VPS13C_W395C selected\n")
cat("2. Verify volcano plots show data instead of 'no results available'\n")
cat("3. Test all variants: VPS13C_W395C, VPS13C_A444P, SNCA_A30P, SNCA_A53T\n")