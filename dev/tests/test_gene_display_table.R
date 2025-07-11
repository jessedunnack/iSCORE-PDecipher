# Test Gene Display in Data Table
# Run this after restarting the Shiny app

cat("=== Testing Gene Display in Plot Details ===\n")
cat("1. Launch the Shiny app\n")
cat("2. Go to Functional Enrichment Results\n")
cat("3. Click on the Plot Details tab\n")
cat("4. Check for Associated_Genes column in the table\n")
cat("5. Verify genes are displayed for each term\n")
cat("\nExpected behavior:\n")
cat("- Table should have an Associated_Genes column\n")
cat("- Each term shows its DE genes (up to 20 + count of remaining)\n")
cat("- Gene text is blue and slightly smaller\n")
cat("- Table is horizontally scrollable\n")

