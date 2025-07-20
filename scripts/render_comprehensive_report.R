#!/usr/bin/env Rscript
# Render the comprehensive PD signature report

cat("Rendering comprehensive PD signature report...\n")

# Check if pandoc is available
if (Sys.which("pandoc") == "") {
  cat("WARNING: pandoc not found. Trying to render anyway...\n")
}

# Set working directory
setwd("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/scripts")

# Try to render
tryCatch({
  rmarkdown::render(
    "pd_signature_comprehensive_report.Rmd",
    output_file = "../results/pd_signatures/comprehensive/pd_signature_comprehensive_report.html",
    output_format = "html_document"
  )
  cat("Report successfully rendered!\n")
  cat("Location: /mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/comprehensive/pd_signature_comprehensive_report.html\n")
}, error = function(e) {
  cat("Error rendering report:\n")
  print(e)
  cat("\nThe report Rmd file has been created successfully.\n")
  cat("You can render it manually in RStudio if needed.\n")
})