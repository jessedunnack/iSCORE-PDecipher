#!/usr/bin/env Rscript

# Reinstall iSCORE.PDecipher package with latest v0.2.6 changes

cat("=== REINSTALLING iSCORE.PDecipher v0.2.6 ===\n\n")

# Remove old package if installed
if ("iSCORE.PDecipher" %in% installed.packages()[,"Package"]) {
  cat("Removing old package version...\n")
  remove.packages("iSCORE.PDecipher")
}

# Unload package if loaded
if ("package:iSCORE.PDecipher" %in% search()) {
  cat("Detaching loaded package...\n")
  detach("package:iSCORE.PDecipher", unload = TRUE)
}

# Install remotes if needed
if (!requireNamespace("remotes", quietly = TRUE)) {
  cat("Installing remotes package...\n")
  install.packages("remotes", repos = "https://cran.rstudio.com/")
}

# Install latest package from GitHub
cat("Installing latest iSCORE.PDecipher from GitHub...\n")
remotes::install_github("jessedunnack/iSCORE-PDecipher", force = TRUE)

# Test installation
cat("\nTesting installation...\n")
library(iSCORE.PDecipher)

# Check if new function parameter is available
if ("use_enhanced_analysis" %in% names(formals(discover_top_signatures))) {
  cat("✅ SUCCESS: Enhanced analysis parameter detected\n")
  cat("✅ Package installation completed successfully\n")
  cat("\nYou can now run:\n")
  cat("  library(iSCORE.PDecipher)\n")
  cat("  launch_iscore_app()\n")
} else {
  cat("❌ ERROR: Enhanced analysis parameter not found\n")
  cat("❌ Package may not have installed correctly\n")
}

cat("\n=== INSTALLATION COMPLETE ===\n")