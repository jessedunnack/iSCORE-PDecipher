#!/usr/bin/env Rscript

# Install missing packages
missing_pkgs <- c("pathview", "heatmaply", "colourpicker")

for (pkg in missing_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

# Also check BiocManager packages
if (!requireNamespace("pathview", quietly = TRUE)) {
  BiocManager::install("pathview")
}

cat("All packages installed\n")