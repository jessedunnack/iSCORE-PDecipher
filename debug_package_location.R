#!/usr/bin/env Rscript

# Debug where the package is actually loading files from

library(iSCORE.PDecipher)

cat("=== PACKAGE LOCATION DEBUGGING ===\n\n")

cat("1. PACKAGE INSTALLATION LOCATION:\n")
pkg_path <- system.file(package = "iSCORE.PDecipher")
cat("   Package root:", pkg_path, "\n")

cat("\n2. UMAP FILE LOCATIONS:\n")
umap_file_30pc <- system.file("extdata", "umap_data", "iSCORE_PD_CRISPRi_umap_data_30pc.rds", 
                              package = "iSCORE.PDecipher")
cat("   30pc file path:", umap_file_30pc, "\n")
cat("   30pc file exists:", file.exists(umap_file_30pc), "\n")

if (file.exists(umap_file_30pc)) {
  cat("   30pc file size:", file.size(umap_file_30pc), "bytes\n")
  cat("   30pc file modified:", file.mtime(umap_file_30pc), "\n")
}

cat("\n3. SHINY APP LOCATION:\n")
shiny_path <- system.file("shiny", package = "iSCORE.PDecipher")
cat("   Shiny app path:", shiny_path, "\n")
cat("   Shiny app exists:", dir.exists(shiny_path), "\n")

cat("\n4. LOCAL DEVELOPMENT FILES:\n")
local_umap <- "inst/extdata/umap_data/iSCORE_PD_CRISPRi_umap_data_30pc.rds"
cat("   Local 30pc file exists:", file.exists(local_umap), "\n")
if (file.exists(local_umap)) {
  cat("   Local 30pc file size:", file.size(local_umap), "bytes\n")
  cat("   Local 30pc file modified:", file.mtime(local_umap), "\n")
}

cat("\nDone.\n")