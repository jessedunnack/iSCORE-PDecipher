# Install Latest Version with Signature Nomination Module
# Run this script to update your package installation

cat("=== Installing Latest iSCORE-PDecipher with Signature Nomination ===\n")

# First, remove the old version if it exists
if ("iSCORE.PDecipher" %in% rownames(installed.packages())) {
  cat("Removing previous version...\n")
  remove.packages("iSCORE.PDecipher")
}

# Clear any loaded package namespaces
if ("iSCORE.PDecipher" %in% loadedNamespaces()) {
  cat("Unloading package namespace...\n")
  try(unloadNamespace("iSCORE.PDecipher"), silent = TRUE)
}

# Install fresh from the current development directory
cat("Installing package from current directory...\n")
install.packages(".", repos = NULL, type = "source")

cat("=== Installation Complete ===\n")
cat("You can now run: library(iSCORE.PDecipher); launch_iscore_app()\n")
cat("The 'Signature Nomination' tab should now appear!\n")