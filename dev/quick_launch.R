# Quick launch script to bypass installation issues
cat("Loading iSCORE-PDecipher functions directly...\n")

# Source the main launch function
source("R/launch_app.R")
source("R/config_manager.R")

cat("Functions loaded. You can now run:\n")
cat("launch_iscore_app()\n\n")

# Try to launch directly
cat("Launching app...\n")
launch_iscore_app()