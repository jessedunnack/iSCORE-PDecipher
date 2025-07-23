#!/usr/bin/env Rscript

# Script to help push to GitHub using HTTPS
# You'll need a GitHub Personal Access Token (PAT)

cat("=== GitHub Push Helper ===\n\n")

cat("To push your changes to GitHub, you need a Personal Access Token (PAT).\n")
cat("If you don't have one, create it at: https://github.com/settings/tokens\n")
cat("Make sure to give it 'repo' permissions.\n\n")

# Check if using RStudio
if (Sys.getenv("RSTUDIO") == "1") {
  cat("Detected RStudio environment.\n\n")
  cat("Option 1: Use RStudio's Git pane\n")
  cat("- Click the Git tab in RStudio\n")
  cat("- Click the 'Push' button (up arrow)\n")
  cat("- Enter your GitHub username and PAT when prompted\n\n")
  
  cat("Option 2: Use credentials package\n")
  cat("Run these commands in the R console:\n\n")
  cat("install.packages('credentials')\n")
  cat("credentials::set_github_pat()\n")
  cat("# Enter your PAT when prompted\n")
  cat("system('git push origin main')\n\n")
  
} else {
  cat("Option 1: Set credentials temporarily (one-time push):\n")
  cat("Run this in R, replacing YOUR_USERNAME and YOUR_PAT:\n\n")
  cat('Sys.setenv(GIT_ASKPASS = "echo")\n')
  cat('system("git push https://YOUR_USERNAME:YOUR_PAT@github.com/jessedunnack/iSCORE-PDecipher.git main")\n\n')
  
  cat("Option 2: Use git credential helper:\n")
  cat('system("git config --global credential.helper cache")\n')
  cat('system("git push origin main")\n')
  cat("# Enter username and PAT when prompted\n\n")
}

cat("Option 3: Use gitcreds package (recommended for long-term use):\n")
cat("install.packages('gitcreds')\n")
cat("gitcreds::gitcreds_set()\n")
cat("# Follow the prompts to set your PAT\n\n")

cat("After setting credentials, you can push with:\n")
cat("system('git push origin main')\n")