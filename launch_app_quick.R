#!/usr/bin/env Rscript

# Quick Launch for iSCORE-PDecipher with Dataset Switcher
# Streamlined version that goes straight to the dataset switcher app

cat("
╔══════════════════════════════════════════════════════════════╗
║           iSCORE-PDecipher Analysis Platform                ║
║     Integrated Analysis of PD Mutations & Perturbations     ║
╚══════════════════════════════════════════════════════════════╝
")

# Check datasets
base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/"
ready_count <- 0

cat("\n📊 Checking Datasets:\n")
cat("─────────────────────────────────────────────────────────\n")

# Check iSCORE-PD
iscore_pd_ready <- file.exists(file.path(base_dir, "iSCORE-PD/full_DE_results.rds")) &&
                   file.exists(file.path(base_dir, "iSCORE-PD/all_enrichment_padj005_complete_with_direction.rds"))

if (iscore_pd_ready) {
  cat("  ✅ iSCORE-PD (Mutations Only) - 210K cells\n")
  ready_count <- ready_count + 1
} else {
  cat("  ❌ iSCORE-PD - Missing data files\n")
}

# Check iSCORE-PD + CRISPRi
iscore_crispri_ready <- file.exists(file.path(base_dir, "iSCORE-PD_plus_CRISPRi/full_DE_results.rds")) &&
                        file.exists(file.path(base_dir, "iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"))

if (iscore_crispri_ready) {
  cat("  ✅ iSCORE-PD + CRISPRi - 232K cells\n")
  ready_count <- ready_count + 1
} else {
  cat("  ❌ iSCORE-PD + CRISPRi - Missing data files\n")
}

cat("─────────────────────────────────────────────────────────\n")

if (ready_count == 0) {
  cat("\n❌ No complete datasets found!\n")
  stop("No datasets available for analysis")
}

cat(sprintf("\n✅ %d dataset(s) ready for analysis\n", ready_count))

# Set default dataset (prefer CRISPRi if available)
if (iscore_crispri_ready) {
  default_path <- file.path(base_dir, "iSCORE-PD_plus_CRISPRi")
  cat("📍 Default dataset: iSCORE-PD + CRISPRi\n")
} else {
  default_path <- file.path(base_dir, "iSCORE-PD")
  cat("📍 Default dataset: iSCORE-PD\n")
}

# Set initial environment variables
Sys.setenv(ISCORE_DATA_DIR = default_path)
Sys.setenv(ISCORE_DE_FILE = file.path(default_path, "full_DE_results.rds"))
Sys.setenv(ISCORE_ENRICHMENT_FILE = file.path(default_path, "all_enrichment_padj005_complete_with_direction.rds"))
Sys.setenv(ISCORE_ENRICHMENT_DIR = file.path(default_path, "enrichment_results"))

# Launch the app with dataset switcher
cat("\n🚀 Launching app with dataset switcher...\n")
cat("   • Switch between datasets using the 'Change Dataset' button\n")
cat("   • No restart required when switching\n")
cat("   • Access at: http://localhost:3838\n\n")

# Check if running in RStudio
if (Sys.getenv("RSTUDIO") == "1") {
  cat("📝 Detected RStudio environment\n")
  cat("   The app will open in the Viewer pane\n\n")
}

# Set app directory
app_dir <- "inst/shiny"
app_file <- file.path(app_dir, "app_with_dataset_selector.R")

# Check if enhanced app exists, fallback to standard if not
if (!file.exists(app_file)) {
  cat("⚠️  Dataset switcher app not found, using standard app\n")
  app_file <- file.path(app_dir, "app.R")
}

# Launch the app
tryCatch({
  shiny::runApp(
    app_file,
    host = "0.0.0.0",
    port = 3838,
    launch.browser = interactive() && Sys.getenv("RSTUDIO") == "1"
  )
}, error = function(e) {
  cat("\n❌ Error launching app:", e$message, "\n")
  cat("\nTroubleshooting:\n")
  cat("1. Make sure you're in the iSCORE-PDecipher directory\n")
  cat("2. Check that required packages are installed\n")
  cat("3. Try: setwd('/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher')\n")
})