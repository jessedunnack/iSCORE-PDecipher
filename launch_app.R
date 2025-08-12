#!/usr/bin/env Rscript

# Launch iSCORE-PDecipher Shiny App with Dataset Selection
# This script provides an easy way to launch the app with the desired dataset

cat("
╔══════════════════════════════════════════════════════════════╗
║           iSCORE-PDecipher Analysis Platform                ║
║     Integrated Analysis of PD Mutations & Perturbations     ║
╚══════════════════════════════════════════════════════════════╝
")

# Check for required packages
required_packages <- c("shiny", "shinyjs", "DT", "plotly", "dplyr")
missing_packages <- required_packages[!required_packages %in% installed.packages()[,"Package"]]

if (length(missing_packages) > 0) {
  cat("\n⚠️  Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Installing missing packages...\n")
  install.packages(missing_packages, repos = "https://cran.r-project.org")
}

# Base directory
base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/"

# Check available datasets
datasets <- list()

# Check iSCORE-PD
iscore_pd_dir <- file.path(base_dir, "iSCORE-PD")
if (dir.exists(iscore_pd_dir)) {
  has_de <- file.exists(file.path(iscore_pd_dir, "full_DE_results.rds"))
  has_enrich <- file.exists(file.path(iscore_pd_dir, "all_enrichment_padj005_complete_with_direction.rds"))
  
  if (has_de && has_enrich) {
    datasets$iscore_pd <- list(
      name = "iSCORE-PD (Mutations Only)",
      path = iscore_pd_dir,
      complete = TRUE
    )
  } else {
    datasets$iscore_pd <- list(
      name = "iSCORE-PD (Mutations Only)",
      path = iscore_pd_dir,
      complete = FALSE,
      missing = paste(
        if(!has_de) "DE results" else NULL,
        if(!has_enrich) "Enrichment data" else NULL,
        collapse = ", "
      )
    )
  }
}

# Check iSCORE-PD + CRISPRi
iscore_crispri_dir <- file.path(base_dir, "iSCORE-PD_plus_CRISPRi")
if (dir.exists(iscore_crispri_dir)) {
  has_de <- file.exists(file.path(iscore_crispri_dir, "full_DE_results.rds"))
  has_enrich <- file.exists(file.path(iscore_crispri_dir, "all_enrichment_padj005_complete_with_direction.rds"))
  
  if (has_de && has_enrich) {
    datasets$iscore_pd_plus_crispri <- list(
      name = "iSCORE-PD + CRISPRi",
      path = iscore_crispri_dir,
      complete = TRUE
    )
  } else {
    datasets$iscore_pd_plus_crispri <- list(
      name = "iSCORE-PD + CRISPRi",
      path = iscore_crispri_dir,
      complete = FALSE,
      missing = paste(
        if(!has_de) "DE results" else NULL,
        if(!has_enrich) "Enrichment data" else NULL,
        collapse = ", "
      )
    )
  }
}

# Display available datasets
cat("\n📊 Available Datasets:\n")
cat("─────────────────────────────────────────────────────────\n")

available_count <- 0
for (key in names(datasets)) {
  dataset <- datasets[[key]]
  if (dataset$complete) {
    available_count <- available_count + 1
    cat(sprintf("  %d. %s ✅\n", available_count, dataset$name))
    cat(sprintf("     Path: %s\n", dataset$path))
  } else {
    cat(sprintf("  ❌ %s (Missing: %s)\n", dataset$name, dataset$missing))
  }
}

if (available_count == 0) {
  cat("\n❌ No complete datasets found!\n")
  cat("Please ensure both DE results and enrichment data are available.\n")
  stop("No datasets available for analysis")
}

cat("\n─────────────────────────────────────────────────────────\n")

# Interactive selection if multiple datasets
selected_dataset <- NULL

if (available_count == 1) {
  # Auto-select if only one dataset
  for (key in names(datasets)) {
    if (datasets[[key]]$complete) {
      selected_dataset <- key
      cat("\n✅ Auto-selected:", datasets[[key]]$name, "\n")
    }
  }
} else {
  # Ask user to select
  cat("\nSelect a dataset to analyze (1-", available_count, "): ", sep = "")
  choice <- as.integer(readline())
  
  # Map choice to dataset
  idx <- 0
  for (key in names(datasets)) {
    if (datasets[[key]]$complete) {
      idx <- idx + 1
      if (idx == choice) {
        selected_dataset <- key
        break
      }
    }
  }
  
  if (is.null(selected_dataset)) {
    cat("\n❌ Invalid selection\n")
    stop("Invalid dataset selection")
  }
}

# Set environment variables for selected dataset
dataset <- datasets[[selected_dataset]]
Sys.setenv(ISCORE_DATA_DIR = dataset$path)
Sys.setenv(ISCORE_DE_FILE = file.path(dataset$path, "full_DE_results.rds"))
Sys.setenv(ISCORE_ENRICHMENT_FILE = file.path(dataset$path, "all_enrichment_padj005_complete_with_direction.rds"))
Sys.setenv(ISCORE_ENRICHMENT_DIR = file.path(dataset$path, "enrichment_results"))

cat("\n🚀 Launching app with:", dataset$name, "\n")

# Launch options
cat("\nLaunch Options:\n")
cat("  1. Standard app (current version)\n")
cat("  2. App with dataset switcher (NEW)\n")
cat("  3. Development mode (with browser)\n")
cat("\nSelect option (1-3): ")

launch_choice <- as.integer(readline())

# Set app directory
app_dir <- system.file("shiny", package = "iSCORE.PDecipher")
if (app_dir == "") {
  # Fallback to local directory if package not installed
  app_dir <- "inst/shiny"
}

# Launch based on choice
if (launch_choice == 2) {
  # Launch with dataset selector
  cat("\n🎯 Starting app with dataset selector...\n")
  if (file.exists(file.path(app_dir, "app_with_dataset_selector.R"))) {
    shiny::runApp(file.path(app_dir, "app_with_dataset_selector.R"), 
                 host = "0.0.0.0", 
                 port = 3838)
  } else {
    cat("Dataset selector app not found. Falling back to standard app.\n")
    shiny::runApp(app_dir, host = "0.0.0.0", port = 3838)
  }
} else if (launch_choice == 3) {
  # Development mode
  cat("\n🔧 Starting in development mode...\n")
  shiny::runApp(app_dir, 
               host = "127.0.0.1", 
               port = 3838,
               launch.browser = TRUE)
} else {
  # Standard launch
  cat("\n📊 Starting standard app...\n")
  cat("Access the app at: http://localhost:3838\n\n")
  shiny::runApp(app_dir, host = "0.0.0.0", port = 3838)
}