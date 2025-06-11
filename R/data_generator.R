# Data generation functions for iSCORE-PDecipher
# These functions run the processing scripts to generate missing files

#' Run data generation pipeline for missing files
#'
#' @param data_dir Path to dataset directory
#' @param missing Vector of missing file types
#' @return TRUE if successful, FALSE otherwise
generate_missing_files <- function(data_dir, missing) {
  success <- TRUE
  
  # Get the package directory
  pkg_dir <- system.file(package = "iSCORE.PDecipher")
  if (pkg_dir == "") {
    # If package not installed, use relative path
    pkg_dir <- dirname(dirname(getwd()))
  }
  
  # Path to scripts
  script_dir <- file.path(pkg_dir, "inst", "scripts")
  
  # Ensure script directory exists
  if (!dir.exists(script_dir)) {
    # Try alternative location
    script_dir <- file.path(dirname(dirname(getwd())), "package_materials")
    if (!dir.exists(script_dir)) {
      stop("Cannot find processing scripts. Please ensure package is properly installed.")
    }
  }
  
  # Save current working directory
  old_wd <- getwd()
  on.exit(setwd(old_wd))
  
  # Change to data directory
  setwd(data_dir)
  
  # Step 1: Generate full_DE_results.rds if missing
  if ("full_DE_results" %in% missing) {
    message("\n=== STEP 1: Compiling differential expression results ===")
    
    # Ask for confirmation
    response <- readline("This will compile MAST and MixScale results (~5-10 minutes). Continue? (Y/N): ")
    if (!tolower(response) %in% c("y", "yes")) {
      message("Skipping DE compilation.")
      return(FALSE)
    }
    
    tryCatch({
      # Source and run data import
      source(file.path(script_dir, "data_import_functions.R"))
      
      # Initialize results list
      full_results <- list()
      
      # Import MAST data if available
      mast_dir <- file.path(data_dir, "iSCORE-PD_MAST_analysis")
      if (dir.exists(mast_dir)) {
        message("Importing MAST data...")
        full_results$iSCORE_PD_MAST <- import_mast_data(mast_dir)
        message("  ✓ Imported ", length(names(full_results$iSCORE_PD_MAST)), " mutations")
      }
      
      # Import MixScale data if available
      # Check multiple possible directory patterns
      mixscale_dirs <- list.dirs(data_dir, recursive = FALSE)
      mixscale_pattern <- "PerturbSeq_MixScale_analysis"
      matching_dirs <- mixscale_dirs[grep(mixscale_pattern, basename(mixscale_dirs))]
      
      if (length(matching_dirs) > 0) {
        mixscale_dir <- matching_dirs[1]
        
        # Check for CRISPRi
        if (dir.exists(file.path(mixscale_dir, "CRISPRi"))) {
          message("Importing CRISPRi MixScale data...")
          full_results$CRISPRi_Mixscale <- import_mixscale_data(mixscale_dir, modality = "CRISPRi")
          message("  ✓ Imported ", length(names(full_results$CRISPRi_Mixscale)), " perturbations")
        }
        
        # Check for CRISPRa
        if (dir.exists(file.path(mixscale_dir, "CRISPRa"))) {
          message("Importing CRISPRa MixScale data...")
          full_results$CRISPRa_Mixscale <- import_mixscale_data(mixscale_dir, modality = "CRISPRa")
          message("  ✓ Imported ", length(names(full_results$CRISPRa_Mixscale)), " perturbations")
        }
      }
      
      # Save results
      saveRDS(full_results, "full_DE_results.rds")
      message("✓ Successfully created full_DE_results.rds")
      
    }, error = function(e) {
      message("✗ Error creating full_DE_results.rds: ", e$message)
      success <<- FALSE
    })
    
    if (!success) return(FALSE)
  }
  
  # Step 2: Run enrichment analysis if missing
  if ("enrichment_results" %in% missing) {
    message("\n=== STEP 2: Running functional enrichment analysis ===")
    message("⚠️  WARNING: This process takes SEVERAL HOURS!")
    
    # Ask for confirmation
    response <- readline("Continue with enrichment analysis? (Y/N): ")
    if (!tolower(response) %in% c("y", "yes")) {
      message("Skipping enrichment analysis.")
      return(FALSE)
    }
    
    tryCatch({
      # Source enrichment functions
      source(file.path(script_dir, "enrichment_analysis.R"))
      
      # Run enrichment analysis
      results <- run_enrichment_analysis(
        input_file = "full_DE_results.rds",
        lfc_threshold = 0.25,
        padj_threshold = 0.05,
        output_dir = "./enrichment_results/",
        run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
        min_genes = 10,
        padj_method = "BH"
      )
      
      message("✓ Successfully completed enrichment analysis")
      
    }, error = function(e) {
      message("✗ Error running enrichment analysis: ", e$message)
      success <<- FALSE
    })
    
    if (!success) return(FALSE)
  }
  
  # Step 3: Process enrichment results if missing
  if ("consolidated_enrichment" %in% missing) {
    message("\n=== STEP 3: Consolidating enrichment results ===")
    
    # Ask for confirmation
    response <- readline("This will consolidate and filter enrichment results (~10-30 minutes). Continue? (Y/N): ")
    if (!tolower(response) %in% c("y", "yes")) {
      message("Skipping result consolidation.")
      return(FALSE)
    }
    
    tryCatch({
      # Source processing functions
      source(file.path(script_dir, "process_enrichment_results.R"))
      
      # Process results
      results <- process_enrichment_results(
        root_dir = "enrichment_results",
        batch_size = 100,
        output_file = "all_enrichment_padj005_complete_with_direction.rds"
      )
      
      message("✓ Successfully consolidated enrichment results")
      
    }, error = function(e) {
      message("✗ Error processing enrichment results: ", e$message)
      success <<- FALSE
    })
    
    if (!success) return(FALSE)
  }
  
  return(success)
}

#' Check if required packages are installed
#'
#' @return TRUE if all packages available, FALSE otherwise
check_required_packages <- function() {
  required_packages <- c(
    # For data import
    "dplyr", "tidyr", "readr",
    # For enrichment analysis
    "clusterProfiler", "enrichplot", "org.Hs.eg.db", 
    "DOSE", "ReactomePA", "pathview", "STRINGdb", 
    "msigdbr", "fgsea",
    # For visualization
    "ggplot2", "pheatmap", "plotly"
  )
  
  missing_packages <- character()
  
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }
  
  if (length(missing_packages) > 0) {
    message("Missing required packages:")
    message(paste("  -", missing_packages, collapse = "\n"))
    message("\nInstall with:")
    message("BiocManager::install(c(", paste0('"', missing_packages, '"', collapse = ", "), "))")
    return(FALSE)
  }
  
  return(TRUE)
}