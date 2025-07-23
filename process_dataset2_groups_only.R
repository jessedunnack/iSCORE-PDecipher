#!/usr/bin/env Rscript

# Script to process only group comparison markers for dataset 2
# Skips cluster markers and UMAP extraction

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(optparse)
})

# Parse command line arguments
option_list <- list(
  make_option(c("-o", "--output-dir"), type="character", 
              default="E:/scRNSeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data",
              help="Output directory for results"),
  
  # FindMarkers parameters for group comparisons
  make_option("--group-test", type="character", default="MAST",
              help="Test to use for FindMarkers: MAST, wilcox, bimod, roc, t, negbinom, poisson, LR, DESeq2 (default: MAST)"),
  make_option("--group-logfc", type="numeric", default=0.25,
              help="Log fold change threshold for FindMarkers (default: 0.25)"),
  make_option("--group-min-pct", type="numeric", default=0.1,
              help="Min percent cells expressing for FindMarkers (default: 0.1)"),
  make_option("--group-min-diff-pct", type="numeric", default=-Inf,
              help="Min difference in percent expressed for FindMarkers (default: -Inf)"),
  make_option("--group-only-pos", type="logical", default=FALSE,
              help="Only return positive markers in FindMarkers (default: FALSE)"),
  make_option("--group-max-cells", type="integer", default=500,
              help="Max cells per ident for FindMarkers, use -1 for Inf (default: 500)"),
  
  make_option("--verbose", type="logical", default=TRUE,
              help="Print detailed progress (default: TRUE)")
)

opt_parser <- OptionParser(option_list=option_list,
                          description="Process group comparison markers for iSCORE_PD_plus_CRISPRi")
opt <- parse_args(opt_parser)

# Configuration
DATASET_NAME <- "iSCORE_PD_CRISPRi"
SEURAT_FILE <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"
OUTPUT_DIR <- opt$`output-dir`

cat("\n=== Group Comparison Markers Processing ===\n")
cat("Working directory:", getwd(), "\n")
cat("Seurat file:", SEURAT_FILE, "\n")
cat("Output directory:", OUTPUT_DIR, "\n\n")

cat("FindMarkers parameters:\n")
cat(sprintf("  Test: %s\n", opt$`group-test`))
cat(sprintf("  Log FC threshold: %.2f\n", opt$`group-logfc`))
cat(sprintf("  Min percent: %.2f\n", opt$`group-min-pct`))
cat(sprintf("  Min diff percent: %s\n", 
            ifelse(is.infinite(opt$`group-min-diff-pct`), "-Inf", sprintf("%.2f", opt$`group-min-diff-pct`))))
cat(sprintf("  Only positive: %s\n", opt$`group-only-pos`))
cat(sprintf("  Max cells: %s\n", ifelse(opt$`group-max-cells` == -1, "all", opt$`group-max-cells`)))
cat("\n")

# Main processing
tryCatch({
  # Check if Seurat file exists
  if (!file.exists(SEURAT_FILE)) {
    stop(sprintf("Seurat file not found: %s", SEURAT_FILE))
  }
  
  cat(sprintf("Loading Seurat object from: %s\n", SEURAT_FILE))
  cat("This may take a minute...\n")
  
  seurat_obj <- readRDS(SEURAT_FILE)
  
  cat(sprintf("\nLoaded dataset with %d cells and %d genes\n", 
              ncol(seurat_obj), nrow(seurat_obj)))
  
  # Check metadata columns
  meta_cols <- colnames(seurat_obj@meta.data)
  cat("\nAvailable metadata columns:\n")
  print(meta_cols)
  
  # Function to calculate group markers with custom parameters
  calculate_group_markers_custom <- function(seurat_obj, opt) {
    cat("\n=== Calculating Group Comparison Markers with FindMarkers ===\n")
    
    group_markers <- list()
    
    # 1. iSCORE-PD: eWT vs all mutants (SKIP if mutation_tidy doesn't exist)
    if ("mutation_tidy" %in% colnames(seurat_obj@meta.data)) {
      cat("\n1. iSCORE-PD comparison: eWT vs all mutants\n")
      
      # Check for non-NA values
      non_na_mutation <- !is.na(seurat_obj@meta.data$mutation_tidy)
      if (sum(non_na_mutation) > 0) {
        seurat_obj@meta.data$iscore_group <- ifelse(
          seurat_obj@meta.data$mutation_tidy == "eWT", "eWT", "mutant"
        )
        
        ewt_cells <- sum(seurat_obj@meta.data$iscore_group == "eWT", na.rm = TRUE)
        mut_cells <- sum(seurat_obj@meta.data$iscore_group == "mutant", na.rm = TRUE)
        cat(sprintf("   - eWT cells: %d\n", ewt_cells))
        cat(sprintf("   - Mutant cells: %d\n", mut_cells))
        
        if (ewt_cells > 10 && mut_cells > 10) {
          # Subset to cells with valid iscore_group
          valid_cells <- !is.na(seurat_obj@meta.data$iscore_group)
          seurat_subset <- seurat_obj[, valid_cells]
          Idents(seurat_subset) <- "iscore_group"
          
          markers <- tryCatch({
            FindMarkers(seurat_subset,
                       ident.1 = "mutant",
                       ident.2 = "eWT",
                       test.use = opt$`group-test`,
                       logfc.threshold = opt$`group-logfc`,
                       min.pct = opt$`group-min-pct`,
                       min.diff.pct = opt$`group-min-diff-pct`,
                       only.pos = opt$`group-only-pos`,
                       max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                       verbose = FALSE)
          }, error = function(e) {
            warning(sprintf("Failed: %s", e$message))
            return(NULL)
          })
          
          if (!is.null(markers)) {
            markers$gene <- rownames(markers)
            markers$comparison <- "iSCORE_mutants_vs_eWT"
            group_markers[["iSCORE_mutants_vs_eWT"]] <- markers
            cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
          }
        } else {
          cat("   - Skipping: insufficient cells in one or both groups\n")
        }
      } else {
        cat("   - Skipping: all mutation_tidy values are NA\n")
      }
    } else {
      cat("\n1. Skipping iSCORE-PD comparison (no mutation_tidy column)\n")
    }
    
    # 2. CRISPRi: Non-targeting vs perturbed
    if ("scMAGeCK_gene_assignment" %in% colnames(seurat_obj@meta.data)) {
      cat("\n2. CRISPRi comparison: Non-targeting vs perturbed cells\n")
      
      # Check for non-NA values
      non_na_crispr <- !is.na(seurat_obj@meta.data$scMAGeCK_gene_assignment)
      if (sum(non_na_crispr) > 0) {
        seurat_obj@meta.data$crispr_group <- ifelse(
          seurat_obj@meta.data$scMAGeCK_gene_assignment == "Non-Targeting", 
          "control", "perturbed"
        )
        
        ctrl_cells <- sum(seurat_obj@meta.data$crispr_group == "control", na.rm = TRUE)
        pert_cells <- sum(seurat_obj@meta.data$crispr_group == "perturbed", na.rm = TRUE)
        cat(sprintf("   - Control cells: %d\n", ctrl_cells))
        cat(sprintf("   - Perturbed cells: %d\n", pert_cells))
        
        if (ctrl_cells > 10 && pert_cells > 10) {
          # Subset to cells with valid crispr_group
          valid_cells <- !is.na(seurat_obj@meta.data$crispr_group)
          seurat_subset <- seurat_obj[, valid_cells]
          Idents(seurat_subset) <- "crispr_group"
          
          markers <- tryCatch({
            FindMarkers(seurat_subset,
                       ident.1 = "perturbed",
                       ident.2 = "control",
                       test.use = opt$`group-test`,
                       logfc.threshold = opt$`group-logfc`,
                       min.pct = opt$`group-min-pct`,
                       min.diff.pct = opt$`group-min-diff-pct`,
                       only.pos = opt$`group-only-pos`,
                       max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                       verbose = FALSE)
          }, error = function(e) {
            warning(sprintf("Failed: %s", e$message))
            return(NULL)
          })
          
          if (!is.null(markers)) {
            markers$gene <- rownames(markers)
            markers$comparison <- "CRISPRi_perturbed_vs_control"
            group_markers[["CRISPRi_perturbed_vs_control"]] <- markers
            cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
          }
        } else {
          cat("   - Skipping: insufficient cells in one or both groups\n")
        }
      } else {
        cat("   - Skipping: all scMAGeCK_gene_assignment values are NA\n")
      }
    } else {
      cat("\n2. Skipping CRISPRi comparison (no scMAGeCK_gene_assignment column)\n")
    }
    
    # 3. Cross-method comparison (only if both metadata exist and have non-NA values)
    if ("mutation_tidy" %in% colnames(seurat_obj@meta.data) && 
        "scMAGeCK_gene_assignment" %in% colnames(seurat_obj@meta.data)) {
      cat("\n3. Cross-method comparison: iSCORE-PD vs PerturbSeq cells\n")
      
      # Create method groups, handling NAs carefully
      seurat_obj@meta.data$method_group <- NA
      
      # Assign iSCORE to cells with non-NA mutation_tidy
      iscore_mask <- !is.na(seurat_obj@meta.data$mutation_tidy)
      seurat_obj@meta.data$method_group[iscore_mask] <- "iSCORE"
      
      # Assign PerturbSeq to cells with non-NA scMAGeCK_gene_assignment (overwrites if both exist)
      perturb_mask <- !is.na(seurat_obj@meta.data$scMAGeCK_gene_assignment)
      seurat_obj@meta.data$method_group[perturb_mask] <- "PerturbSeq"
      
      iscore_cells <- sum(seurat_obj@meta.data$method_group == "iSCORE", na.rm = TRUE)
      perturb_cells <- sum(seurat_obj@meta.data$method_group == "PerturbSeq", na.rm = TRUE)
      cat(sprintf("   - iSCORE cells: %d\n", iscore_cells))
      cat(sprintf("   - PerturbSeq cells: %d\n", perturb_cells))
      
      if (iscore_cells > 10 && perturb_cells > 10) {
        valid_cells <- !is.na(seurat_obj@meta.data$method_group)
        seurat_subset <- seurat_obj[, valid_cells]
        Idents(seurat_subset) <- "method_group"
        
        markers <- tryCatch({
          FindMarkers(seurat_subset,
                     ident.1 = "PerturbSeq",
                     ident.2 = "iSCORE",
                     test.use = opt$`group-test`,
                     logfc.threshold = opt$`group-logfc`,
                     min.pct = opt$`group-min-pct`,
                     min.diff.pct = opt$`group-min-diff-pct`,
                     only.pos = opt$`group-only-pos`,
                     max.cells.per.ident = ifelse(opt$`group-max-cells` == -1, Inf, opt$`group-max-cells`),
                     verbose = FALSE)
        }, error = function(e) {
          warning(sprintf("Failed: %s", e$message))
          return(NULL)
        })
        
        if (!is.null(markers)) {
          markers$gene <- rownames(markers)
          markers$comparison <- "PerturbSeq_vs_iSCORE"
          group_markers[["PerturbSeq_vs_iSCORE"]] <- markers
          cat(sprintf("   - Found %d DE genes\n", nrow(markers)))
        }
      } else {
        cat("   - Skipping: insufficient cells in one or both groups\n")
      }
    }
    
    if (length(group_markers) > 0) {
      return(do.call(rbind, group_markers))
    } else {
      return(NULL)
    }
  }
  
  # Calculate group markers
  group_markers <- calculate_group_markers_custom(seurat_obj, opt)
  
  if (!is.null(group_markers)) {
    group_file <- file.path(OUTPUT_DIR, sprintf("%s_group_comparison_markers.rds", DATASET_NAME))
    saveRDS(group_markers, group_file)
    cat(sprintf("\nSaved group markers to: %s\n", basename(group_file)))
  } else {
    cat("\nNo group markers were successfully calculated.\n")
  }
  
  cat("\n=== PROCESSING COMPLETE ===\n")
  
}, error = function(e) {
  cat(sprintf("\nERROR: %s\n", e$message))
  quit(status = 1)
})