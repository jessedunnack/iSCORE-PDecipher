#' Run Complete iSCORE-PDecipher Analysis Pipeline
#'
#' This function executes the complete analysis pipeline from differential expression
#' results through functional enrichment analysis and result consolidation.
#'
#' @param mast_directory Path to directory containing MAST analysis results
#' @param mixscale_directory Path to directory containing MixScale analysis results
#' @param output_directory Path to directory for saving output files
#' @param lfc_threshold Log2 fold change threshold for enrichment analysis (default: 0.25)
#' @param padj_threshold Adjusted p-value threshold for enrichment analysis (default: 0.05)
#' @param run_methods Vector of enrichment methods to run (default: all methods)
#' @param launch_app Whether to launch Shiny app after completion (default: FALSE)
#'
#' @return List containing analysis results and file paths
#' @export
#'
#' @examples
#' \dontrun{
#' # Run complete pipeline
#' results <- run_complete_pipeline(
#'   mast_directory = "./iSCORE-PD_MAST_analysis/",
#'   mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
#'   output_directory = "./analysis_output/",
#'   launch_app = TRUE
#' )
#' }
run_complete_pipeline <- function(mast_directory = "./iSCORE-PD_MAST_analysis/",
                                 mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
                                 output_directory = "./analysis_output/",
                                 lfc_threshold = 0.25,
                                 padj_threshold = 0.05,
                                 run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
                                 launch_app = FALSE) {
  
  cat("=== FULL ANALYSIS PIPELINE ===\n")
  cat("Starting at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_directory)) {
    dir.create(output_directory, recursive = TRUE)
  }
  
  # 1. Import DE results
  cat("Step 1: Importing DE results...\n")

# 2. Import DE results
cat("\nStep 2: Importing DE results...\n")
mast_data <- import_mast_data("./iSCORE-PD_MAST_analysis/")
cat("  ✓ MAST:", length(names(mast_data)[names(mast_data) != "metadata"]), "mutations\n")

crispi_data <- import_mixscale_data("./PerturbSeq_MixScale_analysis_full_dataset/", modality = "CRISPRi")
cat("  ✓ CRISPRi:", length(names(crispi_data)), "perturbations\n")

crispa_data <- import_mixscale_data("./PerturbSeq_MixScale_analysis_full_dataset/", modality = "CRISPRa")
cat("  ✓ CRISPRa:", length(names(crispa_data)), "perturbations\n")

# 3. Combine results
cat("\nStep 3: Creating combined dataset...\n")
full_results <- list(
  iSCORE_PD_MAST = mast_data,
  CRISPRi_Mixscale = crispi_data,
  CRISPRa_Mixscale = crispa_data
)
saveRDS(full_results, "full_DE_results.rds")
cat("  ✓ Saved: full_DE_results.rds\n")

# 4. Run enrichment analysis
cat("\nStep 4: Loading enrichment functions and packages...\n")
cat("  This will load many bioinformatics packages - may take a minute\n")
source("enrichment_analysis.R")

cat("\nStep 5: Running functional enrichment analysis...\n")
cat("  This will take several hours for the full dataset!\n")
cat("  Methods: GO, KEGG, Reactome, WikiPathways, STRING, GSEA\n\n")

results <- run_enrichment_analysis(
  input_file = "full_DE_results.rds",
  lfc_threshold = 0.25,
  padj_threshold = 0.05,
  output_dir = "./enrichment_results_full/",
  run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
  min_genes = 10,
  padj_method = "BH"
)

# 6. Summary
cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Ended at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Results saved to:", results$output_directory, "\n")
cat("\nNext steps:\n")
cat("1. Explore results in:", results$output_directory, "\n")
cat("2. Use the Shiny app for visualization\n")
cat("3. Generate reports using report_generation.R\n")
  
  # Launch app if requested
  if (launch_app) {
    launch_iscore_app(consolidated_file)
  }
  
  # Return results summary
  return(list(
    output_directory = output_directory,
    de_results_file = file.path(output_directory, "full_DE_results.rds"),
    enrichment_directory = file.path(output_directory, "enrichment_results"),
    consolidated_file = consolidated_file,
    summary = results
  ))
}