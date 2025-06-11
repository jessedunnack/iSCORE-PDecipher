#' @keywords internal
"_PACKAGE"

#' iSCORE-PDecipher: Integrated Analysis of Parkinson's Disease Mutations and Perturbations
#'
#' iSCORE-PDecipher provides comprehensive analysis tools for Parkinson's disease
#' research, combining iSCORE-PD mutation data with PerturbSeq (CRISPRi/CRISPRa) experiments.
#' It includes differential expression analysis using MAST and MixScale, comprehensive
#' functional enrichment analysis, and interactive Shiny visualization tools.
#'
#' @section Main Functions:
#' \itemize{
#'   \item \code{\link{launch_iscore_app}}: Launch the interactive Shiny application
#'   \item \code{\link{run_complete_pipeline}}: Run the complete analysis pipeline
#'   \item \code{\link{import_mast_data}}: Import MAST differential expression results
#'   \item \code{\link{import_mixscale_data}}: Import MixScale differential expression results
#'   \item \code{\link{run_enrichment_analysis}}: Perform functional enrichment analysis
#' }
#'
#' @section Data Types Supported:
#' \itemize{
#'   \item \strong{iSCORE-PD MAST}: Mutation vs. wild-type comparisons (13 mutations, 14 clusters)
#'   \item \strong{CRISPRi MixScale}: Gene knockdown perturbations (10 genes, 10 clusters)  
#'   \item \strong{CRISPRa MixScale}: Gene activation perturbations (10 genes, 1 cluster)
#' }
#'
#' @section Analysis Pipeline:
#' The package provides a complete workflow:
#' \enumerate{
#'   \item Import differential expression results from MAST and MixScale
#'   \item Perform functional enrichment analysis (GO, KEGG, Reactome, WikiPathways, STRING, GSEA)
#'   \item Consolidate results into unified format
#'   \item Launch interactive Shiny application for exploration and visualization
#' }
#'
#' @section Key Features:
#' \itemize{
#'   \item Process 767,337+ enrichment terms across all conditions
#'   \item Interactive heatmap visualizations for cross-condition comparisons
#'   \item Export capabilities for publication-ready figures
#'   \item Modality comparison between mutations and perturbations
#'   \item Smart file detection and validation
#' }
#'
#' @docType package
#' @name iSCORE.PDecipher-package
#' @aliases iSCORE.PDecipher
#'
#' @author Jesse Dunnack \email{jesse.dunnack@example.com}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/yourusername/iSCORE-PDecipher}
#'   \item Report bugs at \url{https://github.com/yourusername/iSCORE-PDecipher/issues}
#' }
#'
#' @examples
#' \dontrun{
#' # Launch the Shiny app with your data
#' launch_iscore_app("path/to/enrichment_results.rds")
#'
#' # Run complete analysis pipeline
#' results <- run_complete_pipeline(
#'   mast_directory = "./iSCORE-PD_MAST_analysis/",
#'   mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
#'   output_directory = "./analysis_output/",
#'   launch_app = TRUE
#' )
#' }
NULL