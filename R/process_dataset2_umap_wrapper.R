#' Process UMAP and marker data for dataset 2 (iSCORE_PD_plus_CRISPRi)
#'
#' @description
#' Wrapper function to process UMAP data and calculate markers for the 
#' iSCORE_PD_plus_CRISPRi dataset. This function provides an R interface
#' to the process_dataset2_umap.R script.
#'
#' @param pc_count Number of PCs for UMAP (default: 30)
#' @param calculate_markers Logical, whether to calculate cluster markers (default: TRUE)
#' @param calculate_group_markers Logical, whether to calculate group comparison markers (default: FALSE)
#' @param max_cells Maximum cells per cluster for marker calculation. Use -1 for all cells (default: 500)
#' @param force_overwrite Logical, force overwrite existing files without prompting (default: FALSE)
#' @param resolution Clustering resolution if reclustering needed (default: NULL)
#' @param top_markers Number of top markers to display per cluster (default: 10)
#' @param save_all_pc Logical, save UMAP for all PC counts (30, 50, 100) (default: FALSE)
#' @param verbose Logical, show detailed progress messages (default: TRUE)
#'
#' @return Invisible NULL. The function saves files to the inst/extdata/umap_data directory.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Process with default settings
#' process_dataset2_umap()
#' 
#' # Process with all cells (no sampling) and group markers
#' process_dataset2_umap(
#'   max_cells = -1,
#'   calculate_group_markers = TRUE,
#'   force_overwrite = TRUE
#' )
#' 
#' # Generate all PC versions with maximum detail
#' process_dataset2_umap(
#'   save_all_pc = TRUE,
#'   max_cells = -1,
#'   calculate_group_markers = TRUE,
#'   top_markers = 20
#' )
#' }
process_dataset2_umap <- function(
  pc_count = 30,
  calculate_markers = TRUE,
  calculate_group_markers = FALSE,
  max_cells = 500,
  force_overwrite = FALSE,
  resolution = NULL,
  top_markers = 10,
  save_all_pc = FALSE,
  verbose = TRUE
) {
  
  # Find the script
  script_path <- system.file("scripts/process_dataset2_umap.R", package = "iSCORE.PDecipher")
  
  if (!file.exists(script_path)) {
    stop("process_dataset2_umap.R script not found in package installation.")
  }
  
  # Build command
  cmd <- paste("Rscript", shQuote(script_path))
  
  # Add arguments
  cmd <- paste(cmd, "--pc", pc_count)
  cmd <- paste(cmd, "--max-cells", max_cells)
  cmd <- paste(cmd, "--top-markers", top_markers)
  
  if (!calculate_markers) {
    cmd <- paste(cmd, "--markers FALSE")
  }
  
  if (calculate_group_markers) {
    cmd <- paste(cmd, "--group-markers")
  }
  
  if (force_overwrite) {
    cmd <- paste(cmd, "--force")
  }
  
  if (!is.null(resolution)) {
    cmd <- paste(cmd, "--resolution", resolution)
  }
  
  if (save_all_pc) {
    cmd <- paste(cmd, "--save-all-pc")
  }
  
  # Show command if verbose
  if (verbose) {
    cat("Running command:\n")
    cat(cmd, "\n\n")
  }
  
  # Run the script
  result <- system(cmd, intern = !verbose)
  
  if (verbose) {
    cat("\nProcessing complete. Check the output directory for results.\n")
    cat("Default location: inst/extdata/umap_data/\n")
  }
  
  invisible(NULL)
}

#' Process UMAP with all features for high-memory systems
#'
#' @description
#' Convenience function for processing with all features enabled,
#' optimized for systems with high memory (e.g., 96GB RAM).
#' This runs with no cell sampling and calculates all marker types.
#'
#' @param pc_count Number of PCs for UMAP (default: 30)
#' @param save_all_pc Logical, save all PC versions (default: FALSE)
#'
#' @return Invisible NULL
#' @export
#'
#' @examples
#' \dontrun{
#' # Run full analysis
#' process_dataset2_umap_full()
#' 
#' # Run full analysis with all PC versions
#' process_dataset2_umap_full(save_all_pc = TRUE)
#' }
process_dataset2_umap_full <- function(pc_count = 30, save_all_pc = FALSE) {
  cat("Running full UMAP processing with all features...\n")
  cat("This will use all cells (no sampling) and calculate all marker types.\n")
  cat("Expected runtime: 10-20 minutes on a high-memory system.\n\n")
  
  process_dataset2_umap(
    pc_count = pc_count,
    calculate_markers = TRUE,
    calculate_group_markers = TRUE,
    max_cells = -1,  # Use all cells
    force_overwrite = TRUE,
    top_markers = 20,
    save_all_pc = save_all_pc,
    verbose = TRUE
  )
}