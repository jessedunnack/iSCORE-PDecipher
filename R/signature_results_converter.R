#' Convert discover_top_signatures results to format expected by analyze_signature_trends
#'
#' This function bridges the gap between the file-based output of discover_top_signatures
#' and the data frame input expected by analyze_signature_trends.
#'
#' @param signature_results Results from discover_top_signatures()
#' @return List with all_signatures and pan_cluster_signatures data frames
#' @export
convert_signature_results_for_trends <- function(signature_results) {
  
  # Check if we have the expected structure
  if (is.null(signature_results) || !is.list(signature_results)) {
    warning("Invalid signature results provided")
    return(list(
      all_signatures = data.frame(),
      pan_cluster_signatures = data.frame()
    ))
  }
  
  # Debug: Log the structure we received
  cat("[CONVERTER] Signature results structure:", paste(names(signature_results), collapse = ", "), "\n")
  
  # Initialize empty data frames
  all_signatures <- data.frame()
  pan_cluster_signatures <- data.frame()
  
  # Try to load data from files if they exist
  if (!is.null(signature_results$files_generated)) {
    cat("[CONVERTER] Files generated structure:", paste(names(signature_results$files_generated), collapse = ", "), "\n")
    cat("[CONVERTER] Top signatures file:", signature_results$files_generated$top_signatures %||% "NULL", "\n")
    
    # Load top signatures file
    if (!is.null(signature_results$files_generated$top_signatures) && 
        file.exists(signature_results$files_generated$top_signatures)) {
      tryCatch({
        all_signatures <- read.csv(signature_results$files_generated$top_signatures, 
                                   stringsAsFactors = FALSE)
        
        # Ensure required columns exist
        if (!"signature_strength" %in% colnames(all_signatures)) {
          # Calculate signature strength if missing
          if (all(c("gene_overlap_count", "pathway_overlap_count", "gene_fisher_p") %in% colnames(all_signatures))) {
            all_signatures$signature_strength <- 
              (all_signatures$gene_overlap_count * 0.4) +
              (all_signatures$pathway_overlap_count * 0.3) +
              (-log10(pmax(all_signatures$gene_fisher_p, 1e-10)) * 0.3)
          } else {
            all_signatures$signature_strength <- 1.0  # Default value
          }
        }
        
        # Ensure cluster column exists
        if (!"cluster" %in% colnames(all_signatures) && "cluster_name" %in% colnames(all_signatures)) {
          all_signatures$cluster <- all_signatures$cluster_name
        }
        
      }, error = function(e) {
        warning("Failed to load top signatures file: ", e$message)
      })
    }
    
    # Load pan-cluster signatures file
    if (!is.null(signature_results$files_generated$pan_cluster) && 
        file.exists(signature_results$files_generated$pan_cluster)) {
      tryCatch({
        pan_cluster_signatures <- read.csv(signature_results$files_generated$pan_cluster, 
                                          stringsAsFactors = FALSE)
        
        # Ensure required columns exist
        if (!"signature_strength" %in% colnames(pan_cluster_signatures)) {
          # Calculate mean signature strength if missing
          if ("mean_gene_overlap" %in% colnames(pan_cluster_signatures)) {
            pan_cluster_signatures$signature_strength <- pan_cluster_signatures$mean_gene_overlap * 0.1
          } else {
            pan_cluster_signatures$signature_strength <- 1.0  # Default value
          }
        }
        
        # Add cluster column for pan-cluster (set to "pan_cluster")
        if (!"cluster" %in% colnames(pan_cluster_signatures)) {
          pan_cluster_signatures$cluster <- "pan_cluster"
        }
        
      }, error = function(e) {
        warning("Failed to load pan-cluster signatures file: ", e$message)
      })
    }
  }
  
  # If no files were loaded, check if data is directly in the results object
  if (nrow(all_signatures) == 0) {
    cat("[CONVERTER] No files loaded, checking for direct data in results object\n")
    
    # Check if signature data is directly available (not in files)
    if (!is.null(signature_results$all_signatures)) {
      cat("[CONVERTER] Found all_signatures directly in results\n")
      all_signatures <- signature_results$all_signatures
    } else if (!is.null(signature_results$top_signatures)) {
      cat("[CONVERTER] Found top_signatures directly in results\n")
      all_signatures <- signature_results$top_signatures
    } else {
      cat("[CONVERTER] No signature data found in any format\n")
      warning("No signature data could be loaded from files or summary stats")
      all_signatures <- data.frame(
        gene_pair = character(0),
        cluster = character(0),
        signature_strength = numeric(0),
        gene_overlap_count = integer(0),
        pathway_overlap_count = integer(0),
        stringsAsFactors = FALSE
      )
    }
  }
  
  # Same check for pan-cluster signatures
  if (nrow(pan_cluster_signatures) == 0 && !is.null(signature_results$pan_cluster_signatures)) {
    cat("[CONVERTER] Found pan_cluster_signatures directly in results\n")
    pan_cluster_signatures <- signature_results$pan_cluster_signatures
  }
  
  return(list(
    all_signatures = all_signatures,
    pan_cluster_signatures = pan_cluster_signatures,
    summary_stats = signature_results$summary_stats,
    files_generated = signature_results$files_generated
  ))
}