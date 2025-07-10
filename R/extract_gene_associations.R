#' Extract Gene-Term Associations from Enrichment Files
#'
#' This script processes all enrichment result files and extracts gene associations
#' to create a consolidated dataset for the R package deployment.
#'
#' @author Claude Code Assistant
#' @date 2025-01-13

library(dplyr)
library(tibble)
library(purrr)
library(stringr)

#' Extract gene associations from a single enrichment file
#'
#' @param file_path Path to the enrichment .rds file
#' @param metadata Metadata extracted from file path
#' @return Data frame with gene associations
extract_from_single_file <- function(file_path, metadata) {
  tryCatch({
    # Read the enrichment result
    result <- readRDS(file_path)
    
    # Handle different result object types
    if (inherits(result, "enrichResult") || inherits(result, "gseaResult")) {
      # clusterProfiler S4 objects
      df <- result@result
    } else if (is.data.frame(result)) {
      # Already a data frame
      df <- result
    } else if (is.list(result) && !is.null(result$result)) {
      # List with result element
      df <- result$result
    } else {
      warning("Unknown result format in file: ", file_path)
      return(NULL)
    }
    
    # Check for required columns
    if (!all(c("ID", "Description") %in% colnames(df))) {
      warning("Missing required columns in file: ", file_path)
      return(NULL)
    }
    
    # Check for gene ID column (different possible names)
    gene_col <- NULL
    if ("geneID" %in% colnames(df)) {
      gene_col <- "geneID"
    } else if ("core_enrichment" %in% colnames(df)) {
      gene_col <- "core_enrichment"  # GSEA results
    } else if ("genes" %in% colnames(df)) {
      gene_col <- "genes"
    } else {
      warning("No gene ID column found in file: ", file_path)
      return(NULL)
    }
    
    # Extract relevant data
    associations <- df %>%
      filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
      mutate(
        # Add metadata
        analysis_type = metadata$analysis_type,
        gene = metadata$gene,
        cluster = metadata$cluster,
        enrichment_type = metadata$enrichment_type,
        direction = metadata$direction,
        experiment = metadata$experiment,
        
        # Process gene list
        associated_genes = .data[[gene_col]],
        gene_count = str_count(associated_genes, "/") + 1,
        
        # Create composite key for fast lookup
        composite_key = paste(analysis_type, gene, cluster, enrichment_type, 
                             direction, experiment, ID, sep = "|"),
        
        # Keep only essential columns
        term_id = ID,
        description = Description
      ) %>%
      select(
        term_id, description, analysis_type, gene, cluster,
        enrichment_type, direction, experiment, associated_genes,
        gene_count, composite_key
      )
    
    return(associations)
    
  }, error = function(e) {
    warning("Error processing file ", file_path, ": ", e$message)
    return(NULL)
  })
}

#' Parse metadata from file path
#'
#' @param file_path Path to enrichment file
#' @return List with metadata components
parse_file_metadata <- function(file_path) {
  # Split path components
  path_parts <- str_split(file_path, "/")[[1]]
  
  # Find enrichment_results index
  er_idx <- which(path_parts == "enrichment_results")
  
  if (length(er_idx) == 0) {
    stop("Cannot find 'enrichment_results' in path: ", file_path)
  }
  
  # Extract components after enrichment_results
  # Pattern: enrichment_results/[METHOD]/[GENE]/[CLUSTER]/default/[ENRICHMENT_TYPE]/[FILE]
  analysis_type <- path_parts[er_idx + 1]  # MAST or MixScale
  gene <- path_parts[er_idx + 2]           # Gene name
  cluster <- path_parts[er_idx + 3]        # cluster_X
  # skip 'default'
  enrichment_type <- path_parts[er_idx + 5] # GO_BP, KEGG, etc.
  
  # Extract direction and experiment from filename
  filename <- path_parts[length(path_parts)]
  filename_base <- str_remove(filename, "\\.rds$")
  
  # Handle different filename patterns
  if (str_detect(filename_base, "GSEA_")) {
    # GSEA files: GSEA_GO_BP.rds
    direction <- "RANKED"
    experiment <- "default"
  } else {
    # Standard files: GO_BP_UP.rds, KEGG_ALL.rds
    parts <- str_split(filename_base, "_")[[1]]
    if (length(parts) >= 2) {
      direction <- parts[length(parts)]  # Last part is direction
    } else {
      direction <- "ALL"
    }
    experiment <- "default"
  }
  
  # Handle MixScale experiment detection (if needed later)
  if (analysis_type == "MixScale") {
    # For MixScale, experiment info might be embedded differently
    # For now, use default
    experiment <- "default"
  }
  
  return(list(
    analysis_type = analysis_type,
    gene = gene,
    cluster = cluster,
    enrichment_type = enrichment_type,
    direction = direction,
    experiment = experiment
  ))
}

#' Extract gene associations from all enrichment files
#'
#' @param dataset_dirs Vector of dataset directory paths
#' @param batch_size Number of files to process per batch
#' @param output_file Path for output file
#' @return Data frame with all gene associations
extract_all_gene_associations <- function(dataset_dirs, batch_size = 100, output_file = NULL) {
  
  cat("=== Starting Gene Association Extraction ===\n")
  cat("Dataset directories:", paste(dataset_dirs, collapse = ", "), "\n")
  
  # Find all enrichment .rds files using system find command
  all_files <- c()
  for (dir in dataset_dirs) {
    if (!dir.exists(file.path(dir, "enrichment_results"))) {
      cat("Warning: enrichment_results directory not found in", dir, "\n")
      next
    }
    
    # Use system find command for recursive search
    cmd <- paste0("find ", file.path(dir, "enrichment_results"), " -name '*.rds'")
    files <- system(cmd, intern = TRUE)
    
    # Filter to only enrichment result files (not metadata)
    enrich_files <- files[str_detect(files, "(GO_|KEGG|Reactome|STRING|GSEA)")]
    all_files <- c(all_files, enrich_files)
    cat("Found", length(enrich_files), "enrichment files in", dir, "\n")
  }
  
  cat("Total files to process:", length(all_files), "\n")
  
  if (length(all_files) == 0) {
    stop("No enrichment files found!")
  }
  
  # Process files in batches
  all_associations <- list()
  n_batches <- ceiling(length(all_files) / batch_size)
  
  for (batch in 1:n_batches) {
    cat("Processing batch", batch, "of", n_batches, "...\n")
    
    # Get batch file indices
    start_idx <- (batch - 1) * batch_size + 1
    end_idx <- min(batch * batch_size, length(all_files))
    batch_files <- all_files[start_idx:end_idx]
    
    # Process batch
    batch_results <- map_dfr(batch_files, function(file) {
      metadata <- parse_file_metadata(file)
      extract_from_single_file(file, metadata)
    })
    
    if (!is.null(batch_results) && nrow(batch_results) > 0) {
      all_associations[[batch]] <- batch_results
      cat("  Batch", batch, "processed:", nrow(batch_results), "associations\n")
    } else {
      cat("  Batch", batch, "produced no results\n")
    }
    
    # Garbage collection after each batch
    gc()
  }
  
  # Combine all results
  cat("Combining all results...\n")
  final_associations <- bind_rows(all_associations)
  
  cat("=== Extraction Summary ===\n")
  cat("Total associations extracted:", nrow(final_associations), "\n")
  cat("Unique terms:", n_distinct(final_associations$term_id), "\n")
  cat("Unique genes:", n_distinct(final_associations$gene), "\n")
  cat("Analysis types:", paste(unique(final_associations$analysis_type), collapse = ", "), "\n")
  cat("Enrichment types:", paste(unique(final_associations$enrichment_type), collapse = ", "), "\n")
  
  # Save results if output file specified
  if (!is.null(output_file)) {
    cat("Saving to:", output_file, "\n")
    saveRDS(final_associations, output_file, compress = "xz")
    
    # Check file size
    file_size <- file.size(output_file) / (1024^2)  # MB
    cat("Output file size:", round(file_size, 1), "MB\n")
  }
  
  return(final_associations)
}

#' Main execution function
#'
#' @export
main <- function() {
  # Define dataset directories
  dataset_dirs <- c(
    "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD",
    "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi"
  )
  
  # Output file path
  output_file <- "inst/extdata/gene_term_associations.rds"
  
  # Ensure output directory exists
  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
  
  # Extract associations
  associations <- extract_all_gene_associations(
    dataset_dirs = dataset_dirs,
    batch_size = 100,
    output_file = output_file
  )
  
  cat("=== Gene Association Extraction Complete ===\n")
  
  return(associations)
}

# Allow script to be run directly
# Commented out for testing
# if (!interactive()) {
#   main()
# }