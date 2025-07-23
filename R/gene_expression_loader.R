#' Gene Expression Loader Functions
#' 
#' Functions to load gene expression data on-demand without modifying SCE files

#' Get gene expression for a specific gene
#' 
#' @param gene_name Gene symbol to retrieve
#' @param dataset_name Dataset identifier (e.g., "iSCORE_PD_CRISPRi")
#' @param cell_ids Vector of cell IDs to match (from SCE object)
#' @param normalized Use normalized data (default TRUE)
#' @return Named vector of expression values matching cell order
get_gene_expression <- function(gene_name, dataset_name, cell_ids, normalized = TRUE) {
  
  # Map dataset names to Seurat file paths
  seurat_files <- list(
    "iSCORE_PD" = "../../all_iSCORE-PD_v1.rds",
    "iSCORE_PD_CRISPRi" = "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
    "Full_Dataset" = "../../iSCORE-PD_plus_CRISPRi_and_CRISPRa/iSCORE-PD_plus_CRISPRi_and_CRISPRa.rds"
  )
  
  seurat_path <- seurat_files[[dataset_name]]
  if (is.null(seurat_path) || !file.exists(seurat_path)) {
    warning(sprintf("Seurat file not found for dataset: %s", dataset_name))
    return(rep(NA, length(cell_ids)))
  }
  
  # Try to get from cache first
  cache_key <- paste0(dataset_name, "_", gene_name, "_", ifelse(normalized, "norm", "raw"))
  cached_expr <- .expression_cache[[cache_key]]
  
  if (!is.null(cached_expr)) {
    # Return cached expression matched to cell order
    return(cached_expr[cell_ids])
  }
  
  # Load expression data
  tryCatch({
    # For large files, we'll use a more efficient approach
    # Load only the specific gene's expression
    expr_values <- .load_gene_from_seurat(seurat_path, gene_name, normalized)
    
    if (!is.null(expr_values)) {
      # Cache the result
      .expression_cache[[cache_key]] <- expr_values
      
      # Return matched to cell order
      return(expr_values[cell_ids])
    } else {
      return(rep(NA, length(cell_ids)))
    }
    
  }, error = function(e) {
    warning(sprintf("Error loading gene expression: %s", e$message))
    return(rep(NA, length(cell_ids)))
  })
}

#' Internal function to load specific gene from Seurat object
#' @noRd
.load_gene_from_seurat <- function(seurat_path, gene_name, normalized = TRUE) {
  # This is a placeholder - in practice, you'd want to use HDF5 or other
  # efficient methods to load just the needed data
  
  # For now, return NULL to indicate not implemented
  # In the actual implementation, this would:
  # 1. Open the Seurat file
  # 2. Extract just the specified gene's expression
  # 3. Return as a named vector
  
  return(NULL)
}

# Expression cache environment
.expression_cache <- new.env(parent = emptyenv())

#' Clear expression cache
#' @export
clear_expression_cache <- function() {
  rm(list = ls(.expression_cache), envir = .expression_cache)
}

#' Preload expression for common genes
#' 
#' @param dataset_name Dataset identifier
#' @param genes Vector of gene names to preload
#' @export
preload_gene_expression <- function(dataset_name, genes) {
  # Preload expression for efficiency
  # This would be called when the app starts or dataset changes
  
  for (gene in genes) {
    get_gene_expression(gene, dataset_name, character(0))
  }
}

#' Get list of available genes for a dataset
#' 
#' @param dataset_name Dataset identifier
#' @return Character vector of gene names
#' @export
get_available_genes <- function(dataset_name) {
  # Return a pre-computed list of available genes
  # This should be stored with the dataset metadata
  
  # For now, return common genes
  return(c(
    # PD genes
    "LRRK2", "SNCA", "PINK1", "PRKN", "PARK2", "PARK7",
    "GBA", "ATP13A2", "VPS35", "FBXO7", "DNAJC6", "SYNJ1", "VPS13C",
    
    # Neuronal markers
    "TH", "SLC18A2", "DDC", "SLC6A3", "MAP2", "RBFOX3", "SYN1", 
    "TUBB3", "NCAM1", "GAP43", "STMN2",
    
    # Glial markers
    "GFAP", "AQP4", "S100B", "ALDH1L1", "CD68", "AIF1", "CX3CR1"
  ))
}