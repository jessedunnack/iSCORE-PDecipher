# Enhanced Import Functions for iSCORE-PDecipher Multi-Dimensional Data
# Handles all data variations: Resolution, FDR Status, Approach, Direction

#' Import enrichment results from FDR-corrected directory structure
#' 
#' @param base_dir Base directory containing enrichment_FDR_corrected folder
#' @param resolution "COARSE" or "FINE"
#' @param approach "Pooled", "Split", or "Combined" (Split/Combined only for COARSE)
#' @param gene Gene name (e.g., "LRRK2", "SNCA")
#' @param cluster Cluster number (0-15 for COARSE, 0-39 for FINE)
#' @param direction "ALL", "UP", or "DOWN"
#' @param database Database name: "GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING"
#' @param experiment For Split approach: "C12_FPD-23", "C12_FPD-24", or "C18_FPD-23"
#' @return Enrichment results or NULL if not found
#' @export
import_enrichment_fdr <- function(base_dir, 
                                 resolution = "COARSE",
                                 approach = "Pooled", 
                                 gene, 
                                 cluster, 
                                 direction = "ALL",
                                 database = "GO_BP",
                                 experiment = NULL) {
  
  # Validate inputs
  if (!resolution %in% c("COARSE", "FINE")) {
    stop("resolution must be 'COARSE' or 'FINE'")
  }
  
  if (!approach %in% c("Pooled", "Split", "Combined")) {
    stop("approach must be 'Pooled', 'Split', or 'Combined'")
  }
  
  if (!direction %in% c("ALL", "UP", "DOWN")) {
    stop("direction must be 'ALL', 'UP', or 'DOWN'")
  }
  
  # Check approach availability for resolution
  if (resolution == "FINE" && approach != "Pooled") {
    warning("FINE resolution only supports Pooled approach currently")
    return(NULL)
  }
  
  # Build path based on approach
  if (approach == "Split") {
    if (is.null(experiment)) {
      stop("experiment must be specified for Split approach")
    }
    enrichment_path <- file.path(
      base_dir, "enrichment_analysis", "enrichment_FDR_corrected",
      resolution, approach, experiment, gene,
      paste0("cluster_", cluster), direction,
      paste0(database, "_ALL.rds")
    )
  } else {
    enrichment_path <- file.path(
      base_dir, "enrichment_analysis", "enrichment_FDR_corrected",
      resolution, approach, gene,
      paste0("cluster_", cluster), direction,
      paste0(database, "_ALL.rds")
    )
  }
  
  # Check if file exists
  if (!file.exists(enrichment_path)) {
    message(sprintf("Enrichment file not found: %s", basename(enrichment_path)))
    return(NULL)
  }
  
  # Load and return data
  tryCatch({
    data <- readRDS(enrichment_path)
    attr(data, "metadata") <- list(
      resolution = resolution,
      approach = approach,
      gene = gene,
      cluster = cluster,
      direction = direction,
      database = database,
      experiment = experiment,
      path = enrichment_path
    )
    return(data)
  }, error = function(e) {
    warning(sprintf("Error reading enrichment file: %s", e$message))
    return(NULL)
  })
}

#' Import convergence analysis results
#' 
#' @param base_dir Base directory
#' @param resolution "COARSE" or "FINE"
#' @param type "gene" or "pathway" convergence
#' @param fdr_status "WITH" or "WITHOUT" FDR correction
#' @return Convergence results data frame
#' @export
import_convergence_results <- function(base_dir,
                                      resolution = "COARSE",
                                      type = "gene",
                                      fdr_status = "WITHOUT") {
  
  # Build filename based on parameters
  if (resolution == "COARSE") {
    if (type == "gene") {
      if (fdr_status == "WITHOUT") {
        filename <- "coarse_cluster_convergence_NO_FDR.rds"
      } else {
        filename <- "coarse_cluster_convergence_results.rds"
      }
    } else { # pathway
      filename <- "coarse_pathway_convergence_results.rds"
    }
    
    file_path <- file.path(base_dir, "enrichment_analysis", filename)
    
  } else { # FINE
    if (type == "gene") {
      filename <- "fine_cluster_convergence_with_fdr.rds"
    } else { # pathway
      filename <- "fine_cluster_pathway_convergence_results.rds"
    }
    
    file_path <- file.path(base_dir, "enrichment_analysis", filename)
  }
  
  # Check if file exists
  if (!file.exists(file_path)) {
    # Try alternate location in figures directory
    alt_path <- file.path(base_dir, "enrichment_analysis", "figures",
                         paste0(resolution, "_clusters_", 
                               ifelse(resolution == "COARSE", "16", "40")),
                         basename(file_path))
    if (file.exists(alt_path)) {
      file_path <- alt_path
    } else {
      warning(sprintf("Convergence file not found: %s", filename))
      return(NULL)
    }
  }
  
  # Load and return data
  tryCatch({
    readRDS(file_path)
  }, error = function(e) {
    warning(sprintf("Error reading convergence file: %s", e$message))
    return(NULL)
  })
}

#' Import gene lists used for enrichment analysis
#' 
#' @param base_dir Base directory
#' @param resolution "COARSE" or "FINE"
#' @param approach "Pooled", "Split", or "Combined"
#' @param gene Gene name
#' @param cluster Cluster number
#' @param direction "ALL", "UP", or "DOWN"
#' @param experiment For Split approach
#' @return Character vector of gene names
#' @export
import_gene_list <- function(base_dir,
                            resolution = "COARSE",
                            approach = "Pooled",
                            gene,
                            cluster,
                            direction = "ALL",
                            experiment = NULL) {
  
  # Build filename
  if (approach == "Split") {
    if (is.null(experiment)) {
      stop("experiment must be specified for Split approach")
    }
    filename <- sprintf("%s_Split_%s_%s_cluster_%d_%s.txt",
                       resolution, experiment, gene, cluster, direction)
  } else {
    filename <- sprintf("%s_%s_%s_cluster_%d_%s.txt",
                       resolution, approach, gene, cluster, direction)
  }
  
  file_path <- file.path(base_dir, "enrichment_analysis", 
                        "enrichment_FDR_corrected", "gene_lists", filename)
  
  # Check if file exists
  if (!file.exists(file_path)) {
    message(sprintf("Gene list not found: %s", filename))
    return(character(0))
  }
  
  # Read gene list
  tryCatch({
    genes <- readLines(file_path)
    # Remove any empty lines
    genes <- genes[genes != ""]
    return(genes)
  }, error = function(e) {
    warning(sprintf("Error reading gene list: %s", e$message))
    return(character(0))
  })
}

#' Import all enrichment results for a specific gene and cluster
#' 
#' @param base_dir Base directory
#' @param resolution "COARSE" or "FINE"
#' @param approach "Pooled", "Split", or "Combined"
#' @param gene Gene name
#' @param cluster Cluster number
#' @param experiment For Split approach
#' @return List of enrichment results by direction and database
#' @export
import_enrichment_complete <- function(base_dir,
                                      resolution = "COARSE",
                                      approach = "Pooled",
                                      gene,
                                      cluster,
                                      experiment = NULL) {
  
  directions <- c("ALL", "UP", "DOWN")
  databases <- c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING")
  
  results <- list()
  
  for (direction in directions) {
    results[[direction]] <- list()
    
    for (database in databases) {
      enrichment <- import_enrichment_fdr(
        base_dir = base_dir,
        resolution = resolution,
        approach = approach,
        gene = gene,
        cluster = cluster,
        direction = direction,
        database = database,
        experiment = experiment
      )
      
      if (!is.null(enrichment)) {
        results[[direction]][[database]] <- enrichment
      }
    }
  }
  
  # Add metadata
  attr(results, "query") <- list(
    resolution = resolution,
    approach = approach,
    gene = gene,
    cluster = cluster,
    experiment = experiment
  )
  
  return(results)
}

#' Get available data inventory
#' 
#' @param base_dir Base directory
#' @return Data frame with available data files
#' @export
get_data_inventory <- function(base_dir) {
  
  enrichment_dir <- file.path(base_dir, "enrichment_analysis", "enrichment_FDR_corrected")
  
  inventory <- data.frame(
    resolution = character(),
    approach = character(),
    gene = character(),
    cluster = character(),
    direction = character(),
    database = character(),
    experiment = character(),
    file_exists = logical(),
    stringsAsFactors = FALSE
  )
  
  # Check COARSE data
  for (approach in c("Pooled", "Split", "Combined")) {
    approach_dir <- file.path(enrichment_dir, "COARSE", approach)
    
    if (!dir.exists(approach_dir)) next
    
    if (approach == "Split") {
      # Handle Split with experiments
      experiments <- list.dirs(approach_dir, recursive = FALSE, full.names = FALSE)
      
      for (exp in experiments) {
        exp_dir <- file.path(approach_dir, exp)
        genes <- list.dirs(exp_dir, recursive = FALSE, full.names = FALSE)
        
        for (gene in genes) {
          gene_dir <- file.path(exp_dir, gene)
          clusters <- list.dirs(gene_dir, recursive = FALSE, full.names = FALSE)
          
          for (cluster_str in clusters) {
            cluster <- as.numeric(gsub("cluster_", "", cluster_str))
            cluster_dir <- file.path(gene_dir, cluster_str)
            directions <- list.dirs(cluster_dir, recursive = FALSE, full.names = FALSE)
            
            for (direction in directions) {
              dir_path <- file.path(cluster_dir, direction)
              files <- list.files(dir_path, pattern = "\\.rds$")
              
              for (file in files) {
                database <- gsub("_ALL\\.rds$", "", file)
                
                inventory <- rbind(inventory, data.frame(
                  resolution = "COARSE",
                  approach = "Split",
                  gene = gene,
                  cluster = cluster,
                  direction = direction,
                  database = database,
                  experiment = exp,
                  file_exists = TRUE,
                  stringsAsFactors = FALSE
                ))
              }
            }
          }
        }
      }
    } else {
      # Handle Pooled and Combined
      genes <- list.dirs(approach_dir, recursive = FALSE, full.names = FALSE)
      
      for (gene in genes) {
        gene_dir <- file.path(approach_dir, gene)
        clusters <- list.dirs(gene_dir, recursive = FALSE, full.names = FALSE)
        
        for (cluster_str in clusters) {
          cluster <- as.numeric(gsub("cluster_", "", cluster_str))
          cluster_dir <- file.path(gene_dir, cluster_str)
          directions <- list.dirs(cluster_dir, recursive = FALSE, full.names = FALSE)
          
          for (direction in directions) {
            dir_path <- file.path(cluster_dir, direction)
            files <- list.files(dir_path, pattern = "\\.rds$")
            
            for (file in files) {
              database <- gsub("_ALL\\.rds$", "", file)
              
              inventory <- rbind(inventory, data.frame(
                resolution = "COARSE",
                approach = approach,
                gene = gene,
                cluster = cluster,
                direction = direction,
                database = database,
                experiment = NA,
                file_exists = TRUE,
                stringsAsFactors = FALSE
              ))
            }
          }
        }
      }
    }
  }
  
  # Check FINE data (Pooled only)
  fine_dir <- file.path(enrichment_dir, "FINE", "Pooled")
  
  if (dir.exists(fine_dir)) {
    genes <- list.dirs(fine_dir, recursive = FALSE, full.names = FALSE)
    
    for (gene in genes) {
      gene_dir <- file.path(fine_dir, gene)
      clusters <- list.dirs(gene_dir, recursive = FALSE, full.names = FALSE)
      
      for (cluster_str in clusters) {
        cluster <- as.numeric(gsub("cluster_", "", cluster_str))
        cluster_dir <- file.path(gene_dir, cluster_str)
        directions <- list.dirs(cluster_dir, recursive = FALSE, full.names = FALSE)
        
        for (direction in directions) {
          dir_path <- file.path(cluster_dir, direction)
          files <- list.files(dir_path, pattern = "\\.rds$")
          
          for (file in files) {
            database <- gsub("_ALL\\.rds$", "", file)
            
            inventory <- rbind(inventory, data.frame(
              resolution = "FINE",
              approach = "Pooled",
              gene = gene,
              cluster = cluster,
              direction = direction,
              database = database,
              experiment = NA,
              file_exists = TRUE,
              stringsAsFactors = FALSE
            ))
          }
        }
      }
    }
  }
  
  return(inventory)
}

#' Compare enrichment results across approaches
#' 
#' @param base_dir Base directory
#' @param gene Gene name
#' @param cluster Cluster number
#' @param direction Direction
#' @param database Database
#' @param resolution Resolution level
#' @return List of enrichment results by approach
#' @export
compare_approaches <- function(base_dir,
                             gene,
                             cluster,
                             direction = "ALL",
                             database = "GO_BP",
                             resolution = "COARSE") {
  
  results <- list()
  
  # Pooled approach
  results$Pooled <- import_enrichment_fdr(
    base_dir, resolution, "Pooled", gene, cluster, direction, database
  )
  
  if (resolution == "COARSE") {
    # Split approaches
    for (exp in c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")) {
      key <- paste0("Split_", exp)
      results[[key]] <- import_enrichment_fdr(
        base_dir, resolution, "Split", gene, cluster, direction, database, exp
      )
    }
    
    # Combined approach
    results$Combined <- import_enrichment_fdr(
      base_dir, resolution, "Combined", gene, cluster, direction, database
    )
  }
  
  # Remove NULL results
  results <- results[!sapply(results, is.null)]
  
  # Add comparison metadata
  attr(results, "comparison") <- list(
    gene = gene,
    cluster = cluster,
    direction = direction,
    database = database,
    resolution = resolution,
    n_approaches = length(results)
  )
  
  return(results)
}

#' Create a data accessor function for lazy loading
#' 
#' @param base_dir Base directory
#' @return Function for accessing data on demand
#' @export
create_data_accessor <- function(base_dir) {
  
  # Cache for loaded data
  cache <- new.env(parent = emptyenv())
  
  # Return accessor function
  function(data_type, ..., use_cache = TRUE) {
    
    # Generate cache key
    cache_key <- paste(data_type, paste(..., collapse = "_"), sep = "_")
    
    # Check cache
    if (use_cache && exists(cache_key, envir = cache)) {
      message("Using cached data")
      return(get(cache_key, envir = cache))
    }
    
    # Load data based on type
    result <- switch(data_type,
      "enrichment" = import_enrichment_fdr(base_dir, ...),
      "convergence" = import_convergence_results(base_dir, ...),
      "gene_list" = import_gene_list(base_dir, ...),
      "complete" = import_enrichment_complete(base_dir, ...),
      "inventory" = get_data_inventory(base_dir),
      stop("Unknown data type: ", data_type)
    )
    
    # Cache result
    if (use_cache && !is.null(result)) {
      assign(cache_key, result, envir = cache)
    }
    
    return(result)
  }
}

# Example usage:
# accessor <- create_data_accessor("/mnt/e/ASAP/scRNASeq/PerturbSeq/final")
# enrichment <- accessor("enrichment", "COARSE", "Pooled", "LRRK2", 6, "UP", "GO_BP")
# convergence <- accessor("convergence", "COARSE", "pathway", "WITHOUT")