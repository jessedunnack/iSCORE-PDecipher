# Comprehensive functional enrichment analysis for MAST and MixScale results
# This script provides a standardized workflow for performing enrichment analysis
# across both differential expression analysis methods.

# Required packages (checked during function execution)

# Custom function for consistent message formatting
print_header <- function(title, level = 1) {
  if (level == 1) {
    message("\n", paste(rep("=", 80), collapse = ""))
    message("  ", title)
    message(paste(rep("=", 80), collapse = ""), "\n")
  } else if (level == 2) {
    message("\n", paste(rep("-", 60), collapse = ""))
    message("  ", title)
    message(paste(rep("-", 60), collapse = ""), "\n")
  } else {
    prefix <- paste(rep("  ", level-2), collapse = "")
    message(prefix, "• ", title)
  }
}

#' Main function to run comprehensive enrichment analysis
#'
#' @param input_file Path to the processed RDS file containing MAST and MixScale results
#' @param lfc_threshold Log2 fold change threshold for defining differentially expressed genes
#' @param padj_threshold Adjusted p-value threshold for defining differentially expressed genes
#' @param output_dir Base directory for saving results
#' @param run_methods Vector of analysis methods to run
#' @param min_genes Minimum number of genes required for enrichment analysis
#' @param padj_method Method for adjusting p-values in enrichment analysis
#' @return Invisible summary of the analysis
run_enrichment_analysis <- function(
  input_file = "full_DE_results.rds",
  lfc_threshold = 0.5,
  padj_threshold = 0.05,
  output_dir = "./enrichment_results/",
  run_methods = c("GO", "GSEA", "KEGG", "Reactome", "WikiPathways", "STRING"),
  min_genes = 5,
  padj_method = "BH"
) {
  # 1. Load data
  print_header("STARTING ENRICHMENT ANALYSIS")
  message("Loading data from ", input_file)
  data <- readRDS(input_file)
  
  # 2. Create main output directory
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 3. Process MAST data
  if(!is.null(data$iSCORE_PD_MAST)) {
    print_header("PROCESSING MAST DATA")
    mast_data <- data$iSCORE_PD_MAST
    
    # Iterate through mutations
    for(mutation in names(mast_data)) {
      if(mutation == "metadata") next
      
      print_header(paste("Processing mutation:", mutation), 2)
      mutation_data <- mast_data[[mutation]]
      
      # Iterate through clusters
      for(cluster in names(mutation_data)) {
        if(cluster == "metadata") next
        
        print_header(paste("Processing cluster:", cluster), 3)
        cluster_data <- mutation_data[[cluster]]
        
        # Process this entry
        process_mast_entry(
          entry = cluster_data,
          mutation = mutation,
          cluster = cluster,
          output_base_dir = output_dir,
          lfc_threshold = lfc_threshold,
          padj_threshold = padj_threshold,
          min_genes = min_genes,
          run_methods = run_methods,
          padj_method = padj_method
        )
      }
    }
  }
  
  # 4. Process MixScale CRISPRi data
  if(!is.null(data$CRISPRi_Mixscale)) {
    print_header("PROCESSING CRISPR-I MIXSCALE DATA")
    mixscale_data <- data$CRISPRi_Mixscale
    
    # Iterate through perturbations
    for(perturbation in names(mixscale_data)) {
      if(perturbation == "metadata") next
      
      print_header(paste("Processing CRISPRi perturbation:", perturbation), 2)
      perturbation_data <- mixscale_data[[perturbation]]
      
      # Iterate through clusters
      for(cluster in names(perturbation_data)) {
        if(cluster == "metadata") next
        
        print_header(paste("Processing cluster:", cluster), 3)
        cluster_data <- perturbation_data[[cluster]]
        
        # Process this entry
        process_mixscale_entry(
          entry = cluster_data,
          perturbation = perturbation,
          cluster = cluster,
          output_base_dir = output_dir,
          lfc_threshold = lfc_threshold,
          padj_threshold = padj_threshold,
          min_genes = min_genes,
          run_methods = run_methods,
          padj_method = padj_method
        )
      }
    }
  }
  
  # 5. Process MixScale CRISPRa data
  if(!is.null(data$CRISPRa_Mixscale)) {
    print_header("PROCESSING CRISPR-A MIXSCALE DATA")
    mixscale_data <- data$CRISPRa_Mixscale
    
    # Iterate through perturbations
    for(perturbation in names(mixscale_data)) {
      if(perturbation == "metadata") next
      
      print_header(paste("Processing CRISPRa perturbation:", perturbation), 2)
      perturbation_data <- mixscale_data[[perturbation]]
      
      # Iterate through clusters
      for(cluster in names(perturbation_data)) {
        if(cluster == "metadata") next
        
        print_header(paste("Processing cluster:", cluster), 3)
        cluster_data <- perturbation_data[[cluster]]
        
        # Process this entry
        process_mixscale_entry(
          entry = cluster_data,
          perturbation = perturbation,
          cluster = cluster,
          output_base_dir = output_dir,
          lfc_threshold = lfc_threshold,
          padj_threshold = padj_threshold,
          min_genes = min_genes,
          run_methods = run_methods,
          padj_method = padj_method
        )
      }
    }
  }
  
  # 6. Generate summary report
  print_header("GENERATING SUMMARY REPORT")
  
  # Create metadata directory for analysis summaries
  metadata_dir <- file.path(output_dir, "metadata")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  
  summary <- list(
    input_file = input_file,
    lfc_threshold = lfc_threshold,
    padj_threshold = padj_threshold,
    methods_used = run_methods,
    min_genes = min_genes,
    padj_method = padj_method,
    completion_time = Sys.time(),
    output_directory = output_dir
  )
  
  saveRDS(summary, file.path(metadata_dir, "analysis_summary.rds"))
  write.csv(data.frame(
    Parameter = names(summary),
    Value = unlist(lapply(summary, function(x) paste(x, collapse=", ")))
  ), file.path(metadata_dir, "analysis_summary.csv"))
  
  # Also save a copy at the root level for backward compatibility
  saveRDS(summary, file.path(output_dir, "analysis_summary.rds"))
  write.csv(data.frame(
    Parameter = names(summary),
    Value = unlist(lapply(summary, function(x) paste(x, collapse=", ")))
  ), file.path(output_dir, "analysis_summary.csv"))
  
  print_header("ENRICHMENT ANALYSIS COMPLETED SUCCESSFULLY!")
  return(invisible(summary))
}

#' Process a MAST entry for enrichment analysis
#'
#' @param entry The entry data containing results, metadata, and background genes
#' @param mutation The mutation identifier
#' @param cluster The cluster identifier
#' @param output_base_dir Base directory for output
#' @param lfc_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @param min_genes Minimum number of genes for enrichment
#' @param run_methods Vector of analysis methods to run
#' @param padj_method Method for p-value adjustment
process_mast_entry <- function(entry, mutation, cluster, output_base_dir, 
                             lfc_threshold, padj_threshold, min_genes, run_methods,
                             padj_method = "BH") {
  # Create output directory for this entry
  output_dir <- file.path(output_base_dir, "MAST", mutation, cluster)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  exp_dir <- file.path(output_dir, "default")
  dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract data and metadata
  results_df <- entry$results
  background_genes <- entry$background_genes
  metadata <- entry$metadata
  
  # Identify columns for filtering
  lfc_col <- "avg_log2FC"
  padj_col <- "p_val_adj"
  
  # Filter genes based on thresholds
  gene_lists <- filter_genes(results_df, lfc_col, padj_col, lfc_threshold, padj_threshold)
  
  # Create gene_lists directory for all gene files
  gene_lists_dir <- file.path(exp_dir, "gene_lists")
  dir.create(gene_lists_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save gene lists in the gene_lists directory
  write.table(gene_lists$up_genes, file.path(gene_lists_dir, "up_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(gene_lists$down_genes, file.path(gene_lists_dir, "down_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(gene_lists$all_de_genes, file.path(gene_lists_dir, "all_de_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(background_genes, file.path(gene_lists_dir, "background_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
             
  # Also save at the root level for backward compatibility
  write.table(gene_lists$up_genes, file.path(exp_dir, "up_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(gene_lists$down_genes, file.path(exp_dir, "down_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(gene_lists$all_de_genes, file.path(exp_dir, "all_de_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(background_genes, file.path(exp_dir, "background_genes.csv"), 
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # For GSEA, create ranked gene list
  if("GSEA" %in% run_methods) {
    # Create ranked list based on log2FC * -log10(padj)
    ranked_genes <- results_df
    # Handle non-finite values appropriately
    clean_padj <- ifelse(!is.finite(ranked_genes[[padj_col]]) | ranked_genes[[padj_col]] == 0, 
                         1e-300, ranked_genes[[padj_col]])
    clean_lfc <- ifelse(!is.finite(ranked_genes[[lfc_col]]), 0, ranked_genes[[lfc_col]])
    
    ranked_genes$ranking_score <- clean_lfc * -log10(clean_padj)
    # Replace any remaining non-finite values with a reasonable limit
    ranked_genes$ranking_score <- ifelse(!is.finite(ranked_genes$ranking_score), 0, ranked_genes$ranking_score)
    
    ranked_genes <- ranked_genes[order(ranked_genes$ranking_score, decreasing = TRUE), ]
    ranked_list <- ranked_genes$ranking_score
    names(ranked_list) <- rownames(ranked_genes)
    
    # Create ranked_lists directory
    ranked_lists_dir <- file.path(exp_dir, "ranked_lists")
    dir.create(ranked_lists_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save ranked list in ranked_lists directory
    write.table(data.frame(gene = names(ranked_list), score = ranked_list),
               file.path(ranked_lists_dir, "ranked_gene_list.csv"),
               row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
               
    # Also save at the root level for backward compatibility
    write.table(data.frame(gene = names(ranked_list), score = ranked_list),
               file.path(exp_dir, "ranked_gene_list.csv"),
               row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # Run enrichment analyses
  enrichment_results <- run_all_enrichment_analyses(
    gene_lists = gene_lists,
    background_genes = background_genes,
    ranked_list = if("GSEA" %in% run_methods) ranked_list else NULL,
    min_genes = min_genes,
    run_methods = run_methods,
    output_dir = exp_dir,
    padj_method = padj_method
  )
  
  # Create metadata
  analysis_metadata <- list(
    analysis_type = "MAST",
    mutation = mutation,
    cluster = cluster,
    lfc_threshold = lfc_threshold,
    padj_threshold = padj_threshold,
    lfc_column = lfc_col,
    padj_column = padj_col,
    date = Sys.Date(),
    original_metadata = metadata,
    up_genes_count = length(gene_lists$up_genes),
    down_genes_count = length(gene_lists$down_genes),
    all_de_genes_count = length(gene_lists$all_de_genes),
    background_genes_count = length(background_genes)
  )
  
  # Create metadata directory and save metadata
  metadata_dir <- file.path(exp_dir, "metadata")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(analysis_metadata, file.path(metadata_dir, "analysis_metadata.rds"))
  
  # Also save at the root level for backward compatibility
  saveRDS(analysis_metadata, file.path(exp_dir, "analysis_metadata.rds"))
  
  # Create results directory for all results tables
  results_dir <- file.path(exp_dir, "results_tables")
  dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Save full DE table with flags
  full_de_table <- results_df
  full_de_table$passes_threshold <- rownames(full_de_table) %in% gene_lists$all_de_genes
  full_de_table$direction <- ifelse(rownames(full_de_table) %in% gene_lists$up_genes, "UP",
                             ifelse(rownames(full_de_table) %in% gene_lists$down_genes, "DOWN", "NONE"))
  
  # Save in both locations
  write.csv(full_de_table, file.path(results_dir, "full_de_table.csv"))
  write.csv(full_de_table, file.path(exp_dir, "full_de_table.csv"))
}

#' Process a MixScale entry for enrichment analysis
#'
#' @param entry The entry data containing results, metadata, and background genes
#' @param perturbation The perturbation gene
#' @param cluster The cluster identifier
#' @param output_base_dir Base directory for output
#' @param lfc_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @param min_genes Minimum number of genes for enrichment
#' @param run_methods Vector of analysis methods to run
#' @param padj_method Method for p-value adjustment
process_mixscale_entry <- function(entry, perturbation, cluster, output_base_dir, 
                                  lfc_threshold, padj_threshold, min_genes, run_methods,
                                  padj_method = "BH") {
  # Create output directory for this entry
  output_dir <- file.path(output_base_dir, "MixScale", perturbation, cluster)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create a diagnostics directory to save mapping information
  diagnostics_dir <- file.path(output_dir, "diagnostics")
  dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Extract data and metadata
  results_df <- entry$results
  global_background_genes <- entry$background_genes
  metadata <- entry$metadata
  
  # Get experiment information
  experiments <- metadata$experiments
  if(length(experiments) == 0) experiments <- c("default")
  
  # Print all available columns for debugging
  message("  Available columns in results dataframe:")
  message("  ", paste(colnames(results_df), collapse=", "))
  
  # Initialize experiment diagnostics tracking
  experiment_diagnostics <- data.frame(
    experiment_id = character(),
    log2fc_column = character(),
    pvalue_column = character(),
    status = character(),
    reason = character(),
    na_count_lfc = integer(),
    na_percent_lfc = numeric(),
    na_count_pval = integer(),
    na_percent_pval = numeric(),
    background_genes_count = integer(),
    stringsAsFactors = FALSE
  )
  
  # First get p-value columns - they will help us identify the true experiment IDs
  p_value_cols <- grep("^p_cell_type", colnames(results_df), value = TRUE)
  p_value_weight_cols <- grep(":weight$", p_value_cols, value = TRUE)
  
  message("  Detected p-value columns:")
  message("  ", paste(p_value_cols, collapse=", "))
  message("  Weighted p-value columns:")
  message("  ", paste(p_value_weight_cols, collapse=", "))
  
  # Extract experiment IDs from p-value column names
  experiment_ids_from_p <- unique(gsub("^p_cell_type([^:]+).*$", "\\1", p_value_weight_cols))
  
  message("  Experiment IDs extracted from p-value columns:")
  message("  ", paste(experiment_ids_from_p, collapse=", "))
  
  # Get log fold change columns
  lfc_cols <- grep("^log2FC_", colnames(results_df), value = TRUE)
  if(length(lfc_cols) == 0) lfc_cols <- "log2FC"
  
  message("  Detected log2FC columns:")
  message("  ", paste(lfc_cols, collapse=", "))
  
  # Extract experiment IDs from log2FC column names
  experiment_ids_from_lfc <- unique(gsub("^log2FC_(.+)$", "\\1", lfc_cols))
  
  message("  Experiment IDs extracted from log2FC columns:")
  message("  ", paste(experiment_ids_from_lfc, collapse=", "))
  
  # Create a diagnostic function to check NA counts in each log2FC column
  lfc_na_counts <- sapply(lfc_cols, function(col) {
    sum(is.na(results_df[[col]]))
  })
  
  message("  NA counts in log2FC columns:")
  for (col in names(lfc_na_counts)) {
    message(sprintf("  - %s: %d NAs (%.1f%%)", col, lfc_na_counts[col], 
                    100 * lfc_na_counts[col] / nrow(results_df)))
  }
  
  # Get p-value NA counts
  pval_na_counts <- sapply(p_value_weight_cols, function(col) {
    sum(is.na(results_df[[col]]))
  })
  
  message("  NA counts in p-value columns:")
  for (col in names(pval_na_counts)) {
    message(sprintf("  - %s: %d NAs (%.1f%%)", col, pval_na_counts[col],
                   100 * pval_na_counts[col] / nrow(results_df)))
  }
  
  # List of valid experiments with matching log2FC and p-value columns
  valid_experiments <- character(0)
  valid_lfc_cols <- character(0) 
  valid_pval_cols <- character(0)
  experiment_background_genes <- list()
  
  # Process each experiment to find matching columns - ONLY include experiments with BOTH matching log2FC and p-value columns
  for (exp_id in experiment_ids_from_p) {
    message(sprintf("  Processing experiment: %s", exp_id))
    
    # Look for weighted p-value column for this experiment
    pval_col <- grep(paste0("p_cell_type", exp_id, ".*:weight"), colnames(results_df), value = TRUE)
    
    if (length(pval_col) == 0) {
      message(sprintf("  - SKIPPING: No weighted p-value column found for experiment %s", exp_id))
      experiment_diagnostics <- rbind(experiment_diagnostics, data.frame(
        experiment_id = exp_id,
        log2fc_column = NA,
        pvalue_column = NA,
        status = "EXCLUDED",
        reason = "No weighted p-value column found",
        na_count_lfc = NA,
        na_percent_lfc = NA,
        na_count_pval = NA,
        na_percent_pval = NA,
        background_genes_count = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # If multiple weighted p-value columns, use the first one
    if (length(pval_col) > 1) {
      message(sprintf("  - MULTIPLE P-VALUE COLUMNS: Found %d p-value columns, using %s", 
                    length(pval_col), pval_col[1]))
      pval_col <- pval_col[1]
    }
    
    # Try exact match for log2FC column - ONLY use exact matches to prevent mixing experiments
    exact_match <- paste0("log2FC_", exp_id)
    
    if (exact_match %in% lfc_cols) {
      message(sprintf("  - EXACT MATCH: Found exact matching log2FC column %s", exact_match))
      lfc_col <- exact_match
    } else {
      message(sprintf("  - NO EXACT MATCH: Could not find exact matching log2FC column %s", exact_match))
      experiment_diagnostics <- rbind(experiment_diagnostics, data.frame(
        experiment_id = exp_id,
        log2fc_column = NA,
        pvalue_column = pval_col,
        status = "EXCLUDED",
        reason = "No exact matching log2FC column found - partial matches disallowed",
        na_count_lfc = NA,
        na_percent_lfc = NA,
        na_count_pval = pval_na_counts[pval_col],
        na_percent_pval = 100 * pval_na_counts[pval_col] / nrow(results_df),
        background_genes_count = NA,
        stringsAsFactors = FALSE
      ))
      next
    }
    
    # Count NAs in both columns
    lfc_na_count <- sum(is.na(results_df[[lfc_col]]))
    pval_na_count <- sum(is.na(results_df[[pval_col]]))
    
    # Create experiment-specific background gene list (non-NA genes in log2FC column)
    gene_col <- "gene_ID"
    if (!gene_col %in% colnames(results_df)) {
      gene_col <- colnames(results_df)[1]  # Fallback to first column
    }
    
    # Background genes are those without NA in log2FC column
    exp_background_genes <- results_df[[gene_col]][!is.na(results_df[[lfc_col]])]
    
    # Add experiment to valid list if it has both matching columns
    valid_experiments <- c(valid_experiments, exp_id)
    valid_lfc_cols <- c(valid_lfc_cols, lfc_col)
    valid_pval_cols <- c(valid_pval_cols, pval_col)
    experiment_background_genes[[exp_id]] <- exp_background_genes
    
    message(sprintf("  - INCLUDED: Experiment %s with log2FC column %s and p-value column %s", 
                   exp_id, lfc_col, pval_col))
    message(sprintf("  - Background genes: %d (%.1f%% of total genes)", 
                   length(exp_background_genes), 
                   100 * length(exp_background_genes) / nrow(results_df)))
    
    experiment_diagnostics <- rbind(experiment_diagnostics, data.frame(
      experiment_id = exp_id,
      log2fc_column = lfc_col,
      pvalue_column = pval_col,
      status = "INCLUDED",
      reason = "Matching log2FC and p-value columns found",
      na_count_lfc = lfc_na_count,
      na_percent_lfc = 100 * lfc_na_count / nrow(results_df),
      na_count_pval = pval_na_count,
      na_percent_pval = 100 * pval_na_count / nrow(results_df),
      background_genes_count = length(exp_background_genes),
      stringsAsFactors = FALSE
    ))
  }
  
  # Save experiment diagnostics
  write.csv(experiment_diagnostics, file.path(diagnostics_dir, "experiment_diagnostics.csv"), row.names = FALSE)
  
  # Check if we have any valid experiments
  if (length(valid_experiments) == 0) {
    message("  ERROR: No valid experiments with matching log2FC and p-value columns found!")
    message("  No analysis will be performed without valid experiment data")
    
    # Create metadata for failed analysis with enhanced diagnostics
    failed_metadata <- list(
      analysis_type = "MixScale_failed",
      perturbation = perturbation,
      cluster = cluster,
      date = Sys.Date(),
      error = "No valid experiments with matching log2FC and p-value columns",
      experiment_diagnostics = experiment_diagnostics,
      experiments_from_p_value_columns = experiment_ids_from_p,
      experiments_from_lfc_columns = experiment_ids_from_lfc,
      detailed_reason = "No experiments found with both matching p-value and log2FC columns. Check experiment naming consistency."
    )
    
    # Create metadata directory for failed analysis data
    metadata_dir <- file.path(output_dir, "metadata")
    dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Save failed metadata in metadata directory
    saveRDS(failed_metadata, file.path(metadata_dir, "analysis_failed_metadata.rds"))
    write.csv(experiment_diagnostics, file.path(metadata_dir, "excluded_experiments.csv"), row.names = FALSE)
    
    # Also save at original locations for backward compatibility
    saveRDS(failed_metadata, file.path(output_dir, "analysis_failed_metadata.rds"))
    write.csv(experiment_diagnostics, file.path(diagnostics_dir, "excluded_experiments.csv"), row.names = FALSE)
    return()
  }
  
  # Name vectors for easier access
  names(valid_lfc_cols) <- valid_experiments
  names(valid_pval_cols) <- valid_experiments
  
  # Print final mapping of valid experiments
  message("  Final validated experiment mappings:")
  for (i in seq_along(valid_experiments)) {
    exp_id <- valid_experiments[i]
    message(sprintf("  - Experiment %s: Using log2FC column %s with p-value column %s", 
                   exp_id, valid_lfc_cols[exp_id], valid_pval_cols[exp_id]))
  }
  
  # Filter genes for all valid experiments
  all_gene_lists <- filter_mixscale_genes(results_df, valid_lfc_cols, valid_pval_cols, 
                                         lfc_threshold, padj_threshold)
  
  # Process each experiment individually
  for(i in seq_along(valid_experiments)) {
    exp_id <- valid_experiments[i]
    print_header(paste("Processing experiment:", exp_id), 4)
    
    exp_dir <- file.path(output_dir, exp_id)
    dir.create(exp_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Get experiment-specific background genes
    exp_bg_genes <- experiment_background_genes[[exp_id]]
    message(sprintf("  Using experiment-specific background with %d genes", length(exp_bg_genes)))
    
    # Get gene lists for this experiment
    gene_lists <- all_gene_lists$by_experiment[[i]]
    
    # Save gene lists without headers
    write.table(gene_lists$up_genes, file.path(exp_dir, "up_genes.csv"), 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(gene_lists$down_genes, file.path(exp_dir, "down_genes.csv"), 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(gene_lists$all_de_genes, file.path(exp_dir, "all_de_genes.csv"), 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    write.table(exp_bg_genes, file.path(exp_dir, "background_genes.csv"), 
               row.names = FALSE, col.names = FALSE, quote = FALSE)
    
    # For GSEA, create ranked gene list
    if("GSEA" %in% run_methods) {
      # Get log2FC and p-value columns for this experiment
      lfc_col <- valid_lfc_cols[exp_id]
      padj_col <- valid_pval_cols[exp_id]
      
      # Create ranked list based on log2FC * -log10(padj)
      ranked_genes <- results_df
      
      # Filter out NA values for either log2FC or p-value
      ranked_genes <- ranked_genes[!is.na(ranked_genes[[lfc_col]]) & !is.na(ranked_genes[[padj_col]]), ]
      
      # Handle non-finite values appropriately
      clean_padj <- ifelse(!is.finite(ranked_genes[[padj_col]]) | ranked_genes[[padj_col]] <= 0, 
                           1e-300, ranked_genes[[padj_col]])
      clean_lfc <- ifelse(!is.finite(ranked_genes[[lfc_col]]), 0, ranked_genes[[lfc_col]])
      
      ranked_genes$ranking_score <- clean_lfc * -log10(clean_padj)
      # Replace any remaining non-finite values with a reasonable limit
      ranked_genes$ranking_score <- ifelse(!is.finite(ranked_genes$ranking_score), 0, ranked_genes$ranking_score)
      
      # Remove rows with zero ranking score
      ranked_genes <- ranked_genes[ranked_genes$ranking_score != 0, ]
      
      ranked_genes <- ranked_genes[order(ranked_genes$ranking_score, decreasing = TRUE), ]
      
      # Extract gene ids from a column for MixScale data
      gene_col <- "gene_ID"
      if (!gene_col %in% colnames(ranked_genes)) {
        gene_col <- colnames(ranked_genes)[1]  # Fallback to first column
      }
      
      ranked_list <- ranked_genes$ranking_score
      names(ranked_list) <- ranked_genes[[gene_col]]
      
      # Save ranked list
      write.table(data.frame(gene = names(ranked_list), score = ranked_list),
                 file.path(exp_dir, "ranked_gene_list.csv"),
                 row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
    }
    
    # Run enrichment analyses with experiment-specific background genes
    enrichment_results <- run_all_enrichment_analyses(
      gene_lists = gene_lists,
      background_genes = exp_bg_genes,
      ranked_list = if("GSEA" %in% run_methods) ranked_list else NULL,
      min_genes = min_genes,
      run_methods = run_methods,
      output_dir = exp_dir,
      padj_method = padj_method
    )
    
    # Create metadata
    is_weighted <- grepl("weight", valid_pval_cols[exp_id], fixed = TRUE)
    analysis_metadata <- list(
      analysis_type = "MixScale",
      perturbation = perturbation,
      cluster = cluster,
      experiment = exp_id,
      lfc_threshold = lfc_threshold,
      padj_threshold = padj_threshold,
      lfc_column = valid_lfc_cols[exp_id],
      padj_column = valid_pval_cols[exp_id],
      is_weighted_pvalue = is_weighted,
      date = Sys.Date(),
      original_metadata = metadata,
      up_genes_count = length(gene_lists$up_genes),
      down_genes_count = length(gene_lists$down_genes),
      all_de_genes_count = length(gene_lists$all_de_genes),
      background_genes_count = length(exp_bg_genes),
      uses_experiment_specific_background = TRUE
    )
    
    # Save metadata
    saveRDS(analysis_metadata, file.path(exp_dir, "analysis_metadata.rds"))
    
    # Save full DE table with flags for this experiment
    full_de_table <- results_df
    full_de_table$passes_threshold <- full_de_table[[gene_col]] %in% gene_lists$all_de_genes
    full_de_table$direction <- ifelse(full_de_table[[gene_col]] %in% gene_lists$up_genes, "UP",
                               ifelse(full_de_table[[gene_col]] %in% gene_lists$down_genes, "DOWN", "NONE"))
    write.csv(full_de_table, file.path(exp_dir, "full_de_table.csv"))
  }
  
  # Process combined analysis across all experiments
  print_header("Processing combined analysis across all experiments", 4)
  combined_dir <- file.path(output_dir, "combined")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Create combined background genes as intersection of all experiment backgrounds
  # This ensures we only use genes that have values in ALL experiments
  combined_bg_genes <- Reduce(intersect, experiment_background_genes)
  message(sprintf("  Combined background has %d genes (intersection of all experiment backgrounds)", 
                 length(combined_bg_genes)))
  
  # If combined background is too small, fall back to union of backgrounds
  if (length(combined_bg_genes) < 1000) {
    combined_bg_genes <- unique(unlist(experiment_background_genes))
    message(sprintf("  WARNING: Intersection background too small. Using union background with %d genes", 
                   length(combined_bg_genes)))
  }
  
  # Save combined background genes
  write.table(combined_bg_genes, file.path(combined_dir, "background_genes.csv"),
             row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Get combined gene lists
  combined_gene_lists <- all_gene_lists$combined
  
  # Save gene tracking information
  write.csv(all_gene_lists$gene_tracking, 
            file.path(combined_dir, "gene_tracking_across_experiments.csv"),
            row.names = FALSE)
  
  # Save lists of genes with inconsistent directionality
  write.table(combined_gene_lists$inconsistent_genes,
            file.path(combined_dir, "inconsistent_direction_genes.csv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Save other gene lists without headers
  write.table(combined_gene_lists$up_genes, 
            file.path(combined_dir, "up_genes.csv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(combined_gene_lists$down_genes, 
            file.path(combined_dir, "down_genes.csv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
  write.table(combined_gene_lists$all_de_genes, 
            file.path(combined_dir, "all_de_genes.csv"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # For GSEA, create a combined ranked gene list
  if("GSEA" %in% run_methods) {
    # Use the mean LFC and minimum p-value from gene tracking
    gene_tracking <- all_gene_lists$gene_tracking
    
    # Restrict to genes in combined background
    gene_tracking <- gene_tracking[gene_tracking$gene %in% combined_bg_genes, ]
    
    # Handle non-finite values
    gene_tracking$clean_lfc <- ifelse(!is.finite(gene_tracking$mean_lfc), 0, gene_tracking$mean_lfc)
    gene_tracking$clean_padj <- ifelse(!is.finite(gene_tracking$min_padj) | gene_tracking$min_padj <= 0, 
                                       1e-300, gene_tracking$min_padj)
    
    gene_tracking$ranking_score <- gene_tracking$clean_lfc * -log10(gene_tracking$clean_padj)
    
    # Remove rows with zero ranking score
    gene_tracking <- gene_tracking[gene_tracking$ranking_score != 0, ]
    
    gene_tracking <- gene_tracking[order(gene_tracking$ranking_score, decreasing = TRUE), ]
    ranked_list <- gene_tracking$ranking_score
    names(ranked_list) <- gene_tracking$gene
    
    # Save ranked list
    write.table(data.frame(gene = names(ranked_list), score = ranked_list),
               file.path(combined_dir, "ranked_gene_list.csv"),
               row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  }
  
  # Run enrichment analyses for combined results with combined background
  enrichment_results <- run_all_enrichment_analyses(
    gene_lists = list(
      up_genes = combined_gene_lists$up_genes,
      down_genes = combined_gene_lists$down_genes,
      all_de_genes = combined_gene_lists$all_de_genes
    ),
    background_genes = combined_bg_genes,
    ranked_list = if("GSEA" %in% run_methods) ranked_list else NULL,
    min_genes = min_genes,
    run_methods = run_methods,
    output_dir = combined_dir,
    padj_method = padj_method
  )
  
  # Create a detailed mapping summary for diagnostics
  mapping_summary <- data.frame(
    experiment_id = valid_experiments,
    p_value_column = valid_pval_cols,
    log2fc_column = valid_lfc_cols,
    is_exact_match = valid_lfc_cols == paste0("log2FC_", valid_experiments),
    background_genes_count = sapply(valid_experiments, function(exp) length(experiment_background_genes[[exp]])),
    stringsAsFactors = FALSE
  )
  
  # Save mapping summary to the diagnostics directory
  write.csv(mapping_summary, file.path(diagnostics_dir, "valid_experiment_mapping.csv"))
  
  # Create metadata for combined analysis with enhanced diagnostic information
  analysis_metadata <- list(
    analysis_type = "MixScale_combined",
    perturbation = perturbation,
    cluster = cluster,
    lfc_threshold = lfc_threshold,
    padj_threshold = padj_threshold,
    lfc_columns = paste(valid_lfc_cols, collapse=", "),
    padj_columns = paste(valid_pval_cols, collapse=", "),
    date = Sys.Date(),
    original_metadata = metadata,
    experiments = valid_experiments,
    experiments_included = valid_experiments,
    experiments_excluded = setdiff(experiment_ids_from_p, valid_experiments),
    experiment_diagnostics = experiment_diagnostics,
    experiments_from_p_value_columns = experiment_ids_from_p,
    experiments_from_lfc_columns = experiment_ids_from_lfc,
    column_mapping = mapping_summary,
    experiment_background_genes = experiment_background_genes,  # Store all background gene lists in metadata
    background_genes_strategy = ifelse(length(combined_bg_genes) == length(unique(unlist(experiment_background_genes))),
                                       "union", "intersection"),
    inconsistent_genes_count = length(combined_gene_lists$inconsistent_genes),
    up_genes_count = length(combined_gene_lists$up_genes),
    down_genes_count = length(combined_gene_lists$down_genes),
    all_de_genes_count = length(combined_gene_lists$all_de_genes),
    background_genes_count = length(combined_bg_genes),
    uses_experiment_specific_background = TRUE,
    processing_notes = "Only experiments with exact matching p-value and log2FC columns were included in the analysis."
  )
  
  # Create metadata directory and save detailed metadata
  metadata_dir <- file.path(combined_dir, "metadata")
  dir.create(metadata_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(analysis_metadata, file.path(metadata_dir, "analysis_metadata.rds"))
  
  # Also save at the root level for backward compatibility
  saveRDS(analysis_metadata, file.path(combined_dir, "analysis_metadata.rds"))
  
  # Also save a more accessible CSV version of the metadata
  metadata_df <- data.frame(
    parameter = c(
      "Analysis Type", 
      "Perturbation", 
      "Cluster", 
      "Log2FC Threshold", 
      "P-value Threshold",
      "Valid Experiments",
      "Excluded Experiments",
      "Log2FC Columns",
      "P-value Columns",
      "Background Strategy",
      "Experiment-Specific Background Used",
      "Inconsistent Genes Count",
      "Up-regulated Genes Count",
      "Down-regulated Genes Count",
      "Total DE Genes Count",
      "Background Genes Count",
      "Match Criteria",
      "Analysis Date"
    ),
    value = c(
      "MixScale_combined",
      perturbation,
      cluster,
      as.character(lfc_threshold),
      as.character(padj_threshold),
      paste(valid_experiments, collapse=", "),
      paste(setdiff(experiment_ids_from_p, valid_experiments), collapse=", "),
      paste(valid_lfc_cols, collapse=", "),
      paste(valid_pval_cols, collapse=", "),
      analysis_metadata$background_genes_strategy,
      "TRUE",
      as.character(length(combined_gene_lists$inconsistent_genes)),
      as.character(length(combined_gene_lists$up_genes)),
      as.character(length(combined_gene_lists$down_genes)),
      as.character(length(combined_gene_lists$all_de_genes)),
      as.character(length(combined_bg_genes)),
      "Exact match between p-value and log2FC columns required",
      as.character(Sys.Date())
    ),
    stringsAsFactors = FALSE
  )
  
  metadata_dir <- file.path(combined_dir, "metadata")
  write.csv(metadata_df, file.path(metadata_dir, "analysis_metadata_summary.csv"), row.names = FALSE)
  
  # Also save at the root level for backward compatibility
  write.csv(metadata_df, file.path(combined_dir, "analysis_metadata_summary.csv"), row.names = FALSE)
}

#' Run all requested enrichment analyses
#'
#' @param gene_lists List containing up_genes, down_genes, and all_de_genes
#' @param background_genes Vector of background genes
#' @param ranked_list Named vector of gene ranking scores for GSEA
#' @param min_genes Minimum number of genes for enrichment
#' @param run_methods Vector of analysis methods to run
#' @param output_dir Directory to save results
#' @param padj_method Method for p-value adjustment
#' @return List of enrichment results for each analysis type
run_all_enrichment_analyses <- function(gene_lists, background_genes, ranked_list = NULL,
                                      min_genes, run_methods, output_dir,
                                      padj_method = "BH") {
  # Store all results
  enrichment_results <- list()
  
  # Run GO enrichment if requested
  if("GO" %in% run_methods) {
    print_header("Gene Ontology (GO) Enrichment Analysis", 4)
    
    # Process each GO ontology separately - we now directly use "ALL" as a valid ontology
    ontologies <- c("BP", "CC", "MF", "ALL")
    
    for(ontology in ontologies) {
      go_dir <- file.path(output_dir, paste0("GO_", ontology))
      dir.create(go_dir, recursive = TRUE, showWarnings = FALSE)
      
      print_header(paste("GO", ontology, "Analysis"), 5)
      
      # Up-regulated genes
      if(length(gene_lists$up_genes) >= min_genes) {
        gene_type <- "up-regulated genes"
        print_header(sprintf("%s - %s (%d genes)", ontology, gene_type, length(gene_lists$up_genes)), 6)
        go_result <- run_go_enrichment(
          gene_list = gene_lists$up_genes,
          background_genes = background_genes,
          ontology = ontology,
          min_genes = min_genes,
          padj_method = padj_method
        )
        
        # Save result separately
        rds_name <- paste0("GO_", ontology, "_UP.rds")
        saveRDS(go_result, file.path(go_dir, rds_name))
      }
      
      # Down-regulated genes
      if(length(gene_lists$down_genes) >= min_genes) {
        gene_type <- "down-regulated genes"
        print_header(sprintf("%s - %s (%d genes)", ontology, gene_type, length(gene_lists$down_genes)), 6)
        go_result <- run_go_enrichment(
          gene_list = gene_lists$down_genes,
          background_genes = background_genes,
          ontology = ontology,
          min_genes = min_genes,
          padj_method = padj_method
        )
        
        # Save result separately
        rds_name <- paste0("GO_", ontology, "_DOWN.rds")
        saveRDS(go_result, file.path(go_dir, rds_name))
      }
      
      # All DE genes
      if(length(gene_lists$all_de_genes) >= min_genes) {
        gene_type <- "all DE genes"
        print_header(sprintf("%s - %s (%d genes)", ontology, gene_type, length(gene_lists$all_de_genes)), 6)
        go_result <- run_go_enrichment(
          gene_list = gene_lists$all_de_genes,
          background_genes = background_genes,
          ontology = ontology,
          min_genes = min_genes,
          padj_method = padj_method
        )
        
        # Save result separately
        rds_name <- paste0("GO_", ontology, "_ALL.rds")
        saveRDS(go_result, file.path(go_dir, rds_name))
      }
    }
  }
  
  # Run KEGG enrichment if requested

    # Create directory for KEGG results
    kegg_dir <- file.path(output_dir, "KEGG")
    dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
  if("KEGG" %in% run_methods) {
    print_header("KEGG Pathway Enrichment Analysis", 4)
    
    # Up-regulated genes
    if(length(gene_lists$up_genes) >= min_genes) {
      print_header(sprintf("KEGG - up-regulated genes (%d genes)", length(gene_lists$up_genes)), 5)
      kegg_result <- run_kegg_enrichment(
        gene_list = gene_lists$up_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Create KEGG directory and save result
      kegg_dir <- file.path(output_dir, "KEGG")
      dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(kegg_result, file.path(kegg_dir, "KEGG_UP.rds"))
    }
    
    # Down-regulated genes
    if(length(gene_lists$down_genes) >= min_genes) {
      print_header(sprintf("KEGG - down-regulated genes (%d genes)", length(gene_lists$down_genes)), 5)
      kegg_result <- run_kegg_enrichment(
        gene_list = gene_lists$down_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in KEGG directory
      kegg_dir <- file.path(output_dir, "KEGG")
      dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(kegg_result, file.path(kegg_dir, "KEGG_DOWN.rds"))
    }
    
    # All DE genes
    if(length(gene_lists$all_de_genes) >= min_genes) {
      print_header(sprintf("KEGG - all DE genes (%d genes)", length(gene_lists$all_de_genes)), 5)
      kegg_result <- run_kegg_enrichment(
        gene_list = gene_lists$all_de_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in KEGG directory
      kegg_dir <- file.path(output_dir, "KEGG")
      dir.create(kegg_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(kegg_result, file.path(kegg_dir, "KEGG_ALL.rds"))
    }
  }
  
  # Run Reactome enrichment if requested

    # Create directory for Reactome results
    reactome_dir <- file.path(output_dir, "Reactome")
    dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)
  if("Reactome" %in% run_methods) {
    print_header("Reactome Pathway Enrichment Analysis", 4)
    
    # Up-regulated genes
    if(length(gene_lists$up_genes) >= min_genes) {
      print_header(sprintf("Reactome - up-regulated genes (%d genes)", length(gene_lists$up_genes)), 5)
      reactome_result <- run_reactome_enrichment(
        gene_list = gene_lists$up_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Create Reactome directory and save result
      reactome_dir <- file.path(output_dir, "Reactome")
      dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(reactome_result, file.path(reactome_dir, "Reactome_UP.rds"))
    }
    
    # Down-regulated genes
    if(length(gene_lists$down_genes) >= min_genes) {
      print_header(sprintf("Reactome - down-regulated genes (%d genes)", length(gene_lists$down_genes)), 5)
      reactome_result <- run_reactome_enrichment(
        gene_list = gene_lists$down_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in Reactome directory
      reactome_dir <- file.path(output_dir, "Reactome")
      dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(reactome_result, file.path(reactome_dir, "Reactome_DOWN.rds"))
    }
    
    # All DE genes
    if(length(gene_lists$all_de_genes) >= min_genes) {
      print_header(sprintf("Reactome - all DE genes (%d genes)", length(gene_lists$all_de_genes)), 5)
      reactome_result <- run_reactome_enrichment(
        gene_list = gene_lists$all_de_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in Reactome directory
      reactome_dir <- file.path(output_dir, "Reactome")
      dir.create(reactome_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(reactome_result, file.path(reactome_dir, "Reactome_ALL.rds"))
    }
  }
  
  # Run WikiPathways enrichment if requested

    # Create directory for WikiPathways results
    wiki_dir <- file.path(output_dir, "WikiPathways")
    dir.create(wiki_dir, recursive = TRUE, showWarnings = FALSE)
  if("WikiPathways" %in% run_methods) {
    print_header("WikiPathways Enrichment Analysis", 4)
    
    # Up-regulated genes
    if(length(gene_lists$up_genes) >= min_genes) {
      print_header(sprintf("WikiPathways - up-regulated genes (%d genes)", length(gene_lists$up_genes)), 5)
      wiki_result <- run_wikipathways_enrichment(
        gene_list = gene_lists$up_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Create WikiPathways directory and save result
      wiki_dir <- file.path(output_dir, "WikiPathways")
      dir.create(wiki_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(wiki_result, file.path(wiki_dir, "WikiPathways_UP.rds"))
    }
    
    # Down-regulated genes
    if(length(gene_lists$down_genes) >= min_genes) {
      print_header(sprintf("WikiPathways - down-regulated genes (%d genes)", length(gene_lists$down_genes)), 5)
      wiki_result <- run_wikipathways_enrichment(
        gene_list = gene_lists$down_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in WikiPathways directory
      wiki_dir <- file.path(output_dir, "WikiPathways")
      dir.create(wiki_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(wiki_result, file.path(wiki_dir, "WikiPathways_DOWN.rds"))
    }
    
    # All DE genes
    if(length(gene_lists$all_de_genes) >= min_genes) {
      print_header(sprintf("WikiPathways - all DE genes (%d genes)", length(gene_lists$all_de_genes)), 5)
      wiki_result <- run_wikipathways_enrichment(
        gene_list = gene_lists$all_de_genes,
        background_genes = background_genes,
        min_genes = min_genes,
        padj_method = padj_method
      )
      
      # Save result in WikiPathways directory
      wiki_dir <- file.path(output_dir, "WikiPathways")
      dir.create(wiki_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(wiki_result, file.path(wiki_dir, "WikiPathways_ALL.rds"))
    }
  }
  
  # Run STRING PPI network analysis if requested
  if("STRING" %in% run_methods) {
    print_header("STRING Protein-Protein Interaction Analysis", 4)
    
    # Up-regulated genes
    if(length(gene_lists$up_genes) >= min_genes) {
      print_header(sprintf("STRING - up-regulated genes (%d genes)", length(gene_lists$up_genes)), 5)
      string_result <- run_string_ppi(
        gene_list = gene_lists$up_genes,
        background_genes = background_genes,
        min_genes = min_genes
      )
      
      # Create STRING directory and save result
      string_dir <- file.path(output_dir, "STRING")
      dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(string_result, file.path(string_dir, "STRING_UP.rds"))
    }
    
    # Down-regulated genes
    if(length(gene_lists$down_genes) >= min_genes) {
      print_header(sprintf("STRING - down-regulated genes (%d genes)", length(gene_lists$down_genes)), 5)
      string_result <- run_string_ppi(
        gene_list = gene_lists$down_genes,
        background_genes = background_genes,
        min_genes = min_genes
      )
      
      # Save result in STRING directory
      string_dir <- file.path(output_dir, "STRING")
      dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(string_result, file.path(string_dir, "STRING_DOWN.rds"))
    }
    
    # All DE genes
    if(length(gene_lists$all_de_genes) >= min_genes) {
      print_header(sprintf("STRING - all DE genes (%d genes)", length(gene_lists$all_de_genes)), 5)
      string_result <- run_string_ppi(
        gene_list = gene_lists$all_de_genes,
        background_genes = background_genes,
        min_genes = min_genes
      )
      
      # Save result in STRING directory
      string_dir <- file.path(output_dir, "STRING")
      dir.create(string_dir, recursive = TRUE, showWarnings = FALSE)
      saveRDS(string_result, file.path(string_dir, "STRING_ALL.rds"))
    }
  }
  
  # Run GSEA if requested and ranked list is provided
  if("GSEA" %in% run_methods && !is.null(ranked_list)) {
    print_header("Gene Set Enrichment Analysis (GSEA)", 4)
    
    # First, ensure the ranked list is clean
    print_header("Validating and cleaning ranked gene list", 5)
    
    # Check for non-finite values
    if(any(!is.finite(ranked_list))) {
      print_header(sprintf("Found %d non-finite values in ranked list, removing", sum(!is.finite(ranked_list))), 5)
      ranked_list <- ranked_list[is.finite(ranked_list)]
    }
    
    # Check for NULL or NA names
    if(any(is.na(names(ranked_list)) | names(ranked_list) == "")) {
      print_header(sprintf("Found %d genes with NA or empty names, removing", sum(is.na(names(ranked_list)) | names(ranked_list) == "")), 5)
      ranked_list <- ranked_list[!is.na(names(ranked_list)) & names(ranked_list) != ""]
    }
    
    # Check for duplicate names
    if(any(duplicated(names(ranked_list)))) {
      print_header(sprintf("Found %d duplicate gene names, keeping the one with highest absolute value", sum(duplicated(names(ranked_list)))), 5)
      # Create data frame for handling duplicates
      gene_df <- data.frame(
        gene = names(ranked_list),
        score = ranked_list,
        abs_score = abs(ranked_list),
        stringsAsFactors = FALSE
      )
      
      gene_df <- gene_df %>%
        group_by(gene) %>%
        arrange(desc(abs_score)) %>%
        slice(1) %>%
        ungroup()
      
      ranked_list <- setNames(gene_df$score, gene_df$gene)
    }
    
    print_header(sprintf("Final ranked list has %d genes with range [%.2f, %.2f]", 
                       length(ranked_list), min(ranked_list), max(ranked_list)), 5)
    
    # Run GSEA with the cleaned list
    gsea_results <- run_gsea(
      ranked_list = ranked_list,
      min_genes = min_genes,
      padj_method = padj_method
    )
    
    # Create GSEA directory and save results by category
    gsea_dir <- file.path(output_dir, "GSEA")
    dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
    for(category in names(gsea_results)) {
      if(!is.null(gsea_results[[category]])) {
        # Create subdirectory for each GSEA category
        category_dir <- file.path(gsea_dir, category)
        dir.create(category_dir, recursive = TRUE, showWarnings = FALSE)
        rds_name <- paste0("GSEA_", category, ".rds")
        saveRDS(gsea_results[[category]], file.path(category_dir, rds_name))
      }
    }
  }
  
  return(enrichment_results)
}

#' Filter genes by LFC and p-value thresholds (for MAST data)
#'
#' @param results_df Data frame of differential expression results
#' @param lfc_col Column name for log2 fold change
#' @param padj_col Column name for adjusted p-value
#' @param lfc_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @return List of up, down, and all DE genes
filter_genes <- function(results_df, lfc_col, padj_col, lfc_threshold, padj_threshold) {
  # Filter for up-regulated genes
  up_mask <- results_df[[lfc_col]] >= lfc_threshold & results_df[[padj_col]] <= padj_threshold
  up_genes <- rownames(results_df)[up_mask]
  
  # Filter for down-regulated genes
  down_mask <- results_df[[lfc_col]] <= -lfc_threshold & results_df[[padj_col]] <= padj_threshold
  down_genes <- rownames(results_df)[down_mask]
  
  # Combine for all DE genes
  all_de_genes <- c(up_genes, down_genes)
  
  print_header(sprintf("Gene Filtering Results: %d up-regulated, %d down-regulated, %d total DE genes", 
                 length(up_genes), length(down_genes), length(all_de_genes)), 4)
  
  return(list(
    up_genes = up_genes,
    down_genes = down_genes,
    all_de_genes = all_de_genes
  ))
}

#' Filter genes from MixScale data, handling multiple experiments
#'
#' @param results_df Data frame of differential expression results
#' @param lfc_cols Named vector of log2 fold change columns
#' @param padj_cols Named vector of adjusted p-value columns
#' @param lfc_threshold Log2 fold change threshold
#' @param padj_threshold Adjusted p-value threshold
#' @return List containing by_experiment results, combined results, and gene tracking information
filter_mixscale_genes <- function(results_df, lfc_cols, padj_cols, lfc_threshold, padj_threshold) {
  print_header("Filtering genes across multiple experiments", 4)
  
  # Get gene ID column
  gene_col <- "gene_ID"
  if(!gene_col %in% colnames(results_df)) {
    gene_col <- colnames(results_df)[1]  # Fallback to first column
  }
  
  # Initialize tracking data frame
  gene_tracking <- data.frame(
    gene = results_df[[gene_col]],
    passes_any = FALSE,
    direction = "NONE",
    consistent_direction = TRUE,
    experiments_passed = "",
    mean_lfc = NA,
    min_padj = NA,
    stringsAsFactors = FALSE
  )
  
  # Process each experiment
  exp_results <- list()
  
  for(i in seq_along(lfc_cols)) {
    lfc_col <- lfc_cols[i]
    padj_col <- padj_cols[i]
    
    # Ensure the columns exist
    if(!lfc_col %in% colnames(results_df) || !padj_col %in% colnames(results_df)) {
      message(sprintf("  Warning: Columns %s or %s not found in results_df", lfc_col, padj_col))
      next
    }
    
    # Get genes passing thresholds for this experiment
    up_mask <- results_df[[lfc_col]] >= lfc_threshold & results_df[[padj_col]] <= padj_threshold
    down_mask <- results_df[[lfc_col]] <= -lfc_threshold & results_df[[padj_col]] <= padj_threshold
    
    # Handle NA values
    up_mask[is.na(up_mask)] <- FALSE
    down_mask[is.na(down_mask)] <- FALSE
    
    up_genes <- results_df[[gene_col]][up_mask]
    down_genes <- results_df[[gene_col]][down_mask]
    all_de_genes <- c(up_genes, down_genes)
    
    message(sprintf("  Experiment %s: %d up-regulated, %d down-regulated, %d total DE genes", 
                   names(lfc_cols)[i], length(up_genes), length(down_genes), length(all_de_genes)))
    
    # Store by experiment
    exp_results[[i]] <- list(
      up_genes = up_genes,
      down_genes = down_genes,
      all_de_genes = all_de_genes
    )
    
    # Update tracking
    for(g in up_genes) {
      idx <- which(gene_tracking$gene == g)
      if(length(idx) == 0) next
      
      gene_tracking$passes_any[idx] <- TRUE
      
      # Check for direction consistency
      if(gene_tracking$direction[idx] == "NONE") {
        gene_tracking$direction[idx] <- "UP"
      } else if(gene_tracking$direction[idx] == "DOWN") {
        gene_tracking$consistent_direction[idx] <- FALSE
      }
      
      # Record experiment
      if(gene_tracking$experiments_passed[idx] == "") {
        gene_tracking$experiments_passed[idx] <- names(lfc_cols)[i]
      } else {
        gene_tracking$experiments_passed[idx] <- paste(gene_tracking$experiments_passed[idx],
                                                      names(lfc_cols)[i], sep=",")
      }
    }
    
    for(g in down_genes) {
      idx <- which(gene_tracking$gene == g)
      if(length(idx) == 0) next
      
      gene_tracking$passes_any[idx] <- TRUE
      
      # Check for direction consistency
      if(gene_tracking$direction[idx] == "NONE") {
        gene_tracking$direction[idx] <- "DOWN"
      } else if(gene_tracking$direction[idx] == "UP") {
        gene_tracking$consistent_direction[idx] <- FALSE
      }
      
      # Record experiment
      if(gene_tracking$experiments_passed[idx] == "") {
        gene_tracking$experiments_passed[idx] <- names(lfc_cols)[i]
      } else {
        gene_tracking$experiments_passed[idx] <- paste(gene_tracking$experiments_passed[idx],
                                                      names(lfc_cols)[i], sep=",")
      }
    }
  }
  
  # Calculate mean LFC and minimum p-value across experiments
  for(i in 1:nrow(gene_tracking)) {
    g <- gene_tracking$gene[i]
    idx <- which(results_df[[gene_col]] == g)
    if(length(idx) > 0) {
      gene_tracking$mean_lfc[i] <- mean(sapply(lfc_cols, function(col) {
        val <- results_df[idx, col]
        if(is.na(val) || !is.finite(val)) return(0)
        return(val)
      }), na.rm=TRUE)
      
      gene_tracking$min_padj[i] <- min(sapply(padj_cols, function(col) {
        val <- results_df[idx, col]
        if(is.na(val) || !is.finite(val)) return(1)
        return(val)
      }), na.rm=TRUE)
    }
  }
  
  # Get combined gene lists
  combined_up_genes <- gene_tracking$gene[gene_tracking$passes_any & gene_tracking$direction == "UP" & 
                                           gene_tracking$consistent_direction]
  combined_down_genes <- gene_tracking$gene[gene_tracking$passes_any & gene_tracking$direction == "DOWN" & 
                                             gene_tracking$consistent_direction]
  combined_inconsistent_genes <- gene_tracking$gene[gene_tracking$passes_any & 
                                                     !gene_tracking$consistent_direction]
  combined_all_de_genes <- gene_tracking$gene[gene_tracking$passes_any]
  
  message(sprintf("  Combined across experiments: %d up-regulated, %d down-regulated, %d inconsistent, %d total DE genes", 
                 length(combined_up_genes), length(combined_down_genes), 
                 length(combined_inconsistent_genes), length(combined_all_de_genes)))
  
  # Return both individual experiment results and combined results
  return(list(
    by_experiment = exp_results,
    combined = list(
      up_genes = combined_up_genes,
      down_genes = combined_down_genes,
      inconsistent_genes = combined_inconsistent_genes,
      all_de_genes = combined_all_de_genes
    ),
    gene_tracking = gene_tracking
  ))
}

#' Function to identify and prioritize weighted p-values in MixScale data
#'
#' @param results_df Data frame of differential expression results
#' @param experiment_ids Vector of experiment identifiers
#' @return Named vector of p-value column names
get_weighted_pvalue_cols <- function(results_df, experiment_ids) {
  print_header("Identifying p-value columns for experiments", 4)
  
  padj_cols <- character(length(experiment_ids))
  names(padj_cols) <- experiment_ids
  
  for(i in seq_along(experiment_ids)) {
    exp_id <- experiment_ids[i]
    
    # First priority: experiment-specific weighted p-value
    weighted_cols <- grep(paste0("p_.*", exp_id, ".*:weight"), colnames(results_df), value = TRUE)
    if(length(weighted_cols) > 0) {
      padj_cols[i] <- weighted_cols[1]
      message(sprintf("  Experiment %s: Using weighted p-value: %s", exp_id, weighted_cols[1]))
      next
    }
    
    # Second priority: global weighted p-value
    if("p_weight" %in% colnames(results_df)) {
      padj_cols[i] <- "p_weight"
      message(sprintf("  Experiment %s: Using global weighted p-value", exp_id))
      next
    }
    
    # Third priority: experiment-specific regular p-value
    p_cols <- grep(paste0("p_.*", exp_id), colnames(results_df), value = TRUE)
    if(length(p_cols) > 0) {
      padj_cols[i] <- p_cols[1]
      message(sprintf("  Experiment %s: Using regular p-value: %s", exp_id, p_cols[1]))
      next
    }
    
    # Fallback: any p-value column
    fallback_col <- grep("^p", colnames(results_df), value = TRUE)[1]
    if(length(fallback_col) > 0) {
      padj_cols[i] <- fallback_col
      message(sprintf("  Experiment %s: Using fallback p-value: %s", exp_id, fallback_col))
    } else {
      message(sprintf("  Experiment %s: No p-value column found", exp_id))
      padj_cols[i] <- NA
    }
  }
  
  return(padj_cols)
}

#' Run GO enrichment analysis 
#'
#' @param gene_list Vector of genes to test
#' @param background_genes Vector of background genes
#' @param ontology GO ontology to use ("BP", "CC", "MF", or "ALL")
#' @param min_genes Minimum number of genes required
#' @param padj_method Method for p-value adjustment
#' @return GO enrichment results
run_go_enrichment <- function(gene_list, background_genes, ontology = "ALL", 
                             min_genes = 5, padj_method = "BH") {
  if(length(gene_list) < min_genes) {
    message("    Not enough genes for GO enrichment analysis")
    return(NULL)
  }
  
  # Use ont="ALL" directly with enrichGO (it's a valid parameter)
  print_header(sprintf("Running GO enrichment with %s ontology", ontology), 6)
  message(sprintf("    Gene list: %d genes, background: %d genes", length(gene_list), length(background_genes)))
  
  result <- tryCatch({
    enrichGO(
      gene = gene_list,
      universe = background_genes,
      OrgDb = org.Hs.eg.db,
      keyType = "SYMBOL",
      ont = ontology,  # This works for "BP", "CC", "MF", and "ALL"
      pAdjustMethod = padj_method,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    message("    Error in GO enrichment: ", e$message)
    return(NULL)
  })
  
  # Report findings
  if(!is.null(result) && nrow(result) > 0) {
    top_terms <- head(result, 5)
    message(sprintf("    Found %d enriched terms; top 5 p-values: %s", 
                   nrow(result),
                   paste(signif(top_terms$pvalue, 3), collapse=", ")))
  } else {
    message("    No significant enrichments found")
  }
  
  return(result)
}

#' Run KEGG pathway enrichment analysis
#'
#' @param gene_list Vector of genes to test
#' @param background_genes Vector of background genes
#' @param min_genes Minimum number of genes required
#' @param padj_method Method for p-value adjustment
#' @return KEGG enrichment results
run_kegg_enrichment <- function(gene_list, background_genes, min_genes = 5, 
                               padj_method = "BH") {
  if(length(gene_list) < min_genes) {
    message("    Not enough genes for KEGG enrichment analysis")
    return(NULL)
  }
  
  message("    Converting gene symbols to ENTREZ IDs for KEGG analysis")
  # Convert gene symbols to ENTREZ IDs (required for KEGG)
  gene_list_entrez <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting gene symbols to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  background_entrez <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting background genes to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  # Store a mapping for later reference
  id_mapping <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error creating ID mapping: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  if(nrow(gene_list_entrez) < min_genes) {
    message(sprintf("    Only %d genes could be mapped to ENTREZ IDs (minimum %d required)",
                   nrow(gene_list_entrez), min_genes))
    return(NULL)
  }
  
  # Run KEGG enrichment
  message(sprintf("    Running KEGG pathway enrichment with %d ENTREZ IDs", nrow(gene_list_entrez)))
  kegg_result <- tryCatch({
    enrichKEGG(
      gene = gene_list_entrez$ENTREZID,
      universe = background_entrez$ENTREZID,
      organism = "hsa",
      pAdjustMethod = padj_method,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    message("    Error in KEGG enrichment: ", e$message)
    return(NULL)
  })
  
  # Report findings
  if(!is.null(kegg_result) && nrow(kegg_result) > 0) {
    top_terms <- head(kegg_result, 5)
    message(sprintf("    Found %d enriched pathways; top 5 p-values: %s", 
                   nrow(kegg_result),
                   paste(signif(top_terms$pvalue, 3), collapse=", ")))
  } else {
    message("    No significant KEGG pathways found")
  }
  
  return(list(
    result = kegg_result,
    id_mapping = id_mapping
  ))
}

#' Run Reactome pathway enrichment analysis
#'
#' @param gene_list Vector of genes to test
#' @param background_genes Vector of background genes
#' @param min_genes Minimum number of genes required
#' @param padj_method Method for p-value adjustment
#' @return Reactome enrichment results
run_reactome_enrichment <- function(gene_list, background_genes, min_genes = 5, 
                                   padj_method = "BH") {
  if(length(gene_list) < min_genes) {
    message("    Not enough genes for Reactome enrichment analysis")
    return(NULL)
  }
  
  message("    Converting gene symbols to ENTREZ IDs for Reactome analysis")
  # Convert gene symbols to ENTREZ IDs (required for Reactome)
  gene_list_entrez <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting gene symbols to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  background_entrez <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting background genes to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  # Store a mapping for later reference
  id_mapping <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error creating ID mapping: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  if(nrow(gene_list_entrez) < min_genes) {
    message(sprintf("    Only %d genes could be mapped to ENTREZ IDs (minimum %d required)",
                   nrow(gene_list_entrez), min_genes))
    return(NULL)
  }
  
  # Run Reactome enrichment
  message(sprintf("    Running Reactome pathway enrichment with %d ENTREZ IDs", nrow(gene_list_entrez)))
  reactome_result <- tryCatch({
    enrichPathway(
      gene = gene_list_entrez$ENTREZID,
      universe = background_entrez$ENTREZID,
      organism = "human",
      pAdjustMethod = padj_method,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    message("    Error in Reactome enrichment: ", e$message)
    return(NULL)
  })
  
  # Report findings
  if(!is.null(reactome_result) && nrow(reactome_result) > 0) {
    top_terms <- head(reactome_result, 5)
    message(sprintf("    Found %d enriched pathways; top 5 p-values: %s", 
                   nrow(reactome_result),
                   paste(signif(top_terms$pvalue, 3), collapse=", ")))
  } else {
    message("    No significant Reactome pathways found")
  }
  
  return(list(
    result = reactome_result,
    id_mapping = id_mapping
  ))
}

#' Run WikiPathways enrichment analysis
#'
#' @param gene_list Vector of genes to test
#' @param background_genes Vector of background genes
#' @param min_genes Minimum number of genes required
#' @param padj_method Method for p-value adjustment
#' @return WikiPathways enrichment results
run_wikipathways_enrichment <- function(gene_list, background_genes, min_genes = 5, 
                                       padj_method = "BH") {
  if(length(gene_list) < min_genes) {
    message("    Not enough genes for WikiPathways enrichment analysis")
    return(NULL)
  }
  
  message("    Converting gene symbols to ENTREZ IDs for WikiPathways analysis")
  # Convert gene symbols to ENTREZ IDs
  gene_list_entrez <- tryCatch({
    bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting gene symbols to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  background_entrez <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error converting background genes to ENTREZ IDs: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  # Store a mapping for later reference
  id_mapping <- tryCatch({
    bitr(background_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  }, error = function(e) {
    message("    Error creating ID mapping: ", e$message)
    return(data.frame(SYMBOL = character(0), ENTREZID = character(0)))
  })
  
  if(nrow(gene_list_entrez) < min_genes) {
    message(sprintf("    Only %d genes could be mapped to ENTREZ IDs (minimum %d required)",
                   nrow(gene_list_entrez), min_genes))
    return(NULL)
  }
  
  # Run WikiPathways enrichment
  message(sprintf("    Running WikiPathways enrichment with %d ENTREZ IDs", nrow(gene_list_entrez)))
  wiki_result <- tryCatch({
    enrichWP(
      gene = gene_list_entrez$ENTREZID,
      universe = background_entrez$ENTREZID,
      organism = "Homo sapiens",
      pAdjustMethod = padj_method,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
  }, error = function(e) {
    message("    Error in WikiPathways enrichment: ", e$message)
    return(NULL)
  })
  
  # Report findings
  if(!is.null(wiki_result) && nrow(wiki_result) > 0) {
    top_terms <- head(wiki_result, 5)
    message(sprintf("    Found %d enriched pathways; top 5 p-values: %s", 
                   nrow(wiki_result),
                   paste(signif(top_terms$pvalue, 3), collapse=", ")))
  } else {
    message("    No significant WikiPathways found")
  }
  
  return(list(
    result = wiki_result,
    id_mapping = id_mapping
  ))
}

#' Run STRING PPI network analysis
#'
#' @param gene_list Vector of genes to test
#' @param background_genes Vector of background genes
#' @param min_genes Minimum number of genes required
#' @return STRING PPI network results
run_string_ppi <- function(gene_list, background_genes, min_genes = 5) {
  if(length(gene_list) < min_genes) {
    message("    Not enough genes for STRING PPI analysis")
    return(NULL)
  }
  
  message("    Initializing STRING database connection")
  # Initialize STRING database connection
  string_db <- tryCatch({
    STRINGdb$new(version = "11.5", species = 9606,
                score_threshold = 400, input_directory = tempdir())
  }, error = function(e) {
    message("    Error initializing STRING database: ", e$message)
    return(NULL)
  })
  
  if(is.null(string_db)) {
    return(NULL)
  }
  
  # Map genes to STRING identifiers
  message("    Mapping genes to STRING identifiers")
  mapped_genes <- tryCatch({
    string_db$map(data.frame(gene = gene_list), "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    message("    Error mapping genes to STRING: ", e$message)
    return(data.frame())
  })
  
  if(nrow(mapped_genes) < min_genes) {
    message(sprintf("    Only %d genes could be mapped to STRING IDs (minimum %d required)",
                   nrow(mapped_genes), min_genes))
    return(NULL)
  }
  
  # Get the PPI network
  message(sprintf("    Retrieving protein-protein interactions for %d proteins", nrow(mapped_genes)))
  interactions <- tryCatch({
    string_db$get_interactions(mapped_genes$STRING_id)
  }, error = function(e) {
    message("    Error retrieving interactions: ", e$message)
    return(data.frame())
  })
  
  # Perform enrichment analysis directly with STRING
  message("    Running STRING enrichment analysis")
  enrichment <- tryCatch({
    string_db$get_enrichment(mapped_genes$STRING_id)
  }, error = function(e) {
    message("    Error in STRING enrichment: ", e$message)
    return(data.frame())
  })
  
  # Report findings
  if(!is.null(enrichment) && nrow(enrichment) > 0) {
    top_terms <- head(enrichment, 5)
    message(sprintf("    Found %d enriched terms; top 5 FDR values: %s", 
                   nrow(enrichment),
                   paste(signif(top_terms$fdr, 3), collapse=", ")))
  } else {
    message("    No significant STRING enrichments found")
  }
  
  return(list(
    mapped_genes = mapped_genes,
    interactions = interactions,
    enrichment = enrichment
  ))
}

#' Run Gene Set Enrichment Analysis (GSEA)
#'
#' @param ranked_list Named numeric vector with gene rankings
#' @param min_genes Minimum number of genes required
#' @param padj_method Method for p-value adjustment
#' @return GSEA results by category
run_gsea <- function(ranked_list, min_genes = 5, padj_method = "BH") {
  if(length(ranked_list) < min_genes) {
  # Create directory for GSEA results
  gsea_dir <- file.path(output_dir, "GSEA")
  dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)

    message("    Not enough genes for GSEA")
    return(NULL)
  }
  
  # Initialize results
  gsea_results <- list()
  
  # Make sure we're working with a clean ranked list
  if(any(!is.finite(ranked_list))) {
    message(sprintf("    Removing %d non-finite values from ranked list", sum(!is.finite(ranked_list))))
    ranked_list <- ranked_list[is.finite(ranked_list)]
  }
  
  # Remove genes with missing names
  bad_names <- which(is.na(names(ranked_list)) | names(ranked_list) == "")
  if(length(bad_names) > 0) {
    message(sprintf("    Removing %d genes with missing names", length(bad_names)))
    ranked_list <- ranked_list[-bad_names]
  }
  
  # Handle duplicate names
  if(any(duplicated(names(ranked_list)))) {
    message(sprintf("    Found %d duplicate gene names, keeping highest absolute value", sum(duplicated(names(ranked_list)))))
    
    # Create a data frame for handling duplicates
    gene_df <- data.frame(
      gene = names(ranked_list),
      score = ranked_list,
      abs_score = abs(ranked_list),
      stringsAsFactors = FALSE
    )
    
    gene_df <- gene_df %>%
      group_by(gene) %>%
      arrange(desc(abs_score)) %>%
      slice(1) %>%
      ungroup()
    
    ranked_list <- setNames(gene_df$score, gene_df$gene)
  }
  
  # Get MSigDB gene sets for human
  print_header("Running GSEA with MSigDB gene sets", 5)
  
  # Try different gene set collections, with error handling
  
  # 1. Run with Hallmark gene sets
  print_header("GSEA with Hallmark gene sets", 6)
  h_gene_sets <- tryCatch({
    msigdbr(species = "Homo sapiens", category = "H")
  }, error = function(e) {
    message("    Error getting Hallmark gene sets: ", e$message)
    return(NULL)
  })
  
  if(!is.null(h_gene_sets) && nrow(h_gene_sets) > 0) {
    # Prepare gene sets for FGSEA
    h_pathways <- split(h_gene_sets$gene_symbol, h_gene_sets$gs_name)
    
    # Ensure gene sets have at least minSize genes
    valid_pathways <- lapply(h_pathways, function(genes) {
      intersect(genes, names(ranked_list))
    })
    valid_pathways <- valid_pathways[sapply(valid_pathways, length) >= min_genes]
    
    if(length(valid_pathways) > 0) {
      message(sprintf("    Running fgsea with %d Hallmark gene sets", length(valid_pathways)))
      
      h_gsea <- tryCatch({
        fgseaMultilevel(
          pathways = valid_pathways,
          stats = ranked_list,
          minSize = min_genes,
          maxSize = 500,
          eps = 0,
          scoreType = "std",
          nPermSimple = 1000
        )
      }, error = function(e) {
        message("    Error in Hallmark GSEA: ", e$message)
        return(NULL)
      })
      
      if(!is.null(h_gsea) && nrow(h_gsea) > 0) {
        # Sort by p-value
        h_gsea <- h_gsea[order(h_gsea$pval), ]
        gsea_results[["Hallmark"]] <- h_gsea
        
        message(sprintf("    Found %d enriched Hallmark pathways; top 5 p-values: %s", 
                     nrow(h_gsea),
                     paste(signif(head(h_gsea$pval, 5), 3), collapse=", ")))
      } else {
        message("    No significant Hallmark pathways found")
      }
    } else {
      message("    No valid Hallmark pathways with minimum gene count")
    }
  }
  
  # 2. Run with GO gene sets (using fgsea instead of gseGO for better control)
  print_header("GSEA with GO gene sets", 6)
  go_gene_sets <- tryCatch({
    msigdbr(species = "Homo sapiens", category = "C5")
  }, error = function(e) {
    message("    Error getting GO gene sets: ", e$message)
    return(NULL)
  })
  
  if(!is.null(go_gene_sets) && nrow(go_gene_sets) > 0) {
    # Separate by subcategory (BP, CC, MF)
    subcats <- c("BP" = "GO:BP", "CC" = "GO:CC", "MF" = "GO:MF")
    
    for(sc_name in names(subcats)) {
      sc_code <- subcats[sc_name]
      sc_sets <- go_gene_sets[go_gene_sets$gs_subcat == sc_code, ]
      
      if(nrow(sc_sets) > 0) {
        # Prepare gene sets for FGSEA
        go_pathways <- split(sc_sets$gene_symbol, sc_sets$gs_name)
        
        # Ensure gene sets have at least minSize genes
        valid_pathways <- lapply(go_pathways, function(genes) {
          intersect(genes, names(ranked_list))
        })
        valid_pathways <- valid_pathways[sapply(valid_pathways, length) >= min_genes]
        
        if(length(valid_pathways) > 0) {
          message(sprintf("    Running fgsea with %d GO:%s gene sets", length(valid_pathways), sc_name))
          
          go_gsea <- tryCatch({
            fgseaMultilevel(
              pathways = valid_pathways,
              stats = ranked_list,
              minSize = min_genes,
              maxSize = 500,
              eps = 0,
              scoreType = "std",
              nPermSimple = 1000
            )
          }, error = function(e) {
            message("    Error in GO:", sc_name, " GSEA: ", e$message)
            return(NULL)
          })
          
          if(!is.null(go_gsea) && nrow(go_gsea) > 0) {
            # Sort by p-value
            go_gsea <- go_gsea[order(go_gsea$pval), ]
            gsea_results[[paste0("GO_", sc_name)]] <- go_gsea
            
            message(sprintf("    Found %d enriched GO:%s pathways; top 5 p-values: %s", 
                         nrow(go_gsea), sc_name,
                         paste(signif(head(go_gsea$pval, 5), 3), collapse=", ")))
          } else {
            message("    No significant GO:", sc_name, " pathways found")
          }
        } else {
          message("    No valid GO:", sc_name, " pathways with minimum gene count")
        }
      }
    }
  }
  
  # 3. Run with KEGG pathways
  print_header("GSEA with KEGG pathways", 6)
  kegg_gene_sets <- tryCatch({
    subset(msigdbr(species = "Homo sapiens", category = "C2"), gs_subcat == "CP:KEGG")
  }, error = function(e) {
    message("    Error getting KEGG gene sets: ", e$message)
    return(NULL)
  })
  
  if(!is.null(kegg_gene_sets) && nrow(kegg_gene_sets) > 0) {
    # Prepare gene sets for FGSEA
    kegg_pathways <- split(kegg_gene_sets$gene_symbol, kegg_gene_sets$gs_name)
    
    # Ensure gene sets have at least minSize genes
    valid_pathways <- lapply(kegg_pathways, function(genes) {
      intersect(genes, names(ranked_list))
    })
    valid_pathways <- valid_pathways[sapply(valid_pathways, length) >= min_genes]
    
    if(length(valid_pathways) > 0) {
      message(sprintf("    Running fgsea with %d KEGG pathways", length(valid_pathways)))
      
      kegg_gsea <- tryCatch({
        fgseaMultilevel(
          pathways = valid_pathways,
          stats = ranked_list,
          minSize = min_genes,
          maxSize = 500,
          eps = 0,
          scoreType = "std",
          nPermSimple = 1000
        )
      }, error = function(e) {
        message("    Error in KEGG GSEA: ", e$message)
        return(NULL)
      })
      
      if(!is.null(kegg_gsea) && nrow(kegg_gsea) > 0) {
        # Sort by p-value
        kegg_gsea <- kegg_gsea[order(kegg_gsea$pval), ]
        gsea_results[["KEGG"]] <- kegg_gsea
        
        message(sprintf("    Found %d enriched KEGG pathways; top 5 p-values: %s", 
                     nrow(kegg_gsea),
                     paste(signif(head(kegg_gsea$pval, 5), 3), collapse=", ")))
      } else {
        message("    No significant KEGG pathways found")
      }
    } else {
      message("    No valid KEGG pathways with minimum gene count")
    }
  }
  
  return(gsea_results)
}

# Example usage (commented out for package):
# run_enrichment_analysis(
#   input_file = "full_DE_results.rds",
#   lfc_threshold = 0.25,
#   padj_threshold = 0.05,
#   output_dir = "./enrichment_results/",
#   run_methods = c("GO", "GSEA", "KEGG", "Reactome", "WikiPathways", "STRING")
# )
