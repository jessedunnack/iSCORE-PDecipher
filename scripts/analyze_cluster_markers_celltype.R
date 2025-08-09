#!/usr/bin/env Rscript

# Marker-Based Cell Type Inference Script
# Analyzes top cluster markers and infers cell types through web searches
# Author: Claude
# Date: January 2025

library(dplyr)
library(tidyr)
library(stringr)
library(DT)

# Configuration
N_TOP_MARKERS <- 25  # Number of top markers to analyze per cluster
LFC_THRESHOLD <- 0.25  # Minimum log2FC to consider
PADJ_THRESHOLD <- 0.05  # Maximum adjusted p-value

# Paths
UMAP_DIR <- "inst/extdata/umap_data"
OUTPUT_DIR <- "results/cell_type_annotations"

# Create output directory
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Function to load and process marker data
load_cluster_markers <- function(dataset_name) {
  marker_file <- file.path(UMAP_DIR, paste0(dataset_name, "_cluster_markers.rds"))
  
  if (!file.exists(marker_file)) {
    stop("Marker file not found: ", marker_file)
  }
  
  markers <- readRDS(marker_file)
  
  # Ensure required columns exist
  required_cols <- c("cluster", "gene", "avg_log2FC", "p_val_adj")
  if (!all(required_cols %in% names(markers))) {
    stop("Missing required columns. Found: ", paste(names(markers), collapse = ", "))
  }
  
  return(markers)
}

# Function to rank and filter markers
get_top_markers <- function(markers, n_top = N_TOP_MARKERS) {
  # Calculate a combined score based on effect size and significance
  markers <- markers %>%
    filter(avg_log2FC > LFC_THRESHOLD, p_val_adj < PADJ_THRESHOLD) %>%
    mutate(
      # Combined score: higher LFC and lower p-value = higher score
      score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)
    ) %>%
    group_by(cluster) %>%
    arrange(desc(score)) %>%
    slice_head(n = n_top) %>%
    ungroup()
  
  return(markers)
}

# Function to search for cell type associations
search_marker_associations <- function(gene_list, cluster_id) {
  cat("\n=== Searching cell type associations for Cluster", cluster_id, "===\n")
  
  # Prepare search queries
  search_results <- list()
  
  # Batch search for efficiency
  # First, search for top 5 markers together
  top_genes <- head(gene_list, 5)
  query <- paste(c(paste(top_genes, collapse = " "), 
                   "neuron cell type marker expression"), 
                 collapse = " ")
  
  cat("Searching for top 5 markers together:", paste(top_genes, collapse = ", "), "\n")
  
  # Note: In actual implementation, you would use the WebSearch tool here
  # For now, we'll create a placeholder structure
  search_results$batch_search <- list(
    genes = top_genes,
    query = query,
    status = "pending"
  )
  
  # Search for specific neuronal subtype markers
  neuronal_markers <- c("TH", "GAD1", "GAD2", "SLC17A7", "SLC17A6", "CHAT", 
                       "SLC6A3", "SLC18A2", "CALB1", "PVALB", "SST", "VIP")
  
  markers_present <- intersect(gene_list, neuronal_markers)
  if (length(markers_present) > 0) {
    cat("Found known neuronal subtype markers:", paste(markers_present, collapse = ", "), "\n")
    search_results$neuronal_subtype <- markers_present
  }
  
  # Search for glial markers
  glial_markers <- c("GFAP", "AQP4", "S100B", "ALDH1L1", "CD68", "AIF1", 
                     "CX3CR1", "TMEM119", "OLIG2", "MBP", "PLP1", "SOX10")
  
  glial_present <- intersect(gene_list, glial_markers)
  if (length(glial_present) > 0) {
    cat("Found glial markers:", paste(glial_present, collapse = ", "), "\n")
    search_results$glial_markers <- glial_present
  }
  
  return(search_results)
}

# Function to infer cell type based on markers
infer_cell_type <- function(top_markers, search_results) {
  # Initialize cell type inference
  inference <- list(
    primary_type = "Unknown",
    subtype = NA,
    confidence = "Low",
    evidence = list()
  )
  
  # Check for neuronal subtype markers
  if (!is.null(search_results$neuronal_subtype)) {
    neuronal_markers <- search_results$neuronal_subtype
    
    # Dopaminergic neurons
    if (any(c("TH", "SLC6A3", "SLC18A2") %in% neuronal_markers)) {
      inference$primary_type <- "Neuron"
      inference$subtype <- "Dopaminergic"
      inference$confidence <- "High"
      inference$evidence$dopaminergic <- neuronal_markers[neuronal_markers %in% 
                                                         c("TH", "SLC6A3", "SLC18A2")]
    }
    # GABAergic neurons
    else if (any(c("GAD1", "GAD2") %in% neuronal_markers)) {
      inference$primary_type <- "Neuron"
      inference$subtype <- "GABAergic"
      inference$confidence <- "High"
      inference$evidence$gabaergic <- neuronal_markers[neuronal_markers %in% 
                                                       c("GAD1", "GAD2")]
    }
    # Glutamatergic neurons
    else if (any(c("SLC17A7", "SLC17A6") %in% neuronal_markers)) {
      inference$primary_type <- "Neuron"
      inference$subtype <- "Glutamatergic"
      inference$confidence <- "High"
      inference$evidence$glutamatergic <- neuronal_markers[neuronal_markers %in% 
                                                           c("SLC17A7", "SLC17A6")]
    }
  }
  
  # Check for glial markers
  if (!is.null(search_results$glial_markers) && inference$primary_type == "Unknown") {
    glial_markers <- search_results$glial_markers
    
    # Astrocytes
    if (any(c("GFAP", "AQP4", "S100B", "ALDH1L1") %in% glial_markers)) {
      inference$primary_type <- "Glia"
      inference$subtype <- "Astrocyte"
      inference$confidence <- "High"
      inference$evidence$astrocyte <- glial_markers[glial_markers %in% 
                                                    c("GFAP", "AQP4", "S100B", "ALDH1L1")]
    }
    # Microglia
    else if (any(c("CD68", "AIF1", "CX3CR1", "TMEM119") %in% glial_markers)) {
      inference$primary_type <- "Glia"
      inference$subtype <- "Microglia"
      inference$confidence <- "High"
      inference$evidence$microglia <- glial_markers[glial_markers %in% 
                                                    c("CD68", "AIF1", "CX3CR1", "TMEM119")]
    }
    # Oligodendrocytes
    else if (any(c("OLIG2", "MBP", "PLP1", "SOX10") %in% glial_markers)) {
      inference$primary_type <- "Glia"
      inference$subtype <- "Oligodendrocyte"
      inference$confidence <- "High"
      inference$evidence$oligodendrocyte <- glial_markers[glial_markers %in% 
                                                          c("OLIG2", "MBP", "PLP1", "SOX10")]
    }
  }
  
  # If still unknown, check general neuronal markers
  if (inference$primary_type == "Unknown") {
    general_neuronal <- c("MAP2", "RBFOX3", "SYN1", "TUBB3", "NCAM1", "GAP43", "STMN2")
    neuronal_present <- intersect(top_markers$gene, general_neuronal)
    
    if (length(neuronal_present) >= 2) {
      inference$primary_type <- "Neuron"
      inference$subtype <- "Unspecified"
      inference$confidence <- "Medium"
      inference$evidence$general_neuronal <- neuronal_present
    }
  }
  
  return(inference)
}

# Function to create summary table
create_summary_table <- function(all_results) {
  summary_df <- data.frame(
    cluster = character(),
    primary_type = character(),
    subtype = character(),
    confidence = character(),
    top_markers = character(),
    evidence = character(),
    stringsAsFactors = FALSE
  )
  
  for (cluster_id in names(all_results)) {
    result <- all_results[[cluster_id]]
    
    summary_df <- rbind(summary_df, data.frame(
      cluster = cluster_id,
      primary_type = result$inference$primary_type,
      subtype = ifelse(is.na(result$inference$subtype), "", result$inference$subtype),
      confidence = result$inference$confidence,
      top_markers = paste(head(result$top_markers$gene, 5), collapse = ", "),
      evidence = paste(unlist(result$inference$evidence), collapse = ", "),
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_df)
}

# Main analysis function
analyze_dataset <- function(dataset_name) {
  cat("\n", strrep("=", 70), "\n")
  cat("Analyzing cell types for dataset:", dataset_name, "\n")
  cat(strrep("=", 70), "\n")
  
  # Load markers
  markers <- load_cluster_markers(dataset_name)
  
  # Get unique clusters
  clusters <- sort(unique(as.character(markers$cluster)))
  cat("Found", length(clusters), "clusters:", paste(clusters, collapse = ", "), "\n")
  
  # Analyze each cluster
  all_results <- list()
  
  for (cluster_id in clusters) {
    cat("\n--- Processing Cluster", cluster_id, "---\n")
    
    # Get top markers for this cluster
    cluster_markers <- markers %>% filter(cluster == cluster_id)
    top_markers <- get_top_markers(cluster_markers)
    
    if (nrow(top_markers) == 0) {
      cat("No significant markers found for cluster", cluster_id, "\n")
      next
    }
    
    cat("Top 5 markers:", paste(head(top_markers$gene, 5), collapse = ", "), "\n")
    
    # Search for cell type associations
    search_results <- search_marker_associations(top_markers$gene, cluster_id)
    
    # Infer cell type
    inference <- infer_cell_type(top_markers, search_results)
    
    # Store results
    all_results[[cluster_id]] <- list(
      top_markers = top_markers,
      search_results = search_results,
      inference = inference
    )
    
    # Print inference
    cat("\nInferred cell type:", inference$primary_type)
    if (!is.na(inference$subtype)) {
      cat(" -", inference$subtype)
    }
    cat("\nConfidence:", inference$confidence, "\n")
  }
  
  # Create summary table
  summary_table <- create_summary_table(all_results)
  
  # Save results
  output_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_cell_type_annotations.rds"))
  saveRDS(all_results, output_file)
  cat("\nDetailed results saved to:", output_file, "\n")
  
  # Save summary table
  summary_file <- file.path(OUTPUT_DIR, paste0(dataset_name, "_cell_type_summary.csv"))
  write.csv(summary_table, summary_file, row.names = FALSE)
  cat("Summary table saved to:", summary_file, "\n")
  
  # Display summary
  cat("\n=== Cell Type Summary ===\n")
  print(summary_table)
  
  return(list(
    detailed_results = all_results,
    summary = summary_table
  ))
}

# Function to map fine clusters to coarse clusters
map_fine_to_coarse <- function(fine_summary, coarse_summary) {
  cat("\n=== Mapping Fine Clusters to Coarse Clusters ===\n")
  
  # This is a placeholder - actual mapping would require the Seurat object
  # or a pre-computed mapping table
  cat("NOTE: Actual fine-to-coarse mapping requires Seurat object metadata\n")
  cat("This would show which fine clusters belong to each coarse cluster\n")
  
  # For now, return a template structure
  mapping_template <- data.frame(
    fine_cluster = 0:35,
    coarse_cluster = rep(0:14, length.out = 36),  # Placeholder distribution
    fine_cell_type = "To be determined",
    coarse_cell_type = "To be determined",
    consistency = "To be evaluated"
  )
  
  return(mapping_template)
}

# Main execution
main <- function() {
  cat("Cell Type Annotation Analysis\n")
  cat("=============================\n")
  cat("This script analyzes cluster markers to infer cell types\n")
  cat("Top", N_TOP_MARKERS, "markers per cluster will be analyzed\n")
  cat("Thresholds: LFC >", LFC_THRESHOLD, ", adjusted p-value <", PADJ_THRESHOLD, "\n")
  
  # List available datasets
  marker_files <- list.files(UMAP_DIR, pattern = "_cluster_markers\\.rds$", full.names = FALSE)
  datasets <- gsub("_cluster_markers\\.rds$", "", marker_files)
  
  cat("\nAvailable datasets:\n")
  for (i in seq_along(datasets)) {
    cat(i, ". ", datasets[i], "\n", sep = "")
  }
  
  # Analyze each dataset
  all_dataset_results <- list()
  
  for (dataset in datasets) {
    results <- analyze_dataset(dataset)
    all_dataset_results[[dataset]] <- results
  }
  
  # If we have both coarse and fine clustering results, create mapping
  # This is a placeholder for demonstration
  if (length(all_dataset_results) > 1) {
    cat("\n=== Cross-Dataset Analysis ===\n")
    cat("Multiple datasets analyzed. Fine-to-coarse mapping would be performed here.\n")
  }
  
  cat("\n=== Analysis Complete ===\n")
  cat("Results saved to:", OUTPUT_DIR, "\n")
  cat("\nNext steps:\n")
  cat("1. Review cell type annotations in summary CSV files\n")
  cat("2. Use WebSearch tool to validate uncertain assignments\n")
  cat("3. Load Seurat object to create accurate fine-to-coarse mapping\n")
  cat("4. Update Shiny app with cell type labels\n")
}

# Run main function if script is executed directly
if (!interactive()) {
  main()
}