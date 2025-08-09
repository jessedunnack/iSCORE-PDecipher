#!/usr/bin/env Rscript

# Comprehensive Marker Analysis Across All Fine Clusters
# Extends functional categorization to all clusters and identifies patterns

library(dplyr)
library(tidyr)

# Source the functional categories
source("scripts/functional_gene_categorization.R")

# Function to analyze all clusters comprehensively
analyze_all_clusters <- function(marker_file) {
  markers <- readRDS(marker_file)
  clusters <- sort(unique(as.character(markers$cluster)))
  
  all_results <- list()
  all_top_genes <- list()
  
  cat("Analyzing all", length(clusters), "clusters...\n\n")
  
  for (cluster_id in clusters) {
    # Get significant markers
    cluster_markers <- markers %>%
      filter(cluster == cluster_id, avg_log2FC > 0.25, p_val_adj < 0.05) %>%
      mutate(score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)) %>%
      arrange(desc(score))
    
    # Store top genes
    top_genes <- head(cluster_markers$gene, 100)  # Get more genes for better coverage
    all_top_genes[[cluster_id]] <- data.frame(
      cluster = cluster_id,
      gene = top_genes,
      rank = 1:length(top_genes),
      stringsAsFactors = FALSE
    )
    
    # Categorize genes
    categories <- categorize_genes(top_genes)
    all_results[[cluster_id]] <- categories
  }
  
  return(list(
    results = all_results,
    top_genes = all_top_genes
  ))
}

# Function to identify frequently uncategorized genes
identify_novel_markers <- function(all_results) {
  uncategorized_counts <- list()
  
  for (cluster in names(all_results)) {
    if ("uncategorized" %in% names(all_results[[cluster]])) {
      genes <- all_results[[cluster]]$uncategorized$genes
      for (gene in genes) {
        if (gene %in% names(uncategorized_counts)) {
          uncategorized_counts[[gene]] <- uncategorized_counts[[gene]] + 1
        } else {
          uncategorized_counts[[gene]] <- 1
        }
      }
    }
  }
  
  # Sort by frequency
  uncategorized_df <- data.frame(
    gene = names(uncategorized_counts),
    frequency = unlist(uncategorized_counts),
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(frequency))
  
  return(uncategorized_df)
}

# Function to create cluster functional profiles
create_functional_profiles <- function(all_results, annotations) {
  profiles <- data.frame(
    cluster = character(),
    cell_type = character(),
    top_categories = character(),
    key_genes = character(),
    functional_signature = character(),
    stringsAsFactors = FALSE
  )
  
  for (cluster_id in names(all_results)) {
    cluster_cats <- all_results[[cluster_id]]
    
    # Get cell type annotation
    cluster_num <- as.integer(cluster_id)
    cell_type <- if (cluster_num %in% annotations$fine_cluster) {
      annotations$cell_type[annotations$fine_cluster == cluster_num]
    } else {
      "Unknown"
    }
    
    # Get top 3 categories (excluding uncategorized)
    cat_counts <- sapply(cluster_cats, function(x) x$count)
    cat_counts <- cat_counts[names(cat_counts) != "uncategorized"]
    top_cats <- head(names(sort(cat_counts, decreasing = TRUE)), 3)
    
    # Get key genes from top categories
    key_genes <- character()
    for (cat in top_cats) {
      if (cat %in% names(cluster_cats)) {
        key_genes <- c(key_genes, head(cluster_cats[[cat]]$genes, 2))
      }
    }
    
    # Determine functional signature
    signature <- determine_functional_signature(cluster_cats)
    
    profiles <- rbind(profiles, data.frame(
      cluster = cluster_id,
      cell_type = cell_type,
      top_categories = paste(top_cats, collapse = "; "),
      key_genes = paste(unique(key_genes), collapse = ", "),
      functional_signature = signature,
      stringsAsFactors = FALSE
    ))
  }
  
  return(profiles)
}

# Function to determine functional signature
determine_functional_signature <- function(categories) {
  signatures <- character()
  
  # Check for neuronal signatures
  neuronal_cats <- c("pan_neuronal_mature", "pan_neuronal_immature", "pan_neuronal_structural",
                     "nt_dopaminergic", "nt_gabaergic", "nt_glutamatergic")
  if (any(neuronal_cats %in% names(categories))) {
    signatures <- c(signatures, "Neuronal")
  }
  
  # Check for progenitor signatures
  prog_cats <- c("tf_progenitor", "dev_early_progenitor", "cell_cycle_g2m", "cell_cycle_s")
  if (any(prog_cats %in% names(categories))) {
    signatures <- c(signatures, "Progenitor")
  }
  
  # Check for glial signatures
  glial_cats <- c("oligodendrocyte", "astrocyte", "ependymal", "meningeal")
  if (any(glial_cats %in% names(categories))) {
    signatures <- c(signatures, "Glial")
  }
  
  # Check for vascular signatures
  if ("vascular" %in% names(categories)) {
    signatures <- c(signatures, "Vascular")
  }
  
  # Check for ECM signatures
  ecm_cats <- c("ecm_collagens", "ecm_proteoglycans", "ecm_glycoproteins")
  if (any(ecm_cats %in% names(categories))) {
    signatures <- c(signatures, "ECM-rich")
  }
  
  # Check for stress signatures
  stress_cats <- c("stress_er", "stress_oxidative", "stress_general", "mitochondrial")
  if (any(stress_cats %in% names(categories))) {
    signatures <- c(signatures, "Stressed")
  }
  
  if (length(signatures) == 0) {
    signatures <- "Uncharacterized"
  }
  
  return(paste(signatures, collapse = "/"))
}

# Function to create heatmap data for visualization
create_category_heatmap_data <- function(all_results) {
  # Get all category names
  all_categories <- unique(unlist(lapply(all_results, function(x) names(x))))
  all_categories <- setdiff(all_categories, "uncategorized")
  
  # Create matrix
  clusters <- names(all_results)
  heatmap_matrix <- matrix(0, 
                          nrow = length(all_categories),
                          ncol = length(clusters),
                          dimnames = list(all_categories, clusters))
  
  # Fill matrix
  for (cluster in clusters) {
    cluster_cats <- all_results[[cluster]]
    for (cat in names(cluster_cats)) {
      if (cat != "uncategorized" && cat %in% all_categories) {
        heatmap_matrix[cat, cluster] <- cluster_cats[[cat]]$count
      }
    }
  }
  
  # Order by total gene count
  row_order <- order(rowSums(heatmap_matrix), decreasing = TRUE)
  heatmap_matrix <- heatmap_matrix[row_order, ]
  
  return(heatmap_matrix)
}

# Function to generate plotting recommendations
generate_plot_recommendations <- function(profiles, heatmap_data) {
  recommendations <- list()
  
  # 1. Transcription factor networks
  tf_rows <- grep("^tf_", rownames(heatmap_data))
  if (length(tf_rows) > 0) {
    tf_clusters <- names(which(colSums(heatmap_data[tf_rows, , drop = FALSE]) > 2))
    recommendations$tf_network <- list(
      title = "Transcription Factor Networks",
      clusters = tf_clusters,
      description = "Clusters with rich TF expression for cell fate analysis"
    )
  }
  
  # 2. Neuronal maturation trajectory
  neuronal_clusters <- profiles$cluster[grep("Neuronal", profiles$functional_signature)]
  if (length(neuronal_clusters) > 0) {
    recommendations$neuronal_trajectory <- list(
      title = "Neuronal Maturation Trajectory",
      clusters = neuronal_clusters,
      description = "Track neuronal differentiation stages"
    )
  }
  
  # 3. Cell adhesion/migration
  adhesion_rows <- grep("adhesion|axon", rownames(heatmap_data))
  if (length(adhesion_rows) > 0) {
    adhesion_clusters <- names(which(colSums(heatmap_data[adhesion_rows, , drop = FALSE]) > 2))
    recommendations$adhesion_migration <- list(
      title = "Cell Adhesion & Migration",
      clusters = adhesion_clusters,
      description = "Clusters involved in neural connectivity"
    )
  }
  
  # 4. Regional identity
  regional_rows <- grep("regional|hox", rownames(heatmap_data), ignore.case = TRUE)
  if (length(regional_rows) > 0) {
    regional_clusters <- names(which(heatmap_data[regional_rows, , drop = FALSE] > 0, arr.ind = TRUE)[, "col"])
    recommendations$regional_identity <- list(
      title = "Regional Identity Markers",
      clusters = unique(regional_clusters),
      description = "Rostrocaudal and dorsoventral patterning"
    )
  }
  
  return(recommendations)
}

# Main comprehensive analysis
main <- function() {
  cat("Comprehensive Functional Marker Analysis\n")
  cat("=======================================\n\n")
  
  # Load annotations
  annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv",
                         stringsAsFactors = FALSE)
  
  # Analyze all clusters
  marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
  analysis_results <- analyze_all_clusters(marker_file)
  
  # Identify novel markers
  novel_markers <- identify_novel_markers(analysis_results$results)
  cat("\nTop 20 frequently uncategorized genes:\n")
  print(head(novel_markers, 20))
  
  # Create functional profiles
  profiles <- create_functional_profiles(analysis_results$results, annotations)
  
  # Create heatmap data
  heatmap_data <- create_category_heatmap_data(analysis_results$results)
  
  # Generate recommendations
  recommendations <- generate_plot_recommendations(profiles, heatmap_data)
  
  # Save all results
  dir.create("results/functional_categorization", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(analysis_results, "results/functional_categorization/comprehensive_analysis.rds")
  write.csv(novel_markers, "results/functional_categorization/novel_markers.csv", row.names = FALSE)
  write.csv(profiles, "results/functional_categorization/cluster_functional_profiles.csv", row.names = FALSE)
  write.csv(heatmap_data, "results/functional_categorization/category_heatmap_data.csv")
  saveRDS(recommendations, "results/functional_categorization/plotting_recommendations.rds")
  
  # Print summary
  cat("\n\n=== Functional Signature Summary ===\n")
  sig_table <- table(profiles$functional_signature)
  for (sig in names(sort(sig_table, decreasing = TRUE))) {
    cat(sprintf("%-20s: %d clusters\n", sig, sig_table[sig]))
  }
  
  cat("\n\n=== Plotting Recommendations ===\n")
  for (rec_name in names(recommendations)) {
    rec <- recommendations[[rec_name]]
    cat("\n", rec$title, ":\n", sep="")
    cat("  ", rec$description, "\n", sep="")
    cat("  Clusters:", paste(rec$clusters, collapse=", "), "\n")
  }
  
  # Create gene lists for specific visualizations
  cat("\n\n=== Key Gene Lists for Visualization ===\n")
  
  # Dopaminergic trajectory genes
  da_genes <- unique(c(
    FUNCTIONAL_CATEGORIES$tf_dopaminergic$genes,
    FUNCTIONAL_CATEGORIES$nt_dopaminergic$genes
  ))
  cat("\nDopaminergic trajectory (", length(da_genes), " genes):\n", sep="")
  cat(paste(da_genes, collapse=", "), "\n")
  
  # Maturation genes
  mat_genes <- unique(c(
    FUNCTIONAL_CATEGORIES$dev_early_progenitor$genes,
    FUNCTIONAL_CATEGORIES$dev_intermediate$genes,
    FUNCTIONAL_CATEGORIES$dev_mature$genes
  ))
  cat("\nMaturation markers (", length(mat_genes), " genes):\n", sep="")
  cat(paste(mat_genes, collapse=", "), "\n")
  
  # Axon guidance genes
  axon_genes <- unique(c(
    FUNCTIONAL_CATEGORIES$axon_guidance_receptors$genes,
    FUNCTIONAL_CATEGORIES$axon_guidance_ligands$genes,
    FUNCTIONAL_CATEGORIES$axon_growth$genes
  ))
  cat("\nAxon guidance (", length(axon_genes), " genes):\n", sep="")
  cat(paste(head(axon_genes, 20), collapse=", "), "...\n")
  
  return(list(
    analysis = analysis_results,
    profiles = profiles,
    novel_markers = novel_markers,
    recommendations = recommendations
  ))
}

# Run if executed directly
if (!interactive()) {
  results <- main()
}