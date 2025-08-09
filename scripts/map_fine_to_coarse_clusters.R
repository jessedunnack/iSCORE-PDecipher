#!/usr/bin/env Rscript

# Map Fine Clusters (36) to Coarse Clusters (15)
# Based on hierarchical relationships and cell type annotations

library(dplyr)

# Function to create hypothetical mapping based on cell type similarities
create_hypothetical_mapping <- function() {
  # This is a hypothetical mapping based on:
  # 1. Cell type similarities
  # 2. Assumption that fine clusters are subdivisions of coarse
  # 3. Roughly 2-3 fine clusters per coarse cluster
  
  fine_to_coarse <- list(
    # Dopaminergic neurons group
    "5" = 0,    # DA general
    "18" = 0,   # DA A10
    "28" = 0,   # DA A9
    
    # Immature neurons/neuroblasts
    "0" = 1,    # Neuroblasts
    
    # Floor plate/progenitors
    "1" = 2,    # FP CORIN+
    "3" = 2,    # FP EN2+
    "13" = 2,   # FP SHH+
    
    # Neural progenitors
    "8" = 3,    # PTN+ progenitor
    "10" = 3,   # MSX1+ progenitor
    "26" = 3,   # NES+ progenitor
    
    # Choroid plexus/epithelial
    "4" = 4,    # Choroid plexus
    
    # Unknown neural 1
    "11" = 5,   # PRSS56+
    "12" = 5,   # CRABP1+
    "14" = 5,   # Melanocyte-like
    
    # Stressed/dying cells
    "6" = 6,    # Mitochondrial stress
    "24" = 6,   # ER stress
    
    # Mesenchymal/fibroblasts
    "7" = 7,    # Collagen-producing
    "21" = 7,   # ECM-producing
    "33" = 7,   # CRYAB+
    
    # Unknown neural 2
    "16" = 8,   # RGCC+ progenitor
    "17" = 8,   # GRIA2+
    "23" = 8,   # LINC02539+
    
    # Non-neural support
    "15" = 9,   # Leptomeningeal
    "27" = 9,   # Vascular smooth muscle
    
    # Proliferating
    "22" = 10,  # G2/M phase
    "35" = 10,  # G2/M phase
    
    # Unknown mixed
    "2" = 11,   # Mixed/transitional
    "9" = 11,   # OTX2+
    "25" = 11,  # Metabolic
    
    # Hypothalamic/diencephalic
    "30" = 12,  # POU3F2+ hypothalamic
    "34" = 12,  # H19+/SHH+
    
    # Caudal/spinal
    "29" = 13,  # HOXD10/11+
    "31" = 13,  # RBP4+
    
    # Technical/other
    "19" = 14,  # Oligodendrocytes
    "20" = 14,  # lncRNA artifact
    "32" = 14   # Ciliated
  )
  
  return(fine_to_coarse)
}

# Function to create mapping table with cell type info
create_mapping_table <- function(fine_to_coarse, annotations) {
  mapping_df <- data.frame(
    fine_cluster = integer(),
    coarse_cluster = integer(),
    fine_cell_type = character(),
    fine_subtype = character(),
    confidence = character(),
    stringsAsFactors = FALSE
  )
  
  for (fine in names(fine_to_coarse)) {
    fine_num <- as.integer(fine)
    coarse_num <- fine_to_coarse[[fine]]
    
    # Get annotation info
    ann_row <- annotations[annotations$fine_cluster == fine_num, ]
    
    mapping_df <- rbind(mapping_df, data.frame(
      fine_cluster = fine_num,
      coarse_cluster = coarse_num,
      fine_cell_type = ifelse(nrow(ann_row) > 0, ann_row$cell_type, "Unknown"),
      fine_subtype = ifelse(nrow(ann_row) > 0, ann_row$subtype, ""),
      confidence = ifelse(nrow(ann_row) > 0, ann_row$confidence, "Low"),
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by fine cluster
  mapping_df <- mapping_df[order(mapping_df$fine_cluster), ]
  
  return(mapping_df)
}

# Function to summarize coarse cluster composition
summarize_coarse_clusters <- function(mapping_df) {
  coarse_summary <- mapping_df %>%
    group_by(coarse_cluster) %>%
    summarise(
      n_fine_clusters = n(),
      fine_clusters = paste(fine_cluster, collapse = ", "),
      cell_types = paste(unique(fine_cell_type), collapse = " / "),
      avg_confidence = mean(confidence == "High"),
      .groups = 'drop'
    ) %>%
    arrange(coarse_cluster)
  
  # Add interpretations
  coarse_interpretations <- c(
    "0" = "Dopaminergic Neurons",
    "1" = "Immature Neurons",
    "2" = "Floor Plate Progenitors", 
    "3" = "Neural Progenitors",
    "4" = "Choroid Plexus",
    "5" = "Unknown Neural 1",
    "6" = "Stressed/Dying Cells",
    "7" = "Mesenchymal/Fibroblasts",
    "8" = "Unknown Neural 2",
    "9" = "Non-neural Support",
    "10" = "Proliferating Cells",
    "11" = "Mixed/Unknown",
    "12" = "Hypothalamic/Diencephalic",
    "13" = "Caudal/Spinal",
    "14" = "Technical/Other"
  )
  
  coarse_summary$interpretation <- coarse_interpretations[as.character(coarse_summary$coarse_cluster)]
  
  return(coarse_summary)
}

# Function to create visual mapping
plot_mapping_summary <- function(mapping_df, output_file = NULL) {
  cat("\n=== Fine to Coarse Cluster Mapping ===\n\n")
  
  for (coarse in 0:14) {
    fine_in_coarse <- mapping_df[mapping_df$coarse_cluster == coarse, ]
    
    if (nrow(fine_in_coarse) > 0) {
      cat(sprintf("Coarse Cluster %d:\n", coarse))
      
      for (i in 1:nrow(fine_in_coarse)) {
        cat(sprintf("  └─ Fine %2d: %-30s %-20s [%s]\n",
                    fine_in_coarse$fine_cluster[i],
                    fine_in_coarse$fine_cell_type[i],
                    fine_in_coarse$fine_subtype[i],
                    fine_in_coarse$confidence[i]))
      }
      cat("\n")
    }
  }
}

# Main function
main <- function() {
  cat("Fine to Coarse Cluster Mapping Analysis\n")
  cat("=======================================\n\n")
  
  # Load annotations
  annotations_file <- "results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv"
  if (!file.exists(annotations_file)) {
    stop("Please run comprehensive_cell_type_annotation.R first")
  }
  
  annotations <- read.csv(annotations_file, stringsAsFactors = FALSE)
  
  # Create mapping
  cat("Creating hypothetical fine-to-coarse mapping...\n")
  fine_to_coarse <- create_hypothetical_mapping()
  
  # Create mapping table
  mapping_df <- create_mapping_table(fine_to_coarse, annotations)
  
  # Save mapping
  write.csv(mapping_df, "results/cell_type_annotations/fine_to_coarse_mapping.csv",
            row.names = FALSE)
  
  # Create coarse cluster summary
  coarse_summary <- summarize_coarse_clusters(mapping_df)
  write.csv(coarse_summary, "results/cell_type_annotations/coarse_cluster_summary.csv",
            row.names = FALSE)
  
  # Display mapping
  plot_mapping_summary(mapping_df)
  
  # Print coarse cluster summary
  cat("\n=== Coarse Cluster Summary ===\n")
  print(coarse_summary[, c("coarse_cluster", "interpretation", "n_fine_clusters", "cell_types")])
  
  cat("\n=== Important Notes ===\n")
  cat("1. This mapping is HYPOTHETICAL based on cell type similarities\n")
  cat("2. The actual mapping requires access to the Seurat object metadata\n")
  cat("3. To get the true mapping, run:\n")
  cat("   seurat_obj <- readRDS('../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds')\n")
  cat("   # Extract both clustering resolutions and create mapping\n")
  cat("4. This hypothetical mapping groups similar cell types together\n")
  
  return(list(
    mapping = mapping_df,
    coarse_summary = coarse_summary
  ))
}

# Function to extract true mapping from Seurat object (template)
extract_true_mapping <- function(seurat_path) {
  cat("\n=== Template for Extracting True Mapping ===\n")
  cat("# Load Seurat object\n")
  cat("seurat_obj <- readRDS('", seurat_path, "')\n\n")
  cat("# Check available metadata\n")
  cat("colnames(seurat_obj@meta.data)\n\n")
  cat("# Look for clustering columns (e.g., seurat_clusters, RNA_snn_res.0.X)\n")
  cat("# Assuming coarse is at lower resolution and fine at higher\n\n")
  cat("# Extract both clusterings\n")
  cat("coarse_clusters <- seurat_obj$seurat_clusters_coarse  # or appropriate column\n")
  cat("fine_clusters <- seurat_obj$seurat_clusters_fine      # or appropriate column\n\n")
  cat("# Create mapping\n")
  cat("mapping_true <- data.frame(\n")
  cat("  cell_barcode = colnames(seurat_obj),\n")
  cat("  coarse = coarse_clusters,\n")
  cat("  fine = fine_clusters\n")
  cat(")\n\n")
  cat("# Summarize relationships\n")
  cat("fine_to_coarse_true <- mapping_true %>%\n")
  cat("  group_by(fine) %>%\n")
  cat("  summarise(coarse = names(sort(table(coarse), decreasing = TRUE))[1])\n")
}

# Run if executed directly
if (!interactive()) {
  results <- main()
}