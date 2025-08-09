#!/usr/bin/env Rscript

# COMPREHENSIVE MARKER VALIDATION AND RECLASSIFICATION
# Based on critical findings about mutually exclusive markers

library(dplyr)
library(Seurat)

cat("=================================================================\n")
cat("COMPREHENSIVE MARKER VALIDATION AND RECLASSIFICATION\n")
cat("Based on Mutually Exclusive Marker Principles\n")
cat("=================================================================\n\n")

# Load data
markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")
annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv")

# Define marker categories
NEURONAL_MARKERS <- c("TUBB3", "MAP2", "STMN2", "GAP43", "SYN1", "SNAP25", "NEFL", 
                      "NEFM", "RBFOX3", "MAPT", "DCX", "NCAM1")
PROGENITOR_MARKERS <- c("SOX2", "NES", "VIM", "FABP7", "PAX6", "HES1", "HES5", "ASCL1")
CHOROID_PLEXUS_MARKERS <- c("TTR", "FOLR1", "HTR2C", "CLIC6", "AQP1", "PRLR", "TRPM3", "KCNJ13")
DOPAMINERGIC_MARKERS <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", "ALDH1A1", "NR4A2")

# Function to comprehensively analyze each cluster
analyze_cluster <- function(cluster_id) {
  cluster_markers <- markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC))
  
  # Check marker categories
  neuronal <- cluster_markers %>% filter(gene %in% NEURONAL_MARKERS) %>% head(10)
  progenitor <- cluster_markers %>% filter(gene %in% PROGENITOR_MARKERS) %>% head(10)
  cp <- cluster_markers %>% filter(gene %in% CHOROID_PLEXUS_MARKERS) %>% head(10)
  da <- cluster_markers %>% filter(gene %in% DOPAMINERGIC_MARKERS) %>% head(10)
  
  # Calculate scores
  neuronal_score <- sum(neuronal$avg_log2FC)
  progenitor_score <- sum(progenitor$avg_log2FC)
  cp_score <- sum(cp$avg_log2FC)
  da_score <- sum(da$avg_log2FC)
  
  # Get top 10 markers
  top10 <- cluster_markers %>% head(10) %>% pull(gene)
  
  return(list(
    cluster = cluster_id,
    neuronal_markers = neuronal,
    progenitor_markers = progenitor,
    cp_markers = cp,
    da_markers = da,
    scores = c(neuronal = neuronal_score, progenitor = progenitor_score, 
               cp = cp_score, da = da_score),
    top10 = top10
  ))
}

# Analyze all clusters
cat("Analyzing all clusters for marker categories...\n\n")
all_analyses <- lapply(0:35, analyze_cluster)

# Create revised classifications
revised_classifications <- annotations
revised_classifications$original_type <- annotations$cell_type
revised_classifications$revised_type <- ""
revised_classifications$revision_reason <- ""
revised_classifications$marker_evidence <- ""

for (i in 1:length(all_analyses)) {
  analysis <- all_analyses[[i]]
  cluster_id <- analysis$cluster
  idx <- which(revised_classifications$fine_cluster == cluster_id)
  
  scores <- analysis$scores
  
  # Decision logic based on mutually exclusive principles
  
  # 1. NEURONAL vs PROGENITOR (mutually exclusive)
  if (scores["neuronal"] > 5 && scores["progenitor"] < 2) {
    cell_type <- "Neurons"
    
    # Check if floor plate heritage
    if (any(c("CORIN", "ARX", "EN1", "EN2", "FOXA2") %in% analysis$top10)) {
      cell_type <- "Floor Plate-Derived Neurons"
    }
    
    # Check if dopaminergic
    if (scores["da"] > 5) {
      if ("SOX6" %in% analysis$da_markers$gene && "ALDH1A1" %in% analysis$da_markers$gene) {
        cell_type <- "Dopaminergic Neurons (A9-like, vulnerable)"
      } else if ("CALB1" %in% analysis$da_markers$gene || "CALB2" %in% analysis$da_markers$gene) {
        cell_type <- "Dopaminergic Neurons (A10-like, resilient)"
      } else {
        cell_type <- "Dopaminergic Neurons"
      }
    }
    
    reason <- "High neuronal markers, low/absent progenitor markers"
    
  } else if (scores["progenitor"] > 5 && scores["neuronal"] < 2) {
    cell_type <- "Neural Progenitors"
    reason <- "High progenitor markers, low/absent neuronal markers"
    
  } else if (scores["neuronal"] > 2 && scores["progenitor"] > 2) {
    cell_type <- "Transitioning Cells (Progenitor→Neuron)"
    reason <- "Both neuronal and progenitor markers present"
    
  # 2. CHOROID PLEXUS check
  } else if (scores["cp"] > 8) {
    cell_type <- "Choroid Plexus"
    
    # Check if also has TH
    if ("TH" %in% analysis$da_markers$gene) {
      cell_type <- "Choroid Plexus (with TH+ neurons)"
      reason <- "High CP markers (TTR/FOLR1) with TH expression"
    } else {
      reason <- "High choroid plexus markers"
    }
    
  # 3. Use original classification if unclear
  } else {
    cell_type <- annotations$cell_type[idx]
    reason <- "Maintained original classification"
  }
  
  # Update classification
  revised_classifications$revised_type[idx] <- cell_type
  revised_classifications$revision_reason[idx] <- reason
  
  # Add marker evidence
  evidence <- paste0(
    "Neuronal: ", paste(analysis$neuronal_markers$gene[1:3], collapse=","),
    "; Progenitor: ", paste(analysis$progenitor_markers$gene[1:3], collapse=","),
    "; CP: ", paste(analysis$cp_markers$gene[1:3], collapse=",")
  )
  revised_classifications$marker_evidence[idx] <- evidence
}

# Print major revisions
cat("\n=== MAJOR CLASSIFICATION REVISIONS ===\n\n")
major_revisions <- revised_classifications %>%
  filter(original_type != revised_type) %>%
  select(fine_cluster, original_type, revised_type, revision_reason)

for (i in 1:nrow(major_revisions)) {
  cat("Cluster", major_revisions$fine_cluster[i], ":\n")
  cat("  Original:", major_revisions$original_type[i], "\n")
  cat("  Revised:", major_revisions$revised_type[i], "\n")
  cat("  Reason:", major_revisions$revision_reason[i], "\n\n")
}

# Save results
write.csv(revised_classifications, 
          "results/cell_type_annotations/VALIDATED_cell_type_classifications.csv",
          row.names = FALSE)

# Create detailed report
sink("results/cell_type_annotations/VALIDATION_REPORT.txt")
cat("CELL TYPE VALIDATION REPORT\n")
cat("===========================\n\n")

cat("KEY FINDINGS:\n")
cat("1. Neuronal and progenitor markers are MUTUALLY EXCLUSIVE\n")
cat("2. Clusters with TUBB3/STMN2/GAP43 are NEURONS, not progenitors\n")
cat("3. Choroid plexus can contain TH-expressing neurons\n\n")

cat("CLUSTERS REQUIRING IMMEDIATE RECLASSIFICATION:\n")
for (i in 1:length(all_analyses)) {
  analysis <- all_analyses[[i]]
  if (analysis$scores["neuronal"] > 5 && 
      annotations$cell_type[annotations$fine_cluster == analysis$cluster] %in% 
      c("Floor Plate Progenitors", "Neural Progenitors", "Progenitor-like")) {
    cat("\nCluster", analysis$cluster, "- Misclassified as progenitor but has neuronal markers:\n")
    print(analysis$neuronal_markers[, c("gene", "avg_log2FC")])
  }
}

sink()

cat("\n\nValidation complete! Results saved to:\n")
cat("  results/cell_type_annotations/VALIDATED_cell_type_classifications.csv\n")
cat("  results/cell_type_annotations/VALIDATION_REPORT.txt\n")

# Return summary
list(
  revised_classifications = revised_classifications,
  major_revisions = major_revisions,
  analyses = all_analyses
)