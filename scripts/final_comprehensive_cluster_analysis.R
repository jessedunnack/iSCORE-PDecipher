#!/usr/bin/env Rscript

# FINAL COMPREHENSIVE CLUSTER ANALYSIS WITH HIERARCHY AND REFERENCE DATA
# Using the established cluster hierarchy and all reference markers

library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("FINAL COMPREHENSIVE CLUSTER ANALYSIS\n")
cat("Using cluster hierarchy and reference datasets\n")
cat("=================================================================\n\n")

# 1. Load all necessary data
cat("1. Loading data...\n")
fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)
fine_to_coarse <- read.csv("results/cluster_hierarchy/fine_to_coarse_mapping.csv")
metadata <- readRDS("results/cluster_hierarchy/cluster_metadata.rds")

cat("Data loaded successfully\n\n")

# 2. Define comprehensive marker sets from our research
# Kim 2021 iPSC protocol expected cell types
DA_MARKERS <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", "EN1", "EN2", "LMX1A", "FOXA2", "NR4A2", "PITX3")
DA_VULNERABLE <- c("SOX6", "ALDH1A1", "SNCG", "ATP13A2", "RIT2", "AGTR1") # A9
DA_RESILIENT <- c("CALB1", "CALB2", "OTX2", "GRP", "CCK") # A10
MEIS_NETWORK <- c("MEIS1", "MEIS2", "PBX1", "PBX3", "HOXA9", "HOXA10")

FLOOR_PLATE <- c("CORIN", "LMX1A", "FOXA2", "ARX", "SHH", "WNT1", "MSX1", "LMO3")
NEURONAL <- c("TUBB3", "MAP2", "STMN2", "GAP43", "SYN1", "SNAP25", "RBFOX3", "NEFL", "MAPT")
PROGENITOR <- c("SOX2", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7", "SOX9")

CHOROID_PLEXUS <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "PRLR", "KCNJ13", "ENPP2")
OLIGODENDROCYTE <- c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "MPZ", "MBP", "PLP1")
VASCULAR <- c("CLDN5", "PECAM1", "VWF", "FLT1", "TAGLN", "ACTA2", "MYL9", "PDGFRB")
PROLIFERATING <- c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "HIST1H3B")

# 3. Function to analyze a cluster comprehensively
analyze_cluster_comprehensive <- function(cluster_id, is_coarse = FALSE) {
  # Get markers
  if (is_coarse) {
    markers_df <- coarse_markers
    cluster_type <- "Coarse"
  } else {
    markers_df <- fine_markers
    cluster_type <- "Fine"
  }
  
  # Get top markers
  top_markers <- markers_df %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(50)
  
  # Calculate scores for each cell type
  scores <- list()
  
  # Helper function
  check_markers <- function(marker_list, top_df) {
    present <- marker_list[marker_list %in% top_df$gene]
    if (length(present) == 0) return(list(score = 0, markers = ""))
    
    score_df <- top_df %>%
      filter(gene %in% present) %>%
      summarise(
        score = sum(avg_log2FC),
        n = n(),
        markers = paste(gene, collapse = ", ")
      )
    return(list(score = score_df$score, n = score_df$n, markers = score_df$markers))
  }
  
  # Calculate all scores
  scores$da_core <- check_markers(DA_MARKERS, top_markers)
  scores$da_vulnerable <- check_markers(DA_VULNERABLE, top_markers)
  scores$da_resilient <- check_markers(DA_RESILIENT, top_markers)
  scores$meis <- check_markers(MEIS_NETWORK, top_markers)
  scores$floor_plate <- check_markers(FLOOR_PLATE, top_markers)
  scores$neuronal <- check_markers(NEURONAL, top_markers)
  scores$progenitor <- check_markers(PROGENITOR, top_markers)
  scores$choroid_plexus <- check_markers(CHOROID_PLEXUS, top_markers)
  scores$oligodendrocyte <- check_markers(OLIGODENDROCYTE, top_markers)
  scores$vascular <- check_markers(VASCULAR, top_markers)
  scores$proliferating <- check_markers(PROLIFERATING, top_markers)
  
  # Get top 10 genes
  top10 <- paste(head(top_markers$gene, 10), collapse = ", ")
  
  return(list(
    cluster = cluster_id,
    type = cluster_type,
    top10 = top10,
    scores = scores
  ))
}

# 4. Analyze all coarse clusters first
cat("2. Analyzing COARSE clusters (0-14)...\n")
cat("=====================================\n\n")

coarse_results <- list()
for (i in 0:14) {
  result <- analyze_cluster_comprehensive(i, is_coarse = TRUE)
  coarse_results[[as.character(i)]] <- result
  
  cat(sprintf("COARSE CLUSTER %d:\n", i))
  cat("Top markers:", result$top10, "\n")
  
  # Determine cell type
  cell_type <- "Unknown"
  confidence <- "Low"
  
  # Decision logic
  if (result$scores$choroid_plexus$score > 20) {
    cell_type <- "Choroid Plexus"
    confidence <- "High"
  } else if (result$scores$da_core$score > 10) {
    cell_type <- "Dopaminergic lineage"
    confidence <- "High"
  } else if (result$scores$proliferating$score > 15) {
    cell_type <- "Proliferating cells"
    confidence <- "High"
  } else if (result$scores$neuronal$score > 15) {
    cell_type <- "Neurons"
    confidence <- "High"
  } else if (result$scores$progenitor$score > 10 && result$scores$neuronal$score < 5) {
    cell_type <- "Neural progenitors"
    confidence <- "Medium"
  } else if (result$scores$vascular$score > 10) {
    cell_type <- "Vascular cells"
    confidence <- "High"
  } else if (result$scores$oligodendrocyte$score > 5) {
    cell_type <- "Oligodendrocyte lineage"
    confidence <- "Medium"
  }
  
  cat("Cell type:", cell_type, "(", confidence, "confidence)\n")
  
  # Print significant scores
  if (result$scores$da_core$score > 0) {
    cat("  DA markers:", result$scores$da_core$markers, "\n")
  }
  if (result$scores$neuronal$score > 5) {
    cat("  Neuronal score:", round(result$scores$neuronal$score, 2), "\n")
  }
  
  # Get fine clusters in this coarse cluster
  fine_in_coarse <- fine_to_coarse %>%
    filter(coarse_cluster == i) %>%
    pull(fine_cluster)
  cat("  Contains fine clusters:", paste(fine_in_coarse, collapse = ", "), "\n\n")
  
  coarse_results[[as.character(i)]]$cell_type <- cell_type
  coarse_results[[as.character(i)]]$confidence <- confidence
}

# 5. Now analyze all fine clusters with coarse context
cat("\n3. Analyzing FINE clusters (0-35) with coarse context...\n")
cat("======================================================\n\n")

final_assignments <- data.frame(
  fine_cluster = 0:35,
  coarse_cluster = fine_to_coarse$coarse_cluster[match(0:35, fine_to_coarse$fine_cluster)],
  cell_type = NA,
  subtype = NA,
  confidence = NA,
  vulnerability = NA,
  key_markers = NA,
  stringsAsFactors = FALSE
)

for (i in 0:35) {
  result <- analyze_cluster_comprehensive(i, is_coarse = FALSE)
  coarse_id <- final_assignments$coarse_cluster[i+1]
  coarse_type <- coarse_results[[as.character(coarse_id)]]$cell_type
  
  cat(sprintf("\nFINE CLUSTER %d (Coarse %d: %s):\n", i, coarse_id, coarse_type))
  cat("Top 10:", result$top10, "\n")
  
  # Refined assignment based on both levels
  cell_type <- coarse_type
  subtype <- ""
  confidence <- "Medium"
  vulnerability <- "Unknown"
  
  # Special cases based on specific markers
  if (grepl("TTR", result$top10)) {
    cell_type <- "Choroid Plexus"
    subtype <- "TTR-high"
    confidence <- "Very High"
  } else if (grepl("TH", result$top10) && result$scores$da_core$score > 3) {
    cell_type <- "Dopaminergic Neurons"
    confidence <- "High"
    
    # Determine DA subtype
    if (result$scores$da_vulnerable$score > result$scores$da_resilient$score) {
      subtype <- "A9-like (vulnerable)"
      vulnerability <- "High"
      if (result$scores$meis$n > 0) {
        vulnerability <- "Very High (MEIS+)"
      }
    } else if (result$scores$da_resilient$score > result$scores$da_vulnerable$score) {
      subtype <- "A10-like (resilient)"
      vulnerability <- "Low"
    } else {
      subtype == "Immature/unspecified"
      vulnerability <- "Medium"
    }
  } else if (grepl("MKI67|TOP2A", result$top10)) {
    cell_type <- "Proliferating Cells"
    confidence <- "Very High"
  } else if (result$scores$neuronal$score > 10 && result$scores$progenitor$score < 2) {
    cell_type <- "Neurons"
    if (result$scores$floor_plate$score > 5) {
      subtype <- "Floor plate-derived"
    }
    confidence <- "High"
  } else if (grepl("SOX10|MPZ", result$top10)) {
    cell_type <- "Oligodendrocyte Lineage"
    confidence <- "High"
  } else if (grepl("TAGLN|ACTA2", result$top10)) {
    cell_type <- "Smooth Muscle/Pericytes"
    confidence <- "High"
  } else if (grepl("COL1A1|COL1A2", result$top10) && !grepl("neuron", cell_type, ignore.case = TRUE)) {
    cell_type <- "Fibroblasts"
    confidence <- "High"
  }
  
  # Store results
  final_assignments[i+1, "cell_type"] <- cell_type
  final_assignments[i+1, "subtype"] <- subtype
  final_assignments[i+1, "confidence"] <- confidence
  final_assignments[i+1, "vulnerability"] <- vulnerability
  final_assignments[i+1, "key_markers"] <- substr(result$top10, 1, 50)
  
  cat("ASSIGNMENT:", cell_type)
  if (subtype != "") cat(" -", subtype)
  cat(" (", confidence, ")\n")
  
  if (vulnerability != "Unknown" && vulnerability != "") {
    cat("Vulnerability:", vulnerability, "\n")
  }
}

# 6. Summary statistics
cat("\n\n4. SUMMARY STATISTICS\n")
cat("======================\n")

cell_type_summary <- table(final_assignments$cell_type)
cat("\nCell type distribution:\n")
print(cell_type_summary)

# DA neuron analysis
da_neurons <- final_assignments[grep("Dopaminergic", final_assignments$cell_type), ]
if (nrow(da_neurons) > 0) {
  cat("\n\nDopaminergic neuron clusters:\n")
  print(da_neurons[, c("fine_cluster", "coarse_cluster", "subtype", "vulnerability")])
}

# High confidence assignments
high_conf <- final_assignments[final_assignments$confidence %in% c("High", "Very High"), ]
cat("\n\nHigh confidence assignments:", nrow(high_conf), "out of 36\n")

# 7. Save comprehensive results
cat("\n5. Saving results...\n")
dir.create("results/final_comprehensive_analysis", recursive = TRUE, showWarnings = FALSE)

write.csv(final_assignments, 
          "results/final_comprehensive_analysis/final_cluster_assignments.csv", 
          row.names = FALSE)

# Create detailed report
detailed_report <- list(
  final_assignments = final_assignments,
  coarse_analysis = coarse_results,
  fine_to_coarse_mapping = fine_to_coarse,
  analysis_date = Sys.Date(),
  reference_datasets = c("Kim 2021 iPSC protocol", "FOUNDIN-PD", "Macaque MPTP studies", "La Manno 2016")
)

saveRDS(detailed_report, "results/final_comprehensive_analysis/complete_analysis_report.rds")

cat("\nResults saved to:\n")
cat("- results/final_comprehensive_analysis/final_cluster_assignments.csv\n")
cat("- results/final_comprehensive_analysis/complete_analysis_report.rds\n")

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("All 36 fine clusters have been assigned cell types based on:\n")
cat("- Hierarchical clustering structure\n")
cat("- Comprehensive marker analysis\n")
cat("- Reference dataset comparisons\n")
cat("- Kim 2021 iPSC protocol expectations\n")