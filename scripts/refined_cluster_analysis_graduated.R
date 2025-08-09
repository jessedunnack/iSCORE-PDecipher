#!/usr/bin/env Rscript

# REFINED CLUSTER ANALYSIS WITH GRADUATED SCORING
# More nuanced approach to identify all cell types, especially dopaminergic neurons

library(dplyr)
library(tidyr)
library(ggplot2)
# library(pheatmap) # Commented out if not available

cat("=================================================================\n")
cat("REFINED CLUSTER ANALYSIS WITH GRADUATED SCORING\n")
cat("Investigating low dopaminergic neuron detection\n")
cat("=================================================================\n\n")

# 1. Load data
cat("1. Loading data...\n")
fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)
fine_to_coarse <- read.csv("results/cluster_hierarchy/fine_to_coarse_mapping.csv")

# 2. Define comprehensive marker sets with importance weights
MARKER_SETS <- list(
  # Dopaminergic markers with weights
  DA_ESSENTIAL = list(
    markers = c("TH", "DDC"),
    weight = 3.0,
    description = "Essential DA markers"
  ),
  DA_MATURE = list(
    markers = c("SLC6A3", "SLC18A2", "DRD2", "KCNJ6"),
    weight = 2.0,
    description = "Mature DA markers"
  ),
  DA_TRANSCRIPTION = list(
    markers = c("LMX1A", "FOXA2", "NR4A2", "PITX3", "EN1", "EN2"),
    weight = 1.5,
    description = "DA transcription factors"
  ),
  DA_VULNERABLE = list(
    markers = c("SOX6", "ALDH1A1", "SNCG", "ATP13A2", "RIT2", "AGTR1"),
    weight = 1.0,
    description = "A9 vulnerability markers"
  ),
  DA_RESILIENT = list(
    markers = c("CALB1", "CALB2", "OTX2", "GRP", "CCK"),
    weight = 1.0,
    description = "A10 resilience markers"
  ),
  MEIS_NETWORK = list(
    markers = c("MEIS1", "MEIS2", "PBX1", "PBX3", "HOXA9", "HOXA10"),
    weight = 1.5,
    description = "MEIS vulnerability network"
  ),
  
  # Other neuronal types
  NEURONAL_GENERAL = list(
    markers = c("TUBB3", "MAP2", "STMN2", "GAP43", "SYN1", "SNAP25", "RBFOX3"),
    weight = 2.0,
    description = "General neuronal markers"
  ),
  GLUTAMATERGIC = list(
    markers = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "GRIA1"),
    weight = 2.0,
    description = "Glutamatergic markers"
  ),
  GABAERGIC = list(
    markers = c("GAD1", "GAD2", "SLC32A1", "GABBR1", "GABRG2"),
    weight = 2.0,
    description = "GABAergic markers"
  ),
  SEROTONERGIC = list(
    markers = c("TPH2", "SLC6A4", "FEV", "GATA2", "GATA3"),
    weight = 2.0,
    description = "Serotonergic markers"
  ),
  
  # Progenitor and non-neuronal
  FLOOR_PLATE = list(
    markers = c("CORIN", "LMX1A", "FOXA2", "ARX", "SHH", "WNT1", "MSX1", "LMO3"),
    weight = 2.0,
    description = "Floor plate markers"
  ),
  PROGENITOR = list(
    markers = c("SOX2", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7"),
    weight = 2.0,
    description = "Neural progenitor markers"
  ),
  RADIAL_GLIA = list(
    markers = c("SOX9", "HOPX", "GFAP", "AQP4", "TNC", "BLBP"),
    weight = 1.5,
    description = "Radial glia markers"
  ),
  CHOROID_PLEXUS = list(
    markers = c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "PRLR", "KCNJ13"),
    weight = 3.0,
    description = "Choroid plexus markers"
  ),
  OLIGODENDROCYTE = list(
    markers = c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "MPZ", "MBP", "PLP1"),
    weight = 2.0,
    description = "Oligodendrocyte markers"
  ),
  VASCULAR = list(
    markers = c("CLDN5", "PECAM1", "VWF", "FLT1", "TAGLN", "ACTA2", "MYL9"),
    weight = 2.0,
    description = "Vascular markers"
  ),
  PROLIFERATING = list(
    markers = c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "HIST1H3B"),
    weight = 3.0,
    description = "Cell cycle markers"
  ),
  STRESS_RESPONSE = list(
    markers = c("FOS", "JUN", "EGR1", "ATF3", "HSPA1A", "DDIT3"),
    weight = 1.0,
    description = "Stress response markers"
  )
)

# 3. Function to calculate graduated scores
calculate_graduated_scores <- function(cluster_id, markers_df, top_n = 100) {
  # Get extended marker list
  cluster_markers <- markers_df %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(top_n)
  
  scores <- list()
  
  for (set_name in names(MARKER_SETS)) {
    set_info <- MARKER_SETS[[set_name]]
    markers_present <- cluster_markers %>%
      filter(gene %in% set_info$markers)
    
    if (nrow(markers_present) == 0) {
      scores[[set_name]] <- list(
        score = 0,
        n_markers = 0,
        pct_of_set = 0,
        markers = "",
        avg_fc = 0,
        weighted_score = 0
      )
    } else {
      scores[[set_name]] <- list(
        score = sum(markers_present$avg_log2FC),
        n_markers = nrow(markers_present),
        pct_of_set = round(100 * nrow(markers_present) / length(set_info$markers), 1),
        markers = paste(markers_present$gene, collapse = ", "),
        avg_fc = mean(markers_present$avg_log2FC),
        weighted_score = sum(markers_present$avg_log2FC) * set_info$weight
      )
    }
  }
  
  return(scores)
}

# 4. Analyze each cluster with graduated scoring
cat("\n2. Analyzing all 36 fine clusters with graduated scoring...\n")
cat("=========================================================\n\n")

all_results <- list()

for (i in 0:35) {
  cat(sprintf("\n--- FINE CLUSTER %d ---\n", i))
  
  # Get coarse cluster
  coarse_id <- fine_to_coarse$coarse_cluster[fine_to_coarse$fine_cluster == i]
  
  # Calculate scores
  scores <- calculate_graduated_scores(i, fine_markers, top_n = 100)
  
  # Get top 10 markers for display
  top10 <- fine_markers %>%
    filter(cluster == as.character(i)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(10) %>%
    pull(gene)
  
  cat("Coarse cluster:", coarse_id, "\n")
  cat("Top 10 markers:", paste(top10, collapse = ", "), "\n\n")
  
  # Display significant scores
  cat("Marker set scores (weighted):\n")
  score_summary <- data.frame(
    Set = character(),
    WeightedScore = numeric(),
    NMarkers = integer(),
    PctCoverage = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (set_name in names(scores)) {
    if (scores[[set_name]]$weighted_score > 0.5) {
      cat(sprintf("  %s: %.2f (n=%d, %.1f%% coverage)\n",
                  set_name,
                  scores[[set_name]]$weighted_score,
                  scores[[set_name]]$n_markers,
                  scores[[set_name]]$pct_of_set))
      
      score_summary <- rbind(score_summary, data.frame(
        Set = set_name,
        WeightedScore = scores[[set_name]]$weighted_score,
        NMarkers = scores[[set_name]]$n_markers,
        PctCoverage = scores[[set_name]]$pct_of_set
      ))
    }
  }
  
  # Store results
  all_results[[as.character(i)]] <- list(
    cluster = i,
    coarse_cluster = coarse_id,
    top10 = paste(top10, collapse = ", "),
    scores = scores,
    score_summary = score_summary
  )
}

# 5. Create score matrix for visualization
cat("\n\n3. Creating score matrix for visualization...\n")

# Extract weighted scores into matrix
score_matrix <- matrix(0, nrow = 36, ncol = length(MARKER_SETS))
rownames(score_matrix) <- paste0("Cluster_", 0:35)
colnames(score_matrix) <- names(MARKER_SETS)

for (i in 0:35) {
  for (j in 1:length(MARKER_SETS)) {
    set_name <- names(MARKER_SETS)[j]
    score_matrix[i+1, j] <- all_results[[as.character(i)]]$scores[[set_name]]$weighted_score
  }
}

# 6. Identify dopaminergic clusters more sensitively
cat("\n4. Investigating dopaminergic neuron detection...\n")
cat("================================================\n\n")

# Calculate DA scores combining multiple marker sets
da_combined_scores <- data.frame(
  cluster = 0:35,
  da_essential = numeric(36),
  da_mature = numeric(36),
  da_tf = numeric(36),
  da_total = numeric(36),
  neuronal = numeric(36),
  progenitor = numeric(36),
  has_TH = logical(36),
  has_DDC = logical(36),
  coarse_cluster = fine_to_coarse$coarse_cluster[match(0:35, fine_to_coarse$fine_cluster)]
)

for (i in 0:35) {
  scores <- all_results[[as.character(i)]]$scores
  da_combined_scores[i+1, "da_essential"] <- scores$DA_ESSENTIAL$weighted_score
  da_combined_scores[i+1, "da_mature"] <- scores$DA_MATURE$weighted_score
  da_combined_scores[i+1, "da_tf"] <- scores$DA_TRANSCRIPTION$weighted_score
  da_combined_scores[i+1, "da_total"] <- 
    scores$DA_ESSENTIAL$weighted_score + 
    scores$DA_MATURE$weighted_score + 
    scores$DA_TRANSCRIPTION$weighted_score
  da_combined_scores[i+1, "neuronal"] <- scores$NEURONAL_GENERAL$weighted_score
  da_combined_scores[i+1, "progenitor"] <- scores$PROGENITOR$weighted_score
  
  # Check for specific markers
  cluster_markers <- fine_markers %>%
    filter(cluster == as.character(i)) %>%
    head(100)
  da_combined_scores[i+1, "has_TH"] <- "TH" %in% cluster_markers$gene
  da_combined_scores[i+1, "has_DDC"] <- "DDC" %in% cluster_markers$gene
}

# Identify potential DA clusters with relaxed criteria
potential_da <- da_combined_scores %>%
  filter(da_total > 1 | has_TH | has_DDC) %>%
  arrange(desc(da_total))

cat("Potential dopaminergic clusters:\n")
print(potential_da)

cat("\n\nClusters with ANY dopaminergic markers:\n")
any_da_markers <- da_combined_scores %>%
  filter(da_essential > 0 | da_mature > 0 | da_tf > 0) %>%
  arrange(desc(da_total))
print(any_da_markers[, c("cluster", "da_essential", "da_mature", "da_tf", "da_total", "coarse_cluster")])

# 7. Create comprehensive classification
cat("\n\n5. Creating refined cluster classifications...\n")
cat("===========================================\n")

refined_assignments <- data.frame(
  fine_cluster = 0:35,
  coarse_cluster = fine_to_coarse$coarse_cluster[match(0:35, fine_to_coarse$fine_cluster)],
  cell_type = NA,
  subtype = NA,
  confidence = NA,
  da_score = da_combined_scores$da_total,
  neuronal_score = da_combined_scores$neuronal,
  key_markers = NA,
  stringsAsFactors = FALSE
)

# Classify each cluster
for (i in 0:35) {
  scores <- all_results[[as.character(i)]]$scores
  da_score <- da_combined_scores$da_total[i+1]
  
  # Default values
  cell_type <- "Unknown"
  subtype <- ""
  confidence <- "Low"
  
  # Hierarchical classification logic
  if (scores$CHOROID_PLEXUS$weighted_score > 10) {
    cell_type <- "Choroid Plexus"
    confidence <- "Very High"
  } else if (scores$PROLIFERATING$weighted_score > 10) {
    cell_type <- "Proliferating Cells"
    confidence <- "Very High"
  } else if (da_score > 5 || da_combined_scores$has_TH[i+1]) {
    cell_type <- "Dopaminergic Neurons"
    confidence <- ifelse(da_score > 10, "Very High", "High")
    
    if (scores$DA_VULNERABLE$weighted_score > scores$DA_RESILIENT$weighted_score) {
      subtype <- "A9-like (vulnerable)"
    } else if (scores$DA_RESILIENT$weighted_score > 0) {
      subtype <- "A10-like (resilient)"
    }
  } else if (da_score > 2 && scores$NEURONAL_GENERAL$weighted_score > 5) {
    cell_type <- "Dopaminergic Lineage"
    subtype <- "Immature/Precursor"
    confidence <- "Medium"
  } else if (scores$NEURONAL_GENERAL$weighted_score > 10) {
    cell_type <- "Neurons"
    confidence <- "High"
    
    if (scores$GLUTAMATERGIC$weighted_score > 5) {
      subtype <- "Glutamatergic"
    } else if (scores$GABAERGIC$weighted_score > 5) {
      subtype <- "GABAergic"
    } else if (scores$SEROTONERGIC$weighted_score > 3) {
      subtype <- "Serotonergic"
    }
  } else if (scores$FLOOR_PLATE$weighted_score > 5 && scores$NEURONAL_GENERAL$weighted_score > 3) {
    cell_type <- "Floor Plate-Derived Neurons"
    confidence <- "High"
  } else if (scores$PROGENITOR$weighted_score > 8 && scores$NEURONAL_GENERAL$weighted_score < 3) {
    cell_type <- "Neural Progenitors"
    confidence <- "High"
  } else if (scores$OLIGODENDROCYTE$weighted_score > 5) {
    cell_type <- "Oligodendrocyte Lineage"
    confidence <- "High"
  } else if (scores$VASCULAR$weighted_score > 8) {
    cell_type <- "Vascular Cells"
    confidence <- "High"
  } else if (scores$NEURONAL_GENERAL$weighted_score > 5) {
    cell_type <- "Neurons"
    subtype <- "Unspecified"
    confidence <- "Medium"
  } else if (scores$STRESS_RESPONSE$weighted_score > 5) {
    cell_type <- "Stressed/Damaged Cells"
    confidence <- "Medium"
  }
  
  # Fill in assignments
  refined_assignments[i+1, "cell_type"] <- cell_type
  refined_assignments[i+1, "subtype"] <- subtype
  refined_assignments[i+1, "confidence"] <- confidence
  refined_assignments[i+1, "key_markers"] <- substr(all_results[[as.character(i)]]$top10, 1, 80)
}

# 8. Summary statistics
cat("\n\n6. SUMMARY STATISTICS\n")
cat("=====================\n")

cell_type_summary <- table(refined_assignments$cell_type)
cat("\nCell type distribution:\n")
print(cell_type_summary)

cat("\n\nDopaminergic-related clusters:\n")
da_related <- refined_assignments %>%
  filter(grepl("Dopaminergic", cell_type))
print(da_related[, c("fine_cluster", "coarse_cluster", "cell_type", "subtype", "confidence", "da_score")])

cat("\n\nConfidence distribution:\n")
print(table(refined_assignments$confidence))

# 9. Create visualizations
cat("\n7. Creating visualizations...\n")

# Heatmap of marker scores
# Commented out pheatmap code - requires pheatmap package
# pdf("results/refined_analysis/marker_score_heatmap.pdf", width = 12, height = 10)
# pheatmap(
#   t(score_matrix),
#   cluster_cols = TRUE,
#   cluster_rows = FALSE,
#   scale = "row",
#   color = colorRampPalette(c("white", "yellow", "orange", "red", "darkred"))(100),
#   main = "Marker Set Expression Scores Across Clusters",
#   fontsize = 8,
#   cellwidth = 15,
#   cellheight = 12,
#   annotation_col = data.frame(
#     CellType = refined_assignments$cell_type,
#     row.names = rownames(score_matrix)
#   )
# )
# dev.off()

# DA score plot
p_da <- ggplot(da_combined_scores, aes(x = factor(cluster), y = da_total)) +
  geom_bar(stat = "identity", aes(fill = da_total > 2)) +
  geom_hline(yintercept = c(2, 5, 10), linetype = "dashed", alpha = 0.5) +
  scale_fill_manual(values = c("gray80", "darkred"), guide = FALSE) +
  labs(x = "Fine Cluster", y = "Combined DA Score", 
       title = "Dopaminergic Marker Scores Across Clusters") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/refined_analysis/da_scores_by_cluster.pdf", p_da, width = 10, height = 6)

# 10. Save results
cat("\n8. Saving refined results...\n")
dir.create("results/refined_analysis", recursive = TRUE, showWarnings = FALSE)

write.csv(refined_assignments, 
          "results/refined_analysis/refined_cluster_assignments.csv", 
          row.names = FALSE)

write.csv(da_combined_scores,
          "results/refined_analysis/dopaminergic_scores.csv",
          row.names = FALSE)

saveRDS(list(
  assignments = refined_assignments,
  da_scores = da_combined_scores,
  all_results = all_results,
  score_matrix = score_matrix,
  marker_sets = MARKER_SETS
), "results/refined_analysis/complete_refined_analysis.rds")

cat("\n\n=== REFINED ANALYSIS COMPLETE ===\n")
cat("Key findings:\n")
cat(sprintf("- Dopaminergic neurons/lineage: %d clusters\n", 
            sum(grepl("Dopaminergic", refined_assignments$cell_type))))
cat(sprintf("- Unknown clusters reduced to: %d\n",
            sum(refined_assignments$cell_type == "Unknown")))
cat(sprintf("- High confidence assignments: %d\n",
            sum(refined_assignments$confidence %in% c("High", "Very High"))))

cat("\nResults saved to results/refined_analysis/\n")