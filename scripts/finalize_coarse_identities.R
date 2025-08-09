#!/usr/bin/env Rscript

# FINALIZE COARSE CLUSTER IDENTITIES
# Based on comprehensive marker analysis

library(Seurat)
library(dplyr)

cat("=================================================================\n")
cat("FINALIZING COARSE CLUSTER IDENTITIES\n")
cat("=================================================================\n\n")

# Load the reclustered Seurat object to check TH expression
seurat_obj <- readRDS("results/seurat_obj_reclustered.rds")
DefaultAssay(seurat_obj) <- "SCT"

# Check TH expression in key clusters
cat("1. Verifying TH expression in key clusters:\n")
cat("==========================================\n")

expr_data <- GetAssayData(seurat_obj, slot = "data")
if ("TH" %in% rownames(expr_data)) {
  for (cl in c(0, 11, 12)) {
    cells <- which(seurat_obj$seurat_clusters_coarse == as.character(cl))
    th_expr <- expr_data["TH", cells]
    pct_th <- 100 * sum(th_expr > 0) / length(cells)
    mean_th <- mean(th_expr[th_expr > 0])
    
    cat(sprintf("Cluster %d: TH in %.1f%% cells (mean expr in positive cells: %.2f)\n", 
                cl, pct_th, mean_th))
  }
}

# Create final identity assignments
cat("\n\n2. Final cluster identity assignments:\n")
cat("=====================================\n")

final_identities <- data.frame(
  cluster = 0:14,
  identity = c(
    "Neurons_Dopaminergic",                    # 0 - STMN2+, check TH
    "Progenitors_Intermediate",                # 1 - Unknown, no clear markers
    "Progenitors_PTPRZ1+",                     # 2 - PTPRZ1 high, progenitor-like
    "Mesenchymal_Fibroblasts",                 # 3 - PRRX1/2, COL1A1
    "Progenitors_Uncommitted",                 # 4 - No clear identity
    "Cells_Unidentified",                      # 5 - LINC RNAs, unclear
    "Cells_Stressed",                          # 6 - GDF15, stress response
    "Choroid_Plexus",                          # 7 - TTR, FOLR1, KCNJ13
    "Fibroblasts_ECM",                         # 8 - DCN, LUM
    "Cells_PTGDS+",                            # 9 - PTGDS high, unclear identity
    "Cells_Proliferating",                     # 10 - MKI67, TOP2A
    "Progenitors_CRABP1+",                     # 11 - CRABP1, some TH
    "Cells_Neuroendocrine",                    # 12 - CGA, CALCA, TPH1
    "Cells_RBP4+",                             # 13 - RBP4 high, unclear
    "Neurons_Hypothalamic_HCRT"                # 14 - HCRT/hypocretin neurons
  ),
  cell_type_broad = c(
    "Neurons",                                 # 0
    "Progenitors",                             # 1
    "Progenitors",                             # 2
    "Non-neuronal",                            # 3
    "Progenitors",                             # 4
    "Unknown",                                 # 5
    "Non-neuronal",                            # 6
    "Non-neuronal",                            # 7
    "Non-neuronal",                            # 8
    "Unknown",                                 # 9
    "Progenitors",                             # 10
    "Progenitors",                             # 11
    "Non-neuronal",                            # 12
    "Unknown",                                 # 13
    "Neurons"                                  # 14
  ),
  confidence = c(
    "Medium",  # 0 - Need to verify TH
    "Low",     # 1
    "Medium",  # 2
    "High",    # 3
    "Low",     # 4
    "Low",     # 5
    "High",    # 6
    "High",    # 7
    "High",    # 8
    "Low",     # 9
    "High",    # 10
    "Medium",  # 11
    "High",    # 12
    "Low",     # 13
    "High"     # 14
  ),
  notes = c(
    "STMN2+, verify TH expression",            # 0
    "No clear markers",                        # 1
    "PTPRZ1+ radial glia-like",                # 2
    "Mesenchymal: PRRX1/2, COL1A1",           # 3
    "No specific markers",                     # 4
    "LncRNA enriched",                         # 5
    "Stress response: GDF15, ATF3",            # 6
    "Classic choroid plexus markers",          # 7
    "ECM-producing: DCN, LUM",                 # 8
    "PTGDS+ unclear identity",                 # 9
    "Actively dividing: MKI67+",               # 10
    "CRABP1+, weak TH (12%)",                  # 11
    "Neuroendocrine: CGA, CALCA",             # 12
    "RBP4+ unclear identity",                  # 13
    "Hypocretin/orexin neurons"                # 14
  ),
  stringsAsFactors = FALSE
)

# Print summary
cat("\nCluster identities:\n")
print(final_identities[, c("cluster", "identity", "cell_type_broad", "confidence")])

# Count by broad type
cat("\n\nSummary by broad cell type:\n")
summary_broad <- final_identities %>%
  group_by(cell_type_broad) %>%
  summarise(
    n_clusters = n(),
    clusters = paste(cluster, collapse = ", ")
  )
print(summary_broad)

# Load cell counts
if (file.exists("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")) {
  prev_results <- read.csv("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")
  final_identities <- final_identities %>%
    left_join(prev_results[, c("cluster", "n_cells")], by = "cluster")
  
  # Calculate percentage of cells
  final_identities$pct_cells <- round(100 * final_identities$n_cells / sum(final_identities$n_cells), 1)
  
  cat("\n\nCell distribution:\n")
  cell_summary <- final_identities %>%
    group_by(cell_type_broad) %>%
    summarise(
      total_cells = sum(n_cells),
      pct_cells = round(100 * total_cells / sum(final_identities$n_cells), 1)
    )
  print(cell_summary)
}

# Save final identities
write.csv(final_identities, "results/reclustered_analysis/coarse_cluster_identities_FINAL.csv", 
          row.names = FALSE)

cat("\n\nKey findings:\n")
cat("- Only 2 clusters are definitively neurons (0 and 14)\n")
cat("- Cluster 0 appears to be dopaminergic (pending TH verification)\n")
cat("- Cluster 14 is hypothalamic HCRT+ neurons\n")
cat("- Many clusters are progenitors or non-neuronal support cells\n")
cat("- High stress markers across all clusters suggest culture conditions\n")
cat("\nSaved to: results/reclustered_analysis/coarse_cluster_identities_FINAL.csv\n")

# Next step reminder
cat("\n\nNEXT STEPS:\n")
cat("1. Review and approve these identities\n")
cat("2. Add to Seurat object as seurat_obj$celltypes_coarse\n")
cat("3. Analyze fine clusters within each coarse identity\n")
cat("4. Consider re-running analyses with proper cell type labels\n")