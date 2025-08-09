#!/usr/bin/env Rscript

# COMPREHENSIVE COARSE CLUSTER REVIEW
# Combining top markers + key fate determinants

library(dplyr)

cat("=================================================================\n")
cat("COMPREHENSIVE COARSE CLUSTER REVIEW\n")
cat("Analyzing top markers + key fate determinants\n")
cat("=================================================================\n\n")

# Load coarse cluster markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# Define key fate determinant gene sets from documentation
FATE_DETERMINANTS <- list(
  # Dopaminergic neurons
  DA_ESSENTIAL = c("TH", "DDC", "AADC"),
  DA_MATURE = c("SLC6A3", "DAT", "SLC18A2", "VMAT2", "DRD2", "KCNJ6", "GIRK2"),
  DA_TRANSCRIPTION = c("LMX1A", "LMX1B", "FOXA2", "NR4A2", "NURR1", "PITX3", "EN1", "EN2"),
  DA_VULNERABLE_A9 = c("SOX6", "ALDH1A1", "SNCG", "ATP13A2", "RIT2", "AGTR1"),
  DA_RESILIENT_A10 = c("CALB1", "CALB2", "OTX2", "GRP", "CCK", "VIP"),
  MEIS_VULNERABILITY = c("MEIS1", "MEIS2", "PBX1", "PBX3", "HOXA9", "HOXA10"),
  
  # Floor plate and midbrain
  FLOOR_PLATE = c("CORIN", "LMX1A", "FOXA2", "ARX", "SHH", "WNT1", "MSX1", "LMO3", "FOXA1"),
  MIDBRAIN = c("EN1", "EN2", "PAX2", "PAX5", "WNT1", "FGF8", "OTX2", "GBX2"),
  
  # General neuronal
  PAN_NEURONAL = c("TUBB3", "MAP2", "STMN2", "GAP43", "SYN1", "SNAP25", "RBFOX3", "NEFL", "NEFM", "MAPT", "NCAM1"),
  MATURE_NEURONAL = c("SYT1", "SYP", "VAMP2", "STX1A", "NRXN1", "NLGN1", "CAMK2A"),
  
  # Neurotransmitter types
  GLUTAMATERGIC = c("SLC17A6", "VGLUT2", "SLC17A7", "VGLUT1", "GRIN1", "GRIN2B", "GRIA1"),
  GABAERGIC = c("GAD1", "GAD2", "SLC32A1", "VGAT", "GABBR1", "GABRG2"),
  SEROTONERGIC = c("TPH2", "SLC6A4", "SERT", "FEV", "GATA2", "GATA3", "PET1"),
  CHOLINERGIC = c("CHAT", "SLC18A3", "VACHT", "ACHE"),
  
  # Progenitor and glial
  NEURAL_PROGENITOR = c("SOX2", "SOX1", "NES", "PAX6", "HES1", "HES5", "VIM", "FABP7", "BLBP"),
  RADIAL_GLIA = c("SOX9", "HOPX", "GFAP", "AQP4", "TNC", "BLBP", "GLAST", "SLC1A3"),
  NEUROBLAST = c("DCX", "NEUROD1", "NEUROD2", "NEUROG2", "ASCL1", "DLL3"),
  
  # Non-neuronal
  OLIGODENDROCYTE = c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "NKX2-2", "MPZ", "MBP", "PLP1", "MOG", "MAG"),
  ASTROCYTE = c("GFAP", "S100B", "ALDH1L1", "AQP4", "GLT1", "SLC1A2", "GLAST", "SLC1A3"),
  EPENDYMAL = c("FOXJ1", "PIFO", "SPEF2", "DNAH5", "CCDC39", "RFX3"),
  CHOROID_PLEXUS = c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "PRLR", "KCNJ13", "ENPP2", "PLPP3", "SLC13A4"),
  
  # Vascular and mesenchymal
  ENDOTHELIAL = c("CLDN5", "PECAM1", "CD31", "VWF", "FLT1", "KDR", "CDH5", "TIE1"),
  PERICYTE = c("PDGFRB", "RGS5", "ABCC9", "KCNJ8", "NOTCH3", "DES"),
  SMOOTH_MUSCLE = c("TAGLN", "ACTA2", "MYL9", "CNN1", "MYLK", "MYOCD"),
  FIBROBLAST = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA", "VIM", "FAP"),
  
  # Cell state markers
  PROLIFERATING = c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "HIST1H3B", "HIST1H4C"),
  STRESS_RESPONSE = c("FOS", "JUN", "EGR1", "ATF3", "HSPA1A", "HSPA1B", "DDIT3", "ATF4"),
  HYPOXIA = c("HIF1A", "VEGFA", "SLC2A1", "ENO1", "LDHA", "PGK1")
)

# Function to check for presence of gene sets
check_gene_sets <- function(cluster_markers, gene_sets) {
  results <- list()
  
  for (set_name in names(gene_sets)) {
    genes <- gene_sets[[set_name]]
    present_genes <- genes[genes %in% cluster_markers$gene]
    
    if (length(present_genes) > 0) {
      # Get expression info for present genes
      expr_info <- cluster_markers %>%
        filter(gene %in% present_genes) %>%
        arrange(desc(avg_log2FC))
      
      results[[set_name]] <- list(
        n_present = length(present_genes),
        n_total = length(genes),
        pct_present = round(100 * length(present_genes) / length(genes), 1),
        genes = present_genes,
        avg_fc = round(mean(expr_info$avg_log2FC), 2),
        max_fc = round(max(expr_info$avg_log2FC), 2),
        expr_details = expr_info
      )
    }
  }
  
  return(results)
}

# Comprehensive analysis for each cluster
analyze_cluster_comprehensive <- function(cluster_id) {
  cat(sprintf("\n========== COARSE CLUSTER %d ==========\n", cluster_id))
  
  # Get all markers for this cluster
  cluster_markers <- coarse_markers %>%
    filter(cluster == as.character(cluster_id))
  
  # Get top 30 markers
  top_markers <- cluster_markers %>%
    arrange(desc(avg_log2FC)) %>%
    head(30)
  
  cat("\nTop 15 marker genes:\n")
  for (i in 1:min(15, nrow(top_markers))) {
    cat(sprintf("%2d. %-12s (FC=%.2f, pct.1=%.3f)\n", 
                i, 
                top_markers$gene[i], 
                top_markers$avg_log2FC[i],
                top_markers$pct.1[i]))
  }
  
  # Check for fate determinants
  cat("\n--- FATE DETERMINANT ANALYSIS ---\n")
  gene_set_results <- check_gene_sets(cluster_markers, FATE_DETERMINANTS)
  
  if (length(gene_set_results) == 0) {
    cat("No major fate determinants detected in top markers\n")
  } else {
    # Sort by percentage present
    sorted_results <- gene_set_results[order(sapply(gene_set_results, function(x) x$pct_present), decreasing = TRUE)]
    
    for (set_name in names(sorted_results)) {
      result <- sorted_results[[set_name]]
      cat(sprintf("\n%s: %d/%d genes (%.1f%%)\n", 
                  set_name, result$n_present, result$n_total, result$pct_present))
      cat(sprintf("  Present: %s\n", paste(result$genes, collapse = ", ")))
      cat(sprintf("  Avg FC: %.2f, Max FC: %.2f\n", result$avg_fc, result$max_fc))
      
      # Show top 3 genes with expression
      if (nrow(result$expr_details) > 0) {
        cat("  Top expressed:\n")
        for (j in 1:min(3, nrow(result$expr_details))) {
          gene_info <- result$expr_details[j,]
          cat(sprintf("    - %s: FC=%.2f, pct.1=%.3f\n", 
                      gene_info$gene, gene_info$avg_log2FC, gene_info$pct.1))
        }
      }
    }
  }
  
  # Additional markers at positions 16-30
  if (nrow(top_markers) > 15) {
    cat("\nAdditional markers (16-30):\n")
    additional <- top_markers[16:min(30, nrow(top_markers)), ]
    cat(paste(additional$gene, collapse = ", "), "\n")
  }
  
  cat("\n")
  return(list(
    cluster = cluster_id,
    top_markers = top_markers,
    gene_sets = gene_set_results
  ))
}

# Analyze all coarse clusters
all_results <- list()
for (i in 0:14) {
  all_results[[as.character(i)]] <- analyze_cluster_comprehensive(i)
}

# Save comprehensive results
saveRDS(all_results, "results/coarse_cluster_comprehensive_review.rds")
cat("\n\nResults saved to: results/coarse_cluster_comprehensive_review.rds\n")