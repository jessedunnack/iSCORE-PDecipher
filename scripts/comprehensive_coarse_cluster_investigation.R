#!/usr/bin/env Rscript

# COMPREHENSIVE COARSE CLUSTER INVESTIGATION
# Using direct Seurat object expression data

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("COMPREHENSIVE COARSE CLUSTER INVESTIGATION\n")
cat("=================================================================\n\n")

# Function to check gene expression
check_gene_expression <- function(seurat_obj, genes, cluster_col = "seurat_clusters_coarse") {
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  genes_present <- genes[genes %in% rownames(expr_data)]
  
  if (length(genes_present) == 0) return(NULL)
  
  clusters <- seurat_obj@meta.data[[cluster_col]]
  results <- list()
  
  for (gene in genes_present) {
    gene_expr <- expr_data[gene, ]
    cluster_stats <- data.frame(
      cluster = levels(factor(clusters)),
      stringsAsFactors = FALSE
    )
    
    for (cl in cluster_stats$cluster) {
      cells_in_cluster <- which(clusters == cl)
      expr_in_cluster <- gene_expr[cells_in_cluster]
      
      cluster_stats[cluster_stats$cluster == cl, "mean_expr"] <- mean(expr_in_cluster)
      cluster_stats[cluster_stats$cluster == cl, "pct_expressing"] <- 
        100 * sum(expr_in_cluster > 0) / length(expr_in_cluster)
      cluster_stats[cluster_stats$cluster == cl, "n_cells"] <- length(cells_in_cluster)
    }
    
    cluster_stats$gene <- gene
    results[[gene]] <- cluster_stats
  }
  
  combined <- do.call(rbind, results)
  return(combined)
}

# Load Seurat object
cat("1. Loading Seurat object...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds")
DefaultAssay(seurat_obj) <- "SCT"

# Check coarse clusters
coarse_clusters <- unique(seurat_obj$seurat_clusters_coarse)
cat("\nCoarse clusters found:", paste(sort(as.numeric(as.character(coarse_clusters))), collapse = ", "), "\n")
cat("Total coarse clusters:", length(coarse_clusters), "\n")

# Define comprehensive marker sets
MARKER_SETS <- list(
  # Dopaminergic
  DA_CORE = c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6"),
  DA_TRANSCRIPTION = c("LMX1A", "LMX1B", "FOXA2", "NR4A2", "PITX3", "EN1", "EN2"),
  DA_SUBTYPE = c("ALDH1A1", "SOX6", "CALB1", "CALB2", "OTX2", "SNCG"),
  
  # Other neuronal types
  GLUTAMATERGIC = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "TBR1", "SATB2"),
  GABAERGIC = c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "SST", "PVALB", "VIP"),
  SEROTONERGIC = c("TPH2", "SLC6A4", "FEV", "GATA2", "GATA3"),
  CHOLINERGIC = c("CHAT", "SLC18A3", "ACHE"),
  
  # General neuronal
  PAN_NEURONAL = c("TUBB3", "MAP2", "STMN2", "SYN1", "SNAP25", "RBFOX3"),
  
  # Progenitors
  NEURAL_PROGENITOR = c("SOX2", "SOX1", "NES", "PAX6", "HES1", "HES5", "VIM"),
  FLOOR_PLATE = c("CORIN", "ARX", "SHH", "WNT1", "MSX1", "LMO3"),
  RADIAL_GLIA = c("SOX9", "HOPX", "GFAP", "AQP4", "TNC", "FABP7"),
  
  # Glial
  OLIGODENDROCYTE = c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "MBP", "PLP1", "MOG"),
  ASTROCYTE = c("GFAP", "S100B", "ALDH1L1", "AQP4", "SLC1A2", "SLC1A3"),
  EPENDYMAL = c("FOXJ1", "PIFO", "SPEF2", "DNAH5", "CCDC39"),
  
  # Non-neuronal
  CHOROID_PLEXUS = c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13", "ENPP2"),
  ENDOTHELIAL = c("CLDN5", "PECAM1", "VWF", "FLT1", "CDH5"),
  PERICYTE = c("PDGFRB", "RGS5", "ABCC9", "KCNJ8", "NOTCH3"),
  SMOOTH_MUSCLE = c("TAGLN", "ACTA2", "MYL9", "CNN1", "MYLK"),
  FIBROBLAST = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "VIM"),
  
  # Cell states
  PROLIFERATING = c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "HIST1H3B"),
  STRESS = c("FOS", "JUN", "EGR1", "ATF3", "HSPA1A", "DDIT3"),
  
  # Regional/special
  HYPOTHALAMIC = c("HCRT", "OXT", "AVP", "CRH", "TRH", "GHRH"),
  MIDBRAIN = c("EN1", "EN2", "PAX2", "PAX5", "WNT1", "FGF8")
)

# Analyze each coarse cluster
all_results <- list()

for (cluster_id in sort(as.numeric(as.character(coarse_clusters)))) {
  cat(sprintf("\n\n========== COARSE CLUSTER %d ==========\n", cluster_id))
  
  # Get number of cells
  n_cells <- sum(seurat_obj$seurat_clusters_coarse == as.character(cluster_id))
  cat(sprintf("Number of cells: %d\n", n_cells))
  
  # Check each marker set
  cluster_scores <- list()
  
  for (set_name in names(MARKER_SETS)) {
    genes <- MARKER_SETS[[set_name]]
    expr_data <- check_gene_expression(seurat_obj, genes, "seurat_clusters_coarse")
    
    if (!is.null(expr_data)) {
      cluster_expr <- expr_data %>% 
        filter(cluster == as.character(cluster_id)) %>%
        filter(pct_expressing > 5)  # At least 5% expression
      
      if (nrow(cluster_expr) > 0) {
        cluster_scores[[set_name]] <- list(
          n_genes = nrow(cluster_expr),
          total_genes = length(genes),
          pct_coverage = round(100 * nrow(cluster_expr) / length(genes), 1),
          avg_pct_expressing = round(mean(cluster_expr$pct_expressing), 1),
          genes = paste(cluster_expr$gene, collapse = ", ")
        )
      }
    }
  }
  
  # Display significant marker sets
  if (length(cluster_scores) > 0) {
    cat("\nSignificant marker sets (>5% cells expressing):\n")
    
    # Sort by coverage
    sorted_scores <- cluster_scores[order(sapply(cluster_scores, function(x) x$pct_coverage), decreasing = TRUE)]
    
    for (set_name in names(sorted_scores)) {
      score <- sorted_scores[[set_name]]
      cat(sprintf("  %s: %d/%d genes (%.1f%%), avg %.1f%% cells\n",
                  set_name, score$n_genes, score$total_genes, 
                  score$pct_coverage, score$avg_pct_expressing))
      cat(sprintf("    Genes: %s\n", score$genes))
    }
  } else {
    cat("\nNo significant marker sets found (>5% expression)\n")
  }
  
  # Get top expressed genes from original markers if available
  if (cluster_id %in% 0:14) {  # Only for clusters we have marker data for
    cat("\nChecking original top markers:\n")
    original_markers <- c("LHX9", "C1QL4", "MAB21L1", "STMN2", "C22orf42",  # Cluster 0
                         "AL138899.1", "TMC5", "LINC02838", "MS4A6A", "TMEM190",  # Cluster 1
                         "PTPRZ1", "AC092957.1", "RPRML", "LRP2", "PTN",  # Cluster 2
                         "TMEM119", "PRRX2", "PRRX1", "TWIST2", "COL6A3",  # Cluster 3
                         "SCN1A-AS1", "NRG3", "CNTNAP5", "MIR137HG", "CCSER1",  # Cluster 4
                         "LINC02539", "PRSS56", "RGCC", "NMU", "AC095030.1",  # Cluster 5
                         "GDF15", "PLCG2", "CRYAB", "ATF5", "DDIT3",  # Cluster 6
                         "ANKRD66", "FAM216B", "TTR", "C11orf97", "AL357093.2",  # Cluster 7
                         "DCN", "H19", "FGL2", "PLSCR5", "MATN3",  # Cluster 8
                         "PTGDS", "LCNL1", "CPXM2", "HPD", "PRDM16-DT",  # Cluster 9
                         "MKI67", "UBE2C", "KIF20A", "HIST1H1B", "NDC80",  # Cluster 10
                         "CRABP1", "MIR217HG", "ST8SIA6", "VSTM2B", "AC017050.1",  # Cluster 11
                         "CGA", "CALCA", "GP2", "TPH1", "SMIM22",  # Cluster 12
                         "RBP4", "FFAR4", "AL450352.1", "AC092142.1", "SLC22A8",  # Cluster 13
                         "HCRT", "C19orf85", "MAL", "TMOD1", "CRYM")  # Cluster 14
    
    # Sample a few markers based on cluster
    sample_markers <- switch(as.character(cluster_id),
                           "0" = original_markers[1:5],
                           "1" = original_markers[6:10],
                           "2" = original_markers[11:15],
                           "3" = original_markers[16:20],
                           "4" = original_markers[21:25],
                           "5" = original_markers[26:30],
                           "6" = original_markers[31:35],
                           "7" = original_markers[36:40],
                           "8" = original_markers[41:45],
                           NULL)
    
    if (!is.null(sample_markers)) {
      sample_expr <- check_gene_expression(seurat_obj, sample_markers, "seurat_clusters_coarse")
      if (!is.null(sample_expr)) {
        cluster_sample <- sample_expr %>% 
          filter(cluster == as.character(cluster_id)) %>%
          filter(pct_expressing > 10) %>%
          select(gene, pct_expressing)
        if (nrow(cluster_sample) > 0) {
          print(cluster_sample)
        }
      }
    }
  }
  
  # Store results
  all_results[[as.character(cluster_id)]] <- list(
    n_cells = n_cells,
    scores = cluster_scores
  )
  
  # Make identity assignment
  cat("\n**PRELIMINARY IDENTITY ASSIGNMENT**:\n")
  
  # Logic for assignment based on marker expression
  if ("DA_CORE" %in% names(cluster_scores) && cluster_scores$DA_CORE$n_genes >= 3) {
    cat(">>> DOPAMINERGIC NEURONS <<<\n")
  } else if ("CHOROID_PLEXUS" %in% names(cluster_scores) && cluster_scores$CHOROID_PLEXUS$pct_coverage > 30) {
    cat(">>> CHOROID PLEXUS <<<\n")
  } else if ("PROLIFERATING" %in% names(cluster_scores) && cluster_scores$PROLIFERATING$pct_coverage > 50) {
    cat(">>> PROLIFERATING CELLS <<<\n")
  } else if ("PAN_NEURONAL" %in% names(cluster_scores) && cluster_scores$PAN_NEURONAL$pct_coverage > 30) {
    if ("GLUTAMATERGIC" %in% names(cluster_scores)) {
      cat(">>> GLUTAMATERGIC NEURONS <<<\n")
    } else if ("GABAERGIC" %in% names(cluster_scores)) {
      cat(">>> GABAERGIC NEURONS <<<\n")
    } else {
      cat(">>> NEURONS (subtype unclear) <<<\n")
    }
  } else if ("NEURAL_PROGENITOR" %in% names(cluster_scores) && cluster_scores$NEURAL_PROGENITOR$pct_coverage > 20) {
    cat(">>> NEURAL PROGENITORS <<<\n")
  } else if ("OLIGODENDROCYTE" %in% names(cluster_scores) && cluster_scores$OLIGODENDROCYTE$pct_coverage > 20) {
    cat(">>> OLIGODENDROCYTE LINEAGE <<<\n")
  } else if ("FIBROBLAST" %in% names(cluster_scores) && cluster_scores$FIBROBLAST$pct_coverage > 30) {
    cat(">>> FIBROBLASTS/MESENCHYMAL <<<\n")
  } else if ("STRESS" %in% names(cluster_scores) && cluster_scores$STRESS$avg_pct_expressing > 30) {
    cat(">>> STRESSED/DAMAGED CELLS <<<\n")
  } else {
    cat(">>> UNCERTAIN/MIXED POPULATION <<<\n")
  }
}

# Summary
cat("\n\n=== SUMMARY OF COARSE CLUSTERS ===\n")
cat("===================================\n\n")

# Save comprehensive results
saveRDS(all_results, "results/comprehensive_coarse_cluster_results.rds")
cat("\nResults saved to: results/comprehensive_coarse_cluster_results.rds\n")

# Print final summary table
summary_df <- data.frame(
  cluster = integer(),
  n_cells = integer(),
  top_marker_set = character(),
  preliminary_identity = character(),
  stringsAsFactors = FALSE
)

for (cl in names(all_results)) {
  result <- all_results[[cl]]
  if (length(result$scores) > 0) {
    top_set <- names(result$scores)[1]  # Already sorted by coverage
    
    # Assign identity
    identity <- if (grepl("DA_CORE", top_set)) "Dopaminergic" 
    else if (grepl("CHOROID", top_set)) "Choroid Plexus"
    else if (grepl("PROLIFERAT", top_set)) "Proliferating"
    else if (grepl("NEURONAL", top_set)) "Neurons"
    else if (grepl("PROGENITOR", top_set)) "Progenitors"
    else "Mixed/Uncertain"
    
    summary_df <- rbind(summary_df, data.frame(
      cluster = as.integer(cl),
      n_cells = result$n_cells,
      top_marker_set = top_set,
      preliminary_identity = identity
    ))
  } else {
    summary_df <- rbind(summary_df, data.frame(
      cluster = as.integer(cl),
      n_cells = result$n_cells,
      top_marker_set = "None",
      preliminary_identity = "Unknown"
    ))
  }
}

summary_df <- summary_df %>% arrange(cluster)
cat("\nCoarse Cluster Summary:\n")
print(summary_df)

write.csv(summary_df, "results/coarse_cluster_identity_summary.csv", row.names = FALSE)