#!/usr/bin/env Rscript

# ANALYZE RECLUSTERED EXPRESSION
# Comprehensive expression analysis after reclustering

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("COMPREHENSIVE EXPRESSION ANALYSIS ON RECLUSTERED DATA\n")
cat("=================================================================\n\n")

# Load the reclustered object
cat("1. Loading reclustered Seurat object...\n")
if (file.exists("results/seurat_obj_reclustered.rds")) {
  seurat_obj <- readRDS("results/seurat_obj_reclustered.rds")
  cat("Loaded from results/seurat_obj_reclustered.rds\n")
} else {
  seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds")
  cat("Loaded from original location\n")
}

DefaultAssay(seurat_obj) <- "SCT"

# Check clusters
n_coarse <- length(unique(seurat_obj$seurat_clusters_coarse))
n_fine <- length(unique(seurat_obj$seurat_clusters_fine))
cat(sprintf("\nCoarse clusters: %d\n", n_coarse))
cat(sprintf("Fine clusters: %d\n", n_fine))

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

# Function to analyze expression
analyze_cluster_expression <- function(seurat_obj, cluster_id, cluster_col = "seurat_clusters_coarse") {
  cells_in_cluster <- which(seurat_obj@meta.data[[cluster_col]] == as.character(cluster_id))
  n_cells <- length(cells_in_cluster)
  
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  results <- list()
  for (set_name in names(MARKER_SETS)) {
    genes <- MARKER_SETS[[set_name]]
    genes_present <- genes[genes %in% rownames(expr_data)]
    
    if (length(genes_present) > 0) {
      set_results <- data.frame(
        gene = genes_present,
        mean_expr = numeric(length(genes_present)),
        pct_expressing = numeric(length(genes_present)),
        stringsAsFactors = FALSE
      )
      
      for (i in 1:length(genes_present)) {
        gene_expr <- expr_data[genes_present[i], cells_in_cluster]
        set_results$mean_expr[i] <- mean(gene_expr)
        set_results$pct_expressing[i] <- 100 * sum(gene_expr > 0) / n_cells
      }
      
      # Filter for significant expression
      sig_genes <- set_results %>% filter(pct_expressing > 5)
      
      if (nrow(sig_genes) > 0) {
        results[[set_name]] <- list(
          n_genes = nrow(sig_genes),
          total_genes = length(genes),
          pct_coverage = round(100 * nrow(sig_genes) / length(genes), 1),
          avg_pct = round(mean(sig_genes$pct_expressing), 1),
          genes = paste(sig_genes$gene, collapse = ", "),
          details = sig_genes %>% arrange(desc(pct_expressing))
        )
      }
    }
  }
  
  return(list(n_cells = n_cells, marker_sets = results))
}

# Analyze all coarse clusters
cat("\n2. Analyzing all coarse clusters...\n")
cat("===================================\n")

coarse_results <- list()
coarse_identities <- data.frame(
  cluster = integer(),
  n_cells = integer(),
  identity = character(),
  confidence = character(),
  key_markers = character(),
  stringsAsFactors = FALSE
)

for (cl in sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_coarse))))) {
  cat(sprintf("\n--- COARSE CLUSTER %d ---\n", cl))
  
  result <- analyze_cluster_expression(seurat_obj, cl, "seurat_clusters_coarse")
  coarse_results[[as.character(cl)]] <- result
  
  cat(sprintf("Cells: %d\n", result$n_cells))
  
  # Display top marker sets
  if (length(result$marker_sets) > 0) {
    sorted_sets <- result$marker_sets[order(sapply(result$marker_sets, function(x) x$pct_coverage), decreasing = TRUE)]
    
    cat("\nTop marker sets:\n")
    for (i in 1:min(5, length(sorted_sets))) {
      set_name <- names(sorted_sets)[i]
      set_data <- sorted_sets[[set_name]]
      cat(sprintf("  %s: %d/%d genes (%.1f%%), avg %.1f%% cells\n",
                  set_name, set_data$n_genes, set_data$total_genes,
                  set_data$pct_coverage, set_data$avg_pct))
      cat(sprintf("    Genes: %s\n", set_data$genes))
    }
  }
  
  # Assign identity
  identity <- "Unknown"
  confidence <- "Low"
  key_markers <- ""
  
  if (length(result$marker_sets) > 0) {
    top_set <- names(sorted_sets)[1]
    
    # Identity logic based on marker expression
    if ("DA_CORE" %in% names(result$marker_sets) && result$marker_sets$DA_CORE$n_genes >= 3) {
      identity <- "Dopaminergic neurons"
      confidence <- "High"
      key_markers <- result$marker_sets$DA_CORE$genes
    } else if ("CHOROID_PLEXUS" %in% names(result$marker_sets) && result$marker_sets$CHOROID_PLEXUS$pct_coverage > 30) {
      identity <- "Choroid plexus"
      confidence <- "High"
      key_markers <- result$marker_sets$CHOROID_PLEXUS$genes
    } else if ("PROLIFERATING" %in% names(result$marker_sets) && result$marker_sets$PROLIFERATING$pct_coverage > 50) {
      identity <- "Proliferating cells"
      confidence <- "High"
      key_markers <- result$marker_sets$PROLIFERATING$genes
    } else if ("PAN_NEURONAL" %in% names(result$marker_sets) && result$marker_sets$PAN_NEURONAL$pct_coverage > 30) {
      if ("GLUTAMATERGIC" %in% names(result$marker_sets)) {
        identity <- "Glutamatergic neurons"
        confidence <- "Medium"
      } else if ("GABAERGIC" %in% names(result$marker_sets)) {
        identity <- "GABAergic neurons"
        confidence <- "Medium"
      } else {
        identity <- "Neurons"
        confidence <- "Medium"
      }
      key_markers <- result$marker_sets$PAN_NEURONAL$genes
    } else if ("NEURAL_PROGENITOR" %in% names(result$marker_sets)) {
      identity <- "Neural progenitors"
      confidence <- "Medium"
      key_markers <- result$marker_sets$NEURAL_PROGENITOR$genes
    } else if ("FIBROBLAST" %in% names(result$marker_sets) && result$marker_sets$FIBROBLAST$pct_coverage > 25) {
      identity <- "Fibroblasts"
      confidence <- "Medium"
      key_markers <- result$marker_sets$FIBROBLAST$genes
    }
  }
  
  cat(sprintf("\nIDENTITY: %s (%s confidence)\n", identity, confidence))
  
  coarse_identities <- rbind(coarse_identities, data.frame(
    cluster = cl,
    n_cells = result$n_cells,
    identity = identity,
    confidence = confidence,
    key_markers = substr(key_markers, 1, 100),
    stringsAsFactors = FALSE
  ))
}

# 3. Create coarse-to-fine mapping
cat("\n\n3. Creating coarse-to-fine cluster mapping...\n")
cat("============================================\n")

mapping <- table(
  Fine = seurat_obj$seurat_clusters_fine,
  Coarse = seurat_obj$seurat_clusters_coarse
)

# Find dominant coarse cluster for each fine cluster
fine_to_coarse <- data.frame(
  fine_cluster = integer(),
  coarse_cluster = integer(),
  n_cells = integer(),
  pct_of_fine = numeric(),
  stringsAsFactors = FALSE
)

for (fine_cl in rownames(mapping)) {
  row_data <- mapping[fine_cl, ]
  total_cells <- sum(row_data)
  dominant_coarse <- which.max(row_data)
  
  fine_to_coarse <- rbind(fine_to_coarse, data.frame(
    fine_cluster = as.integer(fine_cl),
    coarse_cluster = as.integer(names(dominant_coarse)),
    n_cells = total_cells,
    pct_of_fine = round(100 * max(row_data) / total_cells, 1),
    stringsAsFactors = FALSE
  ))
}

fine_to_coarse <- fine_to_coarse %>% arrange(fine_cluster)

cat("\nFine-to-coarse mapping (showing dominant coarse cluster):\n")
print(head(fine_to_coarse, 10))
cat(sprintf("... and %d more fine clusters\n", nrow(fine_to_coarse) - 10))

# 4. Save all results
cat("\n4. Saving results...\n")
dir.create("results/reclustered_analysis", recursive = TRUE, showWarnings = FALSE)

saveRDS(coarse_results, "results/reclustered_analysis/coarse_cluster_expression.rds")
write.csv(coarse_identities, "results/reclustered_analysis/coarse_cluster_identities.csv", row.names = FALSE)
write.csv(fine_to_coarse, "results/reclustered_analysis/fine_to_coarse_mapping.csv", row.names = FALSE)

# 5. Summary
cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("=========================\n\n")

cat("Coarse cluster summary:\n")
print(coarse_identities)

cat("\n\nNext steps:\n")
cat("1. Review coarse cluster identities\n")
cat("2. Analyze fine clusters within each coarse cluster context\n")
cat("3. Assign final cell type labels\n")

cat("\nResults saved to:\n")
cat("- results/reclustered_analysis/coarse_cluster_expression.rds\n")
cat("- results/reclustered_analysis/coarse_cluster_identities.csv\n")
cat("- results/reclustered_analysis/fine_to_coarse_mapping.csv\n")