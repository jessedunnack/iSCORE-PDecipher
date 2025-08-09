#!/usr/bin/env Rscript

# ANALYZE RECLUSTERED EXPRESSION WITH COMPREHENSIVE STRESS MARKERS
# Focus on PD-relevant stress pathways

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("COMPREHENSIVE EXPRESSION ANALYSIS WITH STRESS/PD FOCUS\n")
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

# Define comprehensive marker sets with PD/stress focus
MARKER_SETS <- list(
  # Dopaminergic
  DA_CORE = c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6"),
  DA_TRANSCRIPTION = c("LMX1A", "LMX1B", "FOXA2", "NR4A2", "PITX3", "EN1", "EN2"),
  DA_SUBTYPE = c("ALDH1A1", "SOX6", "CALB1", "CALB2", "OTX2", "SNCG"),
  
  # STRESS MARKERS - COMPREHENSIVE
  STRESS_GENERAL = c("FOS", "JUN", "EGR1", "ATF3", "ATF4", "ATF5", "DDIT3", "XBP1", "ERN1"),
  HEAT_SHOCK = c("HSPA1A", "HSPA1B", "HSPA5", "HSP90AA1", "HSP90AB1", "HSPB1", "DNAJB1", "HSPH1", "HSF1"),
  NEURONAL_STRESS = c("TTR", "CLU", "APOE", "B2M", "SERPINA3", "GFAP", "VIM"),
  
  # AUTOPHAGY PATHWAYS
  AUTOPHAGY_CORE = c("MAP1LC3B", "MAP1LC3A", "GABARAP", "GABARAPL1", "GABARAPL2", "SQSTM1", "NBR1"),
  AUTOPHAGY_REGULATION = c("BECN1", "ATG5", "ATG7", "ATG12", "ATG16L1", "ULK1", "MTOR", "TFEB"),
  LYSOSOME = c("LAMP1", "LAMP2", "CTSD", "CTSB", "CTSL", "GBA", "GBA2", "LIPA", "HEXB"),
  
  # MITOCHONDRIAL FUNCTION
  MITO_RESPIRATION = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ND5", "MT-ND6"),
  MITO_DYNAMICS = c("MFN1", "MFN2", "OPA1", "DNM1L", "FIS1", "PINK1", "PRKN", "BNIP3", "BNIP3L"),
  MITO_BIOGENESIS = c("PPARGC1A", "NRF1", "TFAM", "POLG", "SIRT1", "SIRT3", "PPRC1"),
  
  # OXIDATIVE STRESS
  ROS_DEFENSE = c("SOD1", "SOD2", "SOD3", "CAT", "GPX1", "GPX4", "GSR", "PRDX1", "PRDX2", "PRDX3"),
  ANTIOXIDANT_RESPONSE = c("NFE2L2", "KEAP1", "NQO1", "HMOX1", "GCLC", "GCLM", "GSS", "TXNRD1"),
  OXIDATIVE_DAMAGE = c("4HNE", "HIF1A", "VEGFA", "EPO", "SLC2A1", "LDHA"),
  
  # ER STRESS/UPR
  ER_STRESS = c("HSPA5", "HSP90B1", "CALR", "PDIA4", "CANX", "DNAJC3", "HYOU1"),
  UPR_SENSORS = c("ERN1", "EIF2AK3", "ATF6", "XBP1", "DDIT3", "PPP1R15A", "TRIB3"),
  
  # PROTEIN AGGREGATION/CLEARANCE
  PROTEIN_FOLDING = c("HSPA8", "HSPD1", "CCT2", "CCT3", "CCT4", "CCT5", "TCP1", "PPIB"),
  UBIQUITIN_PROTEASOME = c("UBB", "UBC", "UBA52", "PSMD1", "PSMD2", "PSMB5", "USP9X", "UCHL1"),
  AGGREPHAGY = c("SQSTM1", "NBR1", "OPTN", "CALCOCO2", "TAX1BP1", "TOLLIP"),
  
  # NEUROINFLAMMATION
  INFLAMMATION = c("IL1B", "IL6", "TNF", "NFKB1", "RELA", "IKBKB", "TLR2", "TLR4"),
  MICROGLIAL_ACTIVATION = c("AIF1", "CD68", "TREM2", "TYROBP", "CX3CR1", "P2RY12", "TMEM119"),
  ASTROCYTE_REACTIVE = c("GFAP", "VIM", "S100B", "SERPINA3", "C3", "CXCL10", "TIMP1"),
  
  # PD-SPECIFIC GENES
  PD_RISK_GENES = c("SNCA", "LRRK2", "GBA", "VPS35", "PRKN", "PINK1", "DJ1", "ATP13A2", "PLA2G6", "FBXO7"),
  ALPHA_SYNUCLEIN = c("SNCA", "SNCG", "SNCB", "SYNJ1", "DNAJC6", "AUXILIN"),
  
  # CALCIUM HOMEOSTASIS
  CALCIUM_SIGNALING = c("CALB1", "CALB2", "CALM1", "CALM2", "CALM3", "CAMK2A", "CAMK2B", "ATP2A2"),
  
  # IRON METABOLISM
  IRON_METABOLISM = c("FTH1", "FTL", "TF", "TFRC", "SLC11A2", "SLC40A1", "HAMP", "IREB2"),
  
  # Other neuronal types
  GLUTAMATERGIC = c("SLC17A6", "SLC17A7", "GRIN1", "GRIN2B", "TBR1", "SATB2"),
  GABAERGIC = c("GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2", "SST", "PVALB", "VIP"),
  PAN_NEURONAL = c("TUBB3", "MAP2", "STMN2", "SYN1", "SNAP25", "RBFOX3"),
  
  # Progenitors and glia
  NEURAL_PROGENITOR = c("SOX2", "SOX1", "NES", "PAX6", "HES1", "HES5", "VIM"),
  OLIGODENDROCYTE = c("OLIG1", "OLIG2", "SOX10", "PDGFRA", "MBP", "PLP1", "MOG"),
  ASTROCYTE = c("GFAP", "S100B", "ALDH1L1", "AQP4", "SLC1A2", "SLC1A3"),
  
  # Non-neuronal
  CHOROID_PLEXUS = c("FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13", "ENPP2"),  # Removed TTR
  ENDOTHELIAL = c("CLDN5", "PECAM1", "VWF", "FLT1", "CDH5"),
  FIBROBLAST = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "VIM"),
  
  # Cell states
  PROLIFERATING = c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "HIST1H3B"),
  
  # Regional/special
  HYPOTHALAMIC = c("HCRT", "OXT", "AVP", "CRH", "TRH", "GHRH")
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
cat("\n2. Analyzing all coarse clusters with stress focus...\n")
cat("=====================================================\n")

coarse_results <- list()
coarse_identities <- data.frame(
  cluster = integer(),
  n_cells = integer(),
  identity = character(),
  confidence = character(),
  stress_level = character(),
  key_markers = character(),
  stress_markers = character(),
  stringsAsFactors = FALSE
)

for (cl in sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_coarse))))) {
  cat(sprintf("\n--- COARSE CLUSTER %d ---\n", cl))
  
  result <- analyze_cluster_expression(seurat_obj, cl, "seurat_clusters_coarse")
  coarse_results[[as.character(cl)]] <- result
  
  cat(sprintf("Cells: %d\n", result$n_cells))
  
  # Check stress levels
  stress_score <- 0
  stress_markers <- c()
  
  # Calculate stress score
  for (stress_set in c("STRESS_GENERAL", "HEAT_SHOCK", "NEURONAL_STRESS", "ER_STRESS", 
                       "OXIDATIVE_DAMAGE", "INFLAMMATION")) {
    if (stress_set %in% names(result$marker_sets)) {
      stress_score <- stress_score + result$marker_sets[[stress_set]]$avg_pct
      stress_markers <- c(stress_markers, result$marker_sets[[stress_set]]$genes)
    }
  }
  
  # Determine stress level
  stress_level <- if (stress_score > 100) "Very High" 
                  else if (stress_score > 50) "High"
                  else if (stress_score > 25) "Medium"
                  else "Low"
  
  # Display results
  if (length(result$marker_sets) > 0) {
    # Show stress-related markers first
    stress_sets <- intersect(c("STRESS_GENERAL", "HEAT_SHOCK", "NEURONAL_STRESS", 
                              "AUTOPHAGY_CORE", "MITO_DYNAMICS", "OXIDATIVE_DAMAGE",
                              "ER_STRESS", "INFLAMMATION"),
                            names(result$marker_sets))
    
    if (length(stress_sets) > 0) {
      cat("\nSTRESS-RELATED MARKERS:\n")
      for (set_name in stress_sets) {
        set_data <- result$marker_sets[[set_name]]
        cat(sprintf("  %s: %d/%d genes (%.1f%%), avg %.1f%% cells\n",
                    set_name, set_data$n_genes, set_data$total_genes,
                    set_data$pct_coverage, set_data$avg_pct))
        cat(sprintf("    Genes: %s\n", set_data$genes))
      }
    }
    
    # Show other top markers
    other_sets <- setdiff(names(result$marker_sets), stress_sets)
    if (length(other_sets) > 0) {
      sorted_sets <- result$marker_sets[other_sets]
      sorted_sets <- sorted_sets[order(sapply(sorted_sets, function(x) x$pct_coverage), decreasing = TRUE)]
      
      cat("\nOTHER TOP MARKERS:\n")
      for (i in 1:min(3, length(sorted_sets))) {
        set_name <- names(sorted_sets)[i]
        set_data <- sorted_sets[[set_name]]
        cat(sprintf("  %s: %d/%d genes (%.1f%%), avg %.1f%% cells\n",
                    set_name, set_data$n_genes, set_data$total_genes,
                    set_data$pct_coverage, set_data$avg_pct))
      }
    }
  }
  
  # Assign identity
  identity <- "Unknown"
  confidence <- "Low"
  key_markers <- ""
  
  if (length(result$marker_sets) > 0) {
    # Identity logic
    if ("DA_CORE" %in% names(result$marker_sets) && result$marker_sets$DA_CORE$n_genes >= 3) {
      identity <- "Dopaminergic neurons"
      confidence <- "High"
      key_markers <- result$marker_sets$DA_CORE$genes
    } else if ("PROLIFERATING" %in% names(result$marker_sets) && result$marker_sets$PROLIFERATING$pct_coverage > 50) {
      identity <- "Proliferating cells"
      confidence <- "High"
      key_markers <- result$marker_sets$PROLIFERATING$genes
    } else if ("PAN_NEURONAL" %in% names(result$marker_sets) && result$marker_sets$PAN_NEURONAL$pct_coverage > 30) {
      if ("GLUTAMATERGIC" %in% names(result$marker_sets)) {
        identity <- "Glutamatergic neurons"
      } else if ("GABAERGIC" %in% names(result$marker_sets)) {
        identity <- "GABAergic neurons"
      } else {
        identity <- "Neurons"
      }
      confidence <- "Medium"
      key_markers <- result$marker_sets$PAN_NEURONAL$genes
    } else if ("NEURAL_PROGENITOR" %in% names(result$marker_sets)) {
      identity <- "Neural progenitors"
      confidence <- "Medium"
      key_markers <- result$marker_sets$NEURAL_PROGENITOR$genes
    }
    
    # Append stress status to identity
    if (stress_level %in% c("High", "Very High")) {
      identity <- paste0(identity, " (stressed)")
    }
  }
  
  cat(sprintf("\nIDENTITY: %s (%s confidence)\n", identity, confidence))
  cat(sprintf("STRESS LEVEL: %s\n", stress_level))
  
  coarse_identities <- rbind(coarse_identities, data.frame(
    cluster = cl,
    n_cells = result$n_cells,
    identity = identity,
    confidence = confidence,
    stress_level = stress_level,
    key_markers = substr(key_markers, 1, 100),
    stress_markers = substr(paste(unique(stress_markers), collapse = ", "), 1, 100),
    stringsAsFactors = FALSE
  ))
}

# Special analysis: TTR expression
cat("\n\n3. Special Analysis: TTR as stress marker\n")
cat("==========================================\n")
ttr_analysis <- analyze_cluster_expression(seurat_obj, 0, "seurat_clusters_coarse")
# Check TTR specifically
expr_data <- GetAssayData(seurat_obj, slot = "data")
if ("TTR" %in% rownames(expr_data)) {
  for (cl in sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_coarse))))) {
    cells <- which(seurat_obj$seurat_clusters_coarse == as.character(cl))
    ttr_pct <- 100 * sum(expr_data["TTR", cells] > 0) / length(cells)
    cat(sprintf("Cluster %d: TTR in %.1f%% cells\n", cl, ttr_pct))
  }
}

# 4. Create coarse-to-fine mapping
cat("\n\n4. Creating coarse-to-fine cluster mapping...\n")
cat("============================================\n")

mapping <- table(
  Fine = seurat_obj$seurat_clusters_fine,
  Coarse = seurat_obj$seurat_clusters_coarse
)

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

# 5. Save all results
cat("\n5. Saving results...\n")
dir.create("results/reclustered_analysis", recursive = TRUE, showWarnings = FALSE)

saveRDS(coarse_results, "results/reclustered_analysis/coarse_cluster_expression_with_stress.rds")
write.csv(coarse_identities, "results/reclustered_analysis/coarse_cluster_identities_with_stress.csv", row.names = FALSE)
write.csv(fine_to_coarse, "results/reclustered_analysis/fine_to_coarse_mapping.csv", row.names = FALSE)

# Summary report focusing on stress
cat("\n\n=== STRESS ANALYSIS SUMMARY ===\n")
cat("================================\n\n")

stress_summary <- coarse_identities %>%
  group_by(stress_level) %>%
  summarise(
    n_clusters = n(),
    total_cells = sum(n_cells),
    pct_cells = round(100 * total_cells / sum(coarse_identities$n_cells), 1)
  )

cat("Stress level distribution:\n")
print(stress_summary)

cat("\n\nHigh/Very High stress clusters:\n")
stressed_clusters <- coarse_identities %>%
  filter(stress_level %in% c("High", "Very High")) %>%
  select(cluster, identity, stress_level, stress_markers)
print(stressed_clusters)

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("Results saved to:\n")
cat("- results/reclustered_analysis/coarse_cluster_expression_with_stress.rds\n")
cat("- results/reclustered_analysis/coarse_cluster_identities_with_stress.csv\n")
cat("- results/reclustered_analysis/fine_to_coarse_mapping.csv\n")