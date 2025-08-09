#!/usr/bin/env Rscript

# FINE CLUSTER ANALYSIS FRAMEWORK
# Analyzes fine clusters within their coarse cluster context

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("FINE CLUSTER ANALYSIS WITHIN COARSE CONTEXT\n")
cat("=================================================================\n\n")

# 1. Load required data
cat("1. Loading data...\n")
seurat_obj <- readRDS("results/seurat_obj_reclustered.rds")
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_with_stress.csv")
fine_to_coarse <- read.csv("results/reclustered_analysis/fine_to_coarse_mapping.csv")

# We'll need the coarse expression results for reference
coarse_results <- readRDS("results/reclustered_analysis/coarse_cluster_expression_with_stress.rds")

DefaultAssay(seurat_obj) <- "SCT"

# 2. Define subtype markers for each major cell type
SUBTYPE_MARKERS <- list(
  # Dopaminergic subtypes
  DA_A9_SNc = c("ALDH1A1", "SOX6", "KCNJ6", "GIRK2", "SLC18A2", "SNCA", "NR4A2"),
  DA_A10_VTA = c("CALB1", "CALB2", "OTX2", "GRP", "CCK", "VIP", "NEUROD6"),
  DA_A8_RRF = c("CALB1", "LPL", "ADCYAP1", "LINGO1"),
  
  # Glutamatergic subtypes
  GLUT_CORTICAL_UPPER = c("SATB2", "CUX1", "CUX2", "RORB", "POU3F2"),
  GLUT_CORTICAL_DEEP = c("TBR1", "BCL11B", "FOXP2", "TLE4", "SOX5"),
  GLUT_HIPPOCAMPAL = c("PROX1", "C1QL2", "FIBCD1", "ZBTB20"),
  
  # GABAergic subtypes
  GABA_PV = c("PVALB", "SOX6", "TAC1", "ERBB4", "KCNS3"),
  GABA_SST = c("SST", "ELFN1", "NPY", "NOS1", "GRIK1"),
  GABA_VIP = c("VIP", "CALB2", "CCK", "TAC3", "PROX1"),
  GABA_LAMP5 = c("LAMP5", "PAX6", "EBF1", "TOX", "NR2F2"),
  
  # Progenitor states
  PROG_CYCLING = c("MKI67", "TOP2A", "PCNA", "CDC20", "CENPF"),
  PROG_EARLY = c("SOX2", "HES1", "HES5", "ID1", "ID3"),
  PROG_LATE = c("EOMES", "NEUROG2", "NEUROD1", "NEUROD2", "NEUROD4"),
  PROG_GLIOGENIC = c("OLIG1", "OLIG2", "NKX2-2", "PDGFRA", "SOX10"),
  
  # Stress/disease states
  STRESSED_OXIDATIVE = c("NQO1", "HMOX1", "GPX1", "GPX4", "SOD1", "SOD2"),
  STRESSED_ER = c("HSPA5", "XBP1", "DDIT3", "ATF4", "ATF6", "ERN1"),
  STRESSED_INFLAMMATORY = c("IL1B", "IL6", "TNF", "NFKB1", "CCL2", "CXCL10"),
  PD_PATHOLOGY = c("SNCA", "MAPT", "APP", "APOE", "CLU", "TREM2")
)

# 3. Function to analyze fine cluster within coarse context
analyze_fine_cluster <- function(seurat_obj, fine_cluster_id, coarse_cluster_id, 
                                coarse_identity, marker_sets) {
  
  # Get cells in this fine cluster
  cells_in_fine <- which(seurat_obj$seurat_clusters_fine == as.character(fine_cluster_id))
  n_cells <- length(cells_in_fine)
  
  # Get expression data
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Analyze relevant subtype markers based on coarse identity
  relevant_markers <- list()
  
  # Select relevant marker sets based on coarse identity
  if (grepl("Dopaminergic", coarse_identity, ignore.case = TRUE)) {
    relevant_markers <- c(relevant_markers, 
                         marker_sets[grep("^DA_", names(marker_sets))])
  }
  
  if (grepl("Glutamatergic", coarse_identity, ignore.case = TRUE)) {
    relevant_markers <- c(relevant_markers, 
                         marker_sets[grep("^GLUT_", names(marker_sets))])
  }
  
  if (grepl("GABAergic", coarse_identity, ignore.case = TRUE)) {
    relevant_markers <- c(relevant_markers, 
                         marker_sets[grep("^GABA_", names(marker_sets))])
  }
  
  if (grepl("progenitor", coarse_identity, ignore.case = TRUE)) {
    relevant_markers <- c(relevant_markers, 
                         marker_sets[grep("^PROG_", names(marker_sets))])
  }
  
  # Always check stress markers
  relevant_markers <- c(relevant_markers, 
                       marker_sets[grep("^STRESSED_", names(marker_sets))])
  
  # Analyze expression
  results <- list()
  for (set_name in names(relevant_markers)) {
    genes <- relevant_markers[[set_name]]
    genes_present <- genes[genes %in% rownames(expr_data)]
    
    if (length(genes_present) > 0) {
      expr_stats <- sapply(genes_present, function(g) {
        gene_expr <- expr_data[g, cells_in_fine]
        c(mean_expr = mean(gene_expr),
          pct_expr = 100 * sum(gene_expr > 0) / n_cells)
      })
      
      # Calculate summary stats
      pct_expressing <- expr_stats["pct_expr", ]
      sig_genes <- names(pct_expressing)[pct_expressing > 5]
      
      if (length(sig_genes) > 0) {
        results[[set_name]] <- list(
          n_genes = length(sig_genes),
          avg_pct = round(mean(pct_expressing[sig_genes]), 1),
          genes = paste(sig_genes, collapse = ", ")
        )
      }
    }
  }
  
  return(list(
    fine_cluster = fine_cluster_id,
    coarse_cluster = coarse_cluster_id,
    n_cells = n_cells,
    marker_results = results
  ))
}

# 4. Function to assign fine cluster identity
assign_fine_identity <- function(fine_analysis, coarse_identity) {
  
  # Default to coarse identity
  fine_identity <- coarse_identity
  confidence <- "Low"
  subtype <- ""
  
  # Check for specific subtypes based on marker expression
  if (grepl("Dopaminergic", coarse_identity)) {
    if ("DA_A9_SNc" %in% names(fine_analysis$marker_results)) {
      if (fine_analysis$marker_results$DA_A9_SNc$avg_pct > 30) {
        subtype <- "A9-like"
        confidence <- "High"
      }
    }
    if ("DA_A10_VTA" %in% names(fine_analysis$marker_results)) {
      if (fine_analysis$marker_results$DA_A10_VTA$avg_pct > 30) {
        subtype <- "A10-like"
        confidence <- "High"
      }
    }
  }
  
  if (grepl("Glutamatergic", coarse_identity)) {
    if ("GLUT_CORTICAL_UPPER" %in% names(fine_analysis$marker_results)) {
      if (fine_analysis$marker_results$GLUT_CORTICAL_UPPER$avg_pct > 25) {
        subtype <- "Upper-layer"
        confidence <- "Medium"
      }
    }
    if ("GLUT_CORTICAL_DEEP" %in% names(fine_analysis$marker_results)) {
      if (fine_analysis$marker_results$GLUT_CORTICAL_DEEP$avg_pct > 25) {
        subtype <- "Deep-layer"
        confidence <- "Medium"
      }
    }
  }
  
  # Check stress status
  stress_level <- "Normal"
  if ("STRESSED_OXIDATIVE" %in% names(fine_analysis$marker_results) ||
      "STRESSED_ER" %in% names(fine_analysis$marker_results)) {
    stress_markers_pct <- max(
      fine_analysis$marker_results$STRESSED_OXIDATIVE$avg_pct %||% 0,
      fine_analysis$marker_results$STRESSED_ER$avg_pct %||% 0
    )
    if (stress_markers_pct > 40) {
      stress_level <- "High"
    } else if (stress_markers_pct > 20) {
      stress_level <- "Medium"
    }
  }
  
  # Construct final identity
  if (subtype != "") {
    fine_identity <- paste0(coarse_identity, "_", subtype)
  }
  
  if (stress_level != "Normal") {
    fine_identity <- paste0(fine_identity, "_", stress_level, "Stress")
  }
  
  return(list(
    identity = fine_identity,
    confidence = confidence,
    subtype = subtype,
    stress_level = stress_level
  ))
}

# 5. Main analysis loop
cat("\n2. Analyzing fine clusters within coarse context...\n")
cat("==================================================\n")

fine_cluster_results <- data.frame(
  fine_cluster = integer(),
  coarse_cluster = integer(),
  coarse_identity = character(),
  n_cells = integer(),
  fine_identity = character(),
  confidence = character(),
  subtype = character(),
  stress_level = character(),
  key_markers = character(),
  stringsAsFactors = FALSE
)

# Process each fine cluster
for (i in 1:nrow(fine_to_coarse)) {
  fine_cl <- fine_to_coarse$fine_cluster[i]
  coarse_cl <- fine_to_coarse$coarse_cluster[i]
  
  # Get coarse identity
  coarse_info <- coarse_identities[coarse_identities$cluster == coarse_cl, ]
  coarse_identity <- coarse_info$identity[1]
  
  cat(sprintf("\nFine cluster %d (coarse: %d - %s):\n", 
              fine_cl, coarse_cl, coarse_identity))
  
  # Analyze fine cluster
  fine_analysis <- analyze_fine_cluster(
    seurat_obj, fine_cl, coarse_cl, 
    coarse_identity, SUBTYPE_MARKERS
  )
  
  # Assign identity
  identity_info <- assign_fine_identity(fine_analysis, coarse_identity)
  
  # Extract key markers
  key_markers <- ""
  if (length(fine_analysis$marker_results) > 0) {
    # Get top marker set
    top_set <- names(fine_analysis$marker_results)[1]
    key_markers <- fine_analysis$marker_results[[top_set]]$genes
  }
  
  # Add to results
  fine_cluster_results <- rbind(fine_cluster_results, data.frame(
    fine_cluster = fine_cl,
    coarse_cluster = coarse_cl,
    coarse_identity = coarse_identity,
    n_cells = fine_analysis$n_cells,
    fine_identity = identity_info$identity,
    confidence = identity_info$confidence,
    subtype = identity_info$subtype,
    stress_level = identity_info$stress_level,
    key_markers = substr(key_markers, 1, 100),
    stringsAsFactors = FALSE
  ))
  
  cat(sprintf("  Cells: %d\n", fine_analysis$n_cells))
  cat(sprintf("  Identity: %s (%s confidence)\n", 
              identity_info$identity, identity_info$confidence))
  
  # Show top markers
  if (length(fine_analysis$marker_results) > 0) {
    cat("  Top markers:\n")
    for (j in 1:min(3, length(fine_analysis$marker_results))) {
      set_name <- names(fine_analysis$marker_results)[j]
      set_data <- fine_analysis$marker_results[[set_name]]
      cat(sprintf("    %s: %d genes, %.1f%% cells\n", 
                  set_name, set_data$n_genes, set_data$avg_pct))
    }
  }
}

# 6. Save results
cat("\n\n3. Saving fine cluster analysis results...\n")
write.csv(fine_cluster_results, 
          "results/reclustered_analysis/fine_cluster_identities.csv", 
          row.names = FALSE)

# 7. Create summary
cat("\n\n=== FINE CLUSTER ANALYSIS SUMMARY ===\n")
cat("=====================================\n")

# Summary by coarse cluster
summary_by_coarse <- fine_cluster_results %>%
  group_by(coarse_cluster, coarse_identity) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    subtypes = paste(unique(subtype[subtype != ""]), collapse = ", "),
    stress_levels = paste(unique(stress_level[stress_level != "Normal"]), collapse = ", ")
  )

cat("\nFine clusters per coarse cluster:\n")
print(summary_by_coarse)

# Stress summary
stress_summary <- fine_cluster_results %>%
  filter(stress_level != "Normal") %>%
  group_by(coarse_identity, stress_level) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells)
  )

if (nrow(stress_summary) > 0) {
  cat("\n\nStressed fine clusters:\n")
  print(stress_summary)
}

cat("\n\nAnalysis complete!\n")
cat("Results saved to: results/reclustered_analysis/fine_cluster_identities.csv\n")
cat("\nNext steps:\n")
cat("1. Review fine cluster identities\n")
cat("2. Add celltypes_coarse and celltypes_fine to Seurat object\n")
cat("3. Visualize with updated cell type labels\n")