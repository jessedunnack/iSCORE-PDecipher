#!/usr/bin/env Rscript

# COMPREHENSIVE FINE CLUSTER ANNOTATION - FIXED VERSION
# Handles cases where fine clusters have no markers in the marker file

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("COMPREHENSIVE FINE CLUSTER ANNOTATION (FIXED)\n")
cat("=================================================================\n\n")

# 1. Load all necessary data
cat("1. Loading data files...\n")

# Load annotated Seurat object (with coarse cell types)
if (file.exists("results/seurat_obj_annotated_coarse.rds")) {
  seurat_obj <- readRDS("results/seurat_obj_annotated_coarse.rds")
  cat("  - Loaded annotated Seurat object\n")
} else {
  seurat_obj <- readRDS("results/seurat_obj_reclustered.rds")
  cat("  - Loaded reclustered Seurat object (coarse annotations may be missing)\n")
}

# Load fine-to-coarse mapping
fine_to_coarse <- read.csv("results/reclustered_analysis/fine_to_coarse_mapping.csv")
cat("  - Loaded fine-to-coarse mapping\n")

# Load coarse identities
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")
cat("  - Loaded coarse cluster identities\n")

# Load fine cluster markers
fine_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_fine.rds")
cat("  - Loaded fine cluster markers\n")

DefaultAssay(seurat_obj) <- "SCT"

# 2. Define subtype markers for detailed annotation
cat("\n2. Setting up subtype marker sets...\n")

SUBTYPE_MARKERS <- list(
  # Dopaminergic neuron subtypes
  DA_A9_SNc = c("ALDH1A1", "SOX6", "KCNJ6", "GIRK2", "SLC18A2", "SNCA", "NR4A2", "LMO3"),
  DA_A10_VTA = c("CALB1", "CALB2", "OTX2", "GRP", "CCK", "VIP", "NEUROD6", "LPL"),
  DA_Immature = c("TH", "DDC", "FOXA2", "LMX1A", "EN1", "EN2", "PITX3"),
  
  # Progenitor stages
  PROG_Proliferating = c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C", "CENPF"),
  PROG_Early_Neural = c("SOX2", "PAX6", "NES", "HES1", "HES5", "NOTCH1"),
  PROG_Intermediate = c("EOMES", "NEUROG2", "NEUROD1", "BTG2", "INSM1"),
  PROG_Late_Neuronal = c("DCX", "NEUROD2", "NEUROD4", "MAPT", "STMN2"),
  PROG_Gliogenic = c("OLIG1", "OLIG2", "NKX2-2", "SOX10", "PDGFRA"),
  
  # Stress/damage subtypes
  STRESS_Oxidative = c("NQO1", "HMOX1", "GPX1", "GPX4", "SOD1", "SOD2", "GCLC"),
  STRESS_ER_UPR = c("HSPA5", "XBP1", "DDIT3", "ATF4", "ATF6", "ERN1", "EIF2AK3"),
  STRESS_Inflammatory = c("IL1B", "IL6", "TNF", "NFKB1", "CCL2", "CXCL10", "CXCL12"),
  STRESS_Apoptotic = c("CASP3", "CASP7", "BAX", "BCL2", "PARP1", "TP53"),
  
  # Mesenchymal subtypes
  MESEN_Fibroblast = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "VIM", "FN1"),
  MESEN_Pericyte = c("PDGFRB", "RGS5", "ABCC9", "KCNJ8", "NOTCH3", "CSPG4"),
  MESEN_Smooth_Muscle = c("TAGLN", "ACTA2", "MYL9", "CNN1", "MYLK", "MYOCD"),
  
  # Neuroendocrine subtypes
  NEUROENDO_Serotonergic = c("TPH1", "TPH2", "DDC", "SLC6A4", "FEV", "GATA2"),
  NEUROENDO_Peptidergic = c("SST", "VIP", "NPY", "CALCA", "CALCB", "TAC1"),
  
  # Choroid plexus subtypes
  CHOROID_Epithelial = c("TTR", "FOLR1", "CLIC6", "AQP1", "KCNJ13", "SLC4A10"),
  CHOROID_Secretory = c("CP", "ENPP2", "MSX1", "MSX2", "OTX2", "LMX1A"),
  
  # Regional/functional markers
  HYPOTHALAMIC = c("HCRT", "OXT", "AVP", "CRH", "TRH", "GHRH", "POMC", "NPY", "AGRP"),
  MIDBRAIN = c("EN1", "EN2", "PAX2", "PAX5", "WNT1", "FGF8", "LMX1A", "LMX1B"),
  CORTICAL = c("TBR1", "SATB2", "CUX1", "CUX2", "FOXP2", "BCL11B", "RORB"),
  
  # Additional specific markers
  RETINAL = c("RPE65", "BEST1", "RLBP1", "RDH5", "CRALBP", "MITF"),
  IMMUNE = c("PTPRC", "CD74", "HLA-DRA", "AIF1", "CX3CR1", "TMEM119")
)

# 3. Function to find top variable genes when no markers available
find_top_variable_genes <- function(seurat_obj, cells_in_cluster, n_genes = 20) {
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Calculate variance for each gene in this cluster
  cluster_expr <- expr_data[, cells_in_cluster]
  gene_vars <- apply(cluster_expr, 1, var)
  
  # Get top variable genes
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:n_genes]
  
  # Calculate stats for these genes
  gene_stats <- data.frame(
    gene = top_genes,
    mean_expr = apply(cluster_expr[top_genes, ], 1, mean),
    pct_expressing = 100 * apply(cluster_expr[top_genes, ] > 0, 1, mean),
    variance = gene_vars[top_genes]
  )
  
  return(gene_stats)
}

# 4. Fixed analysis function
analyze_fine_cluster_comprehensive <- function(fine_cl, seurat_obj, fine_markers_df, 
                                               coarse_identity, coarse_cluster) {
  
  cat(sprintf("\n=== FINE CLUSTER %d (Coarse: %d - %s) ===\n", 
              fine_cl, coarse_cluster, coarse_identity))
  
  # Get cells in fine cluster
  cells_in_fine <- which(seurat_obj$seurat_clusters_fine == as.character(fine_cl))
  n_cells <- length(cells_in_fine)
  cat(sprintf("Cells: %d\n", n_cells))
  
  # Get fine cluster markers
  cl_markers <- fine_markers_df %>%
    filter(cluster == as.character(fine_cl)) %>%
    arrange(desc(avg_log2FC))
  
  # Handle case where no markers found
  if (nrow(cl_markers) == 0) {
    cat("WARNING: No markers found in marker file. Using variable gene analysis...\n")
    
    # Find top variable genes in this cluster
    var_genes <- find_top_variable_genes(seurat_obj, cells_in_fine, n_genes = 15)
    
    cat("\nTop 15 variable genes:\n")
    for (i in 1:nrow(var_genes)) {
      cat(sprintf("  %2d. %-12s (mean: %.2f, %.1f%% cells)\n", 
                  i, var_genes$gene[i], var_genes$mean_expr[i], var_genes$pct_expressing[i]))
    }
    
    # Use variable genes as pseudo-markers
    cl_markers <- data.frame(
      gene = var_genes$gene,
      avg_log2FC = log2(var_genes$mean_expr + 1),
      pct.1 = var_genes$pct_expressing / 100,
      stringsAsFactors = FALSE
    )
  } else {
    # Display top 15 markers
    cat("\nTop 15 markers:\n")
    top15 <- head(cl_markers, 15)
    for (i in 1:nrow(top15)) {
      cat(sprintf("  %2d. %-12s (%.2f FC, %.1f%%)\n", 
                  i, top15$gene[i], top15$avg_log2FC[i], top15$pct.1[i] * 100))
    }
  }
  
  # Check expression of subtype markers
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  cat("\nSubtype marker analysis:\n")
  subtype_scores <- list()
  
  for (set_name in names(SUBTYPE_MARKERS)) {
    genes <- SUBTYPE_MARKERS[[set_name]]
    genes_present <- genes[genes %in% rownames(expr_data)]
    
    if (length(genes_present) > 0) {
      # Calculate expression in this fine cluster
      expr_stats <- sapply(genes_present, function(g) {
        gene_expr <- expr_data[g, cells_in_fine]
        c(mean = mean(gene_expr), 
          pct = 100 * sum(gene_expr > 0) / n_cells)
      })
      
      # Score based on percentage of cells expressing
      score <- mean(expr_stats["pct", ])
      
      if (score > 10) {  # At least 10% cells expressing on average
        subtype_scores[[set_name]] <- score
        
        # Report significant markers
        sig_genes <- names(which(expr_stats["pct", ] > 20))
        if (length(sig_genes) > 0) {
          cat(sprintf("  %s: %.1f%% score (%s)\n", 
                      set_name, score, paste(sig_genes, collapse=", ")))
        }
      }
    }
  }
  
  # Determine fine cluster identity
  fine_identity <- coarse_identity  # Default to coarse identity
  confidence <- "Low"
  subtype_info <- ""
  
  # Logic for assigning fine identities based on coarse type
  if (grepl("Neurons_Dopaminergic", coarse_identity)) {
    # Check DA subtypes
    if ("DA_A9_SNc" %in% names(subtype_scores) && subtype_scores[["DA_A9_SNc"]] > 25) {
      fine_identity <- "Neurons_DA_A9-like"
      confidence <- "High"
      subtype_info <- "SNc-like"
    } else if ("DA_A10_VTA" %in% names(subtype_scores) && subtype_scores[["DA_A10_VTA"]] > 25) {
      fine_identity <- "Neurons_DA_A10-like"
      confidence <- "High"
      subtype_info <- "VTA-like"
    } else if ("DA_Immature" %in% names(subtype_scores) && subtype_scores[["DA_Immature"]] > 30) {
      fine_identity <- "Neurons_DA_Immature"
      confidence <- "Medium"
      subtype_info <- "Developing"
    } else {
      fine_identity <- "Neurons_DA_Unspecified"
      confidence <- "Medium"
    }
    
  } else if (grepl("Progenitors", coarse_identity)) {
    # Check progenitor stages
    if ("PROG_Proliferating" %in% names(subtype_scores) && subtype_scores[["PROG_Proliferating"]] > 40) {
      fine_identity <- "Progenitors_Proliferating"
      confidence <- "High"
      subtype_info <- "Cycling"
    } else if ("PROG_Early_Neural" %in% names(subtype_scores) && subtype_scores[["PROG_Early_Neural"]] > 30) {
      fine_identity <- "Progenitors_Early_Neural"
      confidence <- "High"
      subtype_info <- "Neural stem"
    } else if ("PROG_Intermediate" %in% names(subtype_scores) && subtype_scores[["PROG_Intermediate"]] > 25) {
      fine_identity <- "Progenitors_Intermediate"
      confidence <- "High"
      subtype_info <- "Neurogenic"
    } else if ("PROG_Late_Neuronal" %in% names(subtype_scores) && subtype_scores[["PROG_Late_Neuronal"]] > 25) {
      fine_identity <- "Progenitors_Late_Neuronal"
      confidence <- "Medium"
      subtype_info <- "Differentiating"
    } else if ("PROG_Gliogenic" %in% names(subtype_scores) && subtype_scores[["PROG_Gliogenic"]] > 20) {
      fine_identity <- "Progenitors_Gliogenic"
      confidence <- "Medium"
      subtype_info <- "Glial fate"
    } else {
      # Use top marker to refine if available
      if (nrow(cl_markers) > 0) {
        top_gene <- cl_markers$gene[1]
        fine_identity <- paste0("Progenitors_", top_gene, "+")
      } else {
        fine_identity <- paste0(coarse_identity, "_subset")
      }
      confidence <- "Low"
    }
    
  } else if (grepl("Mesenchymal", coarse_identity)) {
    # Check mesenchymal subtypes
    if ("MESEN_Pericyte" %in% names(subtype_scores) && subtype_scores[["MESEN_Pericyte"]] > 20) {
      fine_identity <- "Mesenchymal_Pericytes"
      confidence <- "High"
    } else if ("MESEN_Smooth_Muscle" %in% names(subtype_scores) && subtype_scores[["MESEN_Smooth_Muscle"]] > 20) {
      fine_identity <- "Mesenchymal_Smooth_Muscle"
      confidence <- "High"
    } else {
      fine_identity <- "Mesenchymal_Fibroblasts"
      confidence <- "Medium"
    }
    
  } else if (grepl("Stressed", coarse_identity)) {
    # Check stress types
    max_stress <- ""
    max_score <- 0
    
    stress_types <- grep("^STRESS_", names(subtype_scores), value = TRUE)
    if (length(stress_types) > 0) {
      max_idx <- which.max(unlist(subtype_scores[stress_types]))
      max_stress <- stress_types[max_idx]
      max_score <- subtype_scores[[max_stress]]
    }
    
    if (max_score > 20) {
      stress_type <- gsub("STRESS_", "", max_stress)
      fine_identity <- paste0("Stressed_", stress_type)
      confidence <- "High"
    } else {
      fine_identity <- "Stressed_Mixed"
      confidence <- "Medium"
    }
    
  } else if (grepl("Neuroendocrine", coarse_identity)) {
    # Check neuroendocrine subtypes
    if ("NEUROENDO_Serotonergic" %in% names(subtype_scores) && subtype_scores[["NEUROENDO_Serotonergic"]] > 20) {
      fine_identity <- "Neuroendocrine_Serotonergic"
      confidence <- "High"
    } else if ("NEUROENDO_Peptidergic" %in% names(subtype_scores) && subtype_scores[["NEUROENDO_Peptidergic"]] > 20) {
      fine_identity <- "Neuroendocrine_Peptidergic"
      confidence <- "High"
    } else {
      fine_identity <- "Neuroendocrine_Mixed"
      confidence <- "Medium"
    }
    
  } else if (grepl("Choroid", coarse_identity)) {
    # Choroid plexus subtypes
    if ("CHOROID_Secretory" %in% names(subtype_scores) && subtype_scores[["CHOROID_Secretory"]] > 20) {
      fine_identity <- "Choroid_Plexus_Secretory"
      confidence <- "High"
    } else {
      fine_identity <- "Choroid_Plexus_Epithelial"
      confidence <- "High"
    }
    
  } else {
    # For unknown/unidentified clusters, use top markers if available
    if (nrow(cl_markers) > 2) {
      top_genes <- head(cl_markers$gene, 3)
      fine_identity <- paste0(coarse_identity, "_", paste(top_genes, collapse="/"))
    } else {
      fine_identity <- paste0(coarse_identity, "_subset", fine_cl)
    }
    confidence <- "Low"
  }
  
  # Get key markers for this fine cluster
  if (nrow(cl_markers) > 0) {
    key_markers <- paste(head(cl_markers$gene, 5), collapse=", ")
  } else {
    key_markers <- "No specific markers"
  }
  
  cat(sprintf("\nASSIGNED IDENTITY: %s (%s confidence)\n", fine_identity, confidence))
  if (subtype_info != "") {
    cat(sprintf("Subtype info: %s\n", subtype_info))
  }
  
  return(list(
    fine_cluster = fine_cl,
    identity = fine_identity,
    confidence = confidence,
    key_markers = key_markers,
    n_cells = n_cells,
    subtype_scores = subtype_scores
  ))
}

# 5. Main analysis loop
cat("\n\n3. Analyzing all fine clusters...\n")
cat("=====================================\n")

fine_cluster_results <- list()
fine_identities_df <- data.frame(
  fine_cluster = integer(),
  coarse_cluster = integer(),
  coarse_identity = character(),
  fine_identity = character(),
  confidence = character(),
  n_cells = integer(),
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
  
  # Analyze fine cluster
  result <- analyze_fine_cluster_comprehensive(
    fine_cl, seurat_obj, fine_markers, 
    coarse_identity, coarse_cl
  )
  
  fine_cluster_results[[as.character(fine_cl)]] <- result
  
  # Add to dataframe - ensure all fields are present
  fine_identities_df <- rbind(fine_identities_df, data.frame(
    fine_cluster = result$fine_cluster,
    coarse_cluster = coarse_cl,
    coarse_identity = coarse_identity,
    fine_identity = result$identity,
    confidence = result$confidence,
    n_cells = result$n_cells,
    key_markers = substr(result$key_markers, 1, 100),
    stringsAsFactors = FALSE
  ))
}

# 6. Add fine cell types to Seurat object
cat("\n\n4. Adding cell_type_fine to Seurat object...\n")

seurat_obj$cell_type_fine <- plyr::mapvalues(
  seurat_obj$seurat_clusters_fine,
  from = as.character(fine_identities_df$fine_cluster),
  to = fine_identities_df$fine_identity
)

# 7. Create summary statistics
cat("\n5. Creating summary statistics...\n")

# Summary by fine identity
fine_summary <- as.data.frame(table(seurat_obj$cell_type_fine)) %>%
  arrange(desc(Freq))
colnames(fine_summary) <- c("Fine_Cell_Type", "n_cells")
fine_summary$pct <- round(100 * fine_summary$n_cells / ncol(seurat_obj), 1)

cat("\nFine cell type distribution (top 20):\n")
print(head(fine_summary, 20))

# 8. Save everything
cat("\n\n6. Saving results...\n")

# Save fine cluster identities
write.csv(fine_identities_df, "results/reclustered_analysis/fine_cluster_identities_FINAL.csv", 
          row.names = FALSE)

# Save detailed results
saveRDS(fine_cluster_results, "results/reclustered_analysis/fine_cluster_detailed_results.rds")

# Save final annotated Seurat object to main location only
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")

# Also save metadata as separate RDS file
metadata_df <- seurat_obj@meta.data
saveRDS(metadata_df, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_metadata_FINAL.rds")
cat("  - Saved metadata separately as iSCORE-PD_plus_CRISPRi_metadata_FINAL.rds\n")

# Create a summary report
cat("\n\n=== FINAL SUMMARY ===\n")
cat("====================\n")

# Count broad categories in fine clusters
fine_broad_categories <- fine_identities_df %>%
  mutate(broad_category = case_when(
    grepl("^Neurons", fine_identity) ~ "Neurons",
    grepl("^Progenitors", fine_identity) ~ "Progenitors",
    grepl("^Mesenchymal|^Fibroblast|^Choroid|^Neuroendocrine|^Stressed", fine_identity) ~ "Non-neuronal",
    TRUE ~ "Unknown"
  )) %>%
  group_by(broad_category) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    pct_cells = round(100 * total_cells / sum(fine_identities_df$n_cells), 1)
  )

cat("\nFine clusters by broad category:\n")
print(fine_broad_categories)

cat("\n\nHigh-confidence fine cluster assignments:\n")
high_conf <- fine_identities_df %>%
  filter(confidence == "High") %>%
  select(fine_cluster, fine_identity, n_cells)
print(head(high_conf, 20))

# Check for clusters with no markers
no_marker_clusters <- fine_identities_df %>%
  filter(key_markers == "No specific markers")

if (nrow(no_marker_clusters) > 0) {
  cat("\n\nFine clusters analyzed using variable genes (no markers in file):\n")
  print(no_marker_clusters[, c("fine_cluster", "coarse_identity", "fine_identity", "n_cells")])
}

cat("\n\n=== ANALYSIS COMPLETE ===\n")
cat("Key outputs:\n")
cat("- results/reclustered_analysis/fine_cluster_identities_FINAL.csv\n")
cat("- ../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds\n")
cat("- ../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_metadata_FINAL.rds\n")

cat("\nSeurat object now contains:\n")
cat("- celltypes_coarse: Coarse cluster cell types\n")
cat("- cell_type_broad: Broad categories\n")
cat("- cell_type_fine: Fine cluster cell types with subtypes\n")

cat("\nReady for downstream analyses!\n")
cat("\nREMINDER: Don't forget to run 'execute_cleanup.R' to organize files!\n")