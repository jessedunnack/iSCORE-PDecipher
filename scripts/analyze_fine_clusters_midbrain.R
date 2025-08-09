#!/usr/bin/env Rscript

# Fine Cluster Analysis for Midbrain Floorplate-Derived Cells
# Focuses on dopaminergic, midbrain, and hypothalamic cell types

library(dplyr)
library(tidyr)

# Define key marker sets based on midbrain development research
MIDBRAIN_MARKERS <- list(
  # Floor plate progenitors
  floorplate_prog = c("FOXA2", "LMX1A", "LMX1B", "OTX2", "EN1", "EN2", "SHH", "CORIN"),
  
  # Dopaminergic neuron markers
  dopaminergic = c("TH", "DDC", "SLC6A3", "SLC18A2", "KCNJ6", "CALB1", "ALDH1A1", 
                   "NR4A2", "PITX3", "DRD2", "RET", "SOX6", "NEUROD6", "PBX1", "LMO3"),
  
  # Dopaminergic neuron subtypes
  da_a9 = c("KCNJ6", "ALDH1A1", "SOX6", "CALB1"),  # Substantia nigra
  da_a10 = c("CALB2", "SLC17A6", "OTX2", "PAX5", "PAX7"),  # VTA
  da_a8 = c("RET", "NR4A2"),  # Retrorubral field
  
  # Midbrain progenitors
  midbrain_prog = c("SOX2", "NES", "PAX6", "MSX1", "MSX2", "ASCL1", "NGN2", "MASH1"),
  
  # Neuroblast/immature neurons
  neuroblast = c("DCX", "TUBB3", "STMN2", "NCAM1", "GAP43", "MAP2"),
  
  # Mature neuron markers
  mature_neuron = c("SYN1", "SYP", "SNAP25", "SYT1", "RBFOX3", "NEUN", "MAP2"),
  
  # Hypothalamic markers
  hypothalamic = c("OTP", "SIM1", "SIM2", "POU3F2", "NKX2.1", "NKX2.2", "RAX", 
                   "POMC", "NPY", "AGRP", "HCRT", "AVP", "OXT", "CRH", "TRH", "GHRH"),
  
  # Diencephalic markers
  diencephalic = c("PAX6", "TCF7L2", "LHX2", "LHX5", "LHX9", "FOXD1", "FOXG1"),
  
  # Glial markers
  astrocyte = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "SLC1A3", "SOX9"),
  oligodendrocyte = c("OLIG1", "OLIG2", "MBP", "PLP1", "MOG", "MAG", "CNP", "SOX10"),
  opc = c("PDGFRA", "CSPG4", "SOX10", "OLIG2", "NKX2.2"),
  
  # Radial glia
  radial_glia = c("VIM", "HES1", "HES5", "PAX6", "SOX2", "FABP7", "SLC1A3"),
  
  # Cell cycle/proliferation
  proliferation = c("MKI67", "TOP2A", "PCNA", "MCM2", "CDK1", "CCNB1", "CCNB2"),
  
  # Stress/culture artifact
  stress = c("FOS", "JUN", "EGR1", "HSPA1A", "HSPA1B", "HSP90AA1", "DDIT3", "ATF4")
)

# Function to analyze markers in fine clusters
analyze_fine_clusters <- function(marker_file) {
  cat("\n=== Analyzing Fine Clusters (36) for Midbrain Cell Types ===\n")
  
  # Load markers
  markers <- readRDS(marker_file)
  
  # Get all clusters
  clusters <- sort(unique(as.character(markers$cluster)))
  cat("Found", length(clusters), "fine clusters\n\n")
  
  # Initialize results
  cluster_results <- list()
  
  for (cluster_id in clusters) {
    # Get top markers for this cluster
    cluster_markers <- markers %>%
      filter(cluster == cluster_id, avg_log2FC > 0.25, p_val_adj < 0.05) %>%
      mutate(score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)) %>%
      arrange(desc(score)) %>%
      slice_head(n = 50)  # Get more markers for thorough analysis
    
    # Check against all marker sets
    marker_hits <- list()
    
    for (marker_set in names(MIDBRAIN_MARKERS)) {
      hits <- intersect(cluster_markers$gene, MIDBRAIN_MARKERS[[marker_set]])
      if (length(hits) > 0) {
        marker_hits[[marker_set]] <- hits
      }
    }
    
    # Store results
    cluster_results[[cluster_id]] <- list(
      top_markers = cluster_markers,
      marker_hits = marker_hits,
      top_10_genes = head(cluster_markers$gene, 10)
    )
  }
  
  return(cluster_results)
}

# Function to infer cell type based on marker patterns
infer_midbrain_cell_type <- function(marker_hits, top_genes) {
  inference <- list(
    primary_type = "Unknown",
    subtype = NA,
    confidence = "Low",
    evidence = list()
  )
  
  # Priority order for cell type assignment
  
  # Check for dopaminergic neurons first
  if ("dopaminergic" %in% names(marker_hits) && length(marker_hits$dopaminergic) >= 2) {
    inference$primary_type <- "Dopaminergic Neuron"
    inference$evidence$dopaminergic <- marker_hits$dopaminergic
    
    # Check DA subtypes
    if ("da_a9" %in% names(marker_hits) && "KCNJ6" %in% marker_hits$da_a9) {
      inference$subtype <- "A9 (SNc)"
      inference$evidence$a9 <- marker_hits$da_a9
    } else if ("da_a10" %in% names(marker_hits) && any(c("CALB2", "OTX2") %in% marker_hits$da_a10)) {
      inference$subtype <- "A10 (VTA)"
      inference$evidence$a10 <- marker_hits$da_a10
    }
    
    inference$confidence <- ifelse(length(marker_hits$dopaminergic) >= 3, "High", "Medium")
  }
  
  # Check for floor plate progenitors
  else if ("floorplate_prog" %in% names(marker_hits) && 
           any(c("FOXA2", "LMX1A") %in% marker_hits$floorplate_prog)) {
    inference$primary_type <- "Floor Plate Progenitor"
    inference$evidence$floorplate <- marker_hits$floorplate_prog
    inference$confidence <- ifelse(length(marker_hits$floorplate_prog) >= 2, "High", "Medium")
  }
  
  # Check for hypothalamic cells
  else if ("hypothalamic" %in% names(marker_hits) && length(marker_hits$hypothalamic) >= 2) {
    inference$primary_type <- "Hypothalamic"
    inference$evidence$hypothalamic <- marker_hits$hypothalamic
    
    # Check specific hypothalamic subtypes
    if ("POMC" %in% marker_hits$hypothalamic) {
      inference$subtype <- "POMC neurons"
    } else if ("HCRT" %in% marker_hits$hypothalamic) {
      inference$subtype <- "Orexin/Hypocretin neurons"
    } else if ("AVP" %in% marker_hits$hypothalamic || "OXT" %in% marker_hits$hypothalamic) {
      inference$subtype <- "Magnocellular neurons"
    }
    
    inference$confidence <- "High"
  }
  
  # Check for proliferating cells
  else if ("proliferation" %in% names(marker_hits) && length(marker_hits$proliferation) >= 3) {
    inference$primary_type <- "Proliferating"
    inference$evidence$proliferation <- marker_hits$proliferation
    
    # Check if proliferating progenitor
    if ("midbrain_prog" %in% names(marker_hits)) {
      inference$subtype <- "Midbrain Progenitor"
      inference$evidence$progenitor <- marker_hits$midbrain_prog
    }
    
    inference$confidence <- "High"
  }
  
  # Check for immature neurons
  else if ("neuroblast" %in% names(marker_hits) && length(marker_hits$neuroblast) >= 2) {
    inference$primary_type <- "Immature Neuron"
    inference$evidence$neuroblast <- marker_hits$neuroblast
    inference$confidence <- "Medium"
  }
  
  # Check for glial cells
  else if ("oligodendrocyte" %in% names(marker_hits) && length(marker_hits$oligodendrocyte) >= 2) {
    inference$primary_type <- "Oligodendrocyte"
    inference$evidence$oligodendrocyte <- marker_hits$oligodendrocyte
    inference$confidence <- "High"
  }
  else if ("astrocyte" %in% names(marker_hits) && length(marker_hits$astrocyte) >= 2) {
    inference$primary_type <- "Astrocyte"
    inference$evidence$astrocyte <- marker_hits$astrocyte
    inference$confidence <- "High"
  }
  else if ("opc" %in% names(marker_hits) && length(marker_hits$opc) >= 2) {
    inference$primary_type <- "OPC"
    inference$evidence$opc <- marker_hits$opc
    inference$confidence <- "High"
  }
  
  # Check for radial glia
  else if ("radial_glia" %in% names(marker_hits) && length(marker_hits$radial_glia) >= 2) {
    inference$primary_type <- "Radial Glia"
    inference$evidence$radial_glia <- marker_hits$radial_glia
    inference$confidence <- "Medium"
  }
  
  # Check for stress markers
  if ("stress" %in% names(marker_hits) && length(marker_hits$stress) >= 3) {
    inference$evidence$stress <- marker_hits$stress
    if (inference$primary_type == "Unknown") {
      inference$primary_type <- "Stressed/Artifact"
      inference$confidence <- "Medium"
    }
  }
  
  return(inference)
}

# Function to create summary with midbrain focus
create_midbrain_summary <- function(cluster_results) {
  summary_df <- data.frame(
    fine_cluster = character(),
    primary_type = character(),
    subtype = character(),
    confidence = character(),
    top_5_markers = character(),
    key_evidence = character(),
    da_markers = character(),
    fp_markers = character(),
    hypo_markers = character(),
    stringsAsFactors = FALSE
  )
  
  for (cluster_id in names(cluster_results)) {
    result <- cluster_results[[cluster_id]]
    
    # Infer cell type
    inference <- infer_midbrain_cell_type(result$marker_hits, result$top_10_genes)
    
    # Extract specific marker evidence
    da_markers <- if ("dopaminergic" %in% names(result$marker_hits)) {
      paste(result$marker_hits$dopaminergic, collapse = ";")
    } else { "" }
    
    fp_markers <- if ("floorplate_prog" %in% names(result$marker_hits)) {
      paste(result$marker_hits$floorplate_prog, collapse = ";")
    } else { "" }
    
    hypo_markers <- if ("hypothalamic" %in% names(result$marker_hits)) {
      paste(result$marker_hits$hypothalamic, collapse = ";")
    } else { "" }
    
    summary_df <- rbind(summary_df, data.frame(
      fine_cluster = cluster_id,
      primary_type = inference$primary_type,
      subtype = ifelse(is.na(inference$subtype), "", inference$subtype),
      confidence = inference$confidence,
      top_5_markers = paste(head(result$top_10_genes, 5), collapse = ", "),
      key_evidence = paste(names(inference$evidence), collapse = "; "),
      da_markers = da_markers,
      fp_markers = fp_markers,
      hypo_markers = hypo_markers,
      stringsAsFactors = FALSE
    ))
  }
  
  return(summary_df)
}

# Main analysis function
main <- function() {
  cat("Fine Cluster Analysis for Midbrain Floorplate Cells\n")
  cat("==================================================\n\n")
  
  # Analyze fine clusters
  fine_marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
  
  if (!file.exists(fine_marker_file)) {
    stop("Fine cluster marker file not found: ", fine_marker_file)
  }
  
  # Run analysis
  cluster_results <- analyze_fine_clusters(fine_marker_file)
  
  # Create summary
  summary_df <- create_midbrain_summary(cluster_results)
  
  # Save results
  dir.create("results/cell_type_annotations", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(cluster_results, "results/cell_type_annotations/fine_clusters_midbrain_analysis.rds")
  write.csv(summary_df, "results/cell_type_annotations/fine_clusters_midbrain_summary.csv", 
            row.names = FALSE)
  
  # Print summary statistics
  cat("\n=== Cell Type Distribution ===\n")
  type_counts <- table(summary_df$primary_type)
  for (type in names(sort(type_counts, decreasing = TRUE))) {
    cat(sprintf("%-25s: %d clusters\n", type, type_counts[type]))
  }
  
  # Print dopaminergic clusters
  da_clusters <- summary_df[summary_df$primary_type == "Dopaminergic Neuron", ]
  if (nrow(da_clusters) > 0) {
    cat("\n=== Dopaminergic Neuron Clusters ===\n")
    print(da_clusters[, c("fine_cluster", "subtype", "da_markers")])
  }
  
  # Print floor plate clusters
  fp_clusters <- summary_df[summary_df$primary_type == "Floor Plate Progenitor", ]
  if (nrow(fp_clusters) > 0) {
    cat("\n=== Floor Plate Progenitor Clusters ===\n")
    print(fp_clusters[, c("fine_cluster", "fp_markers")])
  }
  
  # Print hypothalamic clusters
  hypo_clusters <- summary_df[summary_df$primary_type == "Hypothalamic", ]
  if (nrow(hypo_clusters) > 0) {
    cat("\n=== Hypothalamic Clusters ===\n")
    print(hypo_clusters[, c("fine_cluster", "subtype", "hypo_markers")])
  }
  
  cat("\nAnalysis complete. Results saved to results/cell_type_annotations/\n")
  
  return(list(
    results = cluster_results,
    summary = summary_df
  ))
}

# Run if called directly
if (!interactive()) {
  main()
}