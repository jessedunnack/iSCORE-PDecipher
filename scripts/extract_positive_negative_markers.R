#!/usr/bin/env Rscript

# Extract both positive (highly expressed) and negative (lowly expressed) markers for each cluster
# This enables positive/negative marker logic for refined cell type classification

library(Seurat)
library(dplyr)
library(tidyr)

# Function to extract positive and negative markers
extract_pos_neg_markers <- function(seurat_obj, 
                                   n_top = 50,
                                   pos_fc_threshold = 0.5,
                                   neg_fc_threshold = -0.5,
                                   pval_threshold = 0.05) {
  
  cat("Extracting positive and negative markers for all clusters...\n")
  
  # Get all clusters
  clusters <- sort(unique(seurat_obj$seurat_clusters_fine))
  
  all_markers <- list()
  
  for (cluster_id in clusters) {
    cat("\nProcessing cluster", cluster_id, "...\n")
    
    # Find all markers for this cluster vs all others
    markers <- FindMarkers(seurat_obj, 
                          ident.1 = cluster_id,
                          ident.2 = NULL,
                          only.pos = FALSE,  # Important: get both positive and negative
                          min.pct = 0.1,
                          logfc.threshold = 0.25,
                          test.use = "wilcox")
    
    # Add gene names and cluster
    markers$gene <- rownames(markers)
    markers$cluster <- cluster_id
    
    # Classify as positive or negative
    markers <- markers %>%
      mutate(
        direction = case_when(
          avg_log2FC > pos_fc_threshold & p_val_adj < pval_threshold ~ "positive",
          avg_log2FC < neg_fc_threshold & p_val_adj < pval_threshold ~ "negative",
          TRUE ~ "neutral"
        ),
        score = abs(avg_log2FC) * -log10(p_val_adj + 1e-300)
      )
    
    # Get top positive markers
    pos_markers <- markers %>%
      filter(direction == "positive") %>%
      arrange(desc(score)) %>%
      head(n_top)
    
    # Get top negative markers (genes NOT expressed in this cluster)
    neg_markers <- markers %>%
      filter(direction == "negative") %>%
      arrange(desc(score)) %>%
      head(n_top)
    
    all_markers[[as.character(cluster_id)]] <- list(
      positive = pos_markers,
      negative = neg_markers,
      summary = data.frame(
        cluster = cluster_id,
        n_pos = nrow(pos_markers),
        n_neg = nrow(neg_markers),
        top_pos = paste(head(pos_markers$gene, 5), collapse = ", "),
        top_neg = paste(head(neg_markers$gene, 5), collapse = ", ")
      )
    )
  }
  
  return(all_markers)
}

# Function to identify vulnerability markers
check_vulnerability_markers <- function(pos_neg_markers) {
  
  # Define vulnerability markers from MPTP studies
  vulnerability_markers <- list(
    # Vulnerable dopaminergic markers
    vulnerable_da = list(
      positive = c("SOX6", "ALDH1A1", "KCNJ6", "GIRK2", "GRM5", "SNCG"),
      negative = c("CALB1", "CALB2", "OTX2", "FOXP2", "SORCS3")
    ),
    
    # Resilient dopaminergic markers
    resilient_da = list(
      positive = c("CALB1", "CALB2", "OTX2", "FOXP2", "SORCS3", "RELN", "TMEFF2"),
      negative = c("SOX6", "KCNJ6", "GIRK2")
    ),
    
    # Vulnerable glutamatergic markers
    vulnerable_glut = list(
      positive = c("MEIS1", "MEIS2", "SLC17A6", "SLC17A7"),
      negative = c("LMX1A", "FOXP2", "SORCS3")
    ),
    
    # Regional markers
    ventral_snc = list(
      positive = c("ALDH1A1", "ANXA1", "LMO3", "ALDH1A7"),
      negative = c("CALB1", "OTX2")
    ),
    
    dorsal_snc = list(
      positive = c("CALB1", "SOX6", "NDNF"),
      negative = c("ALDH1A1")
    )
  )
  
  # Check each cluster for vulnerability patterns
  vulnerability_assessment <- list()
  
  for (cluster_id in names(pos_neg_markers)) {
    pos_genes <- pos_neg_markers[[cluster_id]]$positive$gene
    neg_genes <- pos_neg_markers[[cluster_id]]$negative$gene
    
    scores <- list()
    
    # Check each vulnerability pattern
    for (pattern_name in names(vulnerability_markers)) {
      pattern <- vulnerability_markers[[pattern_name]]
      
      # Count matches
      pos_matches <- sum(pattern$positive %in% pos_genes)
      neg_matches <- sum(pattern$negative %in% neg_genes)
      
      # Calculate score (weighted by importance)
      score <- (pos_matches * 2) + neg_matches
      
      scores[[pattern_name]] <- list(
        score = score,
        pos_matches = pos_matches,
        neg_matches = neg_matches,
        matched_pos = intersect(pattern$positive, pos_genes),
        matched_neg = intersect(pattern$negative, neg_genes)
      )
    }
    
    # Determine best match
    best_match <- names(which.max(sapply(scores, function(x) x$score)))
    
    vulnerability_assessment[[cluster_id]] <- list(
      scores = scores,
      best_match = best_match,
      vulnerability = ifelse(grepl("vulnerable", best_match), "Vulnerable", "Resilient")
    )
  }
  
  return(vulnerability_assessment)
}

# Function to create decision trees for classification
create_classification_rules <- function() {
  
  rules <- list(
    # Dopaminergic neuron subtypes
    da_a9_snc = list(
      name = "A9 SNc Dopaminergic",
      required_pos = c("TH", "DDC", "KCNJ6", "ALDH1A1", "SOX6"),
      required_neg = c("CALB1", "OTX2"),
      optional_pos = c("SLC6A3", "SLC18A2", "SNCG"),
      confidence_boost = 20
    ),
    
    da_a10_vta = list(
      name = "A10 VTA Dopaminergic",
      required_pos = c("TH", "DDC", "CALB1", "OTX2"),
      required_neg = c("SOX6"),
      optional_pos = c("SLC6A3", "SLC18A2", "OTOP1"),
      confidence_boost = 20
    ),
    
    da_a8_rrf = list(
      name = "A8 RRF Dopaminergic",
      required_pos = c("TH", "DDC", "ALDH1A1"),
      required_neg = c("KCNJ6", "CALB1"),
      optional_pos = c("SLC6A3", "SLC18A2"),
      confidence_boost = 15
    ),
    
    # Floor plate stages
    early_floor_plate = list(
      name = "Early Floor Plate",
      required_pos = c("FOXA2", "SHH", "CORIN"),
      required_neg = c("TH", "DDC", "NR4A2"),
      optional_pos = c("LMX1A", "LMX1B", "MSX1"),
      confidence_boost = 15
    ),
    
    late_floor_plate = list(
      name = "Late Floor Plate/Pre-DA",
      required_pos = c("FOXA2", "LMX1A", "NR4A2"),
      required_neg = c("TH", "SHH"),
      optional_pos = c("EN1", "EN2", "PITX3"),
      confidence_boost = 15
    ),
    
    # Glutamatergic subtypes
    glut_meis1_vulnerable = list(
      name = "MEIS1+ Vulnerable Glutamatergic",
      required_pos = c("MEIS1", "SLC17A6"),
      required_neg = c("TH", "GAD1", "GAD2"),
      optional_pos = c("TBR1", "SATB2"),
      confidence_boost = 10
    ),
    
    glut_meis2_vulnerable = list(
      name = "MEIS2+ Vulnerable Glutamatergic",
      required_pos = c("MEIS2", "SLC17A6"),
      required_neg = c("TH", "GAD1", "GAD2"),
      optional_pos = c("TBR1", "SATB2"),
      confidence_boost = 10
    ),
    
    # Regional identity
    hypothalamic = list(
      name = "Hypothalamic",
      required_pos = c("POU3F2", "NKX2.1"),
      required_neg = c("EN1", "EN2", "PAX6"),
      optional_pos = c("OTP", "SIM1", "SIM2", "RAX"),
      confidence_boost = 15
    ),
    
    caudal_hindbrain = list(
      name = "Caudal/Hindbrain",
      required_pos = c("HOXD10", "HOXD11", "GBX2"),
      required_neg = c("OTX2", "EN1"),
      optional_pos = c("HOXA9", "HOXB9"),
      confidence_boost = 15
    )
  )
  
  return(rules)
}

# Main analysis function
main <- function() {
  cat("Positive/Negative Marker Analysis for Classification Refinement\n")
  cat("=============================================================\n\n")
  
  # Load Seurat object
  seurat_path <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_fixed_annotated.rds"
  if (!file.exists(seurat_path)) {
    stop("Please run the clustering and annotation script first")
  }
  
  cat("Loading Seurat object...\n")
  seurat_obj <- readRDS(seurat_path)
  
  # Extract positive and negative markers
  pos_neg_markers <- extract_pos_neg_markers(seurat_obj)
  
  # Check vulnerability patterns
  cat("\n\nChecking vulnerability marker patterns...\n")
  vulnerability <- check_vulnerability_markers(pos_neg_markers)
  
  # Create classification rules
  rules <- create_classification_rules()
  
  # Save results
  dir.create("results/pos_neg_markers", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(pos_neg_markers, "results/pos_neg_markers/all_pos_neg_markers.rds")
  saveRDS(vulnerability, "results/pos_neg_markers/vulnerability_assessment.rds")
  saveRDS(rules, "results/pos_neg_markers/classification_rules.rds")
  
  # Create summary report
  cat("\n\n=== Vulnerability Assessment Summary ===\n")
  for (cluster_id in names(vulnerability)) {
    assessment <- vulnerability[[cluster_id]]
    cat("\nCluster", cluster_id, ":", assessment$best_match, 
        "(", assessment$vulnerability, ")\n")
  }
  
  # Export key markers for visualization
  key_markers <- c(
    # Vulnerability markers
    "SOX6", "ALDH1A1", "KCNJ6", "CALB1", "MEIS1", "MEIS2",
    # DA markers
    "TH", "DDC", "SLC6A3", "SLC18A2", "NR4A2", "FOXA2",
    # Regional markers
    "OTX2", "GBX2", "EN1", "EN2", "HOXD10", "HOXD11",
    # Resilience markers
    "FOXP2", "SORCS3", "RELN", "TMEFF2"
  )
  
  write.csv(data.frame(marker = key_markers),
            "results/pos_neg_markers/key_vulnerability_markers.csv",
            row.names = FALSE)
  
  cat("\n\nAnalysis complete! Results saved to results/pos_neg_markers/\n")
  
  return(list(
    markers = pos_neg_markers,
    vulnerability = vulnerability,
    rules = rules
  ))
}

# Run if not interactive
if (!interactive()) {
  results <- main()
}