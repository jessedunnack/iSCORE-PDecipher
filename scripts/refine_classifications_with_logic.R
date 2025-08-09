#!/usr/bin/env Rscript

# Refine cell type classifications using positive/negative marker logic
# Incorporates vulnerability markers from MPTP studies

library(dplyr)
library(tidyr)

# Function to score cluster against classification rules
score_cluster_against_rules <- function(pos_genes, neg_genes, rules) {
  
  scores <- list()
  
  for (rule_name in names(rules)) {
    rule <- rules[[rule_name]]
    
    # Check required positive markers
    req_pos_score <- sum(rule$required_pos %in% pos_genes) / length(rule$required_pos)
    
    # Check required negative markers (should NOT be in positive list)
    req_neg_score <- sum(rule$required_neg %in% neg_genes) / length(rule$required_neg)
    
    # Check optional positive markers
    opt_pos_score <- if(length(rule$optional_pos) > 0) {
      sum(rule$optional_pos %in% pos_genes) / length(rule$optional_pos)
    } else { 0 }
    
    # Calculate total score
    total_score <- (req_pos_score * 40) + (req_neg_score * 40) + (opt_pos_score * 20)
    
    # Apply confidence boost if all required criteria are met
    if (req_pos_score == 1 && req_neg_score == 1) {
      total_score <- total_score + rule$confidence_boost
    }
    
    scores[[rule_name]] <- list(
      name = rule$name,
      total_score = total_score,
      req_pos_score = req_pos_score,
      req_neg_score = req_neg_score,
      opt_pos_score = opt_pos_score,
      matched_req_pos = intersect(rule$required_pos, pos_genes),
      matched_req_neg = intersect(rule$required_neg, neg_genes),
      matched_opt_pos = intersect(rule$optional_pos, pos_genes)
    )
  }
  
  return(scores)
}

# Function to refine classifications
refine_all_classifications <- function(pos_neg_markers, rules, current_annotations) {
  
  refined_annotations <- list()
  
  for (cluster_id in names(pos_neg_markers)) {
    cat("\nRefining cluster", cluster_id, "...\n")
    
    # Get positive and negative genes
    pos_genes <- pos_neg_markers[[cluster_id]]$positive$gene
    neg_genes <- pos_neg_markers[[cluster_id]]$negative$gene
    
    # Score against all rules
    rule_scores <- score_cluster_against_rules(pos_genes, neg_genes, rules)
    
    # Find best matching rule
    best_rule <- NULL
    best_score <- 0
    
    for (rule_name in names(rule_scores)) {
      if (rule_scores[[rule_name]]$total_score > best_score) {
        best_score <- rule_scores[[rule_name]]$total_score
        best_rule <- rule_scores[[rule_name]]
      }
    }
    
    # Get current annotation
    cluster_num <- as.numeric(cluster_id)
    current <- current_annotations[current_annotations$fine_cluster == cluster_num, ]
    
    # Determine refined classification
    refined_type <- if (!is.null(best_rule) && best_score >= 60) {
      best_rule$name
    } else {
      current$cell_type
    }
    
    # Calculate confidence based on scores
    refined_confidence <- case_when(
      best_score >= 80 ~ "High",
      best_score >= 50 ~ "Medium",
      TRUE ~ "Low"
    )
    
    # Build evidence string
    evidence <- ""
    if (!is.null(best_rule)) {
      evidence <- paste0(
        "Pos: ", paste(best_rule$matched_req_pos, collapse=","),
        "; Neg: ", paste(best_rule$matched_req_neg, collapse=","),
        "; Score: ", round(best_score, 1)
      )
    }
    
    refined_annotations[[cluster_id]] <- data.frame(
      fine_cluster = cluster_num,
      original_type = current$cell_type,
      original_confidence = current$confidence,
      refined_type = refined_type,
      refined_confidence = refined_confidence,
      best_rule_score = best_score,
      evidence = evidence,
      changed = refined_type != current$cell_type,
      stringsAsFactors = FALSE
    )
  }
  
  # Combine results
  refined_df <- do.call(rbind, refined_annotations)
  
  return(refined_df)
}

# Function to create detailed subtype annotations
create_subtype_annotations <- function(pos_neg_markers, vulnerability_assessment) {
  
  subtype_info <- list()
  
  for (cluster_id in names(pos_neg_markers)) {
    pos_genes <- pos_neg_markers[[cluster_id]]$positive$gene
    neg_genes <- pos_neg_markers[[cluster_id]]$negative$gene
    
    # Check for specific subtype markers
    subtypes <- list()
    
    # Dopaminergic subtypes
    if ("TH" %in% pos_genes && "DDC" %in% pos_genes) {
      if ("SOX6" %in% pos_genes && "ALDH1A1" %in% pos_genes && !"CALB1" %in% pos_genes) {
        subtypes$da_subtype <- "A9-like (SNc, vulnerable)"
      } else if ("CALB1" %in% pos_genes && "OTX2" %in% pos_genes && !"SOX6" %in% pos_genes) {
        subtypes$da_subtype <- "A10-like (VTA, resilient)"
      } else if ("ALDH1A1" %in% pos_genes && !"KCNJ6" %in% pos_genes) {
        subtypes$da_subtype <- "A8-like (RRF)"
      } else {
        subtypes$da_subtype <- "Unspecified DA"
      }
    }
    
    # Glutamatergic subtypes
    if ("SLC17A6" %in% pos_genes || "SLC17A7" %in% pos_genes) {
      if ("MEIS1" %in% pos_genes) {
        subtypes$glut_subtype <- "MEIS1+ (vulnerable)"
      } else if ("MEIS2" %in% pos_genes) {
        subtypes$glut_subtype <- "MEIS2+ (vulnerable)"
      } else if ("LMX1A" %in% pos_genes) {
        subtypes$glut_subtype <- "LMX1A+ (resilient)"
      } else {
        subtypes$glut_subtype <- "Unspecified Glut"
      }
    }
    
    # Regional subtypes
    if ("HOXD10" %in% pos_genes || "HOXD11" %in% pos_genes) {
      subtypes$regional <- "Caudal/Spinal"
    } else if ("OTX2" %in% pos_genes && !"GBX2" %in% pos_genes) {
      subtypes$regional <- "Midbrain"
    } else if ("POU3F2" %in% pos_genes || "NKX2.1" %in% pos_genes) {
      subtypes$regional <- "Hypothalamic"
    }
    
    # Developmental stage
    if ("PCNA" %in% pos_genes || "MKI67" %in% pos_genes) {
      subtypes$stage <- "Proliferating"
    } else if ("DCX" %in% pos_genes && !"MAP2" %in% pos_genes) {
      subtypes$stage <- "Neuroblast"
    } else if ("MAP2" %in% pos_genes && "SYN1" %in% pos_genes) {
      subtypes$stage <- "Mature neuron"
    } else if ("SOX2" %in% pos_genes && "NES" %in% pos_genes) {
      subtypes$stage <- "Progenitor"
    }
    
    # Get vulnerability status
    vuln_status <- vulnerability_assessment[[cluster_id]]$vulnerability
    
    subtype_info[[cluster_id]] <- list(
      cluster = as.numeric(cluster_id),
      subtypes = subtypes,
      vulnerability = vuln_status
    )
  }
  
  return(subtype_info)
}

# Main refinement function
main <- function() {
  cat("Refining Cell Type Classifications with Positive/Negative Logic\n")
  cat("=============================================================\n\n")
  
  # Load required data
  pos_neg_markers <- readRDS("results/pos_neg_markers/all_pos_neg_markers.rds")
  vulnerability <- readRDS("results/pos_neg_markers/vulnerability_assessment.rds")
  rules <- readRDS("results/pos_neg_markers/classification_rules.rds")
  current_annotations <- read.csv("results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv")
  
  # Refine classifications
  cat("Applying classification rules...\n")
  refined <- refine_all_classifications(pos_neg_markers, rules, current_annotations)
  
  # Create subtype annotations
  cat("\nGenerating subtype annotations...\n")
  subtypes <- create_subtype_annotations(pos_neg_markers, vulnerability)
  
  # Combine all information
  final_annotations <- refined
  
  # Add subtype information
  for (i in 1:nrow(final_annotations)) {
    cluster_id <- as.character(final_annotations$fine_cluster[i])
    if (cluster_id %in% names(subtypes)) {
      subtype_list <- subtypes[[cluster_id]]$subtypes
      final_annotations$da_subtype[i] <- ifelse(!is.null(subtype_list$da_subtype), 
                                               subtype_list$da_subtype, "")
      final_annotations$glut_subtype[i] <- ifelse(!is.null(subtype_list$glut_subtype), 
                                                 subtype_list$glut_subtype, "")
      final_annotations$regional_identity[i] <- ifelse(!is.null(subtype_list$regional), 
                                                      subtype_list$regional, "")
      final_annotations$developmental_stage[i] <- ifelse(!is.null(subtype_list$stage), 
                                                        subtype_list$stage, "")
      final_annotations$vulnerability_status[i] <- subtypes[[cluster_id]]$vulnerability
    }
  }
  
  # Save results
  write.csv(final_annotations, 
            "results/cell_type_annotations/refined_classifications.csv",
            row.names = FALSE)
  
  # Create summary report
  cat("\n\n=== Classification Refinement Summary ===\n")
  
  # Changed classifications
  changed <- final_annotations[final_annotations$changed, ]
  if (nrow(changed) > 0) {
    cat("\nClusters with changed classifications:\n")
    for (i in 1:nrow(changed)) {
      cat(sprintf("  Cluster %d: %s -> %s (Score: %.1f)\n",
                  changed$fine_cluster[i],
                  changed$original_type[i],
                  changed$refined_type[i],
                  changed$best_rule_score[i]))
    }
  } else {
    cat("\nNo classifications were changed based on pos/neg logic.\n")
  }
  
  # Confidence changes
  conf_improved <- sum(final_annotations$refined_confidence == "High" & 
                      final_annotations$original_confidence != "High")
  cat("\n", conf_improved, "clusters improved to High confidence\n")
  
  # Vulnerability summary
  cat("\n=== Vulnerability Summary ===\n")
  vuln_summary <- table(final_annotations$vulnerability_status)
  print(vuln_summary)
  
  # DA subtype summary
  da_clusters <- final_annotations[final_annotations$da_subtype != "", ]
  if (nrow(da_clusters) > 0) {
    cat("\n=== Dopaminergic Subtypes ===\n")
    for (i in 1:nrow(da_clusters)) {
      cat(sprintf("  Cluster %d: %s\n", 
                  da_clusters$fine_cluster[i],
                  da_clusters$da_subtype[i]))
    }
  }
  
  cat("\n\nRefinement complete! Results saved to:\n")
  cat("  results/cell_type_annotations/refined_classifications.csv\n")
  
  return(final_annotations)
}

# Run if not interactive
if (!interactive()) {
  refined <- main()
}