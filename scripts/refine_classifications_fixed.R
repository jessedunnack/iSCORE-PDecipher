#!/usr/bin/env Rscript

# Fixed version of refine classifications

# Load required functions from original script
source("scripts/refine_classifications_with_logic.R")

# Override the problematic function
refine_all_classifications <- function(pos_neg_markers, rules, current_annotations) {
  
  refined_annotations <- list()
  
  for (cluster_id in names(pos_neg_markers)) {
    cat("\nRefining cluster", cluster_id, "...\n")
    
    # Get current annotation
    cluster_num <- as.numeric(cluster_id)
    current <- current_annotations[current_annotations$fine_cluster == cluster_num, ]
    
    # Skip if no annotation exists (e.g., cluster 36)
    if (nrow(current) == 0) {
      cat("  No annotation found for cluster", cluster_num, "- skipping\n")
      next
    }
    
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

# Run the fixed version
results2 <- main()