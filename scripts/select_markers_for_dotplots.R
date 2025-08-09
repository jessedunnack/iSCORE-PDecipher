#!/usr/bin/env Rscript

# SELECT MOST DEFINITIVE MARKERS FOR DOTPLOTS
# Identifies optimal markers for both coarse and fine cluster dotplots

library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("MARKER SELECTION FOR DOTPLOTS\n")
cat("=================================================================\n\n")

# 1. Load marker data
cat("1. Loading marker data...\n")

# Load coarse markers
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_coarse.rds")
cat("  - Loaded coarse cluster markers\n")

# Load fine markers
fine_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers_fine.rds")
cat("  - Loaded fine cluster markers\n")

# Load cluster identities
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")
fine_identities <- read.csv("results/reclustered_analysis/fine_cluster_identities_FINAL.csv")
fine_to_coarse <- read.csv("results/reclustered_analysis/fine_to_coarse_mapping.csv")

# 2. Calculate specificity scores
cat("\n2. Calculating marker specificity scores...\n")

calculate_specificity <- function(markers_df) {
  markers_df %>%
    mutate(
      # Ensure pct values are between 0 and 1
      pct.1 = ifelse(pct.1 > 1, pct.1/100, pct.1),
      pct.2 = ifelse(pct.2 > 1, pct.2/100, pct.2),
      # Calculate specificity score
      pct_diff = pct.1 - pct.2,
      specificity_score = pct_diff * avg_log2FC,
      # Additional metrics
      exclusivity = pct.1 / (pct.2 + 0.01),  # How exclusive to cluster
      expression_strength = pct.1 * avg_log2FC  # Combined expression metric
    ) %>%
    arrange(desc(specificity_score))
}

# 3. Select coarse cluster markers
cat("\n3. Selecting coarse cluster markers...\n")

select_coarse_markers <- function(markers_df, n_markers = 3) {
  selected_markers <- list()
  
  # Get unique clusters
  clusters <- sort(unique(as.numeric(gsub("cluster_", "", markers_df$cluster))))
  
  for (cl in clusters) {
    cat(sprintf("  Cluster %d: ", cl))
    
    # Get cluster identity
    cluster_info <- coarse_identities[coarse_identities$cluster == cl, ]
    identity <- ifelse(nrow(cluster_info) > 0, cluster_info$identity[1], "Unknown")
    
    # Filter and score markers
    cl_markers <- markers_df %>%
      filter(cluster == paste0("cluster_", cl)) %>%
      calculate_specificity() %>%
      filter(
        pct.1 > 0.3,      # Expressed in >30% of cluster cells
        pct.2 < 0.2,      # Expressed in <20% of other cells
        avg_log2FC > 0.5  # Substantial fold change
      )
    
    # Special handling for known identity markers
    identity_markers <- character()
    
    if (grepl("Dopaminergic", identity)) {
      identity_markers <- c("TH", "SLC18A2", "SLC6A3", "KCNJ6", "CALB1")
    } else if (grepl("Hypothalamic_HCRT", identity)) {
      identity_markers <- c("HCRT", "PDYN", "NPTX2")
    } else if (grepl("Proliferating", identity)) {
      identity_markers <- c("MKI67", "TOP2A", "PCNA")
    } else if (grepl("Choroid_Plexus", identity)) {
      identity_markers <- c("TTR", "FOLR1", "KCNJ13")
    } else if (grepl("Fibroblast", identity) && grepl("Mesenchymal", identity)) {
      identity_markers <- c("PRRX1", "PRRX2", "COL1A1")
    } else if (grepl("ECM", identity)) {
      identity_markers <- c("DCN", "LUM", "SPARCL1")
    } else if (grepl("Stressed", identity)) {
      identity_markers <- c("GDF15", "ATF3", "DDIT3")
    } else if (grepl("Neuroendocrine", identity)) {
      identity_markers <- c("CHGA", "CALCA", "TPH1")
    } else if (grepl("PTPRZ1", identity)) {
      identity_markers <- c("PTPRZ1", "FABP7", "SOX2")
    } else if (grepl("CRABP1", identity)) {
      identity_markers <- c("CRABP1", "DCC", "NRN1")
    } else if (grepl("PTGDS", identity)) {
      identity_markers <- c("PTGDS", "C3", "CP")
    } else if (grepl("RBP4", identity)) {
      identity_markers <- c("RBP4", "RLBP1", "CRABP2")
    }
    
    # Check if identity markers are in the data
    identity_in_data <- cl_markers %>%
      filter(gene %in% identity_markers) %>%
      arrange(desc(specificity_score))
    
    # Combine identity markers with top specific markers
    top_specific <- cl_markers %>%
      filter(!gene %in% identity_markers) %>%
      head(n_markers)
    
    # Take best identity markers plus top specific markers
    if (nrow(identity_in_data) > 0) {
      selected <- bind_rows(
        head(identity_in_data, 2),  # Up to 2 identity markers
        head(top_specific, n_markers - min(2, nrow(identity_in_data)))
      )
    } else {
      selected <- head(cl_markers, n_markers)
    }
    
    selected_markers[[as.character(cl)]] <- selected
    cat(paste(selected$gene, collapse = ", "), "\n")
  }
  
  return(selected_markers)
}

coarse_selected <- select_coarse_markers(coarse_markers, n_markers = 3)

# 4. Select fine cluster markers
cat("\n\n4. Selecting fine cluster markers...\n")

select_fine_markers <- function(fine_markers_df, fine_to_coarse_df, n_specific = 2, n_shared = 1) {
  selected_markers <- list()
  
  # Group fine clusters by coarse parent
  coarse_groups <- fine_to_coarse_df %>%
    group_by(coarse_cluster) %>%
    summarise(fine_clusters = list(fine_cluster))
  
  for (i in 1:nrow(coarse_groups)) {
    coarse_cl <- coarse_groups$coarse_cluster[i]
    fine_cls <- coarse_groups$fine_clusters[[i]]
    
    cat(sprintf("\nCoarse cluster %d fine clusters:\n", coarse_cl))
    
    # Find shared markers across fine clusters in this coarse group
    shared_markers <- NULL
    
    for (fc in fine_cls) {
      fc_markers <- fine_markers_df %>%
        filter(cluster == as.character(fc)) %>%
        calculate_specificity() %>%
        filter(pct.1 > 0.25, avg_log2FC > 0.25) %>%
        pull(gene)
      
      if (is.null(shared_markers)) {
        shared_markers <- fc_markers
      } else {
        shared_markers <- intersect(shared_markers, fc_markers)
      }
    }
    
    # Score shared markers by average expression across fine clusters
    if (length(shared_markers) > 0) {
      shared_scores <- fine_markers_df %>%
        filter(
          cluster %in% as.character(fine_cls),
          gene %in% shared_markers
        ) %>%
        calculate_specificity() %>%
        group_by(gene) %>%
        summarise(
          avg_score = mean(specificity_score, na.rm = TRUE),
          avg_pct1 = mean(pct.1, na.rm = TRUE),
          .groups = 'drop'
        ) %>%
        arrange(desc(avg_score))
      
      top_shared <- head(shared_scores$gene, n_shared)
    } else {
      top_shared <- character()
    }
    
    # Now select specific markers for each fine cluster
    for (fc in fine_cls) {
      fc_info <- fine_identities[fine_identities$fine_cluster == fc, ]
      fc_identity <- ifelse(nrow(fc_info) > 0, fc_info$fine_identity[1], "Unknown")
      
      # Get specific markers (not in shared list)
      fc_specific <- fine_markers_df %>%
        filter(cluster == as.character(fc)) %>%
        calculate_specificity() %>%
        filter(
          !gene %in% top_shared,
          pct.1 > 0.25,
          pct.2 < 0.15,
          avg_log2FC > 0.5
        ) %>%
        head(n_specific)
      
      # Combine shared and specific
      fc_selected <- c(top_shared, fc_specific$gene)
      selected_markers[[as.character(fc)]] <- fc_selected
      
      cat(sprintf("  Fine cluster %d (%s): %s\n", 
                  fc, fc_identity, paste(fc_selected, collapse = ", ")))
    }
  }
  
  return(selected_markers)
}

fine_selected <- select_fine_markers(fine_markers, fine_to_coarse, n_specific = 2, n_shared = 1)

# 5. Create marker gene lists
cat("\n\n5. Creating final marker gene lists...\n")

# Coarse markers - ordered by cluster grouping
coarse_cluster_order <- c(
  # Progenitors
  10, 1, 4, 2, 11,
  # Non-neuronal/Unknown  
  3, 8, 6, 7, 12, 5, 9, 13,
  # Neurons
  0, 14
)

# Compile coarse markers
coarse_genes <- character()
for (cl in coarse_cluster_order) {
  if (as.character(cl) %in% names(coarse_selected)) {
    coarse_genes <- c(coarse_genes, coarse_selected[[as.character(cl)]]$gene)
  }
}
coarse_genes <- unique(coarse_genes)  # Remove duplicates

cat(sprintf("Selected %d unique genes for coarse dotplot\n", length(coarse_genes)))

# Fine clusters - order by coarse parent then fine number
fine_cluster_order <- fine_to_coarse %>%
  arrange(match(coarse_cluster, coarse_cluster_order), fine_cluster) %>%
  pull(fine_cluster)

# Compile fine markers
fine_genes <- character()
for (fc in fine_cluster_order) {
  if (as.character(fc) %in% names(fine_selected)) {
    fine_genes <- c(fine_genes, fine_selected[[as.character(fc)]])
  }
}
fine_genes <- unique(fine_genes)

cat(sprintf("Selected %d unique genes for fine dotplot\n", length(fine_genes)))

# 6. Save results
cat("\n6. Saving marker selections...\n")

# Create results directory if needed
dir.create("results/dotplot_markers", showWarnings = FALSE, recursive = TRUE)

# Save coarse markers with metadata
coarse_df <- do.call(rbind, lapply(names(coarse_selected), function(cl) {
  df <- coarse_selected[[cl]]
  df$cluster_num <- as.integer(cl)
  df$cluster_identity <- coarse_identities$identity[coarse_identities$cluster == as.integer(cl)]
  return(df)
}))

write.csv(coarse_df, "results/dotplot_markers/selected_markers_coarse.csv", row.names = FALSE)

# Save fine markers list
fine_df <- data.frame(
  fine_cluster = rep(as.integer(names(fine_selected)), sapply(fine_selected, length)),
  gene = unlist(fine_selected),
  order = sequence(sapply(fine_selected, length))
)

# Add identities
fine_df <- fine_df %>%
  left_join(fine_identities[, c("fine_cluster", "fine_identity", "coarse_cluster")], 
            by = "fine_cluster")

write.csv(fine_df, "results/dotplot_markers/selected_markers_fine.csv", row.names = FALSE)

# Save gene lists
writeLines(coarse_genes, "results/dotplot_markers/coarse_genes.txt")
writeLines(fine_genes, "results/dotplot_markers/fine_genes.txt")

# Save cluster orders
write.csv(data.frame(cluster = coarse_cluster_order), 
          "results/dotplot_markers/coarse_cluster_order.csv", row.names = FALSE)
write.csv(data.frame(cluster = fine_cluster_order), 
          "results/dotplot_markers/fine_cluster_order.csv", row.names = FALSE)

cat("\n=== MARKER SELECTION COMPLETE ===\n")
cat("Outputs saved to results/dotplot_markers/\n")
cat("\nNext steps:\n")
cat("1. Run generate_coarse_dotplot.R\n")
cat("2. Run generate_fine_dotplot.R\n")