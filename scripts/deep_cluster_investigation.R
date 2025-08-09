#!/usr/bin/env Rscript

# DEEP CLUSTER INVESTIGATION
# Search for specific markers throughout the full DE results

library(dplyr)

cat("=================================================================\n")
cat("DEEP CLUSTER INVESTIGATION - ALL MARKERS\n")
cat("=================================================================\n\n")

# Load marker data
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# Function to search for specific genes in a cluster
search_genes_in_cluster <- function(cluster_id, gene_list, markers_df) {
  cluster_data <- markers_df %>%
    filter(cluster == as.character(cluster_id))
  
  found_genes <- cluster_data %>%
    filter(gene %in% gene_list) %>%
    arrange(desc(avg_log2FC))
  
  return(found_genes)
}

# Function to get all significant markers for a cluster
get_all_significant_markers <- function(cluster_id, markers_df, padj_cutoff = 0.05) {
  sig_markers <- markers_df %>%
    filter(cluster == as.character(cluster_id)) %>%
    filter(p_val_adj < padj_cutoff) %>%
    arrange(desc(avg_log2FC))
  
  return(sig_markers)
}

# Key marker sets to search for
SEARCH_SETS <- list(
  DA_ALL = c("TH", "DDC", "AADC", "SLC6A3", "DAT", "SLC18A2", "VMAT2", 
             "DRD2", "DRD1", "KCNJ6", "GIRK2", "LMX1A", "LMX1B", "FOXA2", 
             "NR4A2", "NURR1", "PITX3", "EN1", "EN2", "ALDH1A1", "SOX6",
             "CALB1", "CALB2", "SNCG", "SNCA", "ATP13A2", "GBA", "LRRK2"),
  
  FLOOR_PLATE_EXTENDED = c("CORIN", "LMX1A", "LMX1B", "FOXA2", "FOXA1", "ARX", 
                          "SHH", "WNT1", "MSX1", "MSX2", "LMO3", "LMO4", "SLIT2", 
                          "NTN1", "BMP7", "SDF1", "CXCL12"),
  
  CHOROID_PLEXUS_ALL = c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "PRLR", 
                        "KCNJ13", "ENPP2", "PLPP3", "SLC13A4", "KL", "TRPM3",
                        "CFTR", "SLC4A10", "SLC12A2", "KCNE2", "FOLR2"),
  
  NEURONAL_SUBTYPES = c("TBR1", "SATB2", "BCL11B", "FOXP2", "CUX1", "CUX2",
                       "RELN", "SST", "PVALB", "VIP", "LAMP5", "SNCG", "CCK",
                       "TAC1", "PENK", "PDYN", "OPRM1", "OPRD1", "OPRK1")
)

# Analyze all clusters (0-14)
cat("Searching for key markers across all coarse clusters...\n\n")

for (cluster_id in 0:14) {
  cat(sprintf("\n========== COARSE CLUSTER %d ==========\n", cluster_id))
  
  # Get all significant markers
  sig_markers <- get_all_significant_markers(cluster_id, coarse_markers)
  total_sig <- nrow(sig_markers)
  
  cat(sprintf("Total significant markers (p_adj < 0.05): %d\n", total_sig))
  
  # Search for specific marker sets
  cat("\n--- SEARCHING KEY MARKER SETS ---\n")
  
  for (set_name in names(SEARCH_SETS)) {
    genes_to_search <- SEARCH_SETS[[set_name]]
    found <- search_genes_in_cluster(cluster_id, genes_to_search, coarse_markers)
    
    if (nrow(found) > 0) {
      cat(sprintf("\n%s: Found %d genes\n", set_name, nrow(found)))
      for (i in 1:nrow(found)) {
        cat(sprintf("  - %s: FC=%.2f, pct.1=%.3f, p_adj=%.2e\n",
                    found$gene[i], found$avg_log2FC[i], 
                    found$pct.1[i], found$p_val_adj[i]))
      }
    }
  }
  
  # Show markers ranked 31-100 to catch deeper signals
  cat("\n--- MARKERS RANKED 31-60 ---\n")
  if (nrow(sig_markers) >= 60) {
    deep_markers <- sig_markers[31:60, ]
    cat(paste(deep_markers$gene, collapse=", "), "\n")
  }
  
  # Check for any DA-related genes in top 200
  cat("\n--- DA-RELATED IN TOP 200 ---\n")
  top200 <- head(sig_markers, 200)
  da_related <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "LMX1A", "FOXA2", 
                  "NR4A2", "PITX3", "EN1", "EN2", "ALDH1A1", "SOX6", "CALB1")
  found_da <- top200 %>% filter(gene %in% da_related)
  if (nrow(found_da) > 0) {
    print(found_da[, c("gene", "avg_log2FC", "pct.1", "p_val_adj")])
  } else {
    cat("None found\n")
  }
}

# Special investigation: Where is TH?
cat("\n\n=== SPECIAL INVESTIGATION: WHERE IS TH? ===\n")
th_clusters <- coarse_markers %>%
  filter(gene == "TH") %>%
  arrange(desc(avg_log2FC))

cat("\nTH expression across all clusters:\n")
print(th_clusters)

# Find clusters with multiple DA markers
cat("\n\n=== CLUSTERS WITH MULTIPLE DA MARKERS ===\n")
da_core <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "LMX1A", "FOXA2", "NR4A2", "PITX3")

for (cluster_id in 0:14) {
  found_da <- search_genes_in_cluster(cluster_id, da_core, coarse_markers)
  if (nrow(found_da) >= 2) {
    cat(sprintf("\nCluster %d has %d DA markers:\n", cluster_id, nrow(found_da)))
    print(found_da[, c("gene", "avg_log2FC", "pct.1")])
  }
}