#!/usr/bin/env Rscript

# RECLUSTER AND INVESTIGATE - FIXED VERSION
# Overwrites existing cluster columns

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("RECLUSTER AND INVESTIGATE - OVERWRITING EXISTING CLUSTERS\n")
cat("=================================================================\n\n")

# Load Seurat object
cat("1. Loading Seurat object...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds")
cat("Loaded successfully. Number of cells:", ncol(seurat_obj), "\n")

# Set default assay
DefaultAssay(seurat_obj) <- "SCT"
cat("Default assay:", DefaultAssay(seurat_obj), "\n\n")

# Re-run clustering
cat("2. Re-running clustering...\n")
cat("This may take a few minutes...\n")

# Coarse clusters (resolution 0.2)
cat("\nFinding coarse clusters (resolution = 0.2)...\n")
seurat_obj <- FindClusters(seurat_obj, resolution = 0.2, verbose = FALSE)
seurat_obj$seurat_clusters_coarse <- Idents(seurat_obj)
n_coarse <- length(unique(seurat_obj$seurat_clusters_coarse))
cat("Found", n_coarse, "coarse clusters\n")

# Fine clusters (default resolution)
cat("\nFinding fine clusters (default resolution)...\n")
seurat_obj <- FindClusters(seurat_obj, verbose = FALSE)  # Uses default resolution
seurat_obj$seurat_clusters_fine <- Idents(seurat_obj)
n_fine <- length(unique(seurat_obj$seurat_clusters_fine))
cat("Found", n_fine, "fine clusters\n")

# Save the reclustered object
cat("\n3. Saving reclustered object...\n")
saveRDS(seurat_obj, "results/seurat_obj_reclustered.rds")
cat("Saved to: results/seurat_obj_reclustered.rds\n")

# Now investigate coarse clusters
cat("\n4. Investigating coarse clusters...\n")
cat("=========================================\n")

# Function to check gene expression
check_gene_expression <- function(seurat_obj, genes, cluster_col = "seurat_clusters_coarse") {
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  genes_present <- genes[genes %in% rownames(expr_data)]
  
  if (length(genes_present) == 0) return(NULL)
  
  clusters <- seurat_obj@meta.data[[cluster_col]]
  results <- list()
  
  for (gene in genes_present) {
    gene_expr <- expr_data[gene, ]
    cluster_stats <- data.frame(
      cluster = levels(factor(clusters)),
      stringsAsFactors = FALSE
    )
    
    for (cl in cluster_stats$cluster) {
      cells_in_cluster <- which(clusters == cl)
      expr_in_cluster <- gene_expr[cells_in_cluster]
      
      cluster_stats[cluster_stats$cluster == cl, "mean_expr"] <- mean(expr_in_cluster)
      cluster_stats[cluster_stats$cluster == cl, "pct_expressing"] <- 
        100 * sum(expr_in_cluster > 0) / length(expr_in_cluster)
      cluster_stats[cluster_stats$cluster == cl, "n_cells"] <- length(cells_in_cluster)
    }
    
    cluster_stats$gene <- gene
    results[[gene]] <- cluster_stats
  }
  
  combined <- do.call(rbind, results)
  return(combined)
}

# Quick check of key markers
cat("\n5. Quick marker check for coarse clusters...\n")

# Check TH expression
th_expr <- check_gene_expression(seurat_obj, "TH", "seurat_clusters_coarse")
if (!is.null(th_expr)) {
  cat("\nTH expression in coarse clusters:\n")
  th_expr <- th_expr %>% 
    filter(pct_expressing > 0) %>%
    arrange(desc(pct_expressing))
  print(th_expr)
}

# Check key identity markers
key_markers <- list(
  DA = c("TH", "DDC", "SLC18A2"),
  Choroid = c("TTR", "FOLR1", "CLIC6"),
  Proliferating = c("MKI67", "TOP2A"),
  Neuronal = c("TUBB3", "MAP2", "STMN2")
)

cat("\n\n6. Summary of key markers in each coarse cluster:\n")
cat("==================================================\n")

coarse_clusters <- sort(unique(as.numeric(as.character(seurat_obj$seurat_clusters_coarse))))

for (cl in coarse_clusters) {
  cat(sprintf("\nCluster %d (n=%d cells):\n", 
              cl, 
              sum(seurat_obj$seurat_clusters_coarse == as.character(cl))))
  
  for (set_name in names(key_markers)) {
    expr_data <- check_gene_expression(seurat_obj, key_markers[[set_name]], "seurat_clusters_coarse")
    if (!is.null(expr_data)) {
      cl_data <- expr_data %>%
        filter(cluster == as.character(cl)) %>%
        filter(pct_expressing > 5)
      
      if (nrow(cl_data) > 0) {
        cat(sprintf("  %s: %s (%.1f%%)\n", 
                    set_name,
                    paste(cl_data$gene, collapse=", "),
                    mean(cl_data$pct_expressing)))
      }
    }
  }
}

cat("\n\n=== RECLUSTERING COMPLETE ===\n")
cat("Coarse clusters:", n_coarse, "(resolution 0.2)\n")
cat("Fine clusters:", n_fine, "(default resolution)\n")
cat("\nColumns overwritten:\n")
cat("- seurat_clusters_coarse\n")
cat("- seurat_clusters_fine\n")
cat("\nNext steps:\n")
cat("1. Run comprehensive marker analysis\n")
cat("2. Establish hierarchy between coarse and fine clusters\n")
cat("3. Assign cell type identities\n")