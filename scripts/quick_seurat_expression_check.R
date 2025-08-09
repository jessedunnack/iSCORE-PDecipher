#!/usr/bin/env Rscript

# QUICK SEURAT EXPRESSION CHECK
# Direct expression analysis without time-consuming FindMarkers

library(Seurat)
library(dplyr)
library(tidyr)

cat("=================================================================\n")
cat("QUICK SEURAT EXPRESSION CHECK\n")
cat("Direct analysis without FindMarkers\n")
cat("=================================================================\n\n")

# Function to check gene expression across clusters
check_gene_expression <- function(seurat_obj, genes, cluster_col = "seurat_clusters_coarse") {
  cat(sprintf("\nChecking expression of %d genes across %s...\n", length(genes), cluster_col))
  
  # Get expression matrix for genes of interest
  expr_data <- GetAssayData(seurat_obj, slot = "data")
  
  # Filter to genes that exist in the data
  genes_present <- genes[genes %in% rownames(expr_data)]
  genes_missing <- setdiff(genes, genes_present)
  
  if (length(genes_missing) > 0) {
    cat("Warning: Following genes not found in data:", paste(genes_missing, collapse = ", "), "\n")
  }
  
  if (length(genes_present) == 0) {
    cat("No requested genes found in data!\n")
    return(NULL)
  }
  
  # Get cluster assignments
  clusters <- seurat_obj@meta.data[[cluster_col]]
  
  # Calculate stats for each gene
  results <- list()
  
  for (gene in genes_present) {
    gene_expr <- expr_data[gene, ]
    
    # Calculate per cluster
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
  
  # Combine results
  combined <- do.call(rbind, results)
  return(combined)
}

# Load Seurat object
cat("1. Loading Seurat object...\n")
seurat_path <- "E:/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"

# Convert path for Windows
seurat_path <- gsub("/", "\\\\", seurat_path)

if (!file.exists(seurat_path)) {
  cat("ERROR: Seurat object not found at:", seurat_path, "\n")
  cat("Trying alternative path...\n")
  seurat_path <- "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds"
}

seurat_obj <- readRDS(seurat_path)
cat("Seurat object loaded successfully\n")
cat("Number of cells:", ncol(seurat_obj), "\n")

# Set default assay
DefaultAssay(seurat_obj) <- "SCT"
cat("Default assay set to:", DefaultAssay(seurat_obj), "\n")

# Check available metadata columns
cat("\nAvailable cluster columns:\n")
cluster_cols <- grep("cluster", names(seurat_obj@meta.data), value = TRUE, ignore.case = TRUE)
cat(paste(cluster_cols, collapse = ", "), "\n")

# 2. Check TH expression specifically
cat("\n\n2. CHECKING TH EXPRESSION IN ALL CLUSTERS\n")
cat("==========================================\n")

th_expr <- check_gene_expression(seurat_obj, "TH", cluster_col = "seurat_clusters_coarse")
if (!is.null(th_expr)) {
  cat("\nTH expression by coarse cluster:\n")
  th_expr <- th_expr %>% arrange(desc(mean_expr))
  print(th_expr)
  
  # Highlight cluster 0
  cluster0_th <- th_expr[th_expr$cluster == "0", ]
  cat(sprintf("\n*** CLUSTER 0 TH EXPRESSION: mean=%.3f, pct_cells=%.1f%% ***\n",
              cluster0_th$mean_expr, cluster0_th$pct_expressing))
}

# 3. Check all DA markers
cat("\n\n3. CHECKING ALL DOPAMINERGIC MARKERS\n")
cat("=====================================\n")

da_markers <- c("TH", "DDC", "SLC6A3", "SLC18A2", "DRD2", "KCNJ6", 
                "LMX1A", "FOXA2", "NR4A2", "PITX3", "EN1", "EN2",
                "ALDH1A1", "SOX6", "CALB1", "SNCG", "ATP13A2")

da_expr <- check_gene_expression(seurat_obj, da_markers, cluster_col = "seurat_clusters_coarse")

if (!is.null(da_expr)) {
  # Show which clusters have multiple DA markers
  cat("\nClusters with DA marker expression:\n")
  da_summary <- da_expr %>%
    filter(pct_expressing > 5) %>%  # At least 5% of cells
    group_by(cluster) %>%
    summarise(
      n_da_markers = n(),
      markers = paste(gene, collapse = ", "),
      avg_pct = mean(pct_expressing)
    ) %>%
    arrange(desc(n_da_markers))
  
  print(da_summary)
}

# 4. Verify known cluster identities
cat("\n\n4. VERIFYING KNOWN CLUSTER IDENTITIES\n")
cat("======================================\n")

# Check choroid plexus markers in cluster 7
cat("\nCluster 7 - Expected: Choroid Plexus\n")
cp_markers <- c("TTR", "FOLR1", "CLIC6", "HTR2C", "AQP1", "KCNJ13")
cp_expr <- check_gene_expression(seurat_obj, cp_markers, cluster_col = "seurat_clusters_coarse")
if (!is.null(cp_expr)) {
  cluster7_cp <- cp_expr %>% filter(cluster == "7")
  print(cluster7_cp)
}

# Check proliferation markers in cluster 10
cat("\n\nCluster 10 - Expected: Proliferating cells\n")
prolif_markers <- c("MKI67", "TOP2A", "PCNA", "CDC20", "UBE2C")
prolif_expr <- check_gene_expression(seurat_obj, prolif_markers, cluster_col = "seurat_clusters_coarse")
if (!is.null(prolif_expr)) {
  cluster10_prolif <- prolif_expr %>% filter(cluster == "10")
  print(cluster10_prolif)
}

# 5. Check cluster 0 identity more thoroughly
cat("\n\n5. DETAILED ANALYSIS OF CLUSTER 0\n")
cat("==================================\n")

cluster0_genes <- c(
  # Neuronal
  "TUBB3", "MAP2", "STMN2", "SYN1", "SNAP25",
  # DA markers
  "TH", "DDC", "SLC6A3", "SLC18A2",
  # Glutamatergic
  "SLC17A6", "SLC17A7",
  # Other markers from my top 50
  "LHX9", "C1QL4", "MAB21L1", "SNCB", "RIT2"
)

cluster0_expr <- check_gene_expression(seurat_obj, cluster0_genes, cluster_col = "seurat_clusters_coarse")
if (!is.null(cluster0_expr)) {
  cluster0_only <- cluster0_expr %>% 
    filter(cluster == "0") %>%
    arrange(desc(pct_expressing))
  
  cat("\nGenes expressed in >10% of cluster 0 cells:\n")
  print(cluster0_only %>% filter(pct_expressing > 10))
}

# Save results
results <- list(
  th_expression = th_expr,
  da_markers = da_expr,
  da_summary = da_summary,
  cluster0_analysis = cluster0_only
)

saveRDS(results, "results/quick_expression_check_results.rds")
cat("\n\nResults saved to: results/quick_expression_check_results.rds\n")