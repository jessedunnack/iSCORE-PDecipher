#!/usr/bin/env Rscript

# ADD CELL TYPE IDENTITIES TO SEURAT OBJECT

library(Seurat)
library(dplyr)

cat("=================================================================\n")
cat("ADDING CELL TYPE IDENTITIES TO SEURAT OBJECT\n")
cat("=================================================================\n\n")

# Load the reclustered Seurat object
cat("1. Loading Seurat object...\n")
seurat_obj <- readRDS("results/seurat_obj_reclustered.rds")

# Load the final identities
cat("2. Loading final cluster identities...\n")
final_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

# Show what we're adding
cat("\nIdentities to be added:\n")
print(final_identities[, c("cluster", "identity", "cell_type_broad")])

# Add celltypes_coarse
cat("\n3. Adding celltypes_coarse...\n")
seurat_obj$celltypes_coarse <- plyr::mapvalues(
  seurat_obj$seurat_clusters_coarse,
  from = as.character(final_identities$cluster),
  to = final_identities$identity
)

# Add cell_type_broad for convenience
seurat_obj$cell_type_broad <- plyr::mapvalues(
  seurat_obj$seurat_clusters_coarse,
  from = as.character(final_identities$cluster),
  to = final_identities$cell_type_broad
)

# Verify the mapping
cat("\n4. Verifying cell type assignments...\n")
ct_table <- table(seurat_obj$celltypes_coarse, seurat_obj$seurat_clusters_coarse)
cat("Coarse cluster to cell type mapping verified.\n")

# Show summary
cat("\nCell type summary:\n")
print(table(seurat_obj$cell_type_broad))

cat("\nDetailed cell type counts:\n")
ct_summary <- as.data.frame(table(seurat_obj$celltypes_coarse)) %>%
  arrange(desc(Freq))
colnames(ct_summary) <- c("Cell_Type", "n_cells")
ct_summary$pct <- round(100 * ct_summary$n_cells / ncol(seurat_obj), 1)
print(ct_summary)

# Save the annotated object
cat("\n5. Saving annotated Seurat object...\n")
saveRDS(seurat_obj, "results/seurat_obj_annotated_coarse.rds")

# Also save to the original location for downstream use
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_annotated.rds")

cat("\nSuccess! Cell types added:\n")
cat("- celltypes_coarse: Detailed cell type identities\n")
cat("- cell_type_broad: Broad categories (Neurons/Progenitors/Non-neuronal/Unknown)\n")

cat("\nSaved to:\n")
cat("- results/seurat_obj_annotated_coarse.rds\n")
cat("- ../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_annotated.rds\n")

cat("\n\nNEXT STEPS:\n")
cat("1. Visualize cell types on UMAP\n")
cat("2. Analyze fine clusters within each coarse cell type\n")
cat("3. Add celltypes_fine annotations\n")
cat("4. Re-run DE analyses stratified by cell type\n")