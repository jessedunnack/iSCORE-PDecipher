#!/usr/bin/env Rscript

# ADD MISSING CELL TYPE COLUMNS TO FINAL ANNOTATED SEURAT OBJECT
# This script adds the missing celltypes_coarse and cell_type_broad columns

library(Seurat)
library(dplyr)

cat("=================================================================\n")
cat("ADDING MISSING CELL TYPE COLUMNS TO FINAL SEURAT OBJECT\n")
cat("=================================================================\n\n")

# 1. Load the current final annotated object
cat("1. Loading final annotated Seurat object...\n")
seurat_obj <- readRDS("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")
cat(sprintf("  - Loaded object with %d cells\n", ncol(seurat_obj)))

# Check current columns
cat("\n2. Checking current metadata columns...\n")
metadata_cols <- colnames(seurat_obj@meta.data)
cat("  Cell type related columns found:\n")
celltype_cols <- grep("cell|type|cluster", metadata_cols, value = TRUE, ignore.case = TRUE)
for(col in celltype_cols) {
  cat(sprintf("    - %s\n", col))
}

# 3. Load the coarse cluster identities
cat("\n3. Loading coarse cluster identities...\n")
coarse_identities <- read.csv("results/reclustered_analysis/coarse_cluster_identities_FINAL.csv")

cat("\nCoarse cluster identities:\n")
print(coarse_identities[, c("cluster", "identity", "cell_type_broad")])

# 4. Add missing columns
cat("\n4. Adding missing cell type columns...\n")

# Check if seurat_clusters_coarse exists
if(!"seurat_clusters_coarse" %in% metadata_cols) {
  cat("  ERROR: seurat_clusters_coarse column not found!\n")
  cat("  This is required for mapping cell types.\n")
  stop("Missing required clustering column")
}

# Add celltypes_coarse
if(!"celltypes_coarse" %in% metadata_cols) {
  cat("  - Adding celltypes_coarse...\n")
  seurat_obj$celltypes_coarse <- plyr::mapvalues(
    seurat_obj$seurat_clusters_coarse,
    from = as.character(coarse_identities$cluster),
    to = coarse_identities$identity
  )
} else {
  cat("  - celltypes_coarse already exists\n")
}

# Add cell_type_broad
if(!"cell_type_broad" %in% metadata_cols) {
  cat("  - Adding cell_type_broad...\n")
  seurat_obj$cell_type_broad <- plyr::mapvalues(
    seurat_obj$seurat_clusters_coarse,
    from = as.character(coarse_identities$cluster),
    to = coarse_identities$cell_type_broad
  )
} else {
  cat("  - cell_type_broad already exists\n")
}

# 5. Verify all columns are present
cat("\n5. Verifying all expected columns...\n")

expected_cols <- c("seurat_clusters_coarse", "seurat_clusters_fine", 
                   "celltypes_coarse", "cell_type_broad", "cell_type_fine")

for(col in expected_cols) {
  if(col %in% colnames(seurat_obj@meta.data)) {
    n_unique <- length(unique(seurat_obj@meta.data[[col]]))
    cat(sprintf("  ✓ %s: %d unique values\n", col, n_unique))
  } else {
    cat(sprintf("  ✗ %s: MISSING\n", col))
  }
}

# 6. Show summary statistics
cat("\n6. Cell type distribution summary:\n")

# Broad categories
cat("\nBroad cell type categories:\n")
broad_table <- table(seurat_obj$cell_type_broad)
broad_df <- data.frame(
  Category = names(broad_table),
  n_cells = as.numeric(broad_table),
  pct = round(100 * as.numeric(broad_table) / ncol(seurat_obj), 1)
)
print(broad_df[order(broad_df$n_cells, decreasing = TRUE), ])

# Coarse cell types (top 10)
cat("\nTop 10 coarse cell types:\n")
coarse_table <- sort(table(seurat_obj$celltypes_coarse), decreasing = TRUE)
coarse_df <- data.frame(
  Cell_Type = names(coarse_table)[1:10],
  n_cells = as.numeric(coarse_table)[1:10],
  pct = round(100 * as.numeric(coarse_table)[1:10] / ncol(seurat_obj), 1)
)
print(coarse_df)

# 7. Save the updated object
cat("\n7. Saving updated Seurat object...\n")

# Backup the original
backup_file <- paste0("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED_backup_", 
                      format(Sys.Date(), "%Y%m%d"), ".rds")
cat(sprintf("  - Creating backup: %s\n", basename(backup_file)))
file.copy("../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds", 
          backup_file, overwrite = FALSE)

# Save updated object
saveRDS(seurat_obj, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_FINAL_ANNOTATED.rds")
cat("  - Saved updated object with all cell type columns\n")

# Also save updated metadata
metadata_df <- seurat_obj@meta.data
saveRDS(metadata_df, "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi_metadata_FINAL.rds")
cat("  - Updated metadata file\n")

cat("\n=== COMPLETE ===\n")
cat("The Seurat object now contains all expected cell type columns:\n")
cat("- seurat_clusters_coarse: 15 coarse clusters\n")
cat("- seurat_clusters_fine: 36 fine clusters\n")
cat("- celltypes_coarse: Detailed coarse cell type names\n")
cat("- cell_type_broad: Broad categories (Neurons/Progenitors/etc)\n")
cat("- cell_type_fine: Detailed fine cluster identities\n")

cat("\nYou can now use these columns for:\n")
cat("- Visualization (color by cell type)\n")
cat("- Stratified analyses\n")
cat("- Cell type-specific DE analyses\n")