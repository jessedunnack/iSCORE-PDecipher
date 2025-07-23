#!/usr/bin/env Rscript

# Process all UMAP files to add gene expression data

# Configuration
UMAP_DIR <- "inst/extdata/umap_data"
SEURAT_FILES <- list(
  "iSCORE_PD" = "../../iSCORE-PD/iSCORE-PD.rds",
  "iSCORE_PD_CRISPRi" = "../../iSCORE-PD_plus_CRISPRi/iSCORE-PD_plus_CRISPRi.rds",
  "Full_Dataset" = "../../iSCORE-PD_plus_CRISPRi_and_CRISPRa/iSCORE-PD_plus_CRISPRi_and_CRISPRa.rds"
)

# Process each dataset
for (dataset_name in names(SEURAT_FILES)) {
  cat("\n=== Processing", dataset_name, "===\n")
  
  seurat_file <- SEURAT_FILES[[dataset_name]]
  
  # Check if Seurat file exists
  if (!file.exists(seurat_file)) {
    cat("Seurat file not found:", seurat_file, "\n")
    cat("Skipping", dataset_name, "\n")
    next
  }
  
  # Find all UMAP files for this dataset
  umap_pattern <- sprintf("%s_umap_data.*\\.rds", dataset_name)
  umap_files <- list.files(UMAP_DIR, pattern = umap_pattern, full.names = TRUE)
  
  # Also check for markers file
  markers_file <- file.path(UMAP_DIR, sprintf("%s_cluster_markers.rds", dataset_name))
  markers_arg <- if (file.exists(markers_file)) {
    sprintf("--markers %s", markers_file)
  } else {
    ""
  }
  
  for (umap_file in umap_files) {
    cat("\nProcessing:", basename(umap_file), "\n")
    
    # Create output filename
    output_file <- sub("\\.rds$", "_with_expr.rds", umap_file)
    
    # Build command
    cmd <- sprintf(
      "Rscript add_expression_to_umap.R --seurat %s --umap %s --output %s %s --top-markers 100",
      seurat_file, umap_file, output_file, markers_arg
    )
    
    # Run command
    cat("Running:", cmd, "\n")
    result <- system(cmd)
    
    if (result == 0) {
      cat("Success! Created:", basename(output_file), "\n")
      
      # Optionally replace original file
      # file.rename(output_file, umap_file)
    } else {
      cat("Error processing", basename(umap_file), "\n")
    }
  }
}

cat("\n=== All processing complete ===\n")
cat("New files created with '_with_expr' suffix\n")
cat("To use these files, either:\n")
cat("1. Rename them to remove '_with_expr' suffix\n")
cat("2. Update the app to look for '_with_expr' files\n")