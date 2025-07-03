#!/usr/bin/env Rscript

# Deep examination of SingleCellExperiment objects to find cluster information

# Function to examine SCE object without loading SingleCellExperiment package
examine_sce <- function(file_path, dataset_name) {
  cat("\n", strrep("=", 70), "\n")
  cat("Examining SCE file:", dataset_name, "\n")
  cat("File:", basename(file_path), "\n")
  cat(strrep("-", 70), "\n")
  
  if (!file.exists(file_path)) {
    cat("ERROR: File not found!\n")
    return(NULL)
  }
  
  tryCatch({
    # Load the object
    sce <- readRDS(file_path)
    
    # Check S4 slots
    if (isS4(sce)) {
      cat("S4 object with slots:\n")
      slots <- slotNames(sce)
      for (slot_name in slots) {
        cat("  -", slot_name)
        slot_obj <- slot(sce, slot_name)
        if (is.list(slot_obj)) {
          cat(" (list, length =", length(slot_obj), ")")
          if (!is.null(names(slot_obj))) {
            cat(" [", paste(head(names(slot_obj), 3), collapse=", "))
            if (length(names(slot_obj)) > 3) cat(", ...")
            cat("]")
          }
        } else if (is.matrix(slot_obj) || is.data.frame(slot_obj)) {
          cat(" (", class(slot_obj)[1], ", dim =", paste(dim(slot_obj), collapse="x"), ")", sep="")
        } else {
          cat(" (", class(slot_obj)[1], ")", sep="")
        }
        cat("\n")
      }
      
      # Try to access colData
      if ("colData" %in% slots) {
        cat("\nExamining colData slot:\n")
        colData <- slot(sce, "colData")
        
        # colData might be a DataFrame (S4) object
        if (isS4(colData) && "listData" %in% slotNames(colData)) {
          # Extract the actual data
          listData <- slot(colData, "listData")
          if (is.list(listData)) {
            cat("  colData columns:", paste(names(listData), collapse=", "), "\n")
            
            # Look for cluster information
            cluster_cols <- grep("cluster", names(listData), ignore.case = TRUE, value = TRUE)
            if (length(cluster_cols) > 0) {
              cat("  Found cluster columns:", paste(cluster_cols, collapse=", "), "\n")
              
              for (col in cluster_cols) {
                clusters <- listData[[col]]
                if (!is.null(clusters)) {
                  cat("\n  Analyzing column '", col, "':\n", sep="")
                  unique_clusters <- sort(unique(as.character(clusters)))
                  cat("    Unique values:", paste(head(unique_clusters, 20), collapse=", "))
                  if (length(unique_clusters) > 20) cat(", ...")
                  cat("\n")
                  cat("    Number of unique values:", length(unique_clusters), "\n")
                  cat("    Total cells:", length(clusters), "\n")
                  
                  # Check if numeric clusters
                  if (all(grepl("^[0-9]+$", unique_clusters))) {
                    numeric_clusters <- as.numeric(unique_clusters)
                    cat("    Numeric range:", min(numeric_clusters), "-", max(numeric_clusters), "\n")
                  }
                }
              }
            } else {
              cat("  No columns with 'cluster' in name found\n")
              cat("  Available columns:", paste(head(names(listData), 10), collapse=", "))
              if (length(names(listData)) > 10) cat(", ...")
              cat("\n")
            }
          }
        } else if (is.data.frame(colData) || is.matrix(colData)) {
          # Regular data frame or matrix
          cat("  colData class:", class(colData)[1], "\n")
          cat("  colData dimensions:", paste(dim(colData), collapse=" x "), "\n")
          if (!is.null(colnames(colData))) {
            cat("  colData columns:", paste(colnames(colData), collapse=", "), "\n")
          }
        }
      }
      
      # Check the int_colData slot (internal column data)
      if ("int_colData" %in% slots) {
        cat("\nChecking int_colData slot:\n")
        int_colData <- slot(sce, "int_colData")
        if (isS4(int_colData) && "listData" %in% slotNames(int_colData)) {
          listData <- slot(int_colData, "listData")
          if (is.list(listData) && length(listData) > 0) {
            cat("  int_colData elements:", paste(names(listData), collapse=", "), "\n")
          }
        }
      }
      
      # Check metadata
      if ("metadata" %in% slots) {
        cat("\nChecking metadata slot:\n")
        metadata <- slot(sce, "metadata")
        if (is.list(metadata) && !is.null(names(metadata))) {
          cat("  Metadata elements:", paste(names(metadata), collapse=", "), "\n")
          if ("cluster_info" %in% names(metadata)) {
            cat("  Found cluster_info in metadata\n")
          }
        }
      }
    }
    
  }, error = function(e) {
    cat("ERROR:", conditionMessage(e), "\n")
  })
}

# Check each file
cat("Deep examination of SingleCellExperiment files\n")
cat("==============================================\n")

# Path to UMAP data
umap_dir <- "inst/extdata/umap_data"

# Check main UMAP files
examine_sce(file.path(umap_dir, "iSCORE_PD_umap_data.rds"), "iSCORE_PD")
examine_sce(file.path(umap_dir, "iSCORE_PD_CRISPRi_umap_data.rds"), "iSCORE_PD_CRISPRi")
examine_sce(file.path(umap_dir, "Full_Dataset_umap_data.rds"), "Full_Dataset")

# Also check the percentage subsampled versions
cat("\n\nChecking subsampled versions:\n")
cat("==============================\n")
examine_sce(file.path(umap_dir, "Full_Dataset_umap_data_30pc.rds"), "Full_Dataset_30pc")