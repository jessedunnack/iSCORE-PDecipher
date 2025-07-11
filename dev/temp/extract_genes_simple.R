# Simple Gene Association Extraction
# More efficient approach with checkpointing

library(dplyr)
library(stringr)

# Simple function to extract genes from one file
extract_genes_from_file <- function(file_path) {
  tryCatch({
    # Parse file path
    parts <- str_split(file_path, "/")[[1]]
    er_idx <- which(parts == "enrichment_results")
    
    analysis_type <- parts[er_idx + 1]
    gene <- parts[er_idx + 2] 
    cluster <- parts[er_idx + 3]
    enrichment_type <- parts[er_idx + 5]
    
    filename <- parts[length(parts)]
    filename_base <- str_remove(filename, "\\.rds$")
    
    if (str_detect(filename_base, "GSEA_")) {
      direction <- "RANKED"
    } else {
      direction_parts <- str_split(filename_base, "_")[[1]]
      direction <- direction_parts[length(direction_parts)]
    }
    
    # Read file
    result <- readRDS(file_path)
    
    # Extract data
    if (inherits(result, "enrichResult") || inherits(result, "gseaResult")) {
      df <- result@result
    } else if (is.data.frame(result)) {
      df <- result
    } else {
      return(NULL)
    }
    
    # Find gene column
    gene_col <- NULL
    if ("geneID" %in% colnames(df)) {
      gene_col <- "geneID"
    } else if ("core_enrichment" %in% colnames(df)) {
      gene_col <- "core_enrichment"
    } else {
      return(NULL)
    }
    
    # Create simple data frame
    simple_df <- df %>%
      filter(!is.na(.data[[gene_col]]) & .data[[gene_col]] != "") %>%
      select(
        term_id = ID,
        description = Description,
        associated_genes = !!sym(gene_col)
      ) %>%
      mutate(
        analysis_type = analysis_type,
        gene = gene,
        cluster = cluster,
        enrichment_type = enrichment_type,
        direction = direction,
        gene_count = str_count(associated_genes, "/") + 1
      )
    
    return(simple_df)
    
  }, error = function(e) {
    return(NULL)
  })
}

# Process one dataset
process_dataset <- function(dataset_dir, output_file) {
  cat("Processing", dataset_dir, "...\n")
  
  # Find files
  cmd <- paste0("find ", file.path(dataset_dir, "enrichment_results"), " -name '*.rds' | grep -E '(GO_|KEGG|Reactome|STRING|GSEA)'")
  files <- system(cmd, intern = TRUE)
  
  cat("Found", length(files), "files\n")
  
  # Process in small batches
  all_results <- list()
  batch_size <- 50
  n_batches <- ceiling(length(files) / batch_size)
  
  for (i in 1:n_batches) {
    cat("Batch", i, "of", n_batches, "...")
    
    start_idx <- (i - 1) * batch_size + 1
    end_idx <- min(i * batch_size, length(files))
    batch_files <- files[start_idx:end_idx]
    
    batch_results <- lapply(batch_files, extract_genes_from_file)
    batch_results <- batch_results[!sapply(batch_results, is.null)]
    
    if (length(batch_results) > 0) {
      batch_df <- bind_rows(batch_results)
      all_results[[i]] <- batch_df
      cat(" extracted", nrow(batch_df), "associations\n")
    } else {
      cat(" no results\n")
    }
    
    # Save progress every 10 batches
    if (i %% 10 == 0) {
      temp_result <- bind_rows(all_results)
      temp_file <- paste0(output_file, "_temp_", i, ".rds")
      saveRDS(temp_result, temp_file, compress = TRUE)
      cat("Saved progress to", temp_file, "\n")
    }
  }
  
  # Combine and save final result
  final_result <- bind_rows(all_results)
  saveRDS(final_result, output_file, compress = "xz")
  
  cat("Saved", nrow(final_result), "associations to", output_file, "\n")
  return(final_result)
}

# Main execution
cat("=== Simple Gene Association Extraction ===\n")

# Process each dataset separately
dataset1 <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD"
dataset2 <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi"

# Create output directory
dir.create("inst/extdata", recursive = TRUE, showWarnings = FALSE)

# Process datasets
result1 <- process_dataset(dataset1, "inst/extdata/gene_associations_iscore.rds")
result2 <- process_dataset(dataset2, "inst/extdata/gene_associations_iscore_crispi.rds")

# Combine results
cat("Combining results...\n")
combined <- bind_rows(result1, result2)

# Save combined result
saveRDS(combined, "inst/extdata/gene_term_associations.rds", compress = "xz")

# Summary
cat("\n=== Final Summary ===\n")
cat("Total associations:", nrow(combined), "\n")
cat("Unique terms:", n_distinct(combined$term_id), "\n")
cat("Unique genes:", n_distinct(combined$gene), "\n")
cat("Analysis types:", paste(unique(combined$analysis_type), collapse = ", "), "\n")

# Check file size
file_size <- file.size("inst/extdata/gene_term_associations.rds") / (1024^2)
cat("Final file size:", round(file_size, 1), "MB\n")

if (file_size > 100) {
  cat("WARNING: File may be too large for GitHub!\n")
} else {
  cat("File size acceptable for GitHub deployment.\n")
}

cat("=== Extraction Complete ===\n")