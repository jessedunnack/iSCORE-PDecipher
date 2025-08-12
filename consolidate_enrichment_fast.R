#!/usr/bin/env Rscript

# Fast consolidation of enrichment results for iSCORE-PD dataset
# Optimized version with better progress tracking

cat("========================================\n")
cat("Fast Consolidation for iSCORE-PD\n")
cat("========================================\n\n")

suppressPackageStartupMessages({
  library(dplyr)
  library(purrr)
})

# Configuration
base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD"
enrichment_dir <- file.path(base_dir, "enrichment_results")
output_file <- file.path(base_dir, "all_enrichment_padj005_complete_with_direction.rds")

# Get all RDS files
cat("Scanning for enrichment files...\n")
enrichment_files <- list.files(
  enrichment_dir, 
  pattern = "\\.rds$", 
  recursive = TRUE, 
  full.names = TRUE
)

# Filter out metadata and summary files
enrichment_files <- enrichment_files[!grepl("analysis_summary|analysis_metadata", enrichment_files)]

cat(sprintf("Found %d enrichment result files\n\n", length(enrichment_files)))

# Process files in batches
process_batch <- function(files) {
  batch_data <- list()
  
  for (file_path in files) {
    tryCatch({
      # Parse path
      rel_path <- gsub(paste0(enrichment_dir, "/"), "", file_path)
      parts <- strsplit(rel_path, "/")[[1]]
      
      if (length(parts) < 5) next
      
      # Extract metadata
      method <- parts[1]
      mutation <- parts[2]
      cluster <- parts[3]
      experiment <- parts[4]
      
      # Handle GSEA vs regular enrichment
      if (parts[5] == "GSEA" && length(parts) >= 7) {
        enrichment_type <- paste("GSEA", parts[6], sep = "_")
        file_name <- parts[7]
      } else {
        enrichment_type <- parts[5]
        file_name <- if(length(parts) > 5) parts[6] else basename(file_path)
      }
      
      # Parse direction from filename
      direction <- if (grepl("_UP\\.rds$", file_name)) {
        "UP"
      } else if (grepl("_DOWN\\.rds$", file_name)) {
        "DOWN"
      } else {
        "ALL"
      }
      
      # For GSEA, direction is always ALL
      if (grepl("^GSEA", enrichment_type)) {
        direction <- "ALL"
      }
      
      # Read data
      data <- readRDS(file_path)
      
      # Convert to data frame if needed
      if (!is.data.frame(data)) {
        if (is.list(data) && length(data) > 0 && is.data.frame(data[[1]])) {
          data <- data[[1]]
        } else {
          next
        }
      }
      
      # Skip if empty
      if (nrow(data) == 0) next
      
      # Add metadata columns
      data$method <- method
      data$mutation_perturbation <- mutation
      data$cluster <- cluster
      data$experiment <- experiment
      data$enrichment_type <- enrichment_type
      data$direction <- direction
      
      batch_data[[length(batch_data) + 1]] <- data
      
    }, error = function(e) {
      # Silently skip failed files
    })
  }
  
  return(batch_data)
}

# Process in chunks for better memory management
batch_size <- 100
n_batches <- ceiling(length(enrichment_files) / batch_size)

cat(sprintf("Processing %d batches of ~%d files each\n", n_batches, batch_size))
cat("Progress: ")

all_data <- list()

for (i in 1:n_batches) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, length(enrichment_files))
  
  batch_files <- enrichment_files[start_idx:end_idx]
  batch_results <- process_batch(batch_files)
  
  if (length(batch_results) > 0) {
    all_data <- c(all_data, batch_results)
  }
  
  cat(".")
  if (i %% 10 == 0) cat(sprintf(" %d%%", round(100 * i / n_batches)))
}
cat(" Done!\n\n")

# Combine all data
cat("Combining all enrichment data...\n")
consolidated <- bind_rows(all_data)

# Standardize column names
col_mapping <- c(
  "ID" = "term_id",
  "Description" = "term_description",
  "pvalue" = "p.value",
  "qvalue" = "p.adjust",
  "GeneRatio" = "gene_ratio",
  "Count" = "count",
  "geneID" = "geneID"
)

for (old_name in names(col_mapping)) {
  if (old_name %in% names(consolidated) && !col_mapping[old_name] %in% names(consolidated)) {
    names(consolidated)[names(consolidated) == old_name] <- col_mapping[old_name]
  }
}

# Add gene column for compatibility
if (!"gene" %in% names(consolidated)) {
  consolidated$gene <- consolidated$mutation_perturbation
}

# Filter for significant results
if ("p.adjust" %in% names(consolidated)) {
  consolidated <- consolidated %>%
    filter(p.adjust <= 0.05)
}

# Summary
cat("\n========================================\n")
cat("Consolidation Summary\n")
cat("========================================\n\n")

cat(sprintf("Total enrichment terms: %s\n", format(nrow(consolidated), big.mark = ",")))
cat(sprintf("Unique mutations: %d\n", length(unique(consolidated$mutation_perturbation))))
cat(sprintf("Unique clusters: %d\n", length(unique(consolidated$cluster))))

# Breakdown by enrichment type
enrich_summary <- consolidated %>%
  count(enrichment_type) %>%
  arrange(desc(n))

cat("\nEnrichment types:\n")
for (i in 1:min(10, nrow(enrich_summary))) {
  cat(sprintf("  %s: %s terms\n", 
              enrich_summary$enrichment_type[i], 
              format(enrich_summary$n[i], big.mark = ",")))
}

# Save
cat(sprintf("\nSaving to: %s\n", output_file))
saveRDS(consolidated, output_file)

file_size <- file.info(output_file)$size / 1024^2
cat(sprintf("Output file size: %.1f MB\n", file_size))

cat("\n✅ Consolidation complete!\n")