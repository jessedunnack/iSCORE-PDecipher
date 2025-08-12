#!/usr/bin/env Rscript

# Consolidate enrichment results for iSCORE-PD dataset
# Creates all_enrichment_padj005_complete_with_direction.rds from individual enrichment files

cat("========================================\n")
cat("Consolidating iSCORE-PD Enrichment Data\n")
cat("========================================\n\n")

library(dplyr)
library(purrr)

# Base directory
base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD"
enrichment_dir <- file.path(base_dir, "enrichment_results")
output_file <- file.path(base_dir, "all_enrichment_padj005_complete_with_direction.rds")

# Check if enrichment directory exists
if (!dir.exists(enrichment_dir)) {
  stop("Enrichment directory not found: ", enrichment_dir)
}

cat("Scanning enrichment directory...\n")

# Get all RDS files recursively
enrichment_files <- list.files(
  enrichment_dir, 
  pattern = "\\.rds$", 
  recursive = TRUE, 
  full.names = TRUE
)

cat(sprintf("Found %d enrichment files\n", length(enrichment_files)))

# Function to extract metadata from file path
extract_metadata <- function(file_path) {
  # Remove base directory and .rds extension
  rel_path <- gsub(paste0(enrichment_dir, "/"), "", file_path)
  rel_path <- gsub("\\.rds$", "", rel_path)
  
  # Split path components
  # Expected formats:
  # MAST/mutation/cluster_n/default/enrichment_type/file.rds
  # MAST/mutation/cluster_n/default/GSEA/database/file.rds
  parts <- strsplit(rel_path, "/")[[1]]
  
  # Skip summary and metadata files
  if (grepl("analysis_summary|analysis_metadata", basename(file_path))) {
    return(NULL)
  }
  
  if (length(parts) >= 5) {
    # Check if this is a GSEA file (has extra directory level)
    if (parts[5] == "GSEA" && length(parts) >= 7) {
      list(
        method = parts[1],
        mutation_perturbation = parts[2],
        cluster = parts[3],
        experiment = parts[4],
        enrichment_type = paste("GSEA", parts[6], sep = "_"),
        file_name = parts[7]
      )
    } else {
      list(
        method = parts[1],
        mutation_perturbation = parts[2],
        cluster = parts[3],
        experiment = parts[4],
        enrichment_type = parts[5],
        file_name = if(length(parts) > 5) parts[6] else parts[5]
      )
    }
  } else {
    NULL
  }
}

# Function to parse enrichment type and direction from filename
parse_enrichment_info <- function(file_name, enrichment_type) {
  # For GSEA files, direction is always ALL
  if (grepl("^GSEA", enrichment_type)) {
    return(list(
      enrichment_type = enrichment_type,
      direction = "ALL"
    ))
  }
  
  # Handle different naming patterns for regular enrichment
  if (grepl("_UP\\.rds$", file_name)) {
    list(
      enrichment_type = enrichment_type,
      direction = "UP"
    )
  } else if (grepl("_DOWN\\.rds$", file_name)) {
    list(
      enrichment_type = enrichment_type,
      direction = "DOWN"
    )
  } else if (grepl("_ALL\\.rds$", file_name)) {
    list(
      enrichment_type = enrichment_type,
      direction = "ALL"
    )
  } else {
    # Default to ALL if no direction specified
    list(
      enrichment_type = enrichment_type,
      direction = "ALL"
    )
  }
}

# Process all files
all_enrichment <- list()
failed_files <- c()

cat("\nProcessing enrichment files...\n")

pb <- txtProgressBar(min = 0, max = length(enrichment_files), style = 3)

for (i in seq_along(enrichment_files)) {
  file_path <- enrichment_files[i]
  
  tryCatch({
    # Extract metadata
    metadata <- extract_metadata(file_path)
    
    if (!is.null(metadata)) {
      # Parse enrichment type and direction
      enrich_info <- parse_enrichment_info(metadata$file_name, metadata$enrichment_type)
      
      # Read enrichment data
      enrich_data <- readRDS(file_path)
      
      # Convert to data frame if necessary
      if (!is.data.frame(enrich_data)) {
        if (is.list(enrich_data) && length(enrich_data) > 0) {
          # Try to extract the first element if it's a list
          if (is.data.frame(enrich_data[[1]])) {
            enrich_data <- enrich_data[[1]]
          } else {
            next
          }
        } else {
          next
        }
      }
      
      # Add metadata columns
      if (nrow(enrich_data) > 0) {
        enrich_data$method <- metadata$method
        enrich_data$mutation_perturbation <- metadata$mutation_perturbation
        enrich_data$cluster <- metadata$cluster
        enrich_data$experiment <- metadata$experiment
        enrich_data$enrichment_type <- enrich_info$enrichment_type
        enrich_data$direction <- enrich_info$direction
        
        # Add to list
        all_enrichment[[length(all_enrichment) + 1]] <- enrich_data
      }
    }
    
  }, error = function(e) {
    failed_files <- c(failed_files, file_path)
  })
  
  setTxtProgressBar(pb, i)
}

close(pb)

cat("\n\nCombining all enrichment data...\n")

# Combine all data frames
if (length(all_enrichment) > 0) {
  consolidated_data <- bind_rows(all_enrichment)
  
  # Standardize column names
  if ("ID" %in% names(consolidated_data) && !"term_id" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(term_id = ID)
  }
  if ("Description" %in% names(consolidated_data) && !"term_description" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(term_description = Description)
  }
  if ("pvalue" %in% names(consolidated_data) && !"p.value" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(p.value = pvalue)
  }
  if ("qvalue" %in% names(consolidated_data) && !"p.adjust" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(p.adjust = qvalue)
  }
  if ("GeneRatio" %in% names(consolidated_data) && !"gene_ratio" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(gene_ratio = GeneRatio)
  }
  if ("Count" %in% names(consolidated_data) && !"count" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>% rename(count = Count)
  }
  
  # Filter for significant results
  if ("p.adjust" %in% names(consolidated_data)) {
    consolidated_data <- consolidated_data %>%
      filter(p.adjust <= 0.05)
  }
  
  # Add gene column if missing (for compatibility)
  if (!"gene" %in% names(consolidated_data) && "mutation_perturbation" %in% names(consolidated_data)) {
    consolidated_data$gene <- consolidated_data$mutation_perturbation
  }
  
  # Summary statistics
  cat("\n========================================\n")
  cat("Consolidation Summary\n")
  cat("========================================\n\n")
  
  cat(sprintf("Total enrichment terms: %d\n", nrow(consolidated_data)))
  cat(sprintf("Unique mutations: %d\n", length(unique(consolidated_data$mutation_perturbation))))
  cat(sprintf("Unique clusters: %d\n", length(unique(consolidated_data$cluster))))
  cat(sprintf("Failed files: %d\n", length(failed_files)))
  
  # Breakdown by enrichment type
  cat("\nEnrichment types:\n")
  enrich_summary <- consolidated_data %>%
    group_by(enrichment_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    arrange(desc(count))
  
  for (i in 1:nrow(enrich_summary)) {
    cat(sprintf("  %s: %d terms\n", enrich_summary$enrichment_type[i], enrich_summary$count[i]))
  }
  
  # Breakdown by direction
  cat("\nDirections:\n")
  dir_summary <- consolidated_data %>%
    group_by(direction) %>%
    summarise(count = n(), .groups = "drop")
  
  for (i in 1:nrow(dir_summary)) {
    cat(sprintf("  %s: %d terms\n", dir_summary$direction[i], dir_summary$count[i]))
  }
  
  # Save consolidated data
  cat(sprintf("\nSaving to: %s\n", output_file))
  saveRDS(consolidated_data, output_file)
  
  # Check file size
  file_size <- file.info(output_file)$size / 1024^2  # MB
  cat(sprintf("Output file size: %.1f MB\n", file_size))
  
  cat("\n✅ Enrichment consolidation complete!\n")
  
} else {
  cat("\n❌ No enrichment data found to consolidate\n")
}

if (length(failed_files) > 0) {
  cat("\n⚠️  Failed to process these files:\n")
  for (f in failed_files) {
    cat("  -", f, "\n")
  }
}