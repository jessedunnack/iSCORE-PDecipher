# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)
mutation <- args[1]

cat(sprintf("Processing mutation: %s\n", mutation))

# Load libraries
library(Seurat)

# Load and prepare data
cat("Loading Seurat object...\n")
all <- readRDS("full_dataset.rds")

# Let's process the "coarse" cell type clusters:
all <- FindClusters(all, resolution = 0.2)

# Determine which batch(es) the mutation is in
mutation_batches <- unique(all$batch[all$mutation_tidy == mutation])
cat(sprintf("Mutation %s is in batch(es): %s\n", mutation, paste(mutation_batches, collapse=", ")))

# Create a subset that includes only cells from the mutation and eWT cells from the same batch(es)
cells_to_keep <- all$mutation_tidy == mutation | (all$mutation_tidy == "eWT" & all$batch %in% mutation_batches)
all_filtered <- subset(all, cells = colnames(all)[cells_to_keep])

# Determine if we should include batch as a latent variable
remaining_batches <- unique(all_filtered$batch)
use_batch_as_latent <- length(remaining_batches) > 1
if (use_batch_as_latent) {
  latent_vars <- c("subclone_ID", "batch")
  cat("Using both 'subclone_ID' and 'batch' as latent variables\n")
} else {
  latent_vars <- "subclone_ID"
  cat("Using only 'subclone_ID' as latent variable (single batch detected)\n")
}

# First determine which clusters exist in the filtered data
existing_clusters <- sort(as.numeric(as.character(unique(Idents(all_filtered)))))
cat(sprintf("Clusters present in filtered data: %s\n", paste(existing_clusters, collapse=", ")))

# Check if we have sufficient data for each existing cluster
valid_clusters <- c()
for (cluster in existing_clusters) {
  cells_in_cluster <- WhichCells(all_filtered, idents = cluster)
  
  has_mutation <- any(all_filtered$mutation_tidy[cells_in_cluster] == mutation)
  has_ewt <- any(all_filtered$mutation_tidy[cells_in_cluster] == "eWT")
  
  if (!has_mutation || !has_ewt) {
    cat(sprintf("Skipping cluster %s: Missing either mutation or eWT cells\n", cluster))
    next
  }
  
  valid_clusters <- c(valid_clusters, cluster)
}

cat(sprintf("Valid clusters for analysis: %s\n", paste(valid_clusters, collapse=", ")))

# Initialize results list for this mutation
mutation_results <- list()

# Create metadata first
mutation_results$metadata <- list(
  date = as.character(Sys.Date()),
  control = "eWT",
  mutation = mutation,
  test = "MAST",
  latent_vars = latent_vars,
  batches_used = mutation_batches,
  existing_clusters = existing_clusters,
  valid_clusters = valid_clusters,
  assay = "RNA"
)

# Create output directory if it doesn't exist
output_dir <- "./iSCORE-PD_MAST_analysis"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat(sprintf("Created output directory: %s\n", output_dir))
}

if (length(valid_clusters) == 0) {
  cat("No valid clusters found for comparison. Saving empty results.\n")
  mutation_results$metadata$error <- "No valid clusters found for comparison"
  
  output_file <- file.path(output_dir, paste0("mutation_", gsub("[^a-zA-Z0-9]", "_", mutation), "_results.rds"))
  saveRDS(mutation_results, output_file)
  cat(sprintf("Empty results saved to %s\n", output_file))
  quit(save = "no", status = 0)
}

# Process each valid cluster
for (cluster in valid_clusters) {
  cat(sprintf("Processing cluster: %s\n", cluster))
  
  # Try-catch to handle potential errors
  result <- tryCatch({
    # Run differential expression
    markers <- FindMarkers(all_filtered, 
                          subset.ident = as.character(cluster),
                          ident.1 = mutation,
                          ident.2 = "eWT", 
                          group.by = "mutation_tidy",
                          latent.vars = latent_vars,
                          assay = "RNA",
			  slot = "data",
                          test.use = "MAST")
    
    # Return the markers
    markers
    
  }, error = function(e) {
    cat(sprintf("Error in cluster: %s - %s\n", cluster, e$message))
    # Return empty dataframe with error message
    data.frame(error = e$message)
  })
  
  # Store the result
  mutation_results[[paste0("cluster_", cluster)]] <- result
  
  # Force garbage collection
  gc()
}

# Save results for this mutation
output_file <- file.path(output_dir, paste0("mutation_", gsub("[^a-zA-Z0-9]", "_", mutation), "_results.rds"))
saveRDS(mutation_results, output_file)
cat(sprintf("Results saved to %s\n", output_file))
