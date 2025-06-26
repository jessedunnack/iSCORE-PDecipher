# Script to explore the enrichment results RDS file structure
# For iSCORE-PDecipher data exploration

# Load required libraries
library(dplyr)

# Set path to your enrichment results file
# Update this path to match your local setup
enrichment_path <- "path/to/your/dataset/all_enrichment_padj005_complete_with_direction.rds"

# Load the enrichment results
enrichment_data <- readRDS(enrichment_path)

# Basic overview
cat("=== Enrichment Results Overview ===\n")
cat("Total enrichment terms:", nrow(enrichment_data), "\n")
cat("Columns available:", paste(colnames(enrichment_data), collapse = ", "), "\n\n")

# Key columns explained:
# - mutation_perturbation: The gene (e.g., LRRK2, PINK1)
# - cluster: Cell cluster (e.g., cluster_0, cluster_1)
# - enrichment_type: Type of analysis (GO_BP, KEGG, Reactome, etc.)
# - direction: Gene regulation direction (UP, DOWN, ALL)
# - Description: Human-readable pathway/term name
# - p.adjust: Adjusted p-value (significance)
# - method: MAST or MixScale
# - experiment: Specific experiment ID for CRISPRi/a

# Explore available methods
cat("=== Methods/Modalities ===\n")
print(table(enrichment_data$method))

# Explore enrichment types
cat("\n=== Enrichment Types ===\n")
print(table(enrichment_data$enrichment_type))

# Explore genes
cat("\n=== Genes Analyzed ===\n")
unique_genes <- unique(enrichment_data$mutation_perturbation)
cat("Total unique genes:", length(unique_genes), "\n")
cat("Examples:", paste(head(unique_genes, 10), collapse = ", "), "\n")

# Explore clusters
cat("\n=== Clusters ===\n")
unique_clusters <- unique(enrichment_data$cluster)
cat("Clusters:", paste(sort(unique_clusters), collapse = ", "), "\n")

# Example: Find top enriched pathways for a specific gene/cluster
find_top_pathways <- function(data, gene, cluster, enrichment_type = "GO_BP", 
                            direction = "ALL", top_n = 10) {
  filtered <- data %>%
    filter(mutation_perturbation == gene,
           cluster == cluster,
           enrichment_type == enrichment_type,
           direction == direction) %>%
    arrange(p.adjust) %>%
    head(top_n) %>%
    select(Description, p.adjust, direction, method)
  
  return(filtered)
}

# Example usage
cat("\n=== Example: Top GO pathways for LRRK2 in cluster_0 ===\n")
top_pathways <- find_top_pathways(enrichment_data, "LRRK2", "cluster_0")
print(top_pathways)

# Compare enrichment between methods
compare_methods <- function(data, gene, enrichment_type = "KEGG") {
  comparison <- data %>%
    filter(mutation_perturbation == gene,
           enrichment_type == enrichment_type) %>%
    group_by(method, direction) %>%
    summarise(
      n_terms = n(),
      avg_pval = mean(p.adjust),
      min_pval = min(p.adjust),
      .groups = "drop"
    )
  
  return(comparison)
}

# Example comparison
cat("\n=== Example: Compare KEGG enrichment across methods ===\n")
# method_comparison <- compare_methods(enrichment_data, "LRRK2")
# print(method_comparison)

# Quick function to find common pathways between genes
find_common_pathways <- function(data, genes, cluster = "cluster_0", 
                               enrichment_type = "GO_BP", top_n = 20) {
  # Get top pathways for each gene
  all_pathways <- list()
  
  for (gene in genes) {
    pathways <- data %>%
      filter(mutation_perturbation == gene,
             cluster == cluster,
             enrichment_type == enrichment_type,
             p.adjust < 0.05) %>%
      arrange(p.adjust) %>%
      head(top_n) %>%
      pull(Description)
    
    all_pathways[[gene]] <- pathways
  }
  
  # Find intersection
  common <- Reduce(intersect, all_pathways)
  return(common)
}

# Example: Find pathways common to multiple PD genes
cat("\n=== Example: Common pathways between PD genes ===\n")
pd_genes <- c("LRRK2", "PINK1", "PRKN")
# common_pathways <- find_common_pathways(enrichment_data, pd_genes)
# cat("Common pathways:", paste(head(common_pathways, 5), collapse = "\n"), "\n")