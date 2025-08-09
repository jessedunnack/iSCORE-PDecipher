#!/usr/bin/env Rscript

# VALIDATE CLUSTER HIERARCHY ASSIGNMENTS
# Prove that the cluster assignments are correctly mapped to cells

library(dplyr)
library(ggplot2)

cat("=================================================================\n")
cat("VALIDATING CLUSTER HIERARCHY ASSIGNMENTS\n")
cat("=================================================================\n\n")

# 1. Load the saved metadata and mapping
cat("1. Loading cluster hierarchy data...\n")
mapping <- read.csv("results/cluster_hierarchy/fine_to_coarse_mapping.csv")
coarse_summary <- read.csv("results/cluster_hierarchy/coarse_cluster_summary.csv")
metadata <- readRDS("results/cluster_hierarchy/cluster_metadata.rds")

cat("Total cells in metadata:", nrow(metadata), "\n")
cat("Fine clusters:", length(unique(metadata$seurat_clusters_fine)), "\n")
cat("Coarse clusters:", length(unique(metadata$seurat_clusters_coarse)), "\n\n")

# 2. Validate fine-to-coarse mapping consistency
cat("2. Validating fine-to-coarse mapping consistency...\n")

# Check that every cell's fine cluster maps to the correct coarse cluster
# Convert to character for joining
metadata$fine_cluster_char <- as.character(metadata$seurat_clusters_fine)
mapping$fine_cluster_char <- as.character(mapping$fine_cluster)

validation_result <- metadata %>%
  left_join(mapping, by = c("fine_cluster_char" = "fine_cluster_char")) %>%
  mutate(mapping_correct = (as.character(seurat_clusters_coarse) == as.character(coarse_cluster))) %>%
  group_by(seurat_clusters_fine) %>%
  summarise(
    n_cells = n(),
    n_correct = sum(mapping_correct),
    pct_correct = round(100 * n_correct / n_cells, 2),
    actual_coarse = first(seurat_clusters_coarse),
    expected_coarse = first(coarse_cluster)
  )

cat("\nFine cluster validation results:\n")
print(validation_result)

# Check if any mismatches
mismatches <- validation_result %>% filter(pct_correct < 100)
if (nrow(mismatches) > 0) {
  cat("\nWARNING: Found mismatches in fine-to-coarse mapping:\n")
  print(mismatches)
} else {
  cat("\n✓ All fine clusters correctly map to their coarse clusters!\n")
}

# 3. Validate that coarse clusters contain the expected fine clusters
cat("\n3. Validating coarse cluster composition...\n")

coarse_validation <- metadata %>%
  group_by(seurat_clusters_coarse) %>%
  summarise(
    n_cells_actual = n(),
    fine_clusters_actual = paste(sort(unique(seurat_clusters_fine)), collapse = ", ")
  ) %>%
  left_join(
    coarse_summary %>% 
      select(seurat_clusters_coarse, n_cells_expected = n_cells, 
             fine_clusters_expected = fine_clusters),
    by = "seurat_clusters_coarse"
  ) %>%
  mutate(
    cells_match = (n_cells_actual == n_cells_expected),
    clusters_match = (fine_clusters_actual == fine_clusters_expected)
  )

cat("\nCoarse cluster validation:\n")
print(select(coarse_validation, seurat_clusters_coarse, n_cells_actual, 
             n_cells_expected, cells_match, clusters_match))

# 4. Load marker data to validate biological sense
cat("\n4. Loading marker data to validate biological groupings...\n")

fine_markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")
coarse_markers <- readRDS("inst/extdata/umap_data/iSCORE_PD_CRISPRi_cluster_markers.rds")
coarse_markers$cluster <- gsub("cluster_", "", coarse_markers$cluster)

# 5. For each coarse cluster, check if fine clusters share similar markers
cat("\n5. Validating biological coherence of cluster groupings...\n\n")

# Function to get top markers
get_top_markers <- function(markers_df, cluster_id, n = 10) {
  markers_df %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n) %>%
    pull(gene)
}

# Analyze coherence for each coarse cluster
for (coarse_id in 0:14) {
  cat(sprintf("\n=== COARSE CLUSTER %d ===\n", coarse_id))
  
  # Get fine clusters in this coarse cluster
  fine_clusters <- mapping %>%
    filter(coarse_cluster == coarse_id) %>%
    pull(fine_cluster) %>%
    as.character()
  
  cat("Contains fine clusters:", paste(fine_clusters, collapse = ", "), "\n")
  
  # Get top markers for the coarse cluster
  coarse_top <- get_top_markers(coarse_markers, coarse_id, 5)
  cat("Coarse cluster top markers:", paste(coarse_top, collapse = ", "), "\n")
  
  # Check if fine clusters express similar markers
  cat("\nFine cluster markers:\n")
  marker_overlap <- list()
  
  for (fine_id in fine_clusters) {
    fine_top <- get_top_markers(fine_markers, fine_id, 10)
    overlap <- intersect(coarse_top, fine_top)
    marker_overlap[[as.character(fine_id)]] <- length(overlap)
    
    cat(sprintf("  Fine cluster %d: %s (overlap: %d)\n", 
                fine_id, 
                paste(head(fine_top, 5), collapse = ", "),
                length(overlap)))
  }
  
  # Calculate coherence score
  avg_overlap <- mean(unlist(marker_overlap))
  cat(sprintf("\nCoherence score: %.2f/5 markers shared on average\n", avg_overlap))
}

# 6. Create visualization of cluster hierarchy
cat("\n\n6. Creating cluster hierarchy visualization...\n")

# Create a summary plot
hierarchy_summary <- mapping %>%
  group_by(coarse_cluster) %>%
  summarise(
    n_fine_clusters = n(),
    total_cells = sum(n_cells),
    fine_clusters_list = paste(fine_cluster, collapse = ",")
  )

p1 <- ggplot(hierarchy_summary, aes(x = factor(coarse_cluster), y = n_fine_clusters)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_text(aes(label = n_fine_clusters), vjust = -0.5) +
  labs(x = "Coarse Cluster", y = "Number of Fine Clusters", 
       title = "Fine Clusters per Coarse Cluster") +
  theme_minimal()

p2 <- ggplot(hierarchy_summary, aes(x = factor(coarse_cluster), y = total_cells)) +
  geom_bar(stat = "identity", fill = "darkgreen") +
  geom_text(aes(label = total_cells), vjust = -0.5, size = 3) +
  labs(x = "Coarse Cluster", y = "Total Cells", 
       title = "Cell Count per Coarse Cluster") +
  theme_minimal() +
  scale_y_continuous(labels = scales::comma)

library(patchwork)
p_combined <- p1 + p2
ggsave("results/cluster_hierarchy/hierarchy_validation_plots.pdf", 
       p_combined, width = 12, height = 6)

# 7. Final validation summary
cat("\n\n=== VALIDATION SUMMARY ===\n")
cat("==========================\n")

# Check total cells
total_cells_metadata <- nrow(metadata)
total_cells_mapping <- sum(mapping$n_cells)
cat(sprintf("\nTotal cells in metadata: %d\n", total_cells_metadata))
cat(sprintf("Total cells in mapping: %d\n", total_cells_mapping))
cat(sprintf("Match: %s\n", ifelse(total_cells_metadata == total_cells_mapping, "✓ YES", "✗ NO")))

# Check unique clusters
n_fine_metadata <- length(unique(metadata$seurat_clusters_fine))
n_fine_mapping <- nrow(mapping)
cat(sprintf("\nFine clusters in metadata: %d\n", n_fine_metadata))
cat(sprintf("Fine clusters in mapping: %d\n", n_fine_mapping))
cat(sprintf("Match: %s\n", ifelse(n_fine_metadata == n_fine_mapping, "✓ YES", "✗ NO")))

# Overall validation
all_correct <- all(validation_result$pct_correct == 100)
if (all_correct) {
  cat("\n✓ VALIDATION PASSED: All cluster assignments are correct!\n")
} else {
  cat("\n✗ VALIDATION FAILED: Some cluster assignments have issues.\n")
}

# Save validation report
cat("\nSaving validation report...\n")
write.csv(validation_result, 
          "results/cluster_hierarchy/validation_report.csv", 
          row.names = FALSE)

cat("\nValidation complete! Check results/cluster_hierarchy/ for detailed reports.\n")