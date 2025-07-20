# Fast PD Signature Discovery Script
# Optimized version for large datasets
# Date: 2025-07-19

library(dplyr)

# Configuration
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
all_data <- readRDS(data_file)
cat("Loaded", nrow(all_data), "rows\n")

# Add gene column if needed
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Quick PD keyword search (simplified)
pd_terms <- c("mitochondr", "lysosom", "autophagy", "synap", "dopamin", 
              "proteasom", "synuclein", "oxidative", "microglia", "calcium",
              "apoptosis", "ER stress", "unfolded protein", "ubiquitin")

# Filter for PD-relevant terms FIRST (much faster)
cat("Filtering for PD-relevant pathways...\n")
pd_pattern <- paste(pd_terms, collapse = "|")
pd_data <- all_data %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(pd_pattern, Description, ignore.case = TRUE))

cat("Found", nrow(pd_data), "PD-relevant pathways\n")

# Separate by method
mast_data <- pd_data %>% filter(method == "MAST")
mixscale_data <- pd_data %>% filter(method == "MixScale")

cat("\nBreakdown:\n")
cat("- MAST:", nrow(mast_data), "enrichments\n")
cat("- MixScale:", nrow(mixscale_data), "enrichments\n")

# Find top signatures by frequency across samples
get_top_by_frequency <- function(data, n_top = 20) {
  data %>%
    group_by(Description, enrichment_type) %>%
    summarise(
      n_genes = n_distinct(gene),
      n_clusters = n_distinct(cluster),
      genes = paste(unique(gene)[1:min(3, length(unique(gene)))], collapse = ", "),
      mean_neg_log_p = mean(-log10(p.adjust)),
      mean_fold = mean(FoldEnrichment, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(n_genes >= 2 | n_clusters >= 2) %>%
    arrange(desc(n_genes), desc(mean_neg_log_p)) %>%
    head(n_top)
}

# Get top signatures
cat("\nIdentifying top signatures...\n")

# 1. MAST-only
mast_terms <- unique(mast_data$Description)
mixscale_terms <- unique(mixscale_data$Description)
mast_only_terms <- setdiff(mast_terms, mixscale_terms)

mast_top <- mast_data %>%
  filter(Description %in% mast_only_terms) %>%
  get_top_by_frequency(30)

# 2. MixScale-only  
mixscale_only_terms <- setdiff(mixscale_terms, mast_terms)

mixscale_top <- mixscale_data %>%
  filter(Description %in% mixscale_only_terms) %>%
  get_top_by_frequency(30)

# 3. Convergent
convergent_terms <- intersect(mast_terms, mixscale_terms)
cat("\nFound", length(convergent_terms), "convergent pathways\n")

convergent_stats <- pd_data %>%
  filter(Description %in% convergent_terms) %>%
  group_by(Description, enrichment_type) %>%
  summarise(
    n_genes_mast = n_distinct(gene[method == "MAST"]),
    n_genes_mixscale = n_distinct(gene[method == "MixScale"]),
    genes_mast = paste(unique(gene[method == "MAST"])[1:min(3, length(unique(gene[method == "MAST"])))], collapse = ", "),
    genes_mixscale = paste(unique(gene[method == "MixScale"])[1:min(3, length(unique(gene[method == "MixScale"])))], collapse = ", "),
    mean_neg_log_p = mean(-log10(p.adjust)),
    .groups = "drop"
  ) %>%
  filter(n_genes_mast > 0, n_genes_mixscale > 0) %>%
  arrange(desc(n_genes_mast + n_genes_mixscale), desc(mean_neg_log_p)) %>%
  head(30)

# Save results
write.csv(mast_top, file.path(output_dir, "mast_top_fast.csv"), row.names = FALSE)
write.csv(mixscale_top, file.path(output_dir, "mixscale_top_fast.csv"), row.names = FALSE)
write.csv(convergent_stats, file.path(output_dir, "convergent_top_fast.csv"), row.names = FALSE)

# Print top 3 signatures
cat("\n=== TOP 3 COMPELLING SIGNATURES ===\n")

if (nrow(mast_top) > 0) {
  cat("\n1. TOP MAST-ONLY:\n")
  cat("   ", mast_top$Description[1], "\n")
  cat("   Found in", mast_top$n_genes[1], "genes:", mast_top$genes[1], "...\n")
  cat("   Mean -log10(p):", round(mast_top$mean_neg_log_p[1], 2), "\n")
}

if (nrow(mixscale_top) > 0) {
  cat("\n2. TOP MIXSCALE-ONLY:\n")
  cat("   ", mixscale_top$Description[1], "\n")
  cat("   Found in", mixscale_top$n_genes[1], "genes:", mixscale_top$genes[1], "...\n")
  cat("   Mean -log10(p):", round(mixscale_top$mean_neg_log_p[1], 2), "\n")
}

if (nrow(convergent_stats) > 0) {
  cat("\n3. TOP CONVERGENT:\n")
  cat("   ", convergent_stats$Description[1], "\n")
  cat("   MAST genes:", convergent_stats$n_genes_mast[1], "-", convergent_stats$genes_mast[1], "...\n")
  cat("   MixScale genes:", convergent_stats$n_genes_mixscale[1], "-", convergent_stats$genes_mixscale[1], "...\n")
  cat("   Mean -log10(p):", round(convergent_stats$mean_neg_log_p[1], 2), "\n")
}

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")