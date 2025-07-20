# PD Signature Discovery Script
# Purpose: Identify the most compelling functional enrichment signatures for PD
# Categories: 1) iSCORE-PD only (MAST), 2) CRISPRi only, 3) Convergent signatures
# Author: iSCORE-PDecipher Analysis Pipeline
# Date: 2025-07-19

# Load required libraries
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
})

# =============================================================================
# CONFIGURATION
# =============================================================================

# Set paths
# Using the dataset with both iSCORE-PD and CRISPRi data (no CRISPRa)
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Define comprehensive PD pathway keywords
pd_keywords <- list(
  mitochondrial = c("mitochondr", "oxidative phosphorylation", "electron transport", 
                    "complex I", "ATP synth", "OXPHOS", "respiratory chain", 
                    "mitophagy", "PINK1", "Parkin", "cristae", "mitochondrial fission",
                    "mitochondrial fusion", "cytochrome c"),
  
  lysosomal = c("lysosom", "autophagy", "endocytosis", "protein degradation", 
                "cathepsin", "glucocerebrosidase", "GBA", "autophagic", "endosom",
                "vacuolar", "lysosomal storage", "macroautophagy", "chaperone-mediated"),
  
  synaptic = c("synap", "dopamin", "neurotransmitter", "synaptic vesicle", 
               "SNARE", "vesicle fusion", "presynaptic", "postsynaptic", 
               "neurotransmission", "exocytosis", "endocytosis", "clathrin",
               "synaptic plasticity", "long-term potentiation"),
  
  protein_aggregation = c("proteasom", "unfolded protein", "alpha-synuclein", 
                         "synuclein", "Lewy body", "aggregat", "inclusion", 
                         "ubiquitin", "misfolded", "protein folding", "chaperone",
                         "heat shock", "HSP", "UPS", "ubiquitin-proteasome"),
  
  neuroinflammation = c("microglia", "neuroinflam", "cytokine", "TNF", 
                       "interleukin", "IL-1", "IL-6", "inflammasome", 
                       "astrocyte activation", "gliosis", "complement",
                       "toll-like receptor", "TLR", "NF-kappa"),
  
  oxidative_stress = c("reactive oxygen", "oxidative stress", "antioxidant", 
                      "glutathione", "peroxidase", "superoxide", "catalase",
                      "nitrosative stress", "lipid peroxidation", "8-oxo",
                      "oxidative damage", "redox", "NAD+", "NADPH"),
  
  er_calcium = c("endoplasmic reticulum stress", "ER stress", "calcium homeostasis",
                "calnexin", "calreticulin", "UPR", "PERK", "ATF6", "IRE1",
                "calcium signaling", "calcium channel", "calcium binding",
                "sarcoplasmic reticulum", "ryanodine receptor"),
  
  cell_death = c("apoptosis", "apoptotic", "caspase", "cell death", "necrosis",
                "ferroptosis", "programmed cell death", "BCL2", "BAX", "cytochrome c",
                "death receptor", "TRAIL", "FAS", "survival")
)

# Flatten keywords for searching
all_pd_keywords <- unlist(pd_keywords)

# =============================================================================
# FUNCTIONS
# =============================================================================

#' Score pathway for PD relevance
#' @param description Character string of pathway description
#' @param keywords Character vector of PD-related keywords
#' @return Numeric score (0-1) indicating PD relevance
score_pd_relevance <- function(description, keywords = all_pd_keywords) {
  description_lower <- tolower(description)
  matches <- sapply(keywords, function(kw) grepl(kw, description_lower, ignore.case = TRUE))
  
  # Calculate relevance score
  n_matches <- sum(matches)
  relevance_score <- min(n_matches / 5, 1)  # Normalize to 0-1, cap at 1
  
  return(relevance_score)
}

#' Categorize pathway into PD biological category
#' @param description Pathway description
#' @return Character vector of matched categories
categorize_pd_pathway <- function(description) {
  description_lower <- tolower(description)
  categories <- character()
  
  for (cat_name in names(pd_keywords)) {
    cat_keywords <- pd_keywords[[cat_name]]
    if (any(sapply(cat_keywords, function(kw) grepl(kw, description_lower)))) {
      categories <- c(categories, cat_name)
    }
  }
  
  if (length(categories) == 0) categories <- "other"
  return(paste(categories, collapse = ";"))
}

#' Calculate signature strength score
#' @param data Data frame with enrichment results
#' @return Data frame with additional strength metrics
calculate_signature_strength <- function(data) {
  data %>%
    mutate(
      # Combined significance and effect size score
      strength_score = -log10(p.adjust) * abs(log2(FoldEnrichment + 1)),
      
      # Normalized scores for comparison
      pval_score = -log10(p.adjust) / max(-log10(p.adjust), na.rm = TRUE),
      fold_score = abs(log2(FoldEnrichment + 1)) / max(abs(log2(FoldEnrichment + 1)), na.rm = TRUE),
      
      # Combined metric
      combined_score = (pval_score + fold_score) / 2
    )
}

# =============================================================================
# LOAD AND PREPARE DATA
# =============================================================================

cat("Loading consolidated enrichment data...\n")

if (!file.exists(data_file)) {
  stop("Enrichment data file not found at: ", data_file)
}

all_data <- readRDS(data_file)
cat("Loaded", nrow(all_data), "enrichment results\n")

# Add gene column if missing
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Filter for significant results
sig_data <- all_data %>%
  filter(p.adjust < 0.05)

cat("Filtered to", nrow(sig_data), "significant results\n")

# =============================================================================
# IDENTIFY PD-RELEVANT SIGNATURES
# =============================================================================

cat("\nScoring pathways for PD relevance...\n")

# Add PD relevance scores and categories
sig_data <- sig_data %>%
  mutate(
    pd_relevance_score = sapply(Description, score_pd_relevance),
    pd_categories = sapply(Description, categorize_pd_pathway)
  )

# Calculate signature strength
sig_data <- calculate_signature_strength(sig_data)

# Filter for PD-relevant pathways
pd_data <- sig_data %>%
  filter(pd_relevance_score > 0)

cat("Identified", nrow(pd_data), "PD-relevant pathways\n")

# =============================================================================
# SEPARATE BY METHOD
# =============================================================================

# MAST-only signatures
mast_data <- pd_data %>%
  filter(method == "MAST")

# MixScale (CRISPRi) signatures
mixscale_data <- pd_data %>%
  filter(method == "MixScale")

# Find convergent signatures (present in both methods)
mast_terms <- unique(mast_data$Description)
mixscale_terms <- unique(mixscale_data$Description)
convergent_terms <- intersect(mast_terms, mixscale_terms)

cat("\nMethod breakdown:\n")
cat("- MAST-only terms:", length(setdiff(mast_terms, mixscale_terms)), "\n")
cat("- MixScale-only terms:", length(setdiff(mixscale_terms, mast_terms)), "\n")
cat("- Convergent terms:", length(convergent_terms), "\n")

# =============================================================================
# IDENTIFY TOP SIGNATURES BY CATEGORY
# =============================================================================

#' Get top signatures for a dataset
#' @param data Enrichment data
#' @param n_top Number of top signatures to return
#' @param min_genes Minimum number of genes a signature must appear in
#' @param min_clusters Minimum number of clusters a signature must appear in
get_top_signatures <- function(data, n_top = 20, min_genes = 2, min_clusters = 1) {
  # Calculate cross-sample frequency
  cross_sample_stats <- data %>%
    group_by(Description, enrichment_type, pd_categories) %>%
    summarise(
      n_genes = n_distinct(gene),
      n_clusters = n_distinct(cluster),
      genes = paste(unique(gene), collapse = ", "),
      clusters = paste(unique(cluster), collapse = ", "),
      mean_pval = mean(-log10(p.adjust)),
      mean_fold = mean(log2(FoldEnrichment + 1)),
      max_strength = max(strength_score),
      mean_pd_relevance = mean(pd_relevance_score),
      directions = paste(unique(direction), collapse = ";"),
      .groups = "drop"
    )
  
  # Filter by cross-sample criteria
  filtered <- cross_sample_stats %>%
    filter(n_genes >= min_genes | n_clusters >= min_clusters)
  
  # Rank by combined criteria
  ranked <- filtered %>%
    mutate(
      # Composite ranking score
      rank_score = (
        0.3 * (mean_pval / max(mean_pval)) +
        0.2 * (mean_fold / max(mean_fold)) +
        0.2 * (n_genes / max(n_genes)) +
        0.2 * (n_clusters / max(n_clusters)) +
        0.1 * mean_pd_relevance
      )
    ) %>%
    arrange(desc(rank_score))
  
  return(head(ranked, n_top))
}

# Get top signatures for each category
cat("\nIdentifying top signatures...\n")

# 1. MAST-only top signatures
mast_only_data <- pd_data %>%
  filter(method == "MAST", !Description %in% convergent_terms)

mast_top <- get_top_signatures(mast_only_data, n_top = 30, min_genes = 2)

# 2. MixScale-only top signatures
mixscale_only_data <- pd_data %>%
  filter(method == "MixScale", !Description %in% convergent_terms)

mixscale_top <- get_top_signatures(mixscale_only_data, n_top = 30, min_genes = 1)

# Analyze by individual CRISPRi experiments
if ("experiment" %in% names(mixscale_data)) {
  experiments <- unique(mixscale_data$experiment)
  cat("\nCRISPRi experiments found:", paste(experiments, collapse = ", "), "\n")
  
  exp_summaries <- list()
  for (exp in experiments) {
    exp_data <- mixscale_only_data %>% filter(experiment == exp)
    if (nrow(exp_data) > 0) {
      exp_summaries[[exp]] <- get_top_signatures(exp_data, n_top = 10, min_genes = 1, min_clusters = 1)
    }
  }
}

# 3. Convergent signatures
convergent_data <- pd_data %>%
  filter(Description %in% convergent_terms)

convergent_stats <- convergent_data %>%
  group_by(Description, enrichment_type, pd_categories) %>%
  summarise(
    n_methods = n_distinct(method),
    n_genes_mast = n_distinct(gene[method == "MAST"]),
    n_genes_mixscale = n_distinct(gene[method == "MixScale"]),
    n_clusters_total = n_distinct(paste(method, cluster)),
    genes_mast = paste(unique(gene[method == "MAST"]), collapse = ", "),
    genes_mixscale = paste(unique(gene[method == "MixScale"]), collapse = ", "),
    mean_pval_mast = mean(-log10(p.adjust[method == "MAST"])),
    mean_pval_mixscale = mean(-log10(p.adjust[method == "MixScale"])),
    mean_fold_mast = mean(log2(FoldEnrichment[method == "MAST"] + 1)),
    mean_fold_mixscale = mean(log2(FoldEnrichment[method == "MixScale"] + 1)),
    mean_pd_relevance = mean(pd_relevance_score),
    .groups = "drop"
  ) %>%
  filter(n_methods == 2) %>%  # Ensure present in both methods
  mutate(
    # Calculate convergence strength
    convergence_score = (
      0.25 * (mean_pval_mast / max(mean_pval_mast)) +
      0.25 * (mean_pval_mixscale / max(mean_pval_mixscale)) +
      0.2 * (n_genes_mast / max(n_genes_mast)) +
      0.2 * (n_genes_mixscale / max(n_genes_mixscale)) +
      0.1 * mean_pd_relevance
    )
  ) %>%
  arrange(desc(convergence_score))

convergent_top <- head(convergent_stats, 30)

# =============================================================================
# SAVE RESULTS
# =============================================================================

cat("\nSaving results...\n")

# Save detailed signature lists
write.csv(mast_top, file.path(output_dir, "mast_only_top_signatures.csv"), row.names = FALSE)
write.csv(mixscale_top, file.path(output_dir, "mixscale_only_top_signatures.csv"), row.names = FALSE)
write.csv(convergent_top, file.path(output_dir, "convergent_top_signatures.csv"), row.names = FALSE)

# Save full PD-relevant data for further analysis
saveRDS(pd_data, file.path(output_dir, "pd_relevant_enrichments.rds"))

# =============================================================================
# GENERATE SUMMARY REPORT
# =============================================================================

cat("\nGenerating summary report...\n")

# Create summary statistics
summary_stats <- data.frame(
  Category = c("MAST-only", "MixScale-only", "Convergent"),
  Total_Terms = c(
    length(setdiff(mast_terms, mixscale_terms)),
    length(setdiff(mixscale_terms, mast_terms)),
    length(convergent_terms)
  ),
  Top_Filtered = c(nrow(mast_top), nrow(mixscale_top), nrow(convergent_top)),
  Avg_Genes_Per_Term = c(
    round(mean(mast_top$n_genes), 2),
    round(mean(mixscale_top$n_genes), 2),
    round(mean(convergent_top$n_genes_mast + convergent_top$n_genes_mixscale), 2)
  ),
  Top_Category = c(
    names(sort(table(unlist(strsplit(mast_top$pd_categories, ";"))), decreasing = TRUE)[1]),
    names(sort(table(unlist(strsplit(mixscale_top$pd_categories, ";"))), decreasing = TRUE)[1]),
    names(sort(table(unlist(strsplit(convergent_top$pd_categories, ";"))), decreasing = TRUE)[1])
  )
)

write.csv(summary_stats, file.path(output_dir, "summary_statistics.csv"), row.names = FALSE)

# =============================================================================
# IDENTIFY TOP 3 COMPELLING SIGNATURES
# =============================================================================

cat("\n=== TOP 3 COMPELLING SIGNATURES FOR THESIS COMMITTEE ===\n")

# 1. Best MAST-only signature
cat("\n1. MAST-ONLY SIGNATURE:\n")
mast_best <- mast_top[1,]
cat("   Term:", mast_best$Description, "\n")
cat("   Category:", mast_best$pd_categories, "\n")
cat("   Present in", mast_best$n_genes, "genes:", mast_best$genes, "\n")
cat("   Mean -log10(p):", round(mast_best$mean_pval, 3), "\n")
cat("   Biological interpretation: Genetic mutations converge on", 
    gsub("_", " ", mast_best$pd_categories), "dysfunction\n")

# 2. Best MixScale-only signature
cat("\n2. MIXSCALE-ONLY SIGNATURE:\n")
mixscale_best <- mixscale_top[1,]
cat("   Term:", mixscale_best$Description, "\n")
cat("   Category:", mixscale_best$pd_categories, "\n")
cat("   Present in", mixscale_best$n_genes, "genes:", mixscale_best$genes, "\n")
cat("   Mean -log10(p):", round(mixscale_best$mean_pval, 3), "\n")
cat("   Biological interpretation: Acute knockdowns reveal compensatory", 
    gsub("_", " ", mixscale_best$pd_categories), "responses\n")

# 3. Best convergent signature
cat("\n3. CONVERGENT SIGNATURE:\n")
conv_best <- convergent_top[1,]
cat("   Term:", conv_best$Description, "\n")
cat("   Category:", conv_best$pd_categories, "\n")
cat("   MAST genes (", conv_best$n_genes_mast, "):", conv_best$genes_mast, "\n")
cat("   MixScale genes (", conv_best$n_genes_mixscale, "):", conv_best$genes_mixscale, "\n")
cat("   Mean -log10(p) MAST:", round(conv_best$mean_pval_mast, 3), 
    "MixScale:", round(conv_best$mean_pval_mixscale, 3), "\n")
cat("   Biological interpretation: Core", gsub("_", " ", conv_best$pd_categories), 
    "axis disrupted across perturbation methods\n")

# Save top 3 for easy access
top3_signatures <- list(
  mast = mast_best,
  mixscale = mixscale_best,
  convergent = conv_best
)

saveRDS(top3_signatures, file.path(output_dir, "top3_signatures.rds"))

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")
cat("\nNext steps:\n")
cat("1. Run pd_signature_visualization.R to generate heatmaps\n")
cat("2. Review top signatures and refine biological interpretations\n")
cat("3. Prepare presentation slides with key findings\n")