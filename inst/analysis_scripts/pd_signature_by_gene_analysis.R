# PD Signature By-Gene Analysis
# Purpose: Analyze each PD gene individually for their pathway signatures
# Date: 2025-07-19

library(dplyr)
library(ggplot2)
library(tidyr)

# Helper function for gene harmonization
harmonize_gene_names <- function(gene_name) {
  gene_map <- c(
    "PRKN" = "PARK2",
    "VPS13C_A444P" = "VPS13C",
    "VPS13C_W395C" = "VPS13C", 
    "SNCA_A30P" = "SNCA",
    "SNCA_A53T" = "SNCA",
    "SNCA_Triplication" = "SNCA"
  )
  if (gene_name %in% names(gene_map)) {
    return(gene_map[gene_name])
  }
  return(gene_name)
}

# Configuration
data_file <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"
output_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/by_gene"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load data
cat("Loading data...\n")
all_data <- readRDS(data_file)

# Add gene column if needed
if (!"gene" %in% names(all_data) && "mutation_perturbation" %in% names(all_data)) {
  all_data$gene <- all_data$mutation_perturbation
}

# Apply gene harmonization
all_data$gene_original <- all_data$gene
all_data$gene_harmonized <- sapply(all_data$gene, harmonize_gene_names)

# PD keyword search
pd_terms <- c("mitochondr", "lysosom", "autophagy", "synap", "dopamin", 
              "proteasom", "synuclein", "oxidative", "microglia", "calcium",
              "apoptosis", "ER stress", "unfolded protein", "ubiquitin")
pd_pattern <- paste(pd_terms, collapse = "|")

# Filter for PD-relevant and significant
pd_data <- all_data %>%
  filter(p.adjust < 0.05) %>%
  filter(grepl(pd_pattern, Description, ignore.case = TRUE))

# Get unique harmonized genes
all_genes <- unique(pd_data$gene_harmonized)
cat("Found", length(all_genes), "harmonized genes to analyze:", paste(all_genes, collapse = ", "), "\n\n")

# Function to analyze a single gene
analyze_gene_signatures <- function(data, gene_name) {
  cat("Analyzing", gene_name, "...\n")
  
  # Filter data for this harmonized gene
  gene_data <- data %>% filter(gene_harmonized == gene_name)
  
  if (nrow(gene_data) == 0) {
    return(list(
      gene = gene_name,
      mast_only = data.frame(),
      mixscale_only = data.frame(),
      convergent = data.frame(),
      summary = data.frame(
        gene = gene_name,
        n_mast_only = 0,
        n_mixscale_only = 0,
        n_convergent = 0,
        n_total = 0
      )
    ))
  }
  
  # Separate by method
  mast_data <- gene_data %>% filter(method == "MAST")
  mixscale_data <- gene_data %>% filter(method == "MixScale")
  
  # Find unique and convergent terms
  mast_terms <- unique(mast_data$Description)
  mixscale_terms <- unique(mixscale_data$Description)
  convergent_terms <- intersect(mast_terms, mixscale_terms)
  
  # Get top pathways for each category
  get_top_pathways <- function(data, terms, n_top = 10) {
    if (length(terms) == 0) return(data.frame())
    
    data %>%
      filter(Description %in% terms) %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        n_clusters = n_distinct(cluster),
        clusters = paste(unique(cluster), collapse = ", "),
        mean_neg_log_p = mean(-log10(p.adjust)),
        mean_fold = mean(FoldEnrichment, na.rm = TRUE),
        directions = paste(unique(direction), collapse = ";"),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_neg_log_p)) %>%
      head(n_top)
  }
  
  # MAST-only pathways
  mast_only_terms <- setdiff(mast_terms, mixscale_terms)
  mast_only_top <- get_top_pathways(mast_data, mast_only_terms)
  
  # MixScale-only pathways
  mixscale_only_terms <- setdiff(mixscale_terms, mast_terms)
  mixscale_only_top <- get_top_pathways(mixscale_data, mixscale_only_terms)
  
  # Convergent pathways
  convergent_top <- data.frame()
  if (length(convergent_terms) > 0) {
    convergent_top <- data %>%
      filter(gene == gene_name, Description %in% convergent_terms) %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        n_clusters_mast = n_distinct(cluster[method == "MAST"]),
        n_clusters_mixscale = n_distinct(cluster[method == "MixScale"]),
        mean_neg_log_p = mean(-log10(p.adjust)),
        methods = paste(unique(method), collapse = "+"),
        .groups = "drop"
      ) %>%
      arrange(desc(mean_neg_log_p)) %>%
      head(10)
  }
  
  # Summary statistics
  summary_stats <- data.frame(
    gene = gene_name,
    n_mast_only = length(mast_only_terms),
    n_mixscale_only = length(mixscale_only_terms),
    n_convergent = length(convergent_terms),
    n_total = length(unique(c(mast_terms, mixscale_terms))),
    top_mast_only = ifelse(nrow(mast_only_top) > 0, mast_only_top$Description[1], "None"),
    top_mixscale_only = ifelse(nrow(mixscale_only_top) > 0, mixscale_only_top$Description[1], "None"),
    top_convergent = ifelse(nrow(convergent_top) > 0, convergent_top$Description[1], "None")
  )
  
  return(list(
    gene = gene_name,
    mast_only = mast_only_top,
    mixscale_only = mixscale_only_top,
    convergent = convergent_top,
    summary = summary_stats
  ))
}

# Analyze all genes
gene_results <- list()
for (gene in all_genes) {
  gene_results[[gene]] <- analyze_gene_signatures(pd_data, gene)
}

# Combine summaries
all_summaries <- bind_rows(lapply(gene_results, function(x) x$summary))

# Save individual gene results
for (gene in names(gene_results)) {
  result <- gene_results[[gene]]
  if (nrow(result$mast_only) > 0 || nrow(result$mixscale_only) > 0 || nrow(result$convergent) > 0) {
    # Save CSVs
    write.csv(result$mast_only, 
              file.path(output_dir, paste0(gene, "_mast_only_pathways.csv")), 
              row.names = FALSE)
    write.csv(result$mixscale_only, 
              file.path(output_dir, paste0(gene, "_mixscale_only_pathways.csv")), 
              row.names = FALSE)
    write.csv(result$convergent, 
              file.path(output_dir, paste0(gene, "_convergent_pathways.csv")), 
              row.names = FALSE)
  }
}

# Save summary
write.csv(all_summaries, file.path(output_dir, "all_genes_summary.csv"), row.names = FALSE)

# Create visualization functions
create_gene_signature_plot <- function(gene_result, save_path = NULL) {
  gene_name <- gene_result$gene
  
  # Prepare data for plotting
  plot_data <- data.frame(
    Category = c("Mutation\niSCORE-PD\nOnly", "CRISPRi\nPerturbation\nOnly", "Convergent"),
    Count = c(
      gene_result$summary$n_mast_only,
      gene_result$summary$n_mixscale_only,
      gene_result$summary$n_convergent
    ),
    stringsAsFactors = FALSE
  )
  
  # Create bar plot
  p <- ggplot(plot_data, aes(x = Category, y = Count, fill = Category)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5, size = 5) +
    scale_fill_manual(values = c(
      "Mutation\niSCORE-PD\nOnly" = "#1f77b4",
      "CRISPRi\nPerturbation\nOnly" = "#ff7f0e",
      "Convergent" = "#2ca02c"
    )) +
    labs(
      title = paste("PD Pathway Signatures for", gene_name),
      subtitle = paste("Total pathways:", gene_result$summary$n_total),
      y = "Number of Pathways"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "none",
      axis.text.x = element_text(size = 12)
    ) +
    ylim(0, max(plot_data$Count) * 1.2)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 8, height = 6)
  }
  
  return(p)
}

# Create plots for all genes
plots_dir <- file.path(output_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE)

for (gene in names(gene_results)) {
  if (gene_results[[gene]]$summary$n_total > 0) {
    create_gene_signature_plot(
      gene_results[[gene]], 
      file.path(plots_dir, paste0(gene, "_signature_summary.pdf"))
    )
  }
}

# Create comparative analysis
cat("\n=== GENE COMPARISON SUMMARY ===\n")

# Find genes with strongest signatures in each category
top_mast_genes <- all_summaries %>%
  arrange(desc(n_mast_only)) %>%
  head(5) %>%
  select(gene, n_mast_only, top_mast_only)

top_mixscale_genes <- all_summaries %>%
  arrange(desc(n_mixscale_only)) %>%
  head(5) %>%
  select(gene, n_mixscale_only, top_mixscale_only)

top_convergent_genes <- all_summaries %>%
  arrange(desc(n_convergent)) %>%
  head(5) %>%
  select(gene, n_convergent, top_convergent)

cat("\nGenes with most Mutation - iSCORE-PD only pathways:\n")
print(top_mast_genes)

cat("\nGenes with most CRISPRi Perturbation only pathways:\n")
print(top_mixscale_genes)

cat("\nGenes with most convergent pathways:\n")
print(top_convergent_genes)

# Create summary heatmap data
heatmap_data <- all_summaries %>%
  select(gene, `Mutation Only` = n_mast_only, 
         `CRISPRi Only` = n_mixscale_only, 
         `Convergent` = n_convergent) %>%
  pivot_longer(cols = -gene, names_to = "Category", values_to = "Count")

# Summary plot
p_summary <- ggplot(heatmap_data, aes(x = gene, y = Category, fill = Count)) +
  geom_tile() +
  geom_text(aes(label = Count), color = "white", size = 4) +
  scale_fill_gradient(low = "lightblue", high = "darkred") +
  labs(
    title = "PD Pathway Distribution Across All Genes",
    x = "Gene",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold", size = 16)
  )

ggsave(file.path(output_dir, "all_genes_pathway_heatmap.pdf"), p_summary, width = 14, height = 6)

cat("\nAnalysis complete! Results saved to:", output_dir, "\n")