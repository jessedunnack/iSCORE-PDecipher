# Thesis Committee Summary Visualization
# Purpose: Create a single, comprehensive figure for thesis committee presentation
# Date: 2025-07-20

library(ggplot2)
library(dplyr)
library(patchwork)
library(gridExtra)

# Helper function for word wrapping
wrap_text <- function(text, width = 40) {
  sapply(text, function(x) {
    if (is.na(x)) return('')
    stringr::str_wrap(x, width = width)
  })
}

# Improved function for bar plot labels - handles STRING database results
truncate_for_bars <- function(text, source_col = NULL) {
  sapply(seq_along(text), function(i) {
    x <- text[i]
    if (is.na(x)) return('')
    
    # Check if it's a STRING result (either by pattern or source column)
    is_string <- FALSE
    if (!is.null(source_col) && length(source_col) >= i) {
      is_string <- source_col[i] == "STRING"
    } else {
      is_string <- grepl("^\\(\\d{4}\\)", x)
    }
    
    if (is_string) {
      # For STRING: Just show year and first key term
      if (grepl("^\\(\\d{4}\\)", x)) {
        year <- gsub("^(\\(\\d{4}\\)).*", "\\1", x)
        rest <- gsub("^\\(\\d{4}\\)\\s*", "", x)
        
        # Look for key terms
        key_terms <- c("mitochondr", "synap", "ribosom", "oxidative", 
                      "ubiquitin", "dopamin", "transport", "vesicle")
        
        for (term in key_terms) {
          if (grepl(term, rest, ignore.case = TRUE)) {
            # Find the word containing this term
            words <- strsplit(rest, " ")[[1]]
            matching_word <- words[grep(term, words, ignore.case = TRUE)[1]]
            if (!is.na(matching_word)) {
              return(paste(year, matching_word))
            }
          }
        }
        
        # If no key term found, just show year + first 2 words
        words <- strsplit(rest, " ")[[1]]
        return(paste(c(year, words[1:min(2, length(words))], "..."), collapse = " "))
      }
    }
    
    # For regular pathways, limit to 45 characters
    if (nchar(x) > 45) {
      return(paste0(substr(x, 1, 42), "..."))
    }
    
    return(x)
  })
}


# Configuration
results_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
output_dir <- file.path(results_dir, "visualizations/comprehensive")

# Load data
mast_top <- read.csv(file.path(results_dir, "mast_top_fast.csv"), stringsAsFactors = FALSE)
mixscale_top <- read.csv(file.path(results_dir, "mixscale_top_fast.csv"), stringsAsFactors = FALSE)
convergent_top <- read.csv(file.path(results_dir, "convergent_top_fast.csv"), stringsAsFactors = FALSE)
gene_summary <- read.csv(file.path(results_dir, "by_gene/all_genes_summary.csv"), stringsAsFactors = FALSE)

# Load actual pathway totals
pathway_totals <- read.csv(file.path(results_dir, "pathway_totals.csv"), stringsAsFactors = FALSE)

# Set theme
theme_set(theme_minimal(base_size = 12))

# Create 4-panel figure for thesis committee
cat("Creating thesis committee summary figure...\n")

# Panel A: Overview
overview_data <- data.frame(
  Category = c("Mutation-\nOnly", "CRISPRi-\nOnly", "Convergent"),
  Count = c(pathway_totals$mast_only_total, pathway_totals$mixscale_only_total, pathway_totals$convergent_total),
  y_pos = c(1.1, 1.1, 1.1)
)

p_overview <- ggplot(overview_data, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = Count), vjust = -0.5, size = 5, fontface = "bold") +
  scale_fill_manual(values = c("Mutation-\nOnly" = "#1f77b4", 
                              "CRISPRi-\nOnly" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  labs(title = "A. PD Pathway Categories",
       subtitle = "Method-specific discoveries",
       y = "Number of Pathways") +
  theme(legend.position = "none",
        plot.title = element_text(face = "bold"),
        axis.text.x = element_text(size = 11)) +
  ylim(0, max(overview_data$Count) * 1.3)

# Panel B: Top Convergent Pathways
top5_conv <- convergent_top %>%
  head(5) %>%
  mutate(
    Description_short = case_when(
      grepl("synapse", Description, ignore.case = TRUE) ~ "Synapse",
      grepl("vesicle", Description, ignore.case = TRUE) ~ "Synaptic vesicle",
      grepl("dopamin", Description, ignore.case = TRUE) ~ "Dopaminergic",
      grepl("mitochond", Description, ignore.case = TRUE) ~ "Mitochondrial",
      TRUE ~ wrap_text(Description, width = 20)
    ),
    total_genes = n_genes_mast + n_genes_mixscale
  )

p_convergent <- ggplot(top5_conv, aes(x = reorder(Description_short, mean_neg_log_p), 
                                      y = mean_neg_log_p)) +
  geom_segment(aes(xend = Description_short, yend = 0), size = 2, color = "gray50") +
  geom_point(aes(size = total_genes), color = "#2ca02c", alpha = 0.8) +
  scale_size_continuous(range = c(6, 12), guide = "none") +
  coord_flip() +
  labs(title = "B. Top Convergent Pathways",
       subtitle = "Validated across both methods",
       x = "",
       y = "-log10(p-value)") +
  theme(plot.title = element_text(face = "bold")) +
  geom_text(aes(label = round(mean_neg_log_p, 1)), 
            hjust = -0.5, size = 4)

# Panel C: Gene Patterns
gene_patterns <- gene_summary %>%
  mutate(
    total = n_mast_only + n_mixscale_only + n_convergent,
    convergent_ratio = n_convergent / total
  ) %>%
  arrange(desc(convergent_ratio)) %>%
  head(8)

gene_pattern_long <- gene_patterns %>%
  select(gene, `Mutation-Only` = n_mast_only, 
         `CRISPRi-Only` = n_mixscale_only, 
         Convergent = n_convergent) %>%
  tidyr::pivot_longer(cols = -gene, names_to = "Type", values_to = "Count")

p_genes <- ggplot(gene_pattern_long, aes(x = gene, y = Count, fill = Type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = c("Mutation-Only" = "#1f77b4",
                              "CRISPRi-Only" = "#ff7f0e",
                              "Convergent" = "#2ca02c")) +
  coord_flip() +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "C. Gene Signature Patterns",
       subtitle = "Top genes by convergence",
       x = "",
       y = "Proportion of Pathways",
       fill = "") +
  theme(plot.title = element_text(face = "bold"),
        legend.position = "bottom",
        legend.box.spacing = unit(0, "pt"))

# Panel D: Key Statistics
key_stats <- data.frame(
  Statistic = c("Genes Analyzed", "Clusters", "Strongest p-value", "Top Gene"),
  Value = c(
    as.character(length(unique(gene_summary$gene))),
    "15",
    paste0("< 10^-", round(max(convergent_top$mean_neg_log_p), 0)),
    gene_summary$gene[which.max(gene_summary$n_convergent)]
  )
)

p_stats <- ggplot(data = NULL) +
  theme_void() +
  annotate("text", x = 0.1, y = 0.9, label = "D. Key Findings", 
           size = 6, fontface = "bold", hjust = 0) +
  annotate("text", x = 0.1, y = 0.7, 
           label = paste0("• ", key_stats$Statistic[1], ": ", key_stats$Value[1]),
           size = 5, hjust = 0) +
  annotate("text", x = 0.1, y = 0.55, 
           label = paste0("• ", key_stats$Statistic[2], ": ", key_stats$Value[2]),
           size = 5, hjust = 0) +
  annotate("text", x = 0.1, y = 0.4, 
           label = paste0("• ", key_stats$Statistic[3], ": ", key_stats$Value[3]),
           size = 5, hjust = 0, color = "#2ca02c") +
  annotate("text", x = 0.1, y = 0.25, 
           label = paste0("• ", key_stats$Statistic[4], ": ", key_stats$Value[4]),
           size = 5, hjust = 0) +
  annotate("text", x = 0.1, y = 0.05, 
           label = "Methods: Mutation (iSCORE-PD) + CRISPRi",
           size = 4, hjust = 0, fontface = "italic") +
  xlim(0, 1) + ylim(0, 1)

# Combine all panels
thesis_fig <- (p_overview | p_convergent) / (p_genes | p_stats) + 
  plot_annotation(
    title = "iSCORE-PDecipher: Convergent Parkinson's Disease Signatures",
    subtitle = "Orthogonal validation through genetic mutations and CRISPRi perturbations",
    theme = theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 14, hjust = 0.5)
    )
  )

# Save as both PDF and PNG
ggsave(file.path(output_dir, "THESIS_COMMITTEE_SUMMARY.pdf"), 
       thesis_fig, width = 16, height = 12)
ggsave(file.path(output_dir, "THESIS_COMMITTEE_SUMMARY.png"), 
       thesis_fig, width = 16, height = 12, dpi = 300)

cat("\nThesis committee summary figure created!\n")
cat("Files saved:\n")
cat("- THESIS_COMMITTEE_SUMMARY.pdf\n")
cat("- THESIS_COMMITTEE_SUMMARY.png\n\n")

# Also create a simple version focusing on the key message
simple_fig <- convergent_top %>%
  head(10) %>%
  mutate(
    Description_clean = case_when(
      grepl("synapse", Description, ignore.case = TRUE) ~ "Synaptic Function",
      grepl("vesicle", Description, ignore.case = TRUE) ~ "Vesicle Transport",
      grepl("dopamin", Description, ignore.case = TRUE) ~ "Dopaminergic Signaling",
      grepl("mitochond", Description, ignore.case = TRUE) ~ "Mitochondrial Function",
      grepl("transport", Description, ignore.case = TRUE) ~ "Cellular Transport",
      grepl("metabol", Description, ignore.case = TRUE) ~ "Metabolism",
      TRUE ~ wrap_text(Description, width = 30)
    )
  ) %>%
  ggplot(aes(x = reorder(Description_clean, mean_neg_log_p), y = mean_neg_log_p)) +
  geom_bar(stat = "identity", fill = "#2ca02c", alpha = 0.8) +
  geom_text(aes(label = paste0("p < 10^-", round(mean_neg_log_p, 0))), 
            hjust = -0.1, size = 4) +
  coord_flip() +
  labs(
    title = "Convergent Parkinson's Disease Pathways",
    subtitle = "Validated through both genetic mutations and CRISPRi knockdowns",
    x = "",
    y = "-log10(adjusted p-value)"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(face = "bold", size = 22),
    plot.subtitle = element_text(size = 16),
    axis.text.y = element_text(size = 14)
  ) +
  xlim(NA, max(convergent_top$mean_neg_log_p[1:10]) * 1.3)

ggsave(file.path(output_dir, "KEY_CONVERGENT_PATHWAYS.pdf"), 
       simple_fig, width = 12, height = 10)

cat("Additional figure created: KEY_CONVERGENT_PATHWAYS.pdf\n")
