# Quick Plot Generation for PD Signatures
# Simplified version for immediate use

library(ggplot2)
library(dplyr)

# Load data
results_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures"
output_dir <- file.path(results_dir, "visualizations")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read CSV files
mast_top <- read.csv(file.path(results_dir, "mast_top_fast.csv"), stringsAsFactors = FALSE)
mixscale_top <- read.csv(file.path(results_dir, "mixscale_top_fast.csv"), stringsAsFactors = FALSE)
convergent_top <- read.csv(file.path(results_dir, "convergent_top_fast.csv"), stringsAsFactors = FALSE)

# Set theme
theme_set(theme_minimal(base_size = 14))

# 1. Top MAST-only pathways
p1 <- mast_top %>%
  head(15) %>%
  mutate(
    Description = as.character(Description),
    Description_short = ifelse(nchar(Description) > 50,
                              paste0(substr(Description, 1, 50), "..."),
                              Description)
  ) %>%
  ggplot(aes(x = reorder(Description_short, mean_neg_log_p), y = mean_neg_log_p)) +
  geom_bar(stat = "identity", fill = "#1f77b4") +
  coord_flip() +
  labs(
    title = "Top Mutation - iSCORE-PD Only Parkinson's Disease Pathways",
    subtitle = paste("Found in genetic mutation data (", nrow(mast_top), " total pathways)"),
    x = "",
    y = "Mean -log10(adjusted p-value)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "mast_only_top15.pdf"), p1, width = 12, height = 8)

# 2. Top MixScale-only pathways
p2 <- mixscale_top %>%
  head(15) %>%
  mutate(
    Description = as.character(Description),
    Description_short = ifelse(nchar(Description) > 50,
                              paste0(substr(Description, 1, 50), "..."),
                              Description)
  ) %>%
  ggplot(aes(x = reorder(Description_short, mean_neg_log_p), y = mean_neg_log_p)) +
  geom_bar(stat = "identity", fill = "#ff7f0e") +
  coord_flip() +
  labs(
    title = "Top CRISPRi Perturbation Only Parkinson's Disease Pathways",
    subtitle = paste("Found in knockdown data (", nrow(mixscale_top), " total pathways)"),
    x = "",
    y = "Mean -log10(adjusted p-value)"
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "mixscale_only_top15.pdf"), p2, width = 12, height = 8)

# 3. Convergent pathways - gene coverage
convergent_plot_data <- convergent_top %>%
  head(15) %>%
  mutate(
    Description = as.character(Description),
    Description_short = ifelse(nchar(Description) > 40,
                              paste0(substr(Description, 1, 40), "..."),
                              Description)
  )

# Create data for stacked bar
mast_bars <- data.frame(
  Description = convergent_plot_data$Description_short,
  Genes = convergent_plot_data$n_genes_mast,
  Method = "Mutation\niSCORE-PD",
  stringsAsFactors = FALSE
)

mixscale_bars <- data.frame(
  Description = convergent_plot_data$Description_short,
  Genes = convergent_plot_data$n_genes_mixscale,
  Method = "CRISPRi\nPerturbation",
  stringsAsFactors = FALSE
)

combined_bars <- rbind(mast_bars, mixscale_bars)

p3 <- ggplot(combined_bars, aes(x = reorder(Description, Genes), y = Genes, fill = Method)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_flip() +
  scale_fill_manual(values = c("Mutation\niSCORE-PD" = "#1f77b4", "CRISPRi\nPerturbation" = "#ff7f0e")) +
  labs(
    title = "Convergent Pathways: Found in Both Methods",
    subtitle = paste("Top pathways from", nrow(convergent_top), "total convergent signatures"),
    x = "",
    y = "Number of Genes"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom"
  )

ggsave(file.path(output_dir, "convergent_comparison.pdf"), p3, width = 12, height = 8)

# 4. Summary plot for presentation
top3_data <- data.frame(
  Category = c("Mutation - iSCORE-PD\nOnly", "CRISPRi Perturbation\nOnly", "Convergent"),
  Pathway = c(
    mast_top$Description[1],
    mixscale_top$Description[1],
    convergent_top$Description[1]
  ),
  Significance = c(
    mast_top$mean_neg_log_p[1],
    mixscale_top$mean_neg_log_p[1],
    convergent_top$mean_neg_log_p[1]
  ),
  Genes = c(
    mast_top$n_genes[1],
    mixscale_top$n_genes[1],
    convergent_top$n_genes_mast[1] + convergent_top$n_genes_mixscale[1]
  ),
  stringsAsFactors = FALSE
)

# Truncate pathway names
top3_data$Pathway_short <- sapply(top3_data$Pathway, function(x) {
  if (nchar(x) > 40) paste0(substr(x, 1, 40), "...") else x
})

p4 <- ggplot(top3_data, aes(x = Category, y = Significance, fill = Category)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Pathway_short), 
            angle = 90, hjust = 0, vjust = 0.5, 
            y = 0.5, size = 4) +
  scale_fill_manual(values = c(
    "Mutation - iSCORE-PD\nOnly" = "#1f77b4",
    "CRISPRi Perturbation\nOnly" = "#ff7f0e",
    "Convergent" = "#2ca02c"
  )) +
  labs(
    title = "Top Parkinson's Disease Signatures by Category",
    subtitle = "Most significant pathway from each analysis type",
    y = "Mean -log10(adjusted p-value)"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    legend.position = "none",
    axis.text.x = element_text(size = 12, face = "bold")
  )

ggsave(file.path(output_dir, "top3_summary.pdf"), p4, width = 10, height = 8)

cat("Plots saved to:", output_dir, "\n\n")
cat("Top signatures summary:\n")
cat("1. Mutation - iSCORE-PD only:", top3_data$Pathway[1], "\n")
cat("2. CRISPRi Perturbation only:", top3_data$Pathway[2], "\n")
cat("3. Convergent:", top3_data$Pathway[3], "\n")