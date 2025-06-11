# Example: Basic Usage of iSCORE-PDecipher
# This script demonstrates how to use the package for basic analysis

library(iSCORE.PDecipher)

# ============================================================================
# Example 1: Launch Shiny App with Data File
# ============================================================================

# If you have a pre-processed enrichment results file:
# launch_iscore_app("path/to/your/enrichment_results.rds")

# If you want to upload a file through the app interface:
# launch_iscore_app()

# ============================================================================
# Example 2: Import Individual Data Types
# ============================================================================

# Import MAST results (mutation analysis)
# mast_data <- import_mast_data("./iSCORE-PD_MAST_analysis/")
# cat("Imported", length(mast_data) - 1, "mutations\n")  # -1 for metadata

# Import CRISPRi results (knockdown perturbations)  
# crispi_data <- import_mixscale_data("./PerturbSeq_data/", modality = "CRISPRi")
# cat("Imported", length(crispi_data), "CRISPRi perturbations\n")

# Import CRISPRa results (activation perturbations)
# crispa_data <- import_mixscale_data("./PerturbSeq_data/", modality = "CRISPRa") 
# cat("Imported", length(crispa_data), "CRISPRa perturbations\n")

# ============================================================================
# Example 3: Run Enrichment Analysis on Existing DE Results  
# ============================================================================

# If you have a combined DE results file:
# enrichment_results <- run_enrichment_analysis(
#     input_file = "full_DE_results.rds",
#     output_dir = "./enrichment_results/",
#     lfc_threshold = 0.25,
#     padj_threshold = 0.05,
#     run_methods = c("GO", "KEGG", "Reactome")
# )

# ============================================================================
# Example 4: Process and Consolidate Enrichment Results
# ============================================================================

# After running enrichment analysis, consolidate results:
# consolidated_results <- process_enrichment_results(
#     root_dir = "enrichment_results",
#     output_file = "consolidated_enrichment.rds",
#     batch_size = 100
# )

# ============================================================================
# Example 5: Quick Analysis of Specific Genes
# ============================================================================

# Example: Focus on key PD genes
# pd_genes <- c("LRRK2", "PINK1", "PARK7", "SNCA")

# If you have the consolidated data:
# data <- readRDS("consolidated_enrichment.rds")
# 
# # Filter to PD genes of interest
# pd_subset <- data[data$mutation_perturbation %in% pd_genes, ]
# 
# # Get top KEGG pathways for these genes
# kegg_results <- pd_subset[pd_subset$enrichment_type == "KEGG" & 
#                           pd_subset$direction == "UP", ]
# 
# # View top results
# head(kegg_results[order(kegg_results$p.adjust), 
#                   c("mutation_perturbation", "cluster", "Description", "p.adjust")])

cat("Example scripts ready!")
cat("\nUncomment the relevant sections and modify paths to match your data.")
cat("\nSee the package documentation for more detailed examples.")