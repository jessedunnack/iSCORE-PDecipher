# Example: Complete Workflow from Raw Data to Visualization
# This script demonstrates the full iSCORE-PDecipher pipeline

library(iSCORE.PDecipher)

# ============================================================================
# Complete Workflow: Raw Data → Analysis → Visualization
# ============================================================================

# OPTION 1: One-step complete pipeline
# ===================================
# This runs everything from import through enrichment to launching the app:

# results <- run_complete_pipeline(
#     mast_directory = "./iSCORE-PD_MAST_analysis/",
#     mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
#     output_directory = "./analysis_output/",
#     lfc_threshold = 0.25,
#     padj_threshold = 0.05,
#     run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
#     launch_app = TRUE  # Automatically launches Shiny app when done
# )

# OPTION 2: Step-by-step workflow 
# ===============================
# For more control over each step:

cat("=== STEP-BY-STEP WORKFLOW ===\n")

# Step 1: Import differential expression results
cat("Step 1: Importing DE results...\n")
# mast_data <- import_mast_data("./iSCORE-PD_MAST_analysis/")
# crispi_data <- import_mixscale_data("./PerturbSeq_data/", modality = "CRISPRi")
# crispa_data <- import_mixscale_data("./PerturbSeq_data/", modality = "CRISPRa")

# Step 2: Combine into single dataset
cat("Step 2: Combining datasets...\n")
# full_results <- list(
#     iSCORE_PD_MAST = mast_data,
#     CRISPRi_Mixscale = crispi_data,
#     CRISPRa_Mixscale = crispa_data
# )
# saveRDS(full_results, "full_DE_results.rds")

# Step 3: Run functional enrichment analysis
cat("Step 3: Running enrichment analysis...\n")
# enrichment_results <- run_enrichment_analysis(
#     input_file = "full_DE_results.rds",
#     output_dir = "./enrichment_results/",
#     lfc_threshold = 0.25,
#     padj_threshold = 0.05,
#     run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA")
# )

# Step 4: Process and consolidate results
cat("Step 4: Consolidating results...\n")
# consolidated_data <- process_enrichment_results(
#     root_dir = "enrichment_results",
#     output_file = "all_enrichment_complete.rds",
#     batch_size = 100
# )

# Step 5: Launch Shiny app for visualization
cat("Step 5: Launching visualization app...\n")
# launch_iscore_app("all_enrichment_complete.rds")

# ============================================================================
# Working with Results
# ============================================================================

# Once you have consolidated results, you can explore them:

# Load consolidated data
# data <- readRDS("all_enrichment_complete.rds")
# 
# # Basic exploration
# cat("\nDataset Overview:\n")
# cat("Total terms:", nrow(data), "\n")
# cat("Unique genes:", length(unique(data$mutation_perturbation)), "\n")
# cat("Enrichment types:", toString(unique(data$enrichment_type)), "\n")
# cat("Methods:", toString(unique(data$method)), "\n")
# 
# # Example analyses:
# 
# # 1. Top KEGG pathways for LRRK2 mutations
# lrrk2_kegg <- data[data$mutation_perturbation == "LRRK2" & 
#                    data$enrichment_type == "KEGG" & 
#                    data$direction == "UP", ]
# cat("\nTop LRRK2 KEGG pathways:\n")
# print(head(lrrk2_kegg[order(lrrk2_kegg$p.adjust), 
#                       c("cluster", "Description", "p.adjust")]))
# 
# # 2. Compare enrichment across clusters for a specific gene
# pink1_go <- data[data$mutation_perturbation == "PINK1" & 
#                  data$enrichment_type == "GO_BP" & 
#                  data$direction == "DOWN", ]
# 
# # 3. Find commonly enriched pathways across multiple genes
# pd_genes <- c("LRRK2", "PINK1", "PARK7", "SNCA")
# common_pathways <- data[data$mutation_perturbation %in% pd_genes & 
#                         data$enrichment_type == "GO_BP", ]

# ============================================================================
# Customizing Analysis Parameters
# ============================================================================

# You can customize various aspects of the analysis:

# # Custom thresholds
# results_strict <- run_enrichment_analysis(
#     input_file = "full_DE_results.rds",
#     output_dir = "./enrichment_strict/",
#     lfc_threshold = 0.5,      # Stricter fold change
#     padj_threshold = 0.01,    # Stricter p-value
#     run_methods = c("GO", "KEGG")  # Subset of methods
# )
# 
# # Focus on specific enrichment types
# results_pathways <- run_enrichment_analysis(
#     input_file = "full_DE_results.rds",
#     output_dir = "./enrichment_pathways/",
#     run_methods = c("KEGG", "Reactome", "WikiPathways")
# )

cat("\nWorkflow examples ready!")
cat("\nUncomment sections and modify paths to run with your data.")
cat("\nSee package documentation for detailed parameter descriptions.")