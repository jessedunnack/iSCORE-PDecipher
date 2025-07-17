# Comprehensive test to prove the gene extraction logic
# This simulates what should happen when there IS overlap

# First, let's check what the Fisher's test is actually testing
cat("=== Understanding the Fisher's Exact Test ===\n")
cat("The test looks at overlap between:\n")
cat("1. DE genes from MAST (mutation)\n")
cat("2. DE genes from CRISPRi (knockdown)\n")
cat("3. Tests if overlap is significant\n\n")

# Simulate the scenario for ATP13A2 with realistic data
selected_pair <- "ATP13A2_vs_ATP13A2"
selected_cluster <- "cluster_0"
genes <- strsplit(selected_pair, "_vs_")[[1]]
mast_gene <- genes[1]
crispri_gene <- genes[2]

cat("=== Test Case: ATP13A2 in cluster_0 ===\n")
cat("Gene pair:", selected_pair, "\n")
cat("Cluster:", selected_cluster, "\n\n")

# Create realistic enrichment data
# Each row is an enrichment term with its associated genes
enrichment_data <- data.frame(
  mutation_perturbation = c(
    # ATP13A2 MAST terms
    "ATP13A2", "ATP13A2", "ATP13A2",
    # ATP13A2 CRISPRi terms
    "ATP13A2", "ATP13A2", "ATP13A2",
    # Other genes for contrast
    "LRRK2", "PARK2"
  ),
  cluster = c(
    "cluster_0", "cluster_0", "cluster_1",
    "cluster_0", "cluster_0", "cluster_1", 
    "cluster_0", "cluster_0"
  ),
  method = c(
    "iSCORE_PD_MAST", "iSCORE_PD_MAST", "iSCORE_PD_MAST",
    "CRISPRi_Mixscale", "CRISPRi_Mixscale", "CRISPRi_Mixscale",
    "iSCORE_PD_MAST", "CRISPRi_Mixscale"
  ),
  Description = c(
    "mitochondrial ATP synthesis", "oxidative phosphorylation", "protein folding",
    "mitochondrial ATP synthesis", "lysosomal transport", "autophagy",
    "dopamine metabolism", "protein degradation"
  ),
  ID = c(
    "GO:0042775", "GO:0006119", "GO:0006457",
    "GO:0042775", "GO:0007041", "GO:0006914",
    "GO:0042417", "GO:0030163"
  ),
  # These are the DE genes driving each enrichment
  geneID = c(
    "ATP5A1/ATP5B/COX4I1/NDUFS1",  # mitochondrial ATP synthesis (MAST)
    "COX4I1/NDUFS1/UQCRC1/SDHA",   # oxidative phosphorylation (MAST)
    "HSPA8/HSPA5/HSP90AA1",         # protein folding (MAST)
    "ATP5A1/ATP5B/COX4I1/ATP6V1A", # mitochondrial ATP synthesis (CRISPRi) - overlaps!
    "ATP6V1A/CTSD/LAMP1/GBA",      # lysosomal transport (CRISPRi)
    "BECN1/ATG5/MAP1LC3B",          # autophagy (CRISPRi)
    "TH/DDC/COMT",                  # dopamine metabolism
    "PSMA1/PSMB5/UBE2D3"           # protein degradation
  ),
  p.adjust = c(0.001, 0.002, 0.01, 0.001, 0.005, 0.02, 0.001, 0.003),
  enrichment_type = "GO_BP",
  stringsAsFactors = FALSE
)

cat("Simulated enrichment data (", nrow(enrichment_data), "rows ):\n")
print(enrichment_data[1:6, c("mutation_perturbation", "cluster", "method", "Description", "geneID")])
cat("\n")

# Run the extraction logic
cat("=== Running Extraction Logic ===\n")

# Get MAST terms
mast_terms <- enrichment_data[
  enrichment_data$mutation_perturbation == mast_gene & 
  enrichment_data$cluster == selected_cluster &
  enrichment_data$method == "iSCORE_PD_MAST", 
]
cat("MAST enrichment terms found:", nrow(mast_terms), "\n")
print(mast_terms[, c("Description", "geneID")])
cat("\n")

# Get CRISPRi terms
crispri_terms <- enrichment_data[
  enrichment_data$mutation_perturbation == crispri_gene & 
  enrichment_data$cluster == selected_cluster &
  enrichment_data$method == "CRISPRi_Mixscale", 
]
cat("CRISPRi enrichment terms found:", nrow(crispri_terms), "\n")
print(crispri_terms[, c("Description", "geneID")])
cat("\n")

# Find shared enrichment terms
shared_descriptions <- intersect(mast_terms$Description, crispri_terms$Description)
cat("Shared enrichment terms:", length(shared_descriptions), "\n")
cat("- ", paste(shared_descriptions, collapse="\n- "), "\n\n")

# Extract genes from shared terms
if (length(shared_descriptions) > 0) {
  cat("=== Extracting Genes from Shared Terms ===\n")
  
  # Get the shared term rows
  shared_mast_terms <- mast_terms[mast_terms$Description %in% shared_descriptions, ]
  shared_crispri_terms <- crispri_terms[crispri_terms$Description %in% shared_descriptions, ]
  
  cat("MAST genes from 'mitochondrial ATP synthesis':\n")
  mast_genes_atp <- unlist(strsplit(shared_mast_terms$geneID[1], "/"))
  cat("  ", paste(mast_genes_atp, collapse=", "), "\n\n")
  
  cat("CRISPRi genes from 'mitochondrial ATP synthesis':\n")
  crispri_genes_atp <- unlist(strsplit(shared_crispri_terms$geneID[1], "/"))
  cat("  ", paste(crispri_genes_atp, collapse=", "), "\n\n")
  
  # Get ALL genes from shared terms
  all_mast_genes <- unique(unlist(strsplit(shared_mast_terms$geneID, "/")))
  all_crispri_genes <- unique(unlist(strsplit(shared_crispri_terms$geneID, "/")))
  
  cat("All unique MAST genes from shared terms:", length(all_mast_genes), "\n")
  cat("All unique CRISPRi genes from shared terms:", length(all_crispri_genes), "\n")
  
  # Find overlapping genes
  shared_genes <- intersect(all_mast_genes, all_crispri_genes)
  cat("\nOverlapping DE genes:", length(shared_genes), "\n")
  cat("  ", paste(sort(shared_genes), collapse=", "), "\n\n")
  
  cat("=== What This Means ===\n")
  cat("1. Both methods found 'mitochondrial ATP synthesis' enriched\n")
  cat("2. MAST found it enriched due to genes: ATP5A1, ATP5B, COX4I1, NDUFS1\n")
  cat("3. CRISPRi found it enriched due to genes: ATP5A1, ATP5B, COX4I1, ATP6V1A\n")
  cat("4. The overlapping genes (ATP5A1, ATP5B, COX4I1) drive this shared signature\n")
  cat("5. These are the genes that would appear in the 'Shared DE Genes' tab\n")
}