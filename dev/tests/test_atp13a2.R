# Test with ATP13A2 - same name in both methods

# Simulate what the extract_shared_genes function does
selected_pair <- "ATP13A2_vs_ATP13A2"
selected_cluster <- "cluster_0"

# Parse genes
genes <- strsplit(selected_pair, "_vs_")[[1]]
mast_gene <- genes[1]
crispri_gene <- genes[2]

cat("Testing with ATP13A2:\n")
cat("MAST gene:", mast_gene, "\n")
cat("CRISPRi gene:", crispri_gene, "\n")
cat("Cluster:", selected_cluster, "\n\n")

# Create test enrichment data
test_data <- data.frame(
  mutation_perturbation = c("ATP13A2", "ATP13A2", "ATP13A2", "ATP13A2", "PRKN", "PARK2"),
  cluster = c("cluster_0", "cluster_0", "cluster_1", "cluster_1", "cluster_0", "cluster_0"),
  method = c("iSCORE_PD_MAST", "CRISPRi_Mixscale", "iSCORE_PD_MAST", "CRISPRi_Mixscale", "iSCORE_PD_MAST", "CRISPRi_Mixscale"),
  Description = c("Term1", "Term1", "Term2", "Term3", "Term4", "Term4"),
  geneID = c("GENE1/GENE2", "GENE1/GENE3", "GENE4/GENE5", "GENE6/GENE7", "GENE8/GENE9", "GENE8/GENE10"),
  ID = c("GO:0001", "GO:0001", "GO:0002", "GO:0003", "GO:0004", "GO:0004"),
  enrichment_type = c("GO_BP", "GO_BP", "GO_BP", "GO_BP", "GO_BP", "GO_BP"),
  stringsAsFactors = FALSE
)

cat("Test enrichment data:\n")
print(test_data)
cat("\n")

# Test MAST filtering for ATP13A2
mast_terms <- test_data[
  test_data$mutation_perturbation == mast_gene & 
  test_data$cluster == selected_cluster &
  test_data$method == "iSCORE_PD_MAST", 
]

cat("MAST terms for ATP13A2 in cluster_0:", nrow(mast_terms), "\n")
if (nrow(mast_terms) > 0) print(mast_terms)
cat("\n")

# Test CRISPRi filtering for ATP13A2
crispri_terms <- test_data[
  test_data$mutation_perturbation == crispri_gene & 
  test_data$cluster == selected_cluster &
  test_data$method == "CRISPRi_Mixscale", 
]

cat("CRISPRi terms for ATP13A2 in cluster_0:", nrow(crispri_terms), "\n")
if (nrow(crispri_terms) > 0) print(crispri_terms)
cat("\n")

# Find shared terms
if (nrow(mast_terms) > 0 && nrow(crispri_terms) > 0) {
  shared_terms <- intersect(mast_terms$Description, crispri_terms$Description)
  cat("Shared terms:", length(shared_terms), "\n")
  cat("Shared descriptions:", paste(shared_terms, collapse=", "), "\n\n")
  
  # Extract genes from shared terms
  shared_mast_terms <- mast_terms[mast_terms$Description %in% shared_terms, ]
  shared_crispri_terms <- crispri_terms[crispri_terms$Description %in% shared_terms, ]
  
  # Get genes
  mast_genes <- unique(unlist(strsplit(as.character(shared_mast_terms$geneID), "/")))
  crispri_genes <- unique(unlist(strsplit(as.character(shared_crispri_terms$geneID), "/")))
  
  cat("Genes from MAST shared terms:", paste(mast_genes, collapse=", "), "\n")
  cat("Genes from CRISPRi shared terms:", paste(crispri_genes, collapse=", "), "\n")
  
  shared_genes <- intersect(mast_genes, crispri_genes)
  cat("Overlapping genes:", paste(shared_genes, collapse=", "), "\n")
}

cat("\n--- Now testing PRKN vs PARK2 ---\n")
# Now test PRKN vs PARK2
selected_pair <- "PRKN_vs_PARK2"
genes <- strsplit(selected_pair, "_vs_")[[1]]
mast_gene <- genes[1]
crispri_gene <- genes[2]

# PRKN in MAST
prkn_mast <- test_data[
  test_data$mutation_perturbation == "PRKN" & 
  test_data$cluster == selected_cluster &
  test_data$method == "iSCORE_PD_MAST", 
]
cat("PRKN MAST terms:", nrow(prkn_mast), "\n")

# PARK2 in CRISPRi
park2_crispri <- test_data[
  test_data$mutation_perturbation == "PARK2" & 
  test_data$cluster == selected_cluster &
  test_data$method == "CRISPRi_Mixscale", 
]
cat("PARK2 CRISPRi terms:", nrow(park2_crispri), "\n")