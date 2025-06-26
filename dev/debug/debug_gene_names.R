# Debug Gene Names in Dataset
# Check what gene names are actually available for MAST vs CRISPRi comparison

data_file <- "E:/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"

cat("=== Debugging Gene Names ===\n")

if (file.exists(data_file)) {
  data <- readRDS(data_file)
  
  # Check MAST genes
  mast_data <- data[data$method == "MAST", ]
  mast_genes <- unique(mast_data$mutation_perturbation)
  cat("MAST genes (", length(mast_genes), "):", paste(mast_genes, collapse = ", "), "\n\n")
  
  # Check CRISPRi genes 
  crispri_data <- data[data$method == "MixScale", ]
  crispri_genes <- unique(crispri_data$mutation_perturbation)
  cat("CRISPRi genes (", length(crispri_genes), "):", paste(crispri_genes, collapse = ", "), "\n\n")
  
  # Test gene harmonization
  cat("=== Testing Gene Harmonization ===\n")
  
  # Source the gene harmonization functions
  source("R/gene_harmonization.R")
  
  # Get gene pairs
  gene_pairs <- get_comparable_gene_pairs(combine_snca_variants = TRUE, 
                                         combine_vps13c_variants = TRUE,
                                         include_mast_only = FALSE)
  
  cat("Gene pairs generated:\n")
  for (i in 1:nrow(gene_pairs)) {
    mast_gene <- gene_pairs$mast_gene[i]
    crispri_gene <- gene_pairs$crispri_gene[i]
    
    # Check if each gene exists in the data
    mast_exists <- mast_gene %in% mast_genes
    crispri_exists <- crispri_gene %in% crispri_genes
    
    cat(sprintf("%d. %s vs %s - MAST: %s, CRISPRi: %s\n", 
                i, mast_gene, crispri_gene, 
                ifelse(mast_exists, "FOUND", "MISSING"),
                ifelse(crispri_exists, "FOUND", "MISSING")))
    
    if (mast_exists && crispri_exists) {
      # Check cluster overlap
      mast_clusters <- unique(data[data$method == "MAST" & data$mutation_perturbation == mast_gene, ]$cluster)
      crispri_clusters <- unique(data[data$method == "MixScale" & data$mutation_perturbation == crispri_gene, ]$cluster)
      
      common_clusters <- intersect(mast_clusters, crispri_clusters)
      cat(sprintf("   Common clusters: %s\n", paste(common_clusters, collapse = ", ")))
    }
  }
  
} else {
  cat("Data file not found:", data_file, "\n")
}

cat("\n=== Debug Complete ===\n")