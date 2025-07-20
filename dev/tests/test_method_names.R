# Test with actual method names from the data
data_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds'
data <- readRDS(data_file)

# Test case 1: ATP13A2 (same name)
cat("=== Test 1: ATP13A2 vs ATP13A2 ===\n")
mast_atp <- data[data$mutation_perturbation == "ATP13A2" & 
                 data$cluster == "cluster_0" & 
                 data$method == "MAST", ]
mixscale_atp <- data[data$mutation_perturbation == "ATP13A2" & 
                     data$cluster == "cluster_0" & 
                     data$method == "MixScale", ]

cat("ATP13A2 MAST terms:", nrow(mast_atp), "\n")
cat("ATP13A2 MixScale terms:", nrow(mixscale_atp), "\n")

# Test case 2: PRKN vs PARK2 (different names, same gene)
cat("\n=== Test 2: PRKN vs PARK2 ===\n")
mast_prkn <- data[data$mutation_perturbation == "PRKN" & 
                  data$cluster == "cluster_0" & 
                  data$method == "MAST", ]
mixscale_park2 <- data[data$mutation_perturbation == "PARK2" & 
                       data$cluster == "cluster_0" & 
                       data$method == "MixScale", ]

cat("PRKN MAST terms:", nrow(mast_prkn), "\n")
cat("PARK2 MixScale terms:", nrow(mixscale_park2), "\n")

# Check if they have overlapping terms
if (nrow(mast_prkn) > 0 && nrow(mixscale_park2) > 0) {
  shared <- intersect(mast_prkn$Description, mixscale_park2$Description)
  cat("Shared terms between PRKN and PARK2:", length(shared), "\n")
  if (length(shared) > 0) {
    cat("First 3 shared terms:\n")
    print(head(shared, 3))
  }
}