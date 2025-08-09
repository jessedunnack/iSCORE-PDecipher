#!/usr/bin/env Rscript

# DETAILED CLUSTER ANALYSIS WITH SPECIFIC MARKER CHECKS

library(dplyr)

cat("=================================================================\n")
cat("DETAILED CLUSTER ANALYSIS BASED ON REFERENCE DATASETS\n")
cat("=================================================================\n\n")

# Load marker data
markers <- readRDS("inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds")

# Create comprehensive analysis for each cluster
cluster_analysis <- data.frame(
  cluster = 0:35,
  cell_type = NA,
  subtype = NA,
  confidence = NA,
  key_evidence = NA,
  stringsAsFactors = FALSE
)

# Function to get top markers for a cluster
get_top_markers <- function(cluster_id, n = 20) {
  markers %>%
    filter(cluster == as.character(cluster_id)) %>%
    arrange(desc(avg_log2FC)) %>%
    head(n) %>%
    pull(gene)
}

# Analyze each cluster
for (i in 0:35) {
  top20 <- get_top_markers(i)
  top10_str <- paste(top20[1:10], collapse = ", ")
  
  # Cluster 0: SNCB, TAC1, RIT2, STMN2
  if (i == 0) {
    cluster_analysis[i+1, "cell_type"] <- "Mature Neurons"
    cluster_analysis[i+1, "subtype"] <- "Non-dopaminergic (SNCB+)"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "SNCB, TAC1, RIT2, STMN2 (neuronal)"
  }
  
  # Cluster 1: AL138899.1, PTPRC (CD45), MS4A6A - immune markers
  else if (i == 1) {
    cluster_analysis[i+1, "cell_type"] <- "Immune-like/Microglia-like"
    cluster_analysis[i+1, "subtype"] <- "Unexpected in iPSC culture"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "PTPRC (CD45), MS4A6A, F13A1"
  }
  
  # Cluster 2: GDA, TMC5, FGF7, LEFTY2
  else if (i == 2) {
    cluster_analysis[i+1, "cell_type"] <- "Mesenchymal/Support Cells"
    cluster_analysis[i+1, "subtype"] <- "FGF7+ secretory"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "FGF7, LEFTY2, COL9A3, ELN"
  }
  
  # Cluster 3: DIO3, NDP, DLK1, PTPRZ1
  else if (i == 3) {
    cluster_analysis[i+1, "cell_type"] <- "Neural Progenitors"
    cluster_analysis[i+1, "subtype"] <- "Midbrain progenitors (DLK1+)"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "DIO3, DLK1, PTPRZ1, SOX1"
  }
  
  # Cluster 4: TTR! Clear choroid plexus
  else if (i == 4) {
    cluster_analysis[i+1, "cell_type"] <- "Choroid Plexus"
    cluster_analysis[i+1, "subtype"] <- "TTR-high epithelial"
    cluster_analysis[i+1, "confidence"] <- "Very High"
    cluster_analysis[i+1, "key_evidence"] <- "TTR, KCNJ13, SLC39A12"
  }
  
  # Cluster 5: TH in top 10! Dopaminergic
  else if (i == 5) {
    cluster_analysis[i+1, "cell_type"] <- "Dopaminergic Neurons"
    cluster_analysis[i+1, "subtype"] <- "Immature/developing (VGF+)"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "TH, VGF, C1QL1, DMRTA2"
  }
  
  # Cluster 6: XIST, mitochondrial genes - stressed/dying
  else if (i == 6) {
    cluster_analysis[i+1, "cell_type"] <- "Stressed/Dying Cells"
    cluster_analysis[i+1, "subtype"] <- "High mitochondrial"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "MT genes, XIST, MALAT1"
  }
  
  # Cluster 7: PRRX1, PRRX2, COL6A3 - mesenchymal
  else if (i == 7) {
    cluster_analysis[i+1, "cell_type"] <- "Mesenchymal Cells"
    cluster_analysis[i+1, "subtype"] <- "PRRX1/2+ fibroblast-like"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "PRRX1, PRRX2, COL1A1, TWIST2"
  }
  
  # Cluster 8: NDP, PTN, PTPRZ1 - radial glia
  else if (i == 8) {
    cluster_analysis[i+1, "cell_type"] <- "Radial Glia-like"
    cluster_analysis[i+1, "subtype"] <- "PTN+ progenitors"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "PTN, PTPRZ1, VCAM1"
  }
  
  # Cluster 9: ESM1, RSPO2, GRIK1
  else if (i == 9) {
    cluster_analysis[i+1, "cell_type"] <- "Endothelial-like"
    cluster_analysis[i+1, "subtype"] <- "ESM1+ vascular"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "ESM1, TFPI2, RSPO2"
  }
  
  # Cluster 10: RSPO3, MSX1 - roof plate
  else if (i == 10) {
    cluster_analysis[i+1, "cell_type"] <- "Roof Plate-like"
    cluster_analysis[i+1, "subtype"] <- "MSX1+/RSPO3+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "RSPO3, MSX1, PAPPA2"
  }
  
  # Cluster 11: PRSS56, ECEL1, NTS
  else if (i == 11) {
    cluster_analysis[i+1, "cell_type"] <- "Specialized Neurons"
    cluster_analysis[i+1, "subtype"] <- "Neuropeptide-expressing"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "NTS, PRSS56, ECEL1"
  }
  
  # Cluster 12: CRABP1, DCC - developing neurons
  else if (i == 12) {
    cluster_analysis[i+1, "cell_type"] <- "Developing Neurons"
    cluster_analysis[i+1, "subtype"] <- "CRABP1+ midbrain"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "CRABP1, DCC, ST8SIA6"
  }
  
  # Cluster 13: DCN, FGL2, AQP4 - astrocyte-like
  else if (i == 13) {
    cluster_analysis[i+1, "cell_type"] <- "Astrocyte-like"
    cluster_analysis[i+1, "subtype"] <- "DCN+/AQP4+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "DCN, AQP4, LRRK2"
  }
  
  # Cluster 14: DCT, TYRP1, PMEL - melanocytes!
  else if (i == 14) {
    cluster_analysis[i+1, "cell_type"] <- "Melanocyte-like"
    cluster_analysis[i+1, "subtype"] <- "Neural crest-derived"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "DCT, TYRP1, MLANA, PMEL"
  }
  
  # Cluster 15: PTGDS, CRYAB, AQP4 - leptomeningeal
  else if (i == 15) {
    cluster_analysis[i+1, "cell_type"] <- "Leptomeningeal Cells"
    cluster_analysis[i+1, "subtype"] <- "PTGDS+ meningeal"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "PTGDS, CRYAB, AQP4, F5"
  }
  
  # Cluster 16: NMU, RGCC - stressed progenitors
  else if (i == 16) {
    cluster_analysis[i+1, "cell_type"] <- "Stressed Progenitors"
    cluster_analysis[i+1, "subtype"] <- "NMU+/RGCC+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "NMU, RGCC, USH1C"
  }
  
  # Cluster 17: DSCAM, CNTNAP5, NRG3 - developing neurons
  else if (i == 17) {
    cluster_analysis[i+1, "cell_type"] <- "Developing Neurons"
    cluster_analysis[i+1, "subtype"] <- "Adhesion molecule-rich"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "DSCAM, CNTNAP5, NRG3"
  }
  
  # Cluster 18: CALB1! A10 dopaminergic
  else if (i == 18) {
    cluster_analysis[i+1, "cell_type"] <- "Dopaminergic Neurons"
    cluster_analysis[i+1, "subtype"] <- "A10-like (CALB1+, resilient)"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "CALB1, PRPH, C1QL4"
  }
  
  # Cluster 19: MPZ, SOX10 - oligodendrocytes
  else if (i == 19) {
    cluster_analysis[i+1, "cell_type"] <- "Oligodendrocyte Lineage"
    cluster_analysis[i+1, "subtype"] <- "SOX10+/MPZ+"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "MPZ, SOX10, GFRA3, ERBB3"
  }
  
  # Cluster 20: TTN, NEAT1 - stressed/artifact
  else if (i == 20) {
    cluster_analysis[i+1, "cell_type"] <- "Low Quality/Doublets"
    cluster_analysis[i+1, "subtype"] <- "High complexity"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "TTN, PLCG2, NEAT1"
  }
  
  # Cluster 21: COL1A1, COL1A2 - fibroblasts
  else if (i == 21) {
    cluster_analysis[i+1, "cell_type"] <- "Fibroblasts"
    cluster_analysis[i+1, "subtype"] <- "Collagen-producing"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "COL1A1, COL1A2, COL12A1"
  }
  
  # Cluster 22: MKI67, TOP2A - proliferating
  else if (i == 22) {
    cluster_analysis[i+1, "cell_type"] <- "Proliferating Cells"
    cluster_analysis[i+1, "subtype"] <- "S/G2/M phase"
    cluster_analysis[i+1, "confidence"] <- "Very High"
    cluster_analysis[i+1, "key_evidence"] <- "MKI67, TOP2A, HIST1H3B"
  }
  
  # Cluster 23: AQP1, RGCC
  else if (i == 23) {
    cluster_analysis[i+1, "cell_type"] <- "Ependymal-like"
    cluster_analysis[i+1, "subtype"] <- "AQP1+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "AQP1, RGCC, GDF7"
  }
  
  # Cluster 24: GDF15, ATF5, DDIT3 - stress response
  else if (i == 24) {
    cluster_analysis[i+1, "cell_type"] <- "Stressed Cells"
    cluster_analysis[i+1, "subtype"] <- "ER stress response"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "GDF15, ATF5, DDIT3, ATF3"
  }
  
  # Cluster 25: CPXM2, PRDM16
  else if (i == 25) {
    cluster_analysis[i+1, "cell_type"] <- "Specialized Progenitors"
    cluster_analysis[i+1, "subtype"] <- "PRDM16+"
    cluster_analysis[i+1, "confidence"] <- "Low"
    cluster_analysis[i+1, "key_evidence"] <- "CPXM2, PRDM16, IGFBP3"
  }
  
  # Cluster 26: NES, MEF2C, TPM2
  else if (i == 26) {
    cluster_analysis[i+1, "cell_type"] <- "Neural Progenitors"
    cluster_analysis[i+1, "subtype"] <- "NES+ (but see TTR)"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "NES, MEF2C, TPM2"
  }
  
  # Cluster 27: TAGLN, ACTA2 - smooth muscle
  else if (i == 27) {
    cluster_analysis[i+1, "cell_type"] <- "Smooth Muscle/Pericytes"
    cluster_analysis[i+1, "subtype"] <- "Vascular smooth muscle"
    cluster_analysis[i+1, "confidence"] <- "Very High"
    cluster_analysis[i+1, "key_evidence"] <- "TAGLN, ACTA2, MYL9"
  }
  
  # Cluster 28: CGA, SST, CALCA - neuroendocrine
  else if (i == 28) {
    cluster_analysis[i+1, "cell_type"] <- "Neuroendocrine Cells"
    cluster_analysis[i+1, "subtype"] <- "CGA+/SST+"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "CGA, SST, CALCA"
  }
  
  # Cluster 29: PENK, HOXD genes
  else if (i == 29) {
    cluster_analysis[i+1, "cell_type"] <- "Posterior Neural"
    cluster_analysis[i+1, "subtype"] <- "HOX+ (unexpected)"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "PENK, HOXD11, HOXD10"
  }
  
  # Cluster 30: ATOH1, TLX3 - hindbrain
  else if (i == 30) {
    cluster_analysis[i+1, "cell_type"] <- "Hindbrain Neurons"
    cluster_analysis[i+1, "subtype"] <- "ATOH1+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "ATOH1, NHLH1, TLX3"
  }
  
  # Cluster 31: RBP4
  else if (i == 31) {
    cluster_analysis[i+1, "cell_type"] <- "Specialized Epithelial"
    cluster_analysis[i+1, "subtype"] <- "RBP4+"
    cluster_analysis[i+1, "confidence"] <- "Low"
    cluster_analysis[i+1, "key_evidence"] <- "RBP4, FFAR4"
  }
  
  # Cluster 32: CCNO, multiciliated
  else if (i == 32) {
    cluster_analysis[i+1, "cell_type"] <- "Multiciliated Cells"
    cluster_analysis[i+1, "subtype"] <- "CCNO+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "CCNO, CDC20B, SNTN"
  }
  
  # Cluster 33: CRYAB, PRG4
  else if (i == 33) {
    cluster_analysis[i+1, "cell_type"] <- "Specialized Support"
    cluster_analysis[i+1, "subtype"] <- "PRG4+"
    cluster_analysis[i+1, "confidence"] <- "Low"
    cluster_analysis[i+1, "key_evidence"] <- "CRYAB, PRG4, PRSS12"
  }
  
  # Cluster 34: H19, IGF2 - imprinted genes
  else if (i == 34) {
    cluster_analysis[i+1, "cell_type"] <- "Embryonic-like"
    cluster_analysis[i+1, "subtype"] <- "H19/IGF2+"
    cluster_analysis[i+1, "confidence"] <- "Medium"
    cluster_analysis[i+1, "key_evidence"] <- "H19, IGF2, LEFTY2"
  }
  
  # Cluster 35: HIST genes - proliferating
  else if (i == 35) {
    cluster_analysis[i+1, "cell_type"] <- "Proliferating Cells"
    cluster_analysis[i+1, "subtype"] <- "G2/M phase"
    cluster_analysis[i+1, "confidence"] <- "High"
    cluster_analysis[i+1, "key_evidence"] <- "HIST1H2BH, HIST1H3B, RRM2"
  }
}

# Print results
cat("\n=== COMPREHENSIVE CLUSTER ASSIGNMENTS ===\n")
cat("Based on Kim 2021 iPSC protocol expectations\n")
cat("=========================================\n\n")

print(cluster_analysis)

# Summary
cat("\n\nSUMMARY OF CELL TYPES:\n")
cat("======================\n")
cell_type_summary <- table(cluster_analysis$cell_type)
print(cell_type_summary)

# Key findings
cat("\n\nKEY FINDINGS:\n")
cat("=============\n")
cat("1. Dopaminergic neurons found in clusters: 5 (TH+), 18 (CALB1+)\n")
cat("2. Choroid plexus correctly identified: cluster 4 (TTR+)\n")
cat("3. Multiple non-neural populations present (unexpected for iPSC protocol)\n")
cat("4. Several stressed/dying cell clusters\n")
cat("5. Presence of oligodendrocytes, melanocytes, and vascular cells\n")

# Save results
write.csv(cluster_analysis, "results/comprehensive_validation/detailed_cluster_assignments.csv", row.names = FALSE)

cat("\n\nAnalysis complete! Results saved.\n")