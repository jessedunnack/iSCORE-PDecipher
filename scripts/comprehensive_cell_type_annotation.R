#!/usr/bin/env Rscript

# Comprehensive Cell Type Annotation for Midbrain Organoid/Culture
# Incorporates web search findings and biological knowledge

library(dplyr)

# Updated cell type assignments based on marker analysis and web searches
CELL_TYPE_ASSIGNMENTS <- list(
  # Neuronal populations
  cluster_5 = list(
    type = "Dopaminergic Neurons",
    subtype = "General/Immature",
    markers = c("TH", "DDC"),
    confidence = "High"
  ),
  cluster_18 = list(
    type = "Dopaminergic Neurons", 
    subtype = "A10 (VTA)",
    markers = c("NR4A2", "CALB1", "SLC17A6", "CALB2", "LHX9"),
    confidence = "High"
  ),
  cluster_28 = list(
    type = "Dopaminergic Neurons",
    subtype = "A9 (SNc) / Neuroendocrine",
    markers = c("SLC18A2", "KCNJ6", "RET", "CHGA", "CHGB", "SST", "ASCL1"),
    confidence = "High"
  ),
  cluster_0 = list(
    type = "Immature Neurons",
    subtype = "Neuroblasts",
    markers = c("STMN2", "TUBB3", "MAPT", "NSG2", "INA"),
    confidence = "High"
  ),
  cluster_30 = list(
    type = "Hypothalamic Neurons",
    subtype = "POU3F2+",
    markers = c("POU3F2", "ATOH1", "INSM1", "NHLH1", "LHX2"),
    confidence = "Medium"
  ),
  
  # Progenitor populations
  cluster_1 = list(
    type = "Floor Plate Progenitors",
    subtype = "CORIN+",
    markers = c("CORIN", "ALCAM", "ARX"),
    confidence = "Medium"
  ),
  cluster_3 = list(
    type = "Floor Plate/Midbrain Progenitors",
    subtype = "EN2+",
    markers = c("EN2", "DLK1", "PTPRZ1", "MDK"),
    confidence = "Medium"
  ),
  cluster_10 = list(
    type = "Midbrain Progenitors",
    subtype = "MSX1+",
    markers = c("MSX1", "RSPO3", "LY6H"),
    confidence = "Medium"
  ),
  cluster_13 = list(
    type = "Floor Plate-like",
    subtype = "SHH+",
    markers = c("SHH", "DCN", "LUM", "NTN1"),
    confidence = "Medium"
  ),
  cluster_26 = list(
    type = "Neural Progenitors",
    subtype = "NES+",
    markers = c("NES", "TTR", "MEF2C"),
    confidence = "Low"
  ),
  
  # Glial populations
  cluster_19 = list(
    type = "Oligodendrocytes",
    subtype = "Mature",
    markers = c("SOX10", "PLP1", "MPZ", "S100B"),
    confidence = "High"
  ),
  
  # Non-neural populations
  cluster_4 = list(
    type = "Choroid Plexus Epithelial Cells",
    subtype = "TTR+/FOLR1+",
    markers = c("TTR", "FOLR1", "HTR2C", "TRPM3", "CLIC6"),
    confidence = "High"
  ),
  cluster_15 = list(
    type = "Leptomeningeal Cells",
    subtype = "PTGDS+",
    markers = c("PTGDS", "ELN", "LCNL1", "SPARCL1"),
    confidence = "High"
  ),
  cluster_27 = list(
    type = "Vascular Smooth Muscle/Pericytes",
    subtype = "ACTA2+",
    markers = c("TAGLN", "ACTA2", "MYL9", "TPM1", "TPM2"),
    confidence = "High"
  ),
  cluster_7 = list(
    type = "Mesenchymal/Fibroblasts",
    subtype = "Collagen-producing",
    markers = c("COL3A1", "COL1A1", "LGALS1", "PRRX1"),
    confidence = "High"
  ),
  cluster_21 = list(
    type = "Mesenchymal/ECM-producing",
    subtype = "COL1A2+",
    markers = c("COL1A2", "COL1A1", "SPARC", "ELN"),
    confidence = "High"
  ),
  
  # Proliferating populations
  cluster_22 = list(
    type = "Proliferating Cells",
    subtype = "G2/M phase",
    markers = c("TOP2A", "MKI67", "CDK1", "CENPF", "UBE2C"),
    confidence = "High"
  ),
  cluster_35 = list(
    type = "Proliferating Cells",
    subtype = "G2/M phase",
    markers = c("NUSAP1", "TOP2A", "CDK1", "BIRC5"),
    confidence = "High"
  ),
  
  # Stressed/artifact populations
  cluster_24 = list(
    type = "Stressed Cells",
    subtype = "ER stress/UPR",
    markers = c("GDF15", "DDIT3", "ATF5", "TRIB3", "CDKN1A"),
    confidence = "High"
  ),
  cluster_6 = list(
    type = "Stressed/Dying Cells",
    subtype = "Mitochondrial stress",
    markers = c("MT-ND4", "MT-CO2", "MT-ATP6", "MT-CO3"),
    confidence = "High"
  ),
  cluster_20 = list(
    type = "Technical Artifact",
    subtype = "lncRNA-enriched",
    markers = c("MALAT1", "NEAT1", "XIST", "FTX"),
    confidence = "Medium"
  ),
  
  # Caudal neural populations
  cluster_29 = list(
    type = "Caudal Spinal Neurons",
    subtype = "HOXD10/11+",
    markers = c("HOXD11", "HOXD10", "CALB1", "PENK", "ARX"),
    confidence = "High"
  ),
  cluster_31 = list(
    type = "Caudal Neural/Spinal",
    subtype = "RBP4+",
    markers = c("RBP4", "HOXA9", "HOTAIRM1"),
    confidence = "Medium"
  ),
  
  # Other/Unknown populations requiring further investigation
  cluster_2 = list(
    type = "Mixed/Transitional",
    subtype = "Unknown",
    markers = c("CMTM8", "S100A10", "ANXA1", "ARX"),
    confidence = "Low"
  ),
  cluster_8 = list(
    type = "Progenitor-like",
    subtype = "PTN+",
    markers = c("PTN", "PTPRZ1", "NKAIN3", "FEZF1"),
    confidence = "Low"
  ),
  cluster_9 = list(
    type = "Unknown",
    subtype = "OTX2+ (Forebrain/Midbrain?)",
    markers = c("TFPI2", "GNG11", "OTX2"),
    confidence = "Low"
  ),
  cluster_11 = list(
    type = "Unknown Neural",
    subtype = "PRSS56+",
    markers = c("PRSS56", "ECEL1", "DLK1", "RTL1"),
    confidence = "Low"
  ),
  cluster_12 = list(
    type = "Unknown Neural",
    subtype = "CRABP1+",
    markers = c("CRABP1", "IGFBPL1", "DCC", "LHX9"),
    confidence = "Low"
  ),
  cluster_14 = list(
    type = "Unknown/Melanocyte-like?",
    subtype = "Pigment genes",
    markers = c("TYRP1", "DCT", "MLANA", "PMEL", "KCNJ6"),
    confidence = "Low"
  ),
  cluster_16 = list(
    type = "Unknown Progenitor",
    subtype = "RGCC+",
    markers = c("RGCC", "NMU", "NR4A3", "TNFAIP3"),
    confidence = "Low"
  ),
  cluster_17 = list(
    type = "Unknown Neural",
    subtype = "GRIA2+",
    markers = c("GRIA2", "GRIA4", "DSCAM"),
    confidence = "Low"
  ),
  cluster_23 = list(
    type = "Unknown Progenitor",
    subtype = "LINC02539+",
    markers = c("LINC02539", "RGCC", "ECEL1", "RTL1"),
    confidence = "Low"
  ),
  cluster_25 = list(
    type = "Unknown",
    subtype = "Metabolic?",
    markers = c("HPD", "SLCO1B3", "CYP1B1", "TRPM3"),
    confidence = "Low"
  ),
  cluster_32 = list(
    type = "Unknown/Ciliated?",
    subtype = "C11orf88+",
    markers = c("C11orf88", "ROPN1L", "FAM183A", "CAPS"),
    confidence = "Low"
  ),
  cluster_33 = list(
    type = "Unknown Mesenchymal",
    subtype = "CRYAB+",
    markers = c("CRYAB", "ANXA2", "TNC", "ANGPTL2"),
    confidence = "Low"
  ),
  cluster_34 = list(
    type = "Unknown Progenitor",
    subtype = "H19+/SHH+",
    markers = c("H19", "IGF2", "SHH", "SMOC1"),
    confidence = "Low"
  )
)

# Function to create comprehensive annotation table
create_annotation_table <- function() {
  annotation_df <- data.frame(
    fine_cluster = integer(),
    cell_type = character(),
    subtype = character(),
    confidence = character(),
    key_markers = character(),
    biological_notes = character(),
    stringsAsFactors = FALSE
  )
  
  for (cluster_name in names(CELL_TYPE_ASSIGNMENTS)) {
    cluster_num <- as.integer(gsub("cluster_", "", cluster_name))
    info <- CELL_TYPE_ASSIGNMENTS[[cluster_name]]
    
    annotation_df <- rbind(annotation_df, data.frame(
      fine_cluster = cluster_num,
      cell_type = info$type,
      subtype = info$subtype,
      confidence = info$confidence,
      key_markers = paste(info$markers, collapse = "; "),
      biological_notes = "",
      stringsAsFactors = FALSE
    ))
  }
  
  # Sort by cluster number
  annotation_df <- annotation_df[order(annotation_df$fine_cluster), ]
  
  # Add biological notes
  notes <- c(
    "0" = "Immature neurons with high neuroblast markers",
    "1" = "Floor plate progenitor markers suggest ventral midbrain origin",
    "3" = "EN2 expression confirms midbrain identity",
    "4" = "Classic choroid plexus markers; likely contamination or co-culture",
    "5" = "TH+/DDC+ dopaminergic neurons, possibly immature",
    "6" = "High mitochondrial gene expression suggests dying/stressed cells",
    "7" = "Mesenchymal contamination or vascular support cells",
    "15" = "PTGDS is specific for leptomeningeal cells",
    "18" = "A10/VTA dopaminergic neurons with CALB2/SLC17A6 co-expression",
    "19" = "Mature oligodendrocytes with myelin genes",
    "20" = "High lncRNA may indicate technical artifact",
    "22" = "Active cell cycle markers indicate proliferation",
    "24" = "ER stress response genes suggest culture stress",
    "27" = "Smooth muscle actin markers indicate vascular cells",
    "28" = "A9/SNc dopaminergic with neuroendocrine features",
    "29" = "HOX gene expression indicates caudalization",
    "30" = "POU3F2 suggests hypothalamic identity"
  )
  
  for (i in 1:nrow(annotation_df)) {
    cluster <- as.character(annotation_df$fine_cluster[i])
    if (cluster %in% names(notes)) {
      annotation_df$biological_notes[i] <- notes[[cluster]]
    }
  }
  
  return(annotation_df)
}

# Function to summarize cell type distribution
summarize_cell_types <- function(annotation_df) {
  cat("\n=== Cell Type Distribution Summary ===\n\n")
  
  # Group by major cell type
  type_summary <- annotation_df %>%
    group_by(cell_type) %>%
    summarise(
      n_clusters = n(),
      clusters = paste(fine_cluster, collapse = ", "),
      .groups = 'drop'
    ) %>%
    arrange(desc(n_clusters))
  
  print(type_summary)
  
  # Confidence distribution
  cat("\n=== Confidence Distribution ===\n")
  conf_table <- table(annotation_df$confidence)
  print(conf_table)
  
  # Key findings
  cat("\n=== Key Findings ===\n")
  cat("- 3 dopaminergic neuron clusters identified (5, 18, 28)\n")
  cat("- Multiple non-neural populations present (choroid plexus, meninges, vasculature)\n")
  cat("- Evidence of caudalization in some clusters (29, 31)\n")
  cat("- Several stressed/proliferating populations\n")
  cat("- Many clusters require further investigation\n")
}

# Main function
main <- function() {
  cat("Comprehensive Cell Type Annotation Summary\n")
  cat("=========================================\n\n")
  
  # Create annotation table
  annotations <- create_annotation_table()
  
  # Save annotations
  dir.create("results/cell_type_annotations", recursive = TRUE, showWarnings = FALSE)
  write.csv(annotations, "results/cell_type_annotations/comprehensive_fine_cluster_annotations.csv",
            row.names = FALSE)
  
  # Print summary
  summarize_cell_types(annotations)
  
  # Save as RDS for integration
  saveRDS(annotations, "results/cell_type_annotations/fine_cluster_annotations.rds")
  
  cat("\nAnnotations saved to results/cell_type_annotations/\n")
  
  return(annotations)
}

# Run if executed directly
if (!interactive()) {
  annotations <- main()
}