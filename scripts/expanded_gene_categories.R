#!/usr/bin/env Rscript

# Expanded Gene Categories Based on Analysis Results
# Includes newly identified marker genes from cluster analysis

# Additional categories to add to FUNCTIONAL_CATEGORIES
EXPANDED_CATEGORIES <- list(
  # Additional transcription factors found
  tf_additional = list(
    genes = c("POU2F2", "ONECUT1", "ONECUT2", "ONECUT3", "SP9", "IRX3", "LMO2", "LMO3", 
              "DACH2", "SOX5", "NRN1", "HOXB6", "ZIC1", "ZIC2", "ZIC3", "ZIC4", "ZIC5"),
    description = "Additional transcription factors"
  ),
  
  # Cell signaling and receptors
  signaling_gpcr = list(
    genes = c("GPR22", "GPR37", "EDNRB", "PTGER3", "DRD1", "DRD2", "VIPR2", "HTR2C", 
              "SSTR1", "SSTR2", "FFAR4", "ADGRL1", "ADGRL2", "ADGRL3"),
    description = "G-protein coupled receptors"
  ),
  
  signaling_secreted = list(
    genes = c("RGCC", "DLK1", "RTL1", "APCDD1", "GDF7", "GDF15", "CCN2", "CCN3", 
              "ANGPTL1", "ANGPTL2", "ANGPTL4", "BMP6", "LEFTY2", "CHRDL1", "FRZB"),
    description = "Secreted signaling molecules"
  ),
  
  # Metabolic and transport
  transport_channels = list(
    genes = c("KCNK1", "KCNK12", "KCNMB2", "KCNIP2", "KCNQ1", "CACNG2", "TMEM35A",
              "CLIC6", "ATP1A2", "ATP11C", "SLC10A4", "SLCO1B3"),
    description = "Ion channels and transporters"
  ),
  
  # Novel neuronal markers
  neuronal_specific = list(
    genes = c("OLFM1", "MYT1L", "ELAVL4", "TAC1", "SNCB", "SNCG", "CBLN1", "CPLX1",
              "RAB3A", "RAB33A", "CAMK2N2", "HPCAL4", "VGF", "STMN4", "GNG3"),
    description = "Neuron-specific proteins"
  ),
  
  # Immune/inflammation markers
  immune_related = list(
    genes = c("PTPRC", "C6", "CFB", "IL18", "IFI44L", "IFITM1", "IFITM2", "IFITM3",
              "MARCH1", "MS4A6A", "LY6H", "CD36"),
    description = "Immune/inflammation markers"
  ),
  
  # Cell surface and membrane proteins
  membrane_proteins = list(
    genes = c("TMEM190", "TMEM158", "TMEM212", "TMEM35A", "CDHR3", "TMC5", "TPBG",
              "EPCAM", "MCAM", "VCAM1", "CLDN11", "MAL"),
    description = "Transmembrane proteins"
  ),
  
  # Extracellular matrix modifiers
  ecm_modifiers = list(
    genes = c("ADAMTS3", "ADAMTS16", "ADAMTS17", "MMP17", "MMP23B", "SULF1", "NDST4",
              "HS3ST3A1", "HS3ST3B1", "HAPLN1"),
    description = "ECM modifying enzymes"
  ),
  
  # Proliferation/cell cycle related
  proliferation_related = list(
    genes = c("RGCC", "CDKN2D", "GADD45A", "PPP1R17", "PPP1R1A", "PPP2R2C"),
    description = "Cell cycle regulation"
  ),
  
  # Lipid metabolism
  lipid_metabolism = list(
    genes = c("FABP7", "PLTP", "APOE", "APOC1", "LPL", "LDLR", "SOAT1"),
    description = "Lipid metabolism"
  ),
  
  # Ciliary/ependymal specific
  ciliary_proteins = list(
    genes = c("CFAP43", "CFAP44", "DNAH5", "RSPH1", "SPAG17", "ROPN1L", "C11orf88",
              "FAM183A", "CAPS"),
    description = "Ciliary/flagellar proteins"
  ),
  
  # Neural development specific
  neural_development = list(
    genes = c("CRABP1", "IGFBPL1", "PRSS56", "CHRND", "NTS", "CNTFR", "SORL1",
              "PKDCC", "SHISA9", "DCC", "ROBO1", "ROBO2"),
    description = "Neural development"
  ),
  
  # Calcium signaling
  calcium_signaling = list(
    genes = c("CALB1", "CALB2", "CALML4", "HPCAL4", "S100B", "S100A10", "ANXA1", "ANXA2"),
    description = "Calcium binding proteins"
  ),
  
  # Neuropeptides and hormones
  neuropeptides = list(
    genes = c("TAC1", "PENK", "NTS", "NPY", "VIP", "SST", "POMC", "AGRP", "HCRT",
              "AVP", "OXT", "CRH", "TRH", "GHRH", "PTH2", "NMU"),
    description = "Neuropeptides and hormones"
  ),
  
  # Long non-coding RNAs
  lncrna = list(
    genes = c("MALAT1", "NEAT1", "XIST", "H19", "MEG3", "KCNQ1OT1", "FTX", "MIAT",
              "LINC00342", "LINC02539", "LINC01088", "LINC00472", "LINC01508"),
    description = "Long non-coding RNAs"
  )
)

# Function to merge categories
merge_categories <- function() {
  # Load original categories
  source("scripts/functional_gene_categorization.R")
  
  # Merge with expanded categories
  all_categories <- c(FUNCTIONAL_CATEGORIES, EXPANDED_CATEGORIES)
  
  return(all_categories)
}

# Function to analyze with expanded categories
analyze_with_expanded <- function(marker_file) {
  all_categories <- merge_categories()
  
  markers <- readRDS(marker_file)
  clusters <- sort(unique(as.character(markers$cluster)))
  
  # Create summary of categorization coverage
  coverage_summary <- data.frame(
    cluster = character(),
    total_genes = integer(),
    categorized = integer(),
    uncategorized = integer(),
    coverage_pct = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (cluster_id in clusters) {
    cluster_markers <- markers %>%
      filter(cluster == cluster_id, avg_log2FC > 0.25, p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC * -log10(p_val_adj + 1e-300)))
    
    top_genes <- head(cluster_markers$gene, 100)
    
    # Categorize with expanded set
    all_categorized <- character()
    for (cat_name in names(all_categories)) {
      cat_genes <- all_categories[[cat_name]]$genes
      matching <- intersect(top_genes, cat_genes)
      all_categorized <- c(all_categorized, matching)
    }
    
    all_categorized <- unique(all_categorized)
    uncategorized <- setdiff(top_genes, all_categorized)
    
    coverage_summary <- rbind(coverage_summary, data.frame(
      cluster = cluster_id,
      total_genes = length(top_genes),
      categorized = length(all_categorized),
      uncategorized = length(uncategorized),
      coverage_pct = round(100 * length(all_categorized) / length(top_genes), 1),
      stringsAsFactors = FALSE
    ))
  }
  
  return(list(
    categories = all_categories,
    coverage = coverage_summary
  ))
}

# Function to create final gene panels for visualization
create_gene_panels <- function() {
  panels <- list(
    # Core dopaminergic differentiation panel
    dopaminergic_core = c("FOXA2", "LMX1A", "LMX1B", "NR4A2", "PITX3", "EN1", "EN2",
                         "TH", "DDC", "SLC6A3", "SLC18A2", "KCNJ6", "ALDH1A1",
                         "CALB1", "CALB2", "SOX6", "DRD2", "RET"),
    
    # Neuronal maturation panel
    neuronal_maturation = c("SOX2", "NES", "DCX", "TUBB3", "STMN2", "MAP2", "NEUN",
                           "RBFOX3", "SYN1", "SNAP25", "SYT1", "CAMK2A"),
    
    # Regional identity panel
    regional_identity = c("OTX2", "EN1", "EN2", "LMX1A", "FOXA2", "SHH", "CORIN",
                         "HOXD9", "HOXD10", "HOXD11", "HOXA9", "LHX2", "LHX9",
                         "POU3F2", "NKX2.1", "PAX6"),
    
    # Cell adhesion/migration panel
    adhesion_migration = c("DCC", "ROBO1", "ROBO2", "SLIT2", "NTN1", "UNC5C",
                          "L1CAM", "NCAM1", "NRCAM", "ALCAM", "CDH2", "CDH6",
                          "CDH8", "CDH10", "CDH11", "DSCAM"),
    
    # Stress/quality control panel
    stress_quality = c("DDIT3", "ATF3", "ATF4", "ATF5", "GDF15", "CDKN1A",
                      "HSPA1A", "HSP90AA1", "SQSTM1", "FOS", "JUN", "EGR1"),
    
    # Non-neuronal markers panel
    non_neuronal = c("GFAP", "AQP4", "S100B", "OLIG2", "SOX10", "MBP", "PLP1",
                    "TAGLN", "ACTA2", "MYL9", "TTR", "FOLR1", "PTGDS"),
    
    # Novel cluster-defining panel
    novel_markers = c("RGCC", "RTL1", "DLK1", "PRSS56", "TMEM190", "IGFBPL1",
                     "CRABP1", "LY6H", "APCDD1", "OLFM1", "MYT1L", "VGF")
  )
  
  return(panels)
}

# Main function
main <- function() {
  cat("Expanded Functional Gene Categorization\n")
  cat("======================================\n\n")
  
  # Analyze with expanded categories
  marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
  results <- analyze_with_expanded(marker_file)
  
  # Print coverage summary
  cat("Categorization Coverage Summary:\n")
  print(results$coverage)
  
  cat("\nAverage coverage: ", round(mean(results$coverage$coverage_pct), 1), "%\n", sep="")
  
  # Create gene panels
  panels <- create_gene_panels()
  
  cat("\n\n=== Gene Panels for Visualization ===\n")
  for (panel_name in names(panels)) {
    cat("\n", gsub("_", " ", toupper(panel_name)), " (", length(panels[[panel_name]]), " genes):\n", sep="")
    cat(paste(panels[[panel_name]], collapse=", "), "\n")
  }
  
  # Save expanded categories
  saveRDS(results$categories, "results/functional_categorization/expanded_categories.rds")
  saveRDS(panels, "results/functional_categorization/visualization_panels.rds")
  
  # Export as CSV
  cat_df <- data.frame(
    category = character(),
    description = character(),
    n_genes = integer(),
    stringsAsFactors = FALSE
  )
  
  for (cat_name in names(results$categories)) {
    cat_info <- results$categories[[cat_name]]
    cat_df <- rbind(cat_df, data.frame(
      category = cat_name,
      description = cat_info$description,
      n_genes = length(cat_info$genes),
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(cat_df, "results/functional_categorization/expanded_category_summary.csv",
            row.names = FALSE)
  
  return(list(
    categories = results$categories,
    coverage = results$coverage,
    panels = panels
  ))
}

# Run if executed directly
if (!interactive()) {
  library(dplyr)
  results <- main()
}