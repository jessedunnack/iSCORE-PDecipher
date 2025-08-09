#!/usr/bin/env Rscript

# Functional Categorization of Marker Genes from Fine Clusters
# Groups genes by biological function for downstream visualization

library(dplyr)
library(tidyr)

# Define comprehensive functional gene categories
FUNCTIONAL_CATEGORIES <- list(
  # Transcription Factors - Cell Fate Determination
  tf_dopaminergic = list(
    genes = c("FOXA2", "LMX1A", "LMX1B", "NR4A2", "NURR1", "PITX3", "EN1", "EN2", "OTX2"),
    description = "Dopaminergic neuron fate TFs"
  ),
  
  tf_neurogenesis = list(
    genes = c("ASCL1", "NGN2", "NEUROG1", "NEUROG2", "NEUROD1", "NEUROD2", "NEUROD4", "NEUROD6",
              "INSM1", "NHLH1", "ATOH1", "HES1", "HES5", "HES6"),
    description = "Neurogenesis and differentiation TFs"
  ),
  
  tf_regional_identity = list(
    genes = c("HOXD9", "HOXD10", "HOXD11", "HOXA9", "LHX2", "LHX5", "LHX9", "NKX2.1", "NKX2.2",
              "PAX6", "PAX7", "POU3F2", "OTP", "SIM1", "SIM2", "FOXD1", "FOXG1"),
    description = "Regional patterning TFs"
  ),
  
  tf_progenitor = list(
    genes = c("SOX2", "SOX4", "SOX9", "SOX10", "PAX6", "MSX1", "MSX2", "ZIC1", "ZIC3", "ZIC4"),
    description = "Neural progenitor TFs"
  ),
  
  tf_other = list(
    genes = c("NR2F1", "NR4A3", "MEF2C", "ONECUT2", "ARX", "DLX1", "DLX2", "DLX5", "DLX6"),
    description = "Other neural TFs"
  ),
  
  # Pan-Neuronal Markers
  pan_neuronal_mature = list(
    genes = c("MAP2", "NEUN", "RBFOX3", "SYN1", "SYP", "SNAP25", "SYT1", "VAMP2", "STX1A"),
    description = "Mature neuron markers"
  ),
  
  pan_neuronal_structural = list(
    genes = c("TUBB3", "TUBA1A", "MAPT", "TAU", "NEFL", "NEFM", "NEFH", "INA", "PRPH"),
    description = "Neuronal cytoskeleton"
  ),
  
  pan_neuronal_immature = list(
    genes = c("DCX", "STMN1", "STMN2", "STMN3", "STMN4", "GAP43", "NCAM1", "NSG1", "NSG2"),
    description = "Immature neuron markers"
  ),
  
  # Neurotransmitter Identity
  nt_dopaminergic = list(
    genes = c("TH", "DDC", "SLC6A3", "DAT", "SLC18A2", "VMAT2", "DRD2", "ALDH1A1", "KCNJ6", "GIRK2"),
    description = "Dopaminergic markers"
  ),
  
  nt_gabaergic = list(
    genes = c("GAD1", "GAD2", "SLC32A1", "VGAT", "CALB1", "CALB2", "PVALB", "SST", "VIP"),
    description = "GABAergic markers"
  ),
  
  nt_glutamatergic = list(
    genes = c("SLC17A6", "VGLUT2", "SLC17A7", "VGLUT1", "SLC17A8", "VGLUT3", "GRIA2", "GRIA4", "GRIN2B"),
    description = "Glutamatergic markers"
  ),
  
  nt_other = list(
    genes = c("CHAT", "SLC18A3", "TPH1", "TPH2", "HTR2C", "CHGA", "CHGB", "PENK", "POMC"),
    description = "Other neurotransmitter systems"
  ),
  
  # Axon Guidance and Outgrowth
  axon_guidance_receptors = list(
    genes = c("DCC", "UNC5A", "UNC5B", "UNC5C", "UNC5D", "ROBO1", "ROBO2", "ROBO3", "ROBO4",
              "EPHA3", "EPHA4", "EPHA5", "EPHA7", "EPHB1", "EPHB2", "PLXNA1", "PLXNA2", "PLXNB1"),
    description = "Axon guidance receptors"
  ),
  
  axon_guidance_ligands = list(
    genes = c("NTN1", "NTN3", "NTN4", "SLIT1", "SLIT2", "SLIT3", "EFNA1", "EFNA5", "EFNB1", 
              "EFNB2", "EFNB3", "SEMA3A", "SEMA3F", "SEMA4D", "SEMA6A"),
    description = "Axon guidance ligands"
  ),
  
  axon_growth = list(
    genes = c("GAP43", "BASP1", "SCG2", "VGF", "RTN1", "RTN4", "NOGO", "L1CAM", "CHL1"),
    description = "Axon growth and regeneration"
  ),
  
  # Cell Adhesion Molecules
  adhesion_ig_superfamily = list(
    genes = c("NCAM1", "NCAM2", "L1CAM", "CHL1", "NRCAM", "CNTN1", "CNTN2", "CNTNAP1", "CNTNAP2", "CNTNAP5"),
    description = "Ig superfamily CAMs"
  ),
  
  adhesion_cadherins = list(
    genes = c("CDH1", "CDH2", "CDH4", "CDH6", "CDH8", "CDH9", "CDH10", "CDH11", "CDH13", "CDH19", "CDH20"),
    description = "Cadherin family"
  ),
  
  adhesion_other = list(
    genes = c("ALCAM", "MCAM", "PECAM1", "VCAM1", "ICAM1", "DSCAM", "CTNNA2", "NLGN1", "NRXN1"),
    description = "Other adhesion molecules"
  ),
  
  # Synaptic Components
  synaptic_vesicle = list(
    genes = c("SYN1", "SYN2", "SYP", "SYT1", "SYT4", "SNAP25", "VAMP2", "STX1A", "CPLX1", "CPLX2"),
    description = "Synaptic vesicle proteins"
  ),
  
  synaptic_structural = list(
    genes = c("BSN", "PCLO", "RIMS1", "UNC13A", "CADPS", "CADPS2", "SYN1", "SYNPO"),
    description = "Synaptic structural proteins"
  ),
  
  synaptic_receptors = list(
    genes = c("GRIA1", "GRIA2", "GRIA3", "GRIA4", "GRIN1", "GRIN2A", "GRIN2B", "GABRA1", "GABRB2"),
    description = "Synaptic receptors"
  ),
  
  # Growth Factors and Signaling
  growth_factors = list(
    genes = c("BDNF", "NGF", "NTF3", "NTF4", "GDNF", "CDNF", "MANF", "IGF1", "IGF2", "FGF1", "FGF2"),
    description = "Neurotrophic factors"
  ),
  
  growth_factor_receptors = list(
    genes = c("NTRK1", "NTRK2", "NTRK3", "RET", "GFRA1", "GFRA2", "GFRA3", "LIFR", "NGFR"),
    description = "Growth factor receptors"
  ),
  
  wnt_signaling = list(
    genes = c("WNT1", "WNT3A", "WNT5A", "WNT7A", "LRP5", "LRP6", "FZD1", "FZD3", "SFRP1", "SFRP2"),
    description = "Wnt signaling pathway"
  ),
  
  # Cell Cycle and Proliferation
  cell_cycle_g2m = list(
    genes = c("MKI67", "TOP2A", "CDK1", "CCNB1", "CCNB2", "CENPF", "UBE2C", "BIRC5", "TPX2", "NUSAP1"),
    description = "G2/M phase markers"
  ),
  
  cell_cycle_s = list(
    genes = c("PCNA", "MCM2", "MCM3", "MCM4", "MCM5", "MCM6", "MCM7", "TYMS", "RRM2"),
    description = "S phase markers"
  ),
  
  # Extracellular Matrix
  ecm_collagens = list(
    genes = c("COL1A1", "COL1A2", "COL3A1", "COL4A1", "COL5A1", "COL5A3", "COL6A1", "COL6A3", "COL21A1", "COL27A1"),
    description = "Collagen family"
  ),
  
  ecm_proteoglycans = list(
    genes = c("DCN", "LUM", "BGN", "VCAN", "ACAN", "HSPG2", "GPC1", "GPC4", "SDC1", "SDC4"),
    description = "Proteoglycans"
  ),
  
  ecm_glycoproteins = list(
    genes = c("FN1", "LAMA1", "LAMB1", "LAMC1", "TNC", "THBS1", "SPARC", "SPARCL1", "POSTN", "FBLN1"),
    description = "ECM glycoproteins"
  ),
  
  # Metabolic and Mitochondrial
  mitochondrial = list(
    genes = c("MT-CO1", "MT-CO2", "MT-CO3", "MT-ND1", "MT-ND2", "MT-ND3", "MT-ND4", "MT-ATP6", "MT-CYB"),
    description = "Mitochondrial genes"
  ),
  
  metabolic_enzymes = list(
    genes = c("ALDH1A1", "ALDH1A2", "HPD", "CYP1B1", "MAOA", "MAOB", "COMT", "DBH", "PAH"),
    description = "Metabolic enzymes"
  ),
  
  # Stress Response
  stress_er = list(
    genes = c("DDIT3", "ATF3", "ATF4", "ATF5", "HSPA5", "HSP90B1", "PDIA3", "CALR", "CANX"),
    description = "ER stress/UPR"
  ),
  
  stress_oxidative = list(
    genes = c("HMOX1", "NQO1", "GCLM", "GCLC", "GSR", "GPX1", "SOD1", "SOD2", "CAT", "PRDX1"),
    description = "Oxidative stress response"
  ),
  
  stress_general = list(
    genes = c("FOS", "JUN", "EGR1", "HSPA1A", "HSPA1B", "HSP90AA1", "HSPB1", "GDF15", "CDKN1A"),
    description = "General stress response"
  ),
  
  # Specific Cell Type Markers
  oligodendrocyte = list(
    genes = c("OLIG1", "OLIG2", "SOX10", "MBP", "PLP1", "MOG", "MAG", "CNP", "CLDN11", "MOBP"),
    description = "Oligodendrocyte markers"
  ),
  
  astrocyte = list(
    genes = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2", "SLC1A3", "GJA1", "FABP7", "SOX9"),
    description = "Astrocyte markers"
  ),
  
  ependymal = list(
    genes = c("TTR", "FOLR1", "HTR2C", "CLIC6", "TMEM212", "CFAP43", "CFAP44", "DNAH5"),
    description = "Ependymal/choroid plexus markers"
  ),
  
  meningeal = list(
    genes = c("PTGDS", "RARRES2", "LUM", "DCN", "ELN", "FBLN1", "FBLN2", "PRELP"),
    description = "Meningeal cell markers"
  ),
  
  vascular = list(
    genes = c("TAGLN", "ACTA2", "MYL9", "TPM1", "TPM2", "CNN1", "CALD1", "PECAM1", "VWF", "CDH5"),
    description = "Vascular cell markers"
  ),
  
  # Developmental Stage Markers
  dev_early_progenitor = list(
    genes = c("NES", "VIM", "SOX2", "PAX6", "HES1", "HES5", "FABP7", "SLC1A3", "BLBP"),
    description = "Early neural progenitor"
  ),
  
  dev_intermediate = list(
    genes = c("EOMES", "TBR2", "NEUROG2", "DCX", "TUBB3", "MAP2", "NCAM1", "STMN2"),
    description = "Intermediate progenitor/neuroblast"
  ),
  
  dev_mature = list(
    genes = c("RBFOX3", "MAP2", "SYN1", "SNAP25", "SYT1", "NEUN", "ENO2", "CAMK2A"),
    description = "Mature neuron markers"
  )
)

# Function to extract all unique genes from categories
get_all_categorized_genes <- function() {
  all_genes <- unique(unlist(lapply(FUNCTIONAL_CATEGORIES, function(x) x$genes)))
  return(all_genes)
}

# Function to categorize a list of genes
categorize_genes <- function(gene_list) {
  categorization <- list()
  
  for (category_name in names(FUNCTIONAL_CATEGORIES)) {
    category <- FUNCTIONAL_CATEGORIES[[category_name]]
    matching_genes <- intersect(gene_list, category$genes)
    
    if (length(matching_genes) > 0) {
      categorization[[category_name]] <- list(
        genes = matching_genes,
        count = length(matching_genes),
        description = category$description
      )
    }
  }
  
  # Find uncategorized genes
  all_categorized <- unique(unlist(lapply(categorization, function(x) x$genes)))
  uncategorized <- setdiff(gene_list, all_categorized)
  
  if (length(uncategorized) > 0) {
    categorization[["uncategorized"]] <- list(
      genes = uncategorized,
      count = length(uncategorized),
      description = "Not yet categorized"
    )
  }
  
  return(categorization)
}

# Function to analyze cluster markers by functional category
analyze_cluster_by_category <- function(cluster_markers, cluster_id) {
  cat("\n=== Cluster", cluster_id, "Functional Analysis ===\n")
  
  # Get top genes
  top_genes <- head(cluster_markers$gene, 50)
  
  # Categorize
  categories <- categorize_genes(top_genes)
  
  # Sort by count
  categories <- categories[order(sapply(categories, function(x) x$count), decreasing = TRUE)]
  
  # Print summary
  for (cat_name in names(categories)) {
    cat_info <- categories[[cat_name]]
    if (cat_info$count > 0) {
      cat(sprintf("%-30s: %2d genes - %s\n", 
                  cat_name, cat_info$count, cat_info$description))
      cat("  Genes:", paste(cat_info$genes, collapse = ", "), "\n")
    }
  }
  
  return(categories)
}

# Function to create category summary across all clusters
create_category_summary <- function(all_cluster_results) {
  # Initialize summary matrix
  category_names <- names(FUNCTIONAL_CATEGORIES)
  cluster_names <- names(all_cluster_results)
  
  summary_matrix <- matrix(0, 
                          nrow = length(category_names), 
                          ncol = length(cluster_names),
                          dimnames = list(category_names, cluster_names))
  
  # Fill matrix
  for (cluster in cluster_names) {
    cluster_cats <- all_cluster_results[[cluster]]
    for (cat in names(cluster_cats)) {
      if (cat %in% category_names) {
        summary_matrix[cat, cluster] <- cluster_cats[[cat]]$count
      }
    }
  }
  
  return(summary_matrix)
}

# Function to suggest plotting groups
suggest_plot_groups <- function() {
  plot_groups <- list(
    transcription_factors = list(
      categories = c("tf_dopaminergic", "tf_neurogenesis", "tf_regional_identity", 
                    "tf_progenitor", "tf_other"),
      title = "Transcription Factors",
      description = "Cell fate determination and differentiation"
    ),
    
    neuronal_identity = list(
      categories = c("pan_neuronal_mature", "pan_neuronal_structural", "pan_neuronal_immature",
                    "nt_dopaminergic", "nt_gabaergic", "nt_glutamatergic"),
      title = "Neuronal Identity",
      description = "Pan-neuronal and neurotransmitter markers"
    ),
    
    axon_guidance_adhesion = list(
      categories = c("axon_guidance_receptors", "axon_guidance_ligands", "axon_growth",
                    "adhesion_ig_superfamily", "adhesion_cadherins"),
      title = "Axon Guidance & Adhesion",
      description = "Neural connectivity and migration"
    ),
    
    synaptic_function = list(
      categories = c("synaptic_vesicle", "synaptic_structural", "synaptic_receptors"),
      title = "Synaptic Components",
      description = "Synapse formation and function"
    ),
    
    developmental_stage = list(
      categories = c("dev_early_progenitor", "dev_intermediate", "dev_mature", 
                    "cell_cycle_g2m", "cell_cycle_s"),
      title = "Developmental Stage",
      description = "Maturation and proliferation markers"
    ),
    
    cell_type_specific = list(
      categories = c("oligodendrocyte", "astrocyte", "ependymal", "meningeal", "vascular"),
      title = "Cell Type Specific",
      description = "Non-neuronal cell markers"
    ),
    
    signaling_metabolism = list(
      categories = c("growth_factors", "growth_factor_receptors", "wnt_signaling",
                    "metabolic_enzymes"),
      title = "Signaling & Metabolism",
      description = "Growth factors and metabolic pathways"
    ),
    
    stress_ecm = list(
      categories = c("stress_er", "stress_oxidative", "stress_general",
                    "ecm_collagens", "ecm_proteoglycans"),
      title = "Stress Response & ECM",
      description = "Cellular stress and matrix components"
    )
  )
  
  return(plot_groups)
}

# Main analysis function
main <- function() {
  cat("Functional Gene Categorization Analysis\n")
  cat("======================================\n\n")
  
  # Load fine cluster markers
  marker_file <- "inst/extdata/umap_data/backup_markers_20250719_224210/Full_Dataset_cluster_markers.rds"
  if (!file.exists(marker_file)) {
    stop("Marker file not found")
  }
  
  markers <- readRDS(marker_file)
  clusters <- sort(unique(as.character(markers$cluster)))
  
  # Analyze each cluster
  all_results <- list()
  
  for (cluster_id in clusters[1:5]) {  # Demo with first 5 clusters
    cluster_markers <- markers %>%
      filter(cluster == cluster_id, avg_log2FC > 0.25, p_val_adj < 0.05) %>%
      arrange(desc(avg_log2FC * -log10(p_val_adj + 1e-300)))
    
    results <- analyze_cluster_by_category(cluster_markers, cluster_id)
    all_results[[cluster_id]] <- results
  }
  
  # Create summary
  cat("\n\n=== Category Summary Matrix ===\n")
  summary_matrix <- create_category_summary(all_results)
  
  # Show categories with most hits
  category_totals <- rowSums(summary_matrix)
  top_categories <- head(sort(category_totals, decreasing = TRUE), 20)
  
  cat("\nTop functional categories across clusters:\n")
  for (cat in names(top_categories)) {
    if (top_categories[cat] > 0) {
      cat(sprintf("%-30s: %3d total genes\n", cat, top_categories[cat]))
    }
  }
  
  # Suggest plot groups
  cat("\n\n=== Suggested Plotting Groups ===\n")
  plot_groups <- suggest_plot_groups()
  
  for (group_name in names(plot_groups)) {
    group <- plot_groups[[group_name]]
    cat("\n", group$title, ":\n", sep="")
    cat("  Description:", group$description, "\n")
    cat("  Categories:", paste(group$categories, collapse=", "), "\n")
  }
  
  # Save results
  dir.create("results/functional_categorization", recursive = TRUE, showWarnings = FALSE)
  
  saveRDS(FUNCTIONAL_CATEGORIES, "results/functional_categorization/gene_categories.rds")
  saveRDS(all_results, "results/functional_categorization/cluster_analysis_results.rds")
  saveRDS(plot_groups, "results/functional_categorization/plot_groups.rds")
  
  # Export category definitions as CSV
  cat_df <- data.frame(
    category = character(),
    description = character(),
    genes = character(),
    stringsAsFactors = FALSE
  )
  
  for (cat_name in names(FUNCTIONAL_CATEGORIES)) {
    cat_info <- FUNCTIONAL_CATEGORIES[[cat_name]]
    cat_df <- rbind(cat_df, data.frame(
      category = cat_name,
      description = cat_info$description,
      genes = paste(cat_info$genes, collapse = "; "),
      stringsAsFactors = FALSE
    ))
  }
  
  write.csv(cat_df, "results/functional_categorization/gene_category_definitions.csv", 
            row.names = FALSE)
  
  cat("\n\nResults saved to results/functional_categorization/\n")
  
  return(list(
    categories = FUNCTIONAL_CATEGORIES,
    results = all_results,
    summary = summary_matrix,
    plot_groups = plot_groups
  ))
}

# Run if executed directly
if (!interactive()) {
  results <- main()
}