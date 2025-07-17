# Analyze PD biological relevance of top Fisher's exact test results
# Link significant overlaps to known PD pathways and biology

# Load enrichment data to examine pathway overlap
data_file <- '/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds'
enrichment_data <- readRDS(data_file)

# Define PD-relevant pathways
pd_pathways <- list(
  mitophagy = c("mitophagy", "mitochondr", "autophagy of mitochondrion"),
  mitochondrial = c("mitochondrial", "oxidative phosphorylation", "respiratory chain", 
                    "ATP synthesis", "electron transport"),
  lysosomal = c("lysosom", "autophagy", "protein degradation", "vacuolar"),
  synaptic = c("synap", "neurotransmitter", "dopamine", "vesicle", "exocytosis"),
  protein_aggregation = c("protein folding", "misfolded", "ubiquitin", "proteasome"),
  neuroinflammation = c("inflamm", "cytokine", "microglial", "immune response"),
  oxidative_stress = c("oxidative stress", "reactive oxygen", "antioxidant", "redox")
)

# Function to check PD relevance
check_pd_relevance <- function(gene_pair, enrichment_data) {
  genes <- strsplit(gene_pair, "_vs_")[[1]]
  
  # Get enrichment terms for this gene pair
  gene_terms <- enrichment_data[
    enrichment_data$mutation_perturbation %in% genes &
    enrichment_data$p.adjust < 0.05,
  ]
  
  # Count PD-relevant pathways
  pd_counts <- list()
  
  for (category in names(pd_pathways)) {
    keywords <- pd_pathways[[category]]
    pattern <- paste0("(", paste(keywords, collapse="|"), ")")
    
    matching_terms <- gene_terms[
      grepl(pattern, gene_terms$Description, ignore.case = TRUE),
    ]
    
    pd_counts[[category]] <- length(unique(matching_terms$Description))
  }
  
  return(pd_counts)
}

# Analyze top gene pairs from Fisher's results
top_pairs <- c("SYNJ1_vs_SYNJ1", "ATP13A2_vs_ATP13A2", "FBXO7_vs_FBXO7", 
               "LRRK2_vs_LRRK2", "PARK7_vs_PARK7", "PRKN_vs_PARK2")

cat("=== PD BIOLOGICAL RELEVANCE ANALYSIS ===\n\n")

pd_relevance_summary <- data.frame()

for (pair in top_pairs) {
  cat("\n", pair, ":\n", sep="")
  cat(rep("-", nchar(pair) + 1), "\n", sep="")
  
  # Get PD pathway counts
  pd_counts <- check_pd_relevance(pair, enrichment_data)
  
  # Print summary
  total_pd_pathways <- sum(unlist(pd_counts))
  cat("Total PD-relevant pathways:", total_pd_pathways, "\n")
  
  for (category in names(pd_counts)) {
    if (pd_counts[[category]] > 0) {
      cat("  ", category, ": ", pd_counts[[category]], " pathways\n", sep="")
    }
  }
  
  # Get specific gene biology
  gene <- strsplit(pair, "_vs_")[[1]][1]
  
  gene_biology <- switch(gene,
    "SYNJ1" = "Synaptic vesicle uncoating, phosphoinositide metabolism",
    "ATP13A2" = "Lysosomal ATPase, mitochondrial function, α-synuclein clearance",
    "FBXO7" = "Mitophagy regulation, proteasome function, PINK1-Parkin pathway",
    "LRRK2" = "Kinase activity, vesicle trafficking, autophagy regulation",
    "PARK7" = "DJ-1 protein, oxidative stress response, mitochondrial protection",
    "PRKN" = "Parkin E3 ubiquitin ligase, mitophagy, protein quality control",
    "Unknown function"
  )
  
  cat("\nKnown PD biology: ", gene_biology, "\n", sep="")
  
  # Get top enriched pathways for this gene
  gene_name <- strsplit(pair, "_vs_")[[1]][1]
  top_terms <- enrichment_data[
    enrichment_data$mutation_perturbation == gene_name &
    enrichment_data$p.adjust < 0.001,
  ]
  
  if (nrow(top_terms) > 0) {
    top_terms <- top_terms[order(top_terms$p.adjust), ]
    cat("\nTop 5 enriched pathways:\n")
    for (i in 1:min(5, nrow(top_terms))) {
      cat("  - ", top_terms$Description[i], " (", top_terms$enrichment_type[i], 
          ", p=", format(top_terms$p.adjust[i], scientific=TRUE), ")\n", sep="")
    }
  }
  
  # Store summary
  pd_relevance_summary <- rbind(pd_relevance_summary, data.frame(
    gene_pair = pair,
    total_pd_pathways = total_pd_pathways,
    mitophagy = pd_counts$mitophagy,
    mitochondrial = pd_counts$mitochondrial,
    lysosomal = pd_counts$lysosomal,
    synaptic = pd_counts$synaptic,
    protein_aggregation = pd_counts$protein_aggregation,
    neuroinflammation = pd_counts$neuroinflammation,
    oxidative_stress = pd_counts$oxidative_stress
  ))
}

cat("\n\n=== SUMMARY TABLE ===\n")
print(pd_relevance_summary)

# Identify convergent PD signatures
cat("\n\n=== CONVERGENT PD SIGNATURES ===\n")
cat("Genes with overlapping PD-relevant pathways across methods:\n\n")

# Focus on genes with strong mitochondrial/mitophagy signatures
mito_genes <- pd_relevance_summary[
  pd_relevance_summary$mitophagy > 5 | pd_relevance_summary$mitochondrial > 10,
]

if (nrow(mito_genes) > 0) {
  cat("Strong mitochondrial dysfunction signatures:\n")
  print(mito_genes[, c("gene_pair", "mitophagy", "mitochondrial")])
}

# Focus on protein quality control
protein_genes <- pd_relevance_summary[
  pd_relevance_summary$protein_aggregation > 5 | pd_relevance_summary$lysosomal > 10,
]

if (nrow(protein_genes) > 0) {
  cat("\nStrong protein quality control signatures:\n")
  print(protein_genes[, c("gene_pair", "protein_aggregation", "lysosomal")])
}

# Save results
write.csv(pd_relevance_summary, "pd_biological_relevance.csv", row.names = FALSE)
cat("\n\nResults saved to pd_biological_relevance.csv\n")