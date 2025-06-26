# Analyze Current Signature Results with PD Biological Interpretation
# This script takes your current signature analysis results and generates
# clear, biologically meaningful summaries focused on Parkinson's disease

# Load required functions
source("R/signature_analysis.R")
source("R/manuscript_signature_discovery.R") 
source("R/gene_harmonization.R")
source("R/pd_signature_interpretation.R")

# Function to run complete PD-focused analysis
run_pd_signature_analysis <- function(data_file = NULL,
                                     output_dir = "pd_signature_analysis_results") {
  
  # Auto-detect data file if not specified
  if (is.null(data_file)) {
    possible_files <- c(
      "all_enrichment_padj005_complete_with_direction.rds",  # Current directory
      "../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds",  # Dataset 2
      "../../iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds"  # Alternative
    )
    
    for (file in possible_files) {
      if (file.exists(file)) {
        data_file <- file
        break
      }
    }
    
    if (is.null(data_file)) {
      stop("Could not find data file. Please specify the path to all_enrichment_padj005_complete_with_direction.rds")
    }
  }
  
  cat("=== PARKINSON'S DISEASE SIGNATURE ANALYSIS ===\n\n")
  
  # Load data
  cat("[STEP 1] Loading enrichment data...\n")
  if (!file.exists(data_file)) {
    stop("Data file not found: ", data_file)
  }
  
  enrichment_data <- readRDS(data_file)
  cat("Loaded", nrow(enrichment_data), "enrichment terms\n\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Run signature discovery (using existing functions)
  cat("[STEP 2] Running signature discovery analysis...\n")
  
  signature_results <- discover_top_signatures(
    enrichment_data = enrichment_data,
    top_n = 15,
    min_cluster_breadth = 6,  # Slightly lower threshold for more results
    combine_variants = TRUE,
    progress_callback = function(msg, value = NULL, detail = NULL) {
      cat("  ", msg, "\n")
    }
  )
  
  cat("Signature discovery complete!\n")
  cat("- Total signatures found:", nrow(signature_results$all_signatures), "\n")
  cat("- Pan-cluster signatures:", nrow(signature_results$pan_cluster_signatures), "\n")
  cat("- Cluster-specific signatures:", length(signature_results$cluster_specific_signatures), "\n\n")
  
  # Run PD-focused analysis
  cat("[STEP 3] Analyzing signatures for PD biological relevance...\n")
  
  pd_analysis <- analyze_pd_signatures(
    signature_results = signature_results,
    enrichment_data = enrichment_data,
    focus_on_pan_cluster = TRUE
  )
  
  cat("PD analysis complete!\n")
  cat("- Enhanced signatures:", length(pd_analysis$enhanced_signatures), "\n")
  cat("- Analysis type:", pd_analysis$analysis_type, "\n\n")
  
  # Generate summary outputs
  cat("[STEP 4] Generating summary outputs...\n")
  
  # 1. Create detailed signature table
  signature_table <- create_detailed_signature_table(pd_analysis$enhanced_signatures)
  write.csv(signature_table, file.path(output_dir, "detailed_signature_analysis.csv"), row.names = FALSE)
  cat("- Saved detailed signature table\n")
  
  # 2. Create PD pathway summary
  pathway_summary <- create_pd_pathway_summary(pd_analysis$enhanced_signatures)
  write.csv(pathway_summary, file.path(output_dir, "pd_pathway_summary.csv"), row.names = FALSE)
  cat("- Saved PD pathway summary\n")
  
  # 3. Save comprehensive text summary
  writeLines(pd_analysis$pd_summary$manuscript_summary, 
             file.path(output_dir, "manuscript_ready_summary.txt"))
  cat("- Saved manuscript-ready summary\n")
  
  # 4. Create individual signature reports
  create_individual_signature_reports(pd_analysis$enhanced_signatures, output_dir)
  cat("- Saved individual signature reports\n")
  
  # 5. Print key findings to console
  cat("\n" , "="*60, "\n")
  cat("KEY FINDINGS SUMMARY\n")
  cat("="*60, "\n\n")
  
  print_key_findings(pd_analysis)
  
  cat("\n", "="*60, "\n")
  cat("Files saved to:", output_dir, "\n")
  cat("="*60, "\n")
  
  return(list(
    signature_results = signature_results,
    pd_analysis = pd_analysis,
    output_dir = output_dir
  ))
}

# Create detailed signature table for export
create_detailed_signature_table <- function(enhanced_signatures) {
  
  table_data <- data.frame()
  
  for (i in seq_along(enhanced_signatures)) {
    sig <- enhanced_signatures[[i]]
    
    # Extract key metrics
    row_data <- data.frame(
      rank = i,
      gene_pair = sig$signature$gene_pair,
      mast_gene = sig$signature$mast_gene,
      crispri_gene = sig$signature$crispri_gene,
      cluster = sig$signature$cluster,
      signature_strength = round(sig$signature$signature_strength, 3),
      gene_overlap_count = sig$signature$gene_overlap_count,
      gene_jaccard_index = round(sig$signature$gene_jaccard, 3),
      pd_relevance_score = round(sig$pd_relevance_score, 2),
      mast_pd_terms = nrow(sig$mast_pd_terms),
      crispri_pd_terms = nrow(sig$crispri_pd_terms),
      shared_pd_pathways = nrow(sig$shared_pd_pathways),
      mitochondrial_pathways = sig$biological_categories$mitochondrial,
      protein_quality_pathways = sig$biological_categories$protein_quality,
      autophagy_pathways = sig$biological_categories$autophagy,
      dopamine_pathways = sig$biological_categories$dopamine,
      synaptic_pathways = sig$biological_categories$synaptic,
      oxidative_stress_pathways = sig$biological_categories$oxidative_stress,
      neuronal_pathways = sig$biological_categories$neuronal,
      stringsAsFactors = FALSE
    )
    
    table_data <- rbind(table_data, row_data)
  }
  
  # Sort by PD relevance score
  table_data <- table_data[order(-table_data$pd_relevance_score, -table_data$signature_strength), ]
  table_data$pd_rank <- seq_len(nrow(table_data))
  
  return(table_data)
}

# Create PD pathway summary across all signatures
create_pd_pathway_summary <- function(enhanced_signatures) {
  
  # Aggregate all PD-relevant pathways
  all_pathways <- data.frame()
  
  for (sig in enhanced_signatures) {
    if (nrow(sig$shared_pd_pathways) > 0) {
      pathway_data <- sig$shared_pd_pathways
      pathway_data$gene_pair <- sig$signature$gene_pair
      pathway_data$cluster <- sig$signature$cluster
      pathway_data$signature_strength <- sig$signature$signature_strength
      all_pathways <- rbind(all_pathways, pathway_data)
    }
  }
  
  if (nrow(all_pathways) > 0) {
    # Create frequency summary
    pathway_frequency <- as.data.frame(table(all_pathways$pathway))
    names(pathway_frequency) <- c("pathway", "frequency")
    pathway_frequency <- pathway_frequency[order(-pathway_frequency$frequency), ]
    
    # Add average significance
    pathway_frequency$avg_mast_pval <- sapply(pathway_frequency$pathway, function(p) {
      pathway_subset <- all_pathways[all_pathways$pathway == p & !is.na(all_pathways$mast_pval), ]
      if (nrow(pathway_subset) > 0) {
        mean(pathway_subset$mast_pval, na.rm = TRUE)
      } else {
        NA
      }
    })
    
    pathway_frequency$avg_crispri_pval <- sapply(pathway_frequency$pathway, function(p) {
      pathway_subset <- all_pathways[all_pathways$pathway == p & !is.na(all_pathways$crispri_pval), ]
      if (nrow(pathway_subset) > 0) {
        mean(pathway_subset$crispri_pval, na.rm = TRUE)
      } else {
        NA
      }
    })
    
    return(pathway_frequency)
  } else {
    return(data.frame(pathway = character(), frequency = integer(), 
                     avg_mast_pval = numeric(), avg_crispri_pval = numeric()))
  }
}

# Create individual signature reports
create_individual_signature_reports <- function(enhanced_signatures, output_dir) {
  
  reports_dir <- file.path(output_dir, "individual_reports")
  if (!dir.exists(reports_dir)) {
    dir.create(reports_dir)
  }
  
  for (i in seq_along(enhanced_signatures)) {
    sig <- enhanced_signatures[[i]]
    
    # Create filename
    filename <- paste0("signature_", i, "_", 
                      gsub("[^A-Za-z0-9]", "_", sig$signature$gene_pair), ".txt")
    filepath <- file.path(reports_dir, filename)
    
    # Write detailed report
    writeLines(sig$interpretation, filepath)
  }
}

# Print key findings to console
print_key_findings <- function(pd_analysis) {
  
  summary_stats <- pd_analysis$pd_summary$summary_stats
  biological_categories <- pd_analysis$pd_summary$biological_categories
  enhanced_signatures <- pd_analysis$enhanced_signatures
  
  cat("ANALYSIS OVERVIEW:\n")
  cat("- Analysis type:", summary_stats$analysis_type, "\n")
  cat("- Signatures analyzed:", summary_stats$total_signatures, "\n") 
  cat("- Mean PD relevance score:", round(summary_stats$mean_pd_relevance, 2), "\n")
  cat("- Most relevant signature:", summary_stats$most_relevant_signature, "\n\n")
  
  cat("TOP BIOLOGICAL CATEGORIES:\n")
  category_ranking <- biological_categories[order(unlist(biological_categories), decreasing = TRUE)]
  top_categories <- head(names(category_ranking)[unlist(category_ranking) > 0], 5)
  
  for (i in seq_along(top_categories)) {
    cat_name <- top_categories[i]
    count <- biological_categories[[cat_name]]
    cat(paste0(i, ". ", str_to_title(gsub("_", " ", cat_name)), ": ", count, " occurrences\n"))
  }
  
  cat("\nTOP 3 SIGNATURES BY PD RELEVANCE:\n")
  pd_scores <- sapply(enhanced_signatures, function(x) x$pd_relevance_score)
  top_indices <- order(pd_scores, decreasing = TRUE)[1:min(3, length(pd_scores))]
  
  for (i in seq_along(top_indices)) {
    idx <- top_indices[i]
    sig <- enhanced_signatures[[idx]]
    cat(paste0(i, ". ", sig$signature$gene_pair, 
              " (Cluster ", sig$signature$cluster, 
              ") - PD Score: ", round(sig$pd_relevance_score, 2), "\n"))
  }
  
  cat("\nBIOLOGICAL INSIGHTS:\n")
  cat(pd_analysis$pd_summary$biological_insights)
}

# Helper function for string manipulation
str_to_title <- function(x) {
  gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x, perl = TRUE)
}

# Main execution
cat("To run the analysis, use:\n")
cat("results <- run_pd_signature_analysis()\n\n")
cat("This will:\n")
cat("1. Load your enrichment data\n")
cat("2. Run signature discovery\n") 
cat("3. Analyze for PD biological relevance\n")
cat("4. Generate interpretable summaries\n")
cat("5. Save results to 'pd_signature_analysis_results/'\n\n")

# Optional: Run immediately if data file exists
if (file.exists("all_enrichment_padj005_complete_with_direction.rds")) {
  cat("Data file found. Running analysis...\n\n")
  results <- run_pd_signature_analysis()
} else {
  cat("Data file not found. Please ensure 'all_enrichment_padj005_complete_with_direction.rds' exists.\n")
}