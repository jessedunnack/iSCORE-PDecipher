# Debug Gene Display Issue
# Check why the Associated_Genes column isn't appearing

cat("=== Debugging Gene Display Issue ===\n")

# 1. Test if gene associations can be loaded
cat("1. Testing gene association loading...\n")
tryCatch({
  source("R/gene_association_lookup.R")
  result <- load_gene_associations()
  if (!is.null(result)) {
    cat("   ✅ Gene associations loaded:", nrow(result), "associations\n")
    available <- gene_associations_available()
    cat("   ✅ Available check:", available, "\n")
  } else {
    cat("   ❌ Failed to load gene associations\n")
  }
}, error = function(e) {
  cat("   ❌ Error loading gene associations:", e$message, "\n")
})

# 2. Test a specific lookup
cat("\n2. Testing specific gene lookup...\n")
if (exists("gene_associations_available") && gene_associations_available()) {
  # Load sample data to get real parameters
  data <- readRDS("inst/extdata/gene_term_associations.rds")
  sample_row <- data[1, ]
  
  cat("   Sample term:", sample_row$term_id, "\n")
  cat("   Sample gene:", sample_row$gene, "\n")
  cat("   Sample cluster:", sample_row$cluster, "\n")
  
  lookup_result <- get_genes_for_term(
    term_id = sample_row$term_id,
    analysis_type = sample_row$analysis_type,
    gene = sample_row$gene,
    cluster = sample_row$cluster,
    enrichment_type = sample_row$enrichment_type,
    direction = sample_row$direction,
    experiment = "default"
  )
  
  if (!is.null(lookup_result$genes)) {
    cat("   ✅ Lookup successful:", length(lookup_result$genes), "genes\n")
    cat("   Genes:", paste(head(lookup_result$genes, 5), collapse = ", "), "\n")
  } else {
    cat("   ❌ Lookup failed:", lookup_result$error, "\n")
  }
}

# 3. Check file paths from Shiny app perspective
cat("\n3. Checking file paths from Shiny app perspective...\n")
shiny_dir <- "inst/shiny"
cat("   Current working directory:", getwd(), "\n")
cat("   Shiny directory exists:", dir.exists(shiny_dir), "\n")

# From Shiny app, the paths would be:
gene_lookup_path1 <- file.path(shiny_dir, "../../R/gene_association_lookup.R")
gene_lookup_path2 <- "R/gene_association_lookup.R"
gene_data_path1 <- file.path(shiny_dir, "../../inst/extdata/gene_term_associations.rds")
gene_data_path2 <- "inst/extdata/gene_term_associations.rds"

cat("   Gene lookup path 1:", gene_lookup_path1, "exists:", file.exists(gene_lookup_path1), "\n")
cat("   Gene lookup path 2:", gene_lookup_path2, "exists:", file.exists(gene_lookup_path2), "\n")
cat("   Gene data path 1:", gene_data_path1, "exists:", file.exists(gene_data_path1), "\n")
cat("   Gene data path 2:", gene_data_path2, "exists:", file.exists(gene_data_path2), "\n")

# 4. Test the visualization module logic
cat("\n4. Testing visualization module logic simulation...\n")

# Simulate what happens in the renderDataTable
simulate_table_rendering <- function() {
  # Load gene associations if available
  gene_data_available <- FALSE
  tryCatch({
    if (!exists(".gene_associations")) {
      source("R/gene_association_lookup.R")
      load_gene_associations()
    }
    gene_data_available <- gene_associations_available()
  }, error = function(e) {
    cat("     Error in table rendering:", e$message, "\n")
  })
  
  cat("   Gene data available in simulation:", gene_data_available, "\n")
  
  if (gene_data_available) {
    # Test the get_gene_list function
    get_gene_list <- function(term_id) {
      result <- get_genes_for_term(
        term_id = term_id,
        analysis_type = "MAST",
        gene = "LRRK2", 
        cluster = "cluster_0",
        enrichment_type = "GO_BP",
        direction = "UP",
        experiment = "default"
      )
      
      if (!is.null(result$genes) && length(result$genes) > 0) {
        genes_display <- head(result$genes, 20)
        if (length(result$genes) > 20) {
          genes_display <- c(genes_display, paste("...", length(result$genes) - 20, "more"))
        }
        return(paste(genes_display, collapse = ", "))
      } else {
        return("No genes found")
      }
    }
    
    # Test with a real term
    test_terms <- head(unique(data$term_id), 3)
    for (term in test_terms) {
      gene_list <- get_gene_list(term)
      cat("     Term", term, "->", substr(gene_list, 1, 80), "...\n")
    }
  }
}

simulate_table_rendering()

cat("\n=== Debugging Complete ===\n")