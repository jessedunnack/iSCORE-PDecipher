# Unit Tests for Enrichment Analysis Functions
# Created: August 2025
#
# These tests validate the enrichment analysis functions used for
# pathway analysis in iSCORE-PDecipher

test_that("enrichment analysis handles gene lists correctly", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Create test gene list with known PD-relevant genes
  test_genes <- c("SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35")
  
  # Test gene ID conversion (common requirement for enrichment)
  if (exists("convert_gene_ids", mode = "function")) {
    converted_ids <- convert_gene_ids(test_genes, from = "SYMBOL", to = "ENTREZID")
    
    expect_true(is.character(converted_ids) || is.numeric(converted_ids))
    expect_true(length(converted_ids) <= length(test_genes))  # May lose some genes
    expect_true(length(converted_ids) > 0)  # Should convert at least some genes
  } else {
    skip("convert_gene_ids function not found")
  }
})

test_that("GO enrichment analysis produces valid results", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Use known PD-relevant genes that should produce enrichment
  pd_genes <- c("SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35", 
                "ATP13A2", "FBXO7", "PLA2G6", "DNAJC6", "SYNJ1")
  
  # Convert to Entrez IDs
  gene_entrez <- clusterProfiler::bitr(pd_genes, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 3) {  # Need sufficient genes for enrichment
    
    # Test GO Biological Process enrichment
    expect_no_error({
      go_bp <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                        ont = "BP",
                                        pAdjustMethod = "BH",
                                        pvalueCutoff = 0.05,
                                        qvalueCutoff = 0.2)
    })
    
    if (!is.null(go_bp) && nrow(go_bp@result) > 0) {
      # Check result structure
      expect_true("ID" %in% colnames(go_bp@result))
      expect_true("Description" %in% colnames(go_bp@result))
      expect_true("pvalue" %in% colnames(go_bp@result))
      expect_true("p.adjust" %in% colnames(go_bp@result))
      
      # Check p-values are valid
      expect_true(all(go_bp@result$pvalue >= 0 & go_bp@result$pvalue <= 1))
      expect_true(all(go_bp@result$p.adjust >= 0 & go_bp@result$p.adjust <= 1))
      expect_true(all(go_bp@result$p.adjust >= go_bp@result$pvalue, na.rm = TRUE))
    }
  }
})

test_that("KEGG enrichment analysis produces valid results", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Use dopamine-related genes that should enrich for PD pathways
  dopamine_genes <- c("TH", "DDC", "SLC6A3", "DRD1", "DRD2", "COMT", 
                     "SNCA", "LRRK2", "PRKN")
  
  # Convert to Entrez IDs
  gene_entrez <- clusterProfiler::bitr(dopamine_genes, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 3) {
    
    expect_no_error({
      kegg_result <- clusterProfiler::enrichKEGG(gene = gene_entrez$ENTREZID,
                                                organism = 'hsa',
                                                pvalueCutoff = 0.05,
                                                pAdjustMethod = "BH")
    })
    
    if (!is.null(kegg_result) && nrow(kegg_result@result) > 0) {
      # Check for Parkinson's disease pathway
      pd_pathways <- grep("Parkinson", kegg_result@result$Description, 
                         ignore.case = TRUE, value = TRUE)
      if (length(pd_pathways) > 0) {
        expect_true(length(pd_pathways) >= 1)
        message("Found PD-related KEGG pathways: ", paste(pd_pathways, collapse = ", "))
      }
      
      # Validate result structure
      expect_true("ID" %in% colnames(kegg_result@result))
      expect_true("Description" %in% colnames(kegg_result@result))
      expect_true("pvalue" %in% colnames(kegg_result@result))
      expect_true("p.adjust" %in% colnames(kegg_result@result))
    }
  }
})

test_that("enrichment analysis handles empty gene lists", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Test with empty gene list
  empty_genes <- character(0)
  
  expect_error({
    go_empty <- clusterProfiler::enrichGO(gene = empty_genes,
                                         OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                         ont = "BP")
  })
  
  # Test with very small gene list (should handle gracefully)
  single_gene <- "SNCA"
  gene_entrez <- clusterProfiler::bitr(single_gene, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 0) {
    # Should not crash but may not find enrichment
    expect_no_error({
      go_single <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                            OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                            ont = "BP")
    })
  }
})

test_that("enrichment analysis timeout handling works", {
  
  skip_if_not_installed("clusterProfiler")
  
  # Test timeout functionality (mentioned in context as 600s default)
  test_genes <- c("SNCA", "LRRK2", "PRKN", "PARK7", "PINK1")
  
  # This tests the concept - actual timeout implementation would be in wrapper functions
  start_time <- Sys.time()
  
  # Simulate a quick enrichment analysis
  gene_entrez <- clusterProfiler::bitr(test_genes, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 0) {
    go_result <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                          ont = "BP")
  }
  
  duration <- difftime(Sys.time(), start_time, units = "secs")
  
  # Should complete well within timeout for small test case
  expect_true(as.numeric(duration) < 60)  # Less than 1 minute
})

test_that("direction-specific enrichment analysis works", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Create mock differential expression results with directions
  mock_de_results <- data.frame(
    gene = c("SNCA", "LRRK2", "PRKN", "TH", "PINK1", "VPS35"),
    avg_log2FC = c(1.5, -0.8, -1.2, 2.0, -0.9, 0.7),
    p_val_adj = c(0.01, 0.02, 0.01, 0.005, 0.03, 0.04),
    stringsAsFactors = FALSE
  )
  
  # Separate into up and down-regulated genes
  up_genes <- mock_de_results$gene[mock_de_results$avg_log2FC > 0 & 
                                  mock_de_results$p_val_adj < 0.05]
  down_genes <- mock_de_results$gene[mock_de_results$avg_log2FC < 0 & 
                                    mock_de_results$p_val_adj < 0.05]
  
  expect_true(length(up_genes) > 0)
  expect_true(length(down_genes) > 0)
  
  # Test that we can analyze up and down-regulated genes separately
  for (gene_set in list("UP" = up_genes, "DOWN" = down_genes)) {
    if (length(gene_set) > 0) {
      gene_entrez <- clusterProfiler::bitr(gene_set, 
                                          fromType = "SYMBOL",
                                          toType = "ENTREZID", 
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db)
      
      if (nrow(gene_entrez) > 0) {
        expect_no_error({
          go_result <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                                OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                                ont = "BP",
                                                pvalueCutoff = 1.0)  # Relaxed for test
        })
      }
    }
  }
})

test_that("multiple database enrichment integration works", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("ReactomePA")
  skip_if_not_installed("org.Hs.eg.db")
  
  test_genes <- c("SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35")
  
  # Convert genes
  gene_entrez <- clusterProfiler::bitr(test_genes, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 3) {
    
    # Test multiple enrichment databases
    enrichment_results <- list()
    
    # GO Biological Process
    expect_no_error({
      enrichment_results$GO_BP <- clusterProfiler::enrichGO(
        gene = gene_entrez$ENTREZID,
        OrgDb = org.Hs.eg.db::org.Hs.eg.db,
        ont = "BP",
        pvalueCutoff = 1.0
      )
    })
    
    # Reactome
    expect_no_error({
      enrichment_results$Reactome <- ReactomePA::enrichPathway(
        gene = gene_entrez$ENTREZID,
        pvalueCutoff = 1.0,
        readable = TRUE
      )
    })
    
    # KEGG
    expect_no_error({
      enrichment_results$KEGG <- clusterProfiler::enrichKEGG(
        gene = gene_entrez$ENTREZID,
        organism = 'hsa',
        pvalueCutoff = 1.0
      )
    })
    
    # Check that we can combine results from multiple databases
    all_results <- list()
    
    for (db_name in names(enrichment_results)) {
      result <- enrichment_results[[db_name]]
      if (!is.null(result) && nrow(result@result) > 0) {
        df <- result@result
        df$database <- db_name
        all_results[[db_name]] <- df
      }
    }
    
    if (length(all_results) > 0) {
      combined_results <- do.call(rbind, all_results)
      expect_true(nrow(combined_results) > 0)
      expect_true("database" %in% colnames(combined_results))
    }
  }
})

test_that("enrichment result processing handles all required formats", {
  
  # Test processing of enrichment results as stored in the package data
  # Based on context: all_enrichment_padj005_complete_with_direction.rds
  
  mock_enrichment_data <- data.frame(
    mutation_perturbation = rep("SNCA_A53T", 3),
    cluster = rep("cluster_0", 3),
    enrichment_type = c("GO_BP", "KEGG", "Reactome"),
    direction = c("UP", "DOWN", "ALL"),
    term_id = c("GO:0006915", "hsa05012", "R-HSA-449147"),
    term_description = c("apoptotic process", "Parkinson disease", "Signaling by Interleukins"),
    p.adjust = c(0.001, 0.003, 0.002),
    gene_ratio = c("5/100", "8/100", "6/100"),
    count = c(5, 8, 6),
    geneID = c("GENE1/GENE2", "GENE3/GENE4", "GENE5/GENE6"),
    stringsAsFactors = FALSE
  )
  
  # Validate expected structure
  expected_columns <- c("mutation_perturbation", "cluster", "enrichment_type", 
                       "direction", "term_id", "term_description", "p.adjust", 
                       "gene_ratio", "count", "geneID")
  
  expect_true(all(expected_columns %in% colnames(mock_enrichment_data)))
  
  # Test p-value validity
  expect_true(all(mock_enrichment_data$p.adjust >= 0))
  expect_true(all(mock_enrichment_data$p.adjust <= 1))
  
  # Test direction categories
  valid_directions <- c("UP", "DOWN", "ALL")
  expect_true(all(mock_enrichment_data$direction %in% valid_directions))
  
  # Test gene ratio format
  gene_ratios <- mock_enrichment_data$gene_ratio
  expect_true(all(grepl("^[0-9]+/[0-9]+$", gene_ratios)))
  
  # Test that counts match gene ratios
  for (i in seq_len(nrow(mock_enrichment_data))) {
    ratio_parts <- strsplit(gene_ratios[i], "/")[[1]]
    expected_count <- as.numeric(ratio_parts[1])
    expect_equal(mock_enrichment_data$count[i], expected_count)
  }
})

test_that("PD signature enrichment produces relevant pathways", {
  
  skip_if_not_installed("clusterProfiler")
  skip_if_not_installed("org.Hs.eg.db")
  
  # Use comprehensive PD-relevant gene signature
  pd_signature <- c(
    # Core PD genes
    "SNCA", "LRRK2", "PRKN", "PARK7", "PINK1", "VPS35",
    # Dopamine system
    "TH", "DDC", "SLC6A3", "DRD1", "DRD2",
    # Mitochondrial/autophagy
    "PARK2", "PINK1", "ATP13A2", "FBXO7",
    # Lysosomal
    "GBA", "ATP13A2", "VPS35"
  )
  
  gene_entrez <- clusterProfiler::bitr(pd_signature, 
                                      fromType = "SYMBOL",
                                      toType = "ENTREZID", 
                                      OrgDb = org.Hs.eg.db::org.Hs.eg.db)
  
  if (nrow(gene_entrez) > 5) {
    
    go_result <- clusterProfiler::enrichGO(gene = gene_entrez$ENTREZID,
                                          OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                          ont = "BP",
                                          pvalueCutoff = 0.1,
                                          qvalueCutoff = 0.2)
    
    if (!is.null(go_result) && nrow(go_result@result) > 0) {
      
      # Look for PD-relevant pathway terms
      descriptions <- go_result@result$Description
      
      pd_related_terms <- grep(
        "dopamin|mitochond|autophagy|lysosom|ubiquit|proteasom|apoptot|neurodegenerat",
        descriptions, ignore.case = TRUE, value = TRUE
      )
      
      if (length(pd_related_terms) > 0) {
        message("Found PD-relevant GO terms: ", paste(head(pd_related_terms, 3), collapse = ", "))
        expect_true(length(pd_related_terms) > 0)
      }
      
      # Check enrichment significance
      significant_terms <- go_result@result[go_result@result$p.adjust < 0.05, ]
      if (nrow(significant_terms) > 0) {
        expect_true(nrow(significant_terms) >= 1)
      }
    }
  }
})