# Core Overlap & Heatmap Functions - iSCORE-PDecipher
**Date:** October 26, 2025
**Purpose:** Identify key functions for DE gene and enrichment overlap analysis with heatmap visualization

---

## 🎯 Overview

This document identifies the core functions for:
1. **DE Gene Overlap Analysis**: Find shared differentially expressed genes across perturbations
2. **Enrichment Overlap Analysis**: Identify common pathways/terms across perturbations
3. **Heatmap Visualization**: Binary presence/absence OR expression value heatmaps

---

## 📊 Part 1: DE Gene Overlap Analysis

### Primary Functions (R/signature_analysis.R)

#### 1. **`calculate_gene_overlap_significance()`** (Line 214)
**Purpose:** Calculate overlap between two gene sets with Fisher's exact test

**Key Features:**
- Fisher's exact test for significance
- Jaccard index for similarity
- Overlap coefficient
- Lists overlapping genes

**Usage:**
```r
overlap_stats <- calculate_gene_overlap_significance(
  mast_genes = c("LRRK2", "ATP13A2", "SNCA"),
  crispri_genes = c("LRRK2", "SNCA", "PARK7"),
  background_genes = all_detected_genes,
  alternative = "greater"
)

# Returns:
# $overlap_count      # Number of shared genes
# $overlap_genes      # Names of shared genes
# $fisher_p          # Fisher's exact test p-value
# $jaccard_index     # Jaccard similarity
# $overlap_coefficient
```

#### 2. **`calculate_gene_overlap_significance_proper()`** (Line 15)
**Purpose:** Enhanced version with better statistical handling

**Key Features:**
- Improved p-value calculations
- Handles edge cases
- More robust for small gene sets

#### 3. **`analyze_gene_pair_signatures()`** (Line 609)
**Purpose:** Complete signature analysis for a mutation-perturbation pair

**Key Features:**
- Gene overlap analysis
- Pathway overlap analysis
- Effect size correlations
- Direction consistency checks
- Composite signature scores

**Usage:**
```r
signature <- analyze_gene_pair_signatures(
  gene_pair = c("LRRK2_mutation", "LRRK2_perturbation"),
  enrichment_data = enrichment_results,
  de_data = full_DE_results,
  clusters = c("cluster_0", "cluster_1"),
  padj_threshold = 0.05,
  lfc_threshold = 0.25
)
```

#### 4. **`calculate_effect_size_correlation()`** (Line 292)
**Purpose:** Correlation of log2FC values between two datasets

**Usage:**
```r
correlation <- calculate_effect_size_correlation(
  mast_data = data.frame(gene = genes, log2FC = fc1),
  crispri_data = data.frame(gene = genes, log2FC = fc2),
  method = "pearson"  # or "spearman"
)
```

#### 5. **`calculate_direction_consistency()`** (Line 372)
**Purpose:** Check if DE direction (up/down) is consistent

**Returns:**
```r
# $consistent_direction  # Genes with same direction
# $inconsistent_direction # Genes with opposite direction
# $consistency_ratio     # Proportion consistent
```

---

### Comprehensive Overlap Analysis (R/comprehensive_fishers_analysis.R)

#### 6. **`run_comprehensive_fishers_analysis()`** (Line 24)
**Purpose:** Complete pairwise overlap analysis across ALL mutations and perturbations

**Key Features:**
- Compares every mutation vs every perturbation
- Across all clusters
- Fisher's exact test for each comparison
- Exports to CSV for easy analysis

**Usage:**
```r
results <- run_comprehensive_fishers_analysis(
  de_data_path = "full_DE_results.rds",
  output_dir = "fisher_overlap_results",
  padj_threshold = 0.05,
  lfc_threshold = 0.25,
  background_genes = 15000
)

# Creates output files:
# - comprehensive_fishers_results.csv
# - summary_statistics.csv
# - significant_overlaps_only.csv
```

**Output Structure:**
```
mutation,perturbation,cluster,
overlap_count,fisher_p,jaccard_index,
mast_deg_count,crispri_deg_count,
overlapping_genes
```

---

## 🔬 Part 2: Enrichment Term Overlap Analysis

### Primary Functions

#### 7. **`calculate_pathway_overlap_by_database()`** (Line 548, R/signature_analysis.R)
**Purpose:** Calculate pathway/term overlap for specific databases (GO, KEGG, Reactome, etc.)

**Key Features:**
- Database-specific analysis (GO_BP, GO_MF, GO_CC, KEGG, Reactome, WikiPathways)
- Fisher's exact test per database
- Jaccard similarity per database

**Usage:**
```r
pathway_overlap <- calculate_pathway_overlap_by_database(
  mast_data = mast_enrichment,  # Filtered to one gene/cluster
  crispri_data = crispri_enrichment,
  enrichment_databases = c("GO_BP", "GO_MF", "KEGG", "Reactome"),
  term_id_col = "term_id",
  pvalue_col = "p.adjust",
  pvalue_threshold = 0.05
)

# Returns list with one element per database:
# $GO_BP$overlap_count
# $GO_BP$fisher_p
# $GO_BP$jaccard_index
# $GO_BP$overlapping_pathways (term IDs)
```

#### 8. **`identify_pd_relevant_enrichments()`** (Line 509, R/signature_analysis.R)
**Purpose:** Filter enrichment results for PD-relevant pathways

**Key Features:**
- Filters for dopamine/mitochondria/lysosome/synapse pathways
- Customizable PD term list
- Priority scoring

**Usage:**
```r
pd_terms <- identify_pd_relevant_enrichments(
  enrichment_data = all_enrichment_results,
  pd_pathway_terms = c("dopamine", "mitochondria", "lysosome",
                       "autophagy", "synapse", "parkinson")
)
```

---

## 📈 Part 3: Heatmap Visualization Functions

### Interactive Heatmaps (R/signature_visualization_functions.R)

#### 9. **`create_interactive_signature_heatmap()`** (Line 112)
**Purpose:** Create interactive heatmap of signature metrics

**Key Features:**
- Plotly-based interactivity
- Multiple metrics (gene overlap, pathway overlap, correlations)
- Color-coded by significance
- Hover details

**Usage:**
```r
heatmap <- create_interactive_signature_heatmap(
  signature_data = all_signatures,
  metric = "gene_overlap",  # or "pathway_overlap", "correlation"
  significance_threshold = 0.05,
  title = "DE Gene Overlap Across Perturbations"
)

# Display
heatmap
```

**Available Metrics:**
- `"gene_overlap"`: -log10(p-value) from Fisher's test
- `"pathway_overlap"`: -log10(p-value) for pathway overlap
- `"correlation"`: Pearson/Spearman correlation
- `"jaccard"`: Jaccard similarity index

#### 10. **`create_interactive_signature_heatmap_enhanced()`** (Line 192)
**Purpose:** Enhanced version with clustering and better visualization

**Key Features:**
- Hierarchical clustering
- Dendrograms
- Better color schemes
- Annotation tracks

---

## 🔨 How to Use These Functions for Your Goals

### Goal 1: DE Gene Overlap Heatmap (Binary Presence/Absence)

**Step-by-step workflow:**

```r
library(iSCORE.PDecipher)

# 1. Load your converted pooled data
full_de <- readRDS("E:/THESIS/scRNASeq/mixscale/FPD_BH_dataset/full_DE_results.rds")

# 2. Extract perturbation data
mixscale_data <- full_de$CRISPRi_Mixscale

# 3. For each perturbation-cluster combination, get significant DEGs
padj_threshold <- 0.05
lfc_threshold <- 0.25

# Create matrix: rows = genes, columns = perturbations
# Value = 1 if DE, 0 if not DE

perturbations <- names(mixscale_data)
clusters <- names(mixscale_data[[1]])

for (cluster in clusters) {
  # Get all genes across all perturbations in this cluster
  all_genes <- unique(unlist(lapply(perturbations, function(p) {
    results <- mixscale_data[[p]][[cluster]]$results
    rownames(results)[results$p_val_adj < padj_threshold &
                     abs(results$avg_log2FC) > lfc_threshold]
  })))

  # Create binary matrix
  binary_matrix <- matrix(0, nrow = length(all_genes), ncol = length(perturbations))
  rownames(binary_matrix) <- all_genes
  colnames(binary_matrix) <- perturbations

  for (i in seq_along(perturbations)) {
    pert <- perturbations[i]
    results <- mixscale_data[[pert]][[cluster]]$results
    sig_genes <- rownames(results)[results$p_val_adj < padj_threshold &
                                   abs(results$avg_log2FC) > lfc_threshold]
    binary_matrix[sig_genes, i] <- 1
  }

  # Create heatmap
  heatmaply::heatmaply(
    binary_matrix,
    main = paste("DE Gene Presence/Absence -", cluster),
    colors = c("white", "darkblue"),
    scale = "none",
    dendrogram = "both"
  )
}
```

### Goal 2: DE Gene Expression Value Heatmap

```r
# Same setup as above, but fill with log2FC values instead of binary

for (cluster in clusters) {
  # Get significant genes
  all_sig_genes <- unique(unlist(lapply(perturbations, function(p) {
    results <- mixscale_data[[p]][[cluster]]$results
    rownames(results)[results$p_val_adj < padj_threshold &
                     abs(results$avg_log2FC) > lfc_threshold]
  })))

  # Create expression matrix (log2FC values)
  expr_matrix <- matrix(NA, nrow = length(all_sig_genes), ncol = length(perturbations))
  rownames(expr_matrix) <- all_sig_genes
  colnames(expr_matrix) <- perturbations

  for (i in seq_along(perturbations)) {
    pert <- perturbations[i]
    results <- mixscale_data[[pert]][[cluster]]$results

    # Fill with log2FC for significant genes (NA for non-significant)
    for (gene in all_sig_genes) {
      if (gene %in% rownames(results)) {
        is_sig <- results[gene, "p_val_adj"] < padj_threshold &
                 abs(results[gene, "avg_log2FC"]) > lfc_threshold
        if (is_sig) {
          expr_matrix[gene, i] <- results[gene, "avg_log2FC"]
        }
      }
    }
  }

  # Create heatmap (centered at 0, red = up, blue = down)
  heatmaply::heatmaply(
    expr_matrix,
    main = paste("DE Gene Expression (log2FC) -", cluster),
    colors = colorRampPalette(c("blue", "white", "red"))(100),
    scale = "none",
    dendrogram = "both",
    na.value = "grey90"
  )
}
```

### Goal 3: Enrichment Term Overlap Analysis

```r
# 1. Load enrichment data
enrichment <- readRDS("E:/THESIS/scRNASeq/mixscale/FPD_BH_dataset/all_enrichment_padj005_complete_with_direction.rds")

# 2. For each database (GO_BP, KEGG, Reactome, etc.)
databases <- c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome", "WikiPathways")

for (db in databases) {
  # Filter to this database
  db_enrichment <- enrichment[enrichment$enrichment_type == db, ]

  # Create binary matrix: rows = terms, columns = perturbations
  perturbations <- unique(db_enrichment$mutation_perturbation)
  terms <- unique(db_enrichment$term_id)

  binary_matrix <- matrix(0, nrow = length(terms), ncol = length(perturbations))
  rownames(binary_matrix) <- terms
  colnames(binary_matrix) <- perturbations

  for (i in seq_along(perturbations)) {
    pert <- perturbations[i]
    sig_terms <- db_enrichment$term_id[db_enrichment$mutation_perturbation == pert &
                                       db_enrichment$p.adjust < 0.05]
    binary_matrix[sig_terms, i] <- 1
  }

  # Create heatmap
  heatmaply::heatmaply(
    binary_matrix,
    main = paste(db, "Term Overlap Across Perturbations"),
    colors = c("white", "darkgreen"),
    scale = "none",
    dendrogram = "both"
  )
}
```

---

## 📝 Key Files to Examine

1. **R/signature_analysis.R** - Core overlap statistics functions
2. **R/comprehensive_fishers_analysis.R** - Automated pairwise overlap analysis
3. **R/signature_visualization_functions.R** - Heatmap and visualization functions
4. **R/manuscript_signature_discovery.R** - High-level signature discovery workflows

---

## 🚀 Quick Start: Run Pre-built Comprehensive Analysis

```r
# If you just want to run the comprehensive overlap analysis:

library(iSCORE.PDecipher)

# Option 1: For your pooled data
results <- run_comprehensive_fishers_analysis(
  de_data_path = "E:/THESIS/scRNASeq/mixscale/FPD_BH_dataset/full_DE_results.rds",
  output_dir = "FPD_BH_overlap_analysis",
  padj_threshold = 0.05,
  lfc_threshold = 0.25
)

# This creates CSV files with all pairwise overlaps!
# Then visualize with:
overlap_data <- read.csv("FPD_BH_overlap_analysis/comprehensive_fishers_results.csv")
library(heatmaply)
# ... create heatmaps from overlap_data
```

---

## 💡 Notes

- All functions handle missing data gracefully (return NA or empty results)
- Fisher's exact test is the gold standard for overlap significance
- Jaccard index is good for similarity when set sizes differ greatly
- For large gene sets, consider using `calculate_gene_overlap_significance_proper()` instead of base version
- Heatmaps work best with 10-100 items per dimension (genes/perturbations)

---

**END OF DOCUMENT**
