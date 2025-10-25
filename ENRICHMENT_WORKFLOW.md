# Enrichment Analysis Workflow for iSCORE-PDecipher

**Date:** October 24, 2025
**Purpose:** Guide for generating enrichment data for the v0.5.0 Perturb-seq module

---

## Overview

The iSCORE-PDecipher package v0.5.0 supports pooled MixScale data with three p-value correction options. To enable enrichment visualization in the Shiny app, you need to run enrichment analysis on the MixScale DE results.

---

## 📁 Required Input Files

### MixScale DE Results (Already Generated ✅)

**FPD Dataset:**
- Location: `../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/`
- Files: 7 RDS files (one per cluster)
- Structure: Pooled MixScale with p_weight, p_weight_BH, p_weight_bonferroni columns

**CRISPRi Dataset:**
- Location: `../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/`
- Files: 6 RDS files (one per cluster)
- Structure: Pooled MixScale with p_weight, p_weight_BH, p_weight_bonferroni columns

---

## 🔧 Enrichment Analysis Script

### Main Script Location

**File:** `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/enrichment_analysis_simple.R`

**Purpose:** Run functional enrichment analysis on pooled MixScale DE results

**Main Function:** `run_mixscale_enrichment_simple()`

### Function Parameters

```r
run_mixscale_enrichment_simple(
  mixscale_dir,        # Directory containing *_mixscale_DEGs.rds files
  output_dir,          # Where to save enrichment results
  pval_column = "p_weight",  # Which p-value column to use
  lfc_threshold = 0.25,      # Log2FC threshold for DEG selection
  padj_threshold = 0.05,     # P-value threshold for significance
  min_genes = 5,             # Minimum genes required for enrichment
  run_methods = c("GO", "KEGG", "Reactome", "WikiPathways", "STRING", "GSEA"),
  padj_method = "BH"         # Method for p-value adjustment
)
```

### Supported Enrichment Methods

1. **GO** - Gene Ontology (BP, CC, MF, ALL)
2. **KEGG** - KEGG pathways
3. **Reactome** - Reactome pathways
4. **WikiPathways** - WikiPathways
5. **STRING** - Protein-protein interaction enrichment
6. **GSEA** - Gene Set Enrichment Analysis

---

## 📊 How to Run Enrichment Analysis

### Step 1: Set Up Environment

```bash
conda activate seuratv4
cd /mnt/e/ASAP/scRNASeq/PerturbSeq/final
R
```

### Step 2: Load the Script

```r
source("enrichment_analysis_simple.R")
```

### Step 3: Run for Each Dataset and P-value Correction

#### **FPD with Uncorrected P-values**

```r
run_mixscale_enrichment_simple(
  mixscale_dir = "final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/",
  output_dir = "enrichment_results_FPD_p_weight",
  pval_column = "p_weight"
)
```

#### **FPD with BH Correction (Recommended)**

```r
run_mixscale_enrichment_simple(
  mixscale_dir = "final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/",
  output_dir = "enrichment_results_FPD_p_weight_BH",
  pval_column = "p_weight_BH"
)
```

#### **FPD with Bonferroni Correction**

```r
run_mixscale_enrichment_simple(
  mixscale_dir = "final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/",
  output_dir = "enrichment_results_FPD_p_weight_bonferroni",
  pval_column = "p_weight_bonferroni"
)
```

#### **CRISPRi with BH Correction**

```r
run_mixscale_enrichment_simple(
  mixscale_dir = "final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/",
  output_dir = "enrichment_results_CRISPRi_p_weight_BH",
  pval_column = "p_weight_BH"
)
```

**Repeat for uncorrected and Bonferroni for CRISPRi...**

---

## 📂 Expected Output Directory Structure

```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/

├── enrichment_results_FPD_p_weight/
│   ├── Cluster0/
│   │   ├── LRRK2/
│   │   │   ├── GO_enrichment.rds
│   │   │   ├── KEGG_enrichment.rds
│   │   │   ├── Reactome_enrichment.rds
│   │   │   ├── WikiPathways_enrichment.rds
│   │   │   ├── STRING_enrichment.rds
│   │   │   └── GSEA_enrichment.rds
│   │   ├── SNCA_A53T/
│   │   └── ...
│   ├── Cluster1/
│   └── ...
│   └── enrichment_summary_log.csv
│
├── enrichment_results_FPD_p_weight_BH/
│   └── (same structure)
│
├── enrichment_results_FPD_p_weight_bonferroni/
│   └── (same structure)
│
├── enrichment_results_CRISPRi_p_weight/
│   └── (same structure, 6 clusters)
│
├── enrichment_results_CRISPRi_p_weight_BH/
│   └── (same structure)
│
└── enrichment_results_CRISPRi_p_weight_bonferroni/
    └── (same structure)
```

---

## ⏱️ Estimated Runtime

**Per dataset:**
- Small clusters (FPD): ~30-60 minutes per p-value correction
- Large clusters (CRISPRi): ~2-4 hours per p-value correction

**Total for all 6 combinations:** ~10-15 hours

**Recommendation:** Run overnight or use HPC if available

---

## 🔍 Monitoring Progress

The script provides detailed progress messages:

```
===============================================
  STARTING MIXSCALE ENRICHMENT ANALYSIS
===============================================

Found 7 MixScale RDS files to process

-----------------------------------------------
  Processing file: all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds
-----------------------------------------------
Cluster: Cluster0

  • Processing perturbation: LRRK2

    Running GO enrichment with ALL ontology
    Gene list: 1234 genes, background: 6347 genes
    ✓ GO enrichment completed (123 terms)

    Running KEGG enrichment...
    ✓ KEGG enrichment completed (45 terms)

    ...
```

---

## 🚀 After Enrichment Completes

### Verify Output

```r
# Check that output directories exist
list.dirs("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/enrichment_results_FPD_p_weight_BH",
          recursive = FALSE)

# Load summary log
summary <- read.csv("enrichment_results_FPD_p_weight_BH/enrichment_summary_log.csv")
table(summary$enrichment_run)  # Should show TRUE for successful runs
```

### Update Shiny App Environment Variables (Optional)

If you want the app to auto-load enrichment:

```r
# In your .Renviron or at startup:
Sys.setenv(ISCORE_ENRICHMENT_BASE = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final")
```

### Launch Perturb-seq Module

```r
library(iSCORE.PDecipher)
run_app(app_type = "perturbseq")
```

The module will automatically detect and load available enrichment results!

---

## 📝 Key Package Functions

### From `R/import_pooled_mixscale_functions.R`

```r
# Import enrichment with specific correction
import_enrichment_with_correction(
  base_dir = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final",
  dataset = "FPD",           # or "CRISPRi"
  pval_correction = "BH"     # "none", "BH", or "bonferroni"
)
```

### From `inst/shiny/R/data_manager.R`

```r
# Get pooled enrichment data (used by Shiny modules)
get_pooled_enrichment_data(
  dataset = "FPD",
  pval_correction = "BH",
  force_reload = FALSE
)
```

---

## 🐛 Troubleshooting

### Issue: "No *_mixscale_DEGs.rds files found"

**Solution:** Check that `mixscale_dir` path is correct and contains RDS files

### Issue: Enrichment timeouts

**Solution:** The script has 10-minute timeout per method. If it times out frequently:
- Increase `timeout_seconds` parameter in `run_with_timeout()` calls
- Or run fewer methods at once

### Issue: Memory errors with large datasets

**Solution:** The script uses `gc()` after each cluster. If still problematic:
- Run one cluster at a time
- Increase available RAM
- Process in batches

---

## 📧 Support

For issues or questions:
1. Check the enrichment summary log: `enrichment_summary_log.csv`
2. Review console output for error messages
3. Verify input RDS files have required columns (gene_ID, log2FC, p_weight variants)

---

**Created:** October 24, 2025
**For:** iSCORE-PDecipher v0.5.0
**Script:** enrichment_analysis_simple.R
