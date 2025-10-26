# Data Conversion Requirements
## Pooled MixScale → launch_app() Format

**Date:** October 25, 2025
**Purpose:** Convert pooled MixScale data to work with existing launch_app() infrastructure

---

## Input Structure (Pooled MixScale)

**Location:** `E:/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/`

**File Pattern:** `all_FPD_no_multiplets_noExptSplit_clust_X_mixscale_DEGs.rds`

**Structure:**
```r
cluster_data = list(
  perturbation1 = dataframe(
    gene_ID, log2FC, p_weight, p_weight_BH, p_weight_bonferroni,
    beta_Intercept, beta_weight, beta_log_ct,
    p_Intercept, p_log_ct, DE_method
  ),
  perturbation2 = dataframe(...),
  ...
)
```

**Key Facts:**
- FPD: 7 clusters (0-6), 41 perturbations per cluster
- CRISPRi: 6 clusters (0-5), 338 perturbations per cluster
- Gene names in `gene_ID` column
- Three p-value corrections available: p_weight, p_weight_BH, p_weight_bonferroni

---

## Output Structure (Required by launch_app())

**Location:** `E:/THESIS/scRNASeq/mixscale/FPD_BH_dataset/full_DE_results.rds`

**Structure:**
```r
full_DE_results = list(
  iSCORE_PD_MAST = NULL,  # No MAST data for pooled
  CRISPRi_Mixscale = list(
    perturbation1 = list(
      cluster_0 = list(
        results = dataframe(
          # Rownames: gene symbols
          # Columns: avg_log2FC, p_val_adj, pct.1, pct.2
        ),
        metadata = list(...),  # Optional
        background_genes = character()  # Optional
      ),
      cluster_1 = list(results = ...),
      ...
    ),
    perturbation2 = list(...)
  )
)
```

**Required Columns in results dataframe:**
- **Rownames:** Gene symbols (from gene_ID column)
- **avg_log2FC:** Log2 fold change (renamed from log2FC)
- **p_val_adj:** Adjusted p-value (from p_weight_BH or chosen correction)
- **pct.1:** Percent expressing in group 1 (can be NA or 0)
- **pct.2:** Percent expressing in group 2 (can be NA or 0)

**Optional columns** (keep for compatibility):
- p_val: Original p-value (from p_weight)
- All other columns from pooled data

---

## Conversion Steps

### Step 1: Invert Nesting
**From:** `cluster → perturbation`
**To:** `perturbation → cluster`

```r
# Read all cluster files
# For each perturbation:
#   For each cluster where perturbation appears:
#     Store in new structure
```

### Step 2: Column Transformations

```r
# Move gene_ID to rownames
rownames(results) <- results$gene_ID

# Rename columns
results$avg_log2FC <- results$log2FC

# Map p-value column based on correction type
results$p_val_adj <- results$p_weight_BH  # For BH correction
# OR: results$p_val_adj <- results$p_weight  # For uncorrected
# OR: results$p_val_adj <- results$p_weight_bonferroni  # For Bonferroni

# Add dummy pct columns
results$pct.1 <- NA  # Or 0
results$pct.2 <- NA  # Or 0

# Keep p_val for compatibility
results$p_val <- results$p_weight
```

### Step 3: Wrap in Required Structure

```r
# Create nested structure
full_DE_results <- list(
  iSCORE_PD_MAST = NULL,
  CRISPRi_Mixscale = list()
)

# For each perturbation:
for (pert in perturbations) {
  full_DE_results$CRISPRi_Mixscale[[pert]] <- list()

  # For each cluster:
  for (cluster in clusters) {
    full_DE_results$CRISPRi_Mixscale[[pert]][[paste0("cluster_", cluster)]] <- list(
      results = transformed_dataframe,
      metadata = list(),  # Optional
      background_genes = character()  # Optional
    )
  }
}
```

### Step 4: Save and Copy Enrichment

```r
# Save converted DE results
saveRDS(full_DE_results, "full_DE_results.rds", compress = TRUE)

# Copy enrichment file
file.copy(
  from = "enrichment_results_FPD_p_weight_BH/all_enrichment_padj005_complete_with_direction.rds",
  to = "FPD_BH_dataset/all_enrichment_padj005_complete_with_direction.rds"
)
```

---

## Dataset Mapping

### P-value Correction → Dataset Directory

| Dataset | P-value Column | Source Enrichment |
|---------|---------------|-------------------|
| FPD_BH_dataset | p_weight_BH | enrichment_results_FPD_p_weight_BH/ |
| FPD_uncorrected_dataset | p_weight | enrichment_results_FPD_p_weight/ |
| FPD_bonferroni_dataset | p_weight_bonferroni | enrichment_results_FPD_p_weight_bonferroni/ |
| CRISPRi_BH_dataset | p_weight_BH | enrichment_results_CRISPRi_p_weight_BH/ |
| CRISPRi_uncorrected_dataset | p_weight | enrichment_results_CRISPRi_p_weight/ |
| CRISPRi_bonferroni_dataset | p_weight_bonferroni | enrichment_results_CRISPRi_p_weight_bonferroni/ |

---

## Module Compatibility

**mod_de_results.R** detects format by checking columns:
- If `avg_log2FC` and `p_val_adj` exist → MAST code path ✓
- If `log2FC_*` pattern exists → Experiment-split MixScale code path
- **Our converted data will use MAST code path** ✓

**mod_heatmap.R** expects:
- Nested structure: `full_DE_results$CRISPRi_Mixscale[[gene]][[cluster]]$results` ✓
- Column checking similar to mod_de_results.R ✓

**mod_enrichment_gene_display_v2.R** expects:
- Consolidated enrichment file (already have!) ✓

---

## Validation Checklist

After conversion, verify:
- [ ] File exists: `full_DE_results.rds`
- [ ] File exists: `all_enrichment_padj005_complete_with_direction.rds`
- [ ] Top-level structure: `list(iSCORE_PD_MAST = NULL, CRISPRi_Mixscale = list(...))`
- [ ] At least one perturbation exists
- [ ] At least one cluster exists for that perturbation
- [ ] Results dataframe has rownames (gene symbols)
- [ ] Results dataframe has columns: `avg_log2FC`, `p_val_adj`
- [ ] Can launch: `launch_app(data_dir = "path/to/dataset")`
- [ ] Volcano plots render without errors
- [ ] Heatmaps work correctly
- [ ] Enrichment displays correctly

---

## Next Steps

1. Create conversion script: `convert_pooled_to_full_de.R`
2. Test with one dataset: FPD_BH_dataset
3. Validate with launch_app()
4. Scale to all 6 datasets
5. Update get_dataset_options()
6. Final testing and documentation
