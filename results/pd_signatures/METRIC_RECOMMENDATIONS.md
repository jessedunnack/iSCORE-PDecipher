# Recommendations for Convergence Strength Metric
**Date**: January 20, 2025

## Current Metric Analysis

### What You're Using Now
- **Formula**: `strength = n_genes × -log10(p-value)`
- **Rationale**: Combines biological relevance (gene count) with statistical significance

### Issues with Current Approach
1. **Non-standard**: This exact metric isn't documented in enrichment analysis literature
2. **Conflates dimensions**: Multiplies two independent aspects that are usually visualized separately
3. **Scale issues**: A pathway with 10 genes and p=0.01 scores the same as 20 genes with p=0.1
4. **Hard to interpret**: What does a "strength" of 50 actually mean biologically?

## Standard Metrics in the Field

### 1. **Gene Ratio** (Recommended)
- **Formula**: `n_genes_in_pathway / total_genes_in_pathway`
- **Range**: 0 to 1 
- **Interpretation**: Proportion of pathway genes that are significant
- **Used by**: clusterProfiler, DAVID, most enrichment tools

### 2. **Enrichment Factor**
- **Formula**: `(n_hits/n_list) / (n_pathway/n_genome)`
- **Interpretation**: Fold enrichment over expected by chance
- **Range**: 0 to infinity (1 = no enrichment)
- **Used by**: GSEA, Enrichr

### 3. **-log10(adjusted p-value)**
- **Standard for significance visualization**
- **Always use adjusted p-values (FDR/Bonferroni)**
- **Higher values = more significant

## Recommended Visualization Approaches

### Option 1: Standard Dot Plot (Best Practice)
```r
# X-axis: Gene Ratio or Enrichment Factor
# Y-axis: Pathway names
# Dot size: Gene count
# Dot color: -log10(adjusted p-value)
```

### Option 2: Scatter Plot Comparing Methods
```r
# X-axis: MAST Gene Ratio
# Y-axis: CRISPRi Gene Ratio
# Size: Average gene count
# Color: Average -log10(p-value)
# Diagonal line: y=x for perfect agreement
```

### Option 3: Composite Score (If You Need One)
```r
# Harmonic mean of normalized metrics:
# score = 2 / (1/gene_ratio + 1/normalized_significance)
# Where normalized_significance = -log10(p) / max(-log10(p))
```

## Implementation Suggestion

Replace the current convergence strength plot with:

```r
# Calculate gene ratios (need total pathway sizes)
conv_strength <- convergent_top %>%
  mutate(
    # Assuming you have pathway sizes
    mast_gene_ratio = n_genes_mast / pathway_size,
    crispri_gene_ratio = n_genes_mixscale / pathway_size,
    avg_significance = mean_neg_log_p,
    total_genes = n_genes_mast + n_genes_mixscale
  )

# Plot with standard metrics
ggplot(conv_strength, aes(x = mast_gene_ratio, y = crispri_gene_ratio)) +
  geom_point(aes(size = total_genes, color = avg_significance)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_color_viridis(name = "-log10(p-value)") +
  labs(
    title = "Convergent Pathway Analysis",
    subtitle = "Comparing gene coverage between methods",
    x = "MAST Gene Ratio",
    y = "CRISPRi Gene Ratio"
  )
```

## Why This Matters

1. **Defensibility**: Standard metrics are accepted in peer review
2. **Interpretability**: Reviewers understand gene ratios immediately
3. **Comparability**: Can compare your results to other studies
4. **Best practices**: Following established visualization guidelines

## Quick Fix vs Full Implementation

### Quick Fix (Keep current, add disclaimer):
Add a note in methods: "Pathway strength calculated as the product of gene count and -log10(p-value) to capture both biological coverage and statistical confidence"

### Better Fix (Recommended):
Switch to gene ratio or enrichment factor visualization, which is standard and more interpretable

### Best Fix:
Calculate proper enrichment factors using the full pathway gene counts and create standard enrichment visualizations

Would you like me to implement one of these standard approaches?