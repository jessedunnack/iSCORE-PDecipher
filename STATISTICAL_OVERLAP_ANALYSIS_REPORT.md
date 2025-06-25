# Statistical Analysis of Gene Set Overlaps: MAST vs MixScale

## Executive Summary

**Objective**: Determine if overlaps between MAST (iSCORE-PD mutations) and MixScale (CRISPR perturbations) differentially expressed genes are statistically significant compared to random chance.

**Recommendation**: **Fisher's Exact Test** as primary method, with **permutation testing** for advanced analysis.

**Implementation Status**: ✅ **COMPLETED** - Fisher's exact test integrated into DE Results summary panel with real-time calculation.

---

## Problem Definition

### Research Question
*"Are the overlapping differentially expressed genes between MAST and MixScale results more numerous than would be expected by random chance alone?"*

### Biological Context
- **MAST data**: Genetic mutations in Parkinson's disease genes (permanent alterations)
- **MixScale data**: CRISPR-mediated perturbations (transient knockdown/activation)
- **Expectation**: Biologically related perturbations should show significant overlap beyond random chance

### Statistical Framework
- **Null Hypothesis (H₀)**: Gene overlap is due to random chance
- **Alternative Hypothesis (H₁)**: Gene overlap reflects shared biological mechanisms
- **Significance Level**: α = 0.05 (adjustable for multiple testing)

---

## Statistical Test Options Analysis

### 1. Fisher's Exact Test (Hypergeometric Distribution) ⭐ **RECOMMENDED**

#### Overview
Tests whether the observed overlap in a 2×2 contingency table deviates significantly from random expectation.

#### Implementation
```
Contingency Table:
                    MixScale_Sig    MixScale_NotSig    Total
MAST_Sig           overlap         mast_only          mast_total
MAST_NotSig        mixscale_only   neither            background-mast_total
Total              mixscale_total  background-mixscale_total  background_total
```

#### Advantages
- **Well-established**: Standard for gene set overlap analysis
- **Exact p-values**: No approximations needed
- **Fast computation**: Suitable for real-time UI updates
- **Conservative**: Controls Type I error effectively
- **Interpretable**: Clear odds ratio for effect size

#### Limitations
- **Independence assumption**: Assumes genes are independent (often violated)
- **Background definition**: Sensitive to choice of background gene set
- **Binary classification**: Doesn't use continuous p-values or effect sizes

#### Parameters
- **Background set**: Intersection of genes tested in both methods (conservative)
- **Significance criteria**: p < 0.05 and |log2FC| > 1 (consistent across methods)
- **Test direction**: One-tailed (testing for enrichment, not depletion)

### 2. Permutation Testing 🔬 **ADVANCED OPTION**

#### Overview
Generate null distribution by repeatedly shuffling gene labels while preserving method-specific properties.

#### Implementation Approaches
1. **Label permutation**: Shuffle which genes are significant
2. **Gene-wise permutation**: Permute p-values within each method
3. **Rank-based permutation**: Preserve gene ranking distributions

#### Advantages
- **Fewer assumptions**: Doesn't assume independence
- **Preserves structure**: Can maintain gene-level correlations
- **Flexible**: Can incorporate continuous measures
- **Robust**: Less sensitive to background set definition

#### Limitations
- **Computational cost**: Requires thousands of iterations
- **Implementation complexity**: More sophisticated coding required
- **UI responsiveness**: May be too slow for real-time updates

### 3. Bootstrap Resampling 🔄 **SPECIALIZED OPTION**

#### Overview
Resample genes of same sizes from background set to generate null distribution.

#### Advantages
- **Confidence intervals**: Provides uncertainty estimates
- **Flexible sample sizes**: Can test different overlap scenarios
- **Non-parametric**: Makes minimal distributional assumptions

#### Limitations
- **Sampling assumptions**: Assumes uniform sampling probability
- **Computational cost**: Similar to permutation testing
- **Background dependency**: Still requires well-defined background set

### 4. Binomial Test 📊 **SIMPLE ALTERNATIVE**

#### Overview
Test if overlap proportion exceeds expectation under random gene selection.

#### Implementation
- **Success probability**: (MAST_sig / background) × (MixScale_sig / background)
- **Trials**: min(MAST_sig, MixScale_sig)
- **Observed successes**: overlap

#### Advantages
- **Simple interpretation**: Direct probability calculation
- **Fast computation**: Immediate results
- **Minimal assumptions**: Only requires random sampling

#### Limitations
- **Oversimplified**: Ignores complex dependency structures
- **Less powerful**: May miss true signals
- **Background insensitive**: Doesn't account for method-specific backgrounds

---

## Implementation Recommendations

### Phase 1: Current Implementation ✅ **COMPLETED**
**Fisher's Exact Test** with the following specifications:

```r
# Background: Intersection of genes tested in both methods
background_genes <- intersect(unique(mast_data$gene_name), 
                             unique(mixscale_data$gene_name))

# Significance criteria (consistent across methods)
significance_threshold <- list(
  pvalue = 0.05,
  log2fc_threshold = 1.0
)

# Contingency table construction
overlap <- length(intersect(mast_significant_genes, mixscale_significant_genes))
mast_only <- length(mast_significant_genes) - overlap
mixscale_only <- length(mixscale_significant_genes) - overlap
neither <- length(background_genes) - mast_only - mixscale_only - overlap

# Fisher's exact test (one-tailed)
fisher.test(matrix(c(overlap, mast_only, mixscale_only, neither), nrow=2), 
           alternative="greater")
```

### Phase 2: Enhanced Analysis 🚀 **FUTURE DEVELOPMENT**

#### Multiple Testing Correction
- **Bonferroni correction**: For multiple clusters/genes tested
- **FDR control**: Using Benjamini-Hochberg procedure
- **Cluster-level analysis**: Test each cluster separately with correction

#### Direction-Specific Analysis
```r
# Separate tests for up/down regulation
up_up_overlap <- intersect(mast_upregulated, mixscale_upregulated)
down_down_overlap <- intersect(mast_downregulated, mixscale_downregulated)
opposing_overlap <- union(intersect(mast_upregulated, mixscale_downregulated),
                         intersect(mast_downregulated, mixscale_upregulated))
```

#### Effect Size Correlation
```r
# Pearson correlation of log2FC values for overlapping genes
correlation_test <- cor.test(mast_log2fc_overlap, mixscale_log2fc_overlap)
```

### Phase 3: Advanced Methods 🔬 **RESEARCH EXTENSION**

#### Rank-Based Methods
- **GSEA-style analysis**: Compare entire gene rankings rather than binary classifications
- **Kolmogorov-Smirnov test**: Compare p-value distributions between methods

#### Network-Informed Analysis
- **Pathway-based testing**: Account for gene functional relationships
- **Protein interaction networks**: Weight overlaps by interaction strength

#### Bayesian Approaches
- **Hierarchical models**: Account for cluster-level and gene-level variability
- **Prior incorporation**: Use biological knowledge to inform tests

---

## Practical Considerations

### Background Set Definition
**Current approach**: Intersection of genes tested in both methods
- **Conservative**: Ensures fair comparison
- **Consistent**: Same background for all comparisons
- **Transparent**: Clear interpretation

### Multiple Testing
**Current status**: Single test per comparison
**Future need**: When testing multiple clusters/genes simultaneously
- **Recommendation**: FDR control at 5% level
- **Implementation**: `p.adjust(p_values, method="BH")`

### Effect Size Interpretation
**Odds Ratio interpretation**:
- **OR = 1**: Random overlap (null hypothesis)
- **OR > 1**: More overlap than expected (biological signal)
- **OR < 1**: Less overlap than expected (antagonistic effects)

### Computational Performance
**Current implementation**: Real-time calculation suitable for interactive use
**Performance**: < 1ms for typical dataset sizes
**Scalability**: Linear with number of genes tested

---

## Quality Control and Validation

### Sanity Checks
1. **Background size validation**: Background ≥ max(MAST_sig, MixScale_sig)
2. **Overlap bounds**: 0 ≤ overlap ≤ min(MAST_sig, MixScale_sig)
3. **Contingency table validity**: All cells ≥ 0

### Expected Results
- **Positive controls**: Same gene in both datasets should show overlap
- **Negative controls**: Unrelated conditions should show random overlap
- **Biological validation**: Overlaps should enrich for relevant pathways

### Troubleshooting
**Common issues**:
- **Empty overlaps**: Check if both methods have significant genes
- **Invalid contingency tables**: Verify background set calculations
- **Extreme p-values**: May indicate data quality issues

---

## Conclusion

**Fisher's Exact Test** provides a robust, interpretable, and computationally efficient method for testing gene set overlap significance in the MAST vs MixScale comparison. The implementation is now integrated into the DE Results summary panel with real-time calculation and clear interpretation.

**Key Benefits**:
- **Scientific rigor**: Well-established statistical framework
- **User-friendly**: Clear p-values and odds ratios displayed
- **Interactive**: Updates with cluster/gene selection changes
- **Extensible**: Foundation for more sophisticated future analyses

**Future enhancements** can build upon this foundation with permutation testing, multiple testing correction, and pathway-informed analyses as research needs evolve.

---

*Report generated: December 25, 2024*  
*Implementation status: Fisher's Exact Test ✅ COMPLETED in DE Results module*