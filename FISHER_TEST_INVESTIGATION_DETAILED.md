# Deep Investigation: Fisher's Exact Test Discrepancy Analysis
**Date:** January 16, 2025  
**Status:** In Progress - Systematic Investigation

## Executive Summary

Investigating the substantial p-value discrepancy between Fisher's exact tests in:
1. **Signature Nomination Module** (Gene-pair analysis page)
2. **DE Results Module** (DE Genes page)

## Key Findings So Far

### 1. Data Source Differences ⚠️

**Signature Nomination Module:**
- **Data source**: `app_data$consolidated_data` (enrichment results)
- **Gene extraction**: `unique(unlist(strsplit(cluster_mast$geneID, "/")))`
- **Gene set**: Only genes that contributed to significantly enriched pathways
- **Background**: Genes from enrichment analysis background

**DE Results Module:**
- **Data source**: `values$de_data_mast` and `values$de_data_mixscale` (raw DE results)
- **Gene extraction**: Direct from DE results with thresholds
- **Gene set**: ALL differentially expressed genes (p < 0.05, |log2FC| > 0.25)
- **Background**: All tested genes in the experiment

### 2. Statistical Corrections 🔍

**Signature Nomination Module:**
- ✅ **Uses FDR correction**: "hierarchical FDR correction" (Benjamini-Hochberg)
- ✅ **Enhanced correction**: `intersection_fisher_p_fdr_enhanced_hierarchical`
- ✅ **Multiple approaches**: Intersection (conservative) vs Union (liberal)
- ✅ **Direction-aware**: Includes biological direction expectations

**DE Results Module:**
- ❌ **No FDR correction**: Uses raw Fisher's exact test p-values
- ❌ **Single approach**: Both intersection and union but no FDR
- ❌ **No direction analysis**: Simple overlap without biological context

### 3. Direction Filtering Analysis 📊

**Signature Nomination Module:**
- Uses consolidated data directly without global direction filtering
- Enrichment data may contain direction information (UP/DOWN/ALL)
- Direction filtering capabilities exist but unclear if applied

**DE Results Module:**
- Uses raw DE data without direction filtering
- Applies simple thresholds (p < 0.05, |log2FC| > 0.25)
- No direction-specific analysis

### 4. Background Gene Set Differences 🎯

**Signature Nomination Module:**
- **Background**: Genes associated with enrichment pathways
- **Approach**: Intersection (conservative) vs Union (liberal)
- **Size**: Likely smaller, more focused gene sets

**DE Results Module:**
- **Background**: All tested genes in the experiment
- **Approach**: Intersection vs Union of all tested genes
- **Size**: Likely larger, includes all experimental genes

## Critical Questions Still Unanswered

### 1. Gene Set Size Comparison
- **Q**: Are the gene sets consistently different in size?
- **Q**: Are there cases where they're similar? If so, are p-values similar?
- **Need**: Actual numerical comparison across multiple gene pairs and clusters

### 2. Direction Considerations
- **Q**: Does the signature nomination module filter by direction?
- **Q**: Does this affect the gene sets being compared?
- **Need**: Check if direction filtering is applied to consolidated data

### 3. FDR Correction Impact
- **Q**: Should FDR correction make p-values higher or lower?
- **Q**: Why are signature nomination p-values often lower despite FDR correction?
- **Need**: Understand the hierarchical FDR correction methodology

### 4. Background Gene Set Validation
- **Q**: What are the actual background gene set sizes?
- **Q**: Are they consistently different or sometimes similar?
- **Need**: Numerical comparison of background sizes

## Hypotheses to Test

### Primary Hypothesis: Gene Set Size Effect
The signature nomination module tests smaller, more focused gene sets (pathway-associated genes) leading to higher overlap ratios and lower p-values, even with FDR correction.

### Secondary Hypothesis: Background Size Effect
The signature nomination module uses smaller background gene sets (pathway-associated genes) while DE results uses larger backgrounds (all tested genes), affecting the Fisher's test calculation.

### Tertiary Hypothesis: FDR Correction Paradox
The hierarchical FDR correction in signature nomination may actually be correcting for a different type of multiple testing, potentially making the remaining significant comparisons more robust.

## Next Investigation Steps

### 1. Numerical Comparison (HIGH PRIORITY)
- [ ] Extract actual gene counts for both methods for 5-10 gene pairs
- [ ] Compare overlap counts between methods
- [ ] Compare background gene set sizes
- [ ] Calculate expected vs actual p-values

### 2. Direction Analysis (MEDIUM PRIORITY)
- [ ] Check if signature nomination applies direction filtering
- [ ] Compare UP vs DOWN vs ALL results
- [ ] Verify if direction affects gene set composition

### 3. FDR Methodology (MEDIUM PRIORITY)
- [ ] Understand hierarchical FDR correction implementation
- [ ] Verify if correction is applied appropriately
- [ ] Check if multiple testing burden is different between methods

### 4. Case Study Analysis (HIGH PRIORITY)
- [ ] Find specific cases where gene counts are similar
- [ ] Compare p-values in these cases
- [ ] Identify patterns in discrepancies

## Code Locations for Investigation

**Signature Nomination Fisher's Test:**
- `R/signature_analysis.R`: Lines 72-74, 113-115 (calculate_fisher_test calls)
- `R/signature_analysis.R`: Lines 169-204 (calculate_fisher_test function)

**DE Results Fisher's Test:**
- `inst/shiny/modules/mod_de_results.R`: Lines 1959, 1979 (direct fisher.test calls)
- `inst/shiny/modules/mod_de_results.R`: Lines 1945-1987 (Fisher's test implementation)

**FDR Correction:**
- `inst/shiny/modules/mod_signature_nomination.R`: Lines 1524, 1536, 1863, 1878 (FDR column usage)
- Need to find where hierarchical FDR correction is actually implemented

## Working Theory

The substantial p-value discrepancy is likely caused by:
1. **Different gene universes**: Pathway-focused vs all DE genes
2. **Different background sizes**: Smaller pathway backgrounds vs larger experiment backgrounds  
3. **Statistical correction differences**: FDR-corrected vs raw p-values
4. **Biological focus**: Pathway-enriched genes vs all significant genes

The signature nomination module is testing more biologically meaningful overlaps (pathway-focused) with appropriate statistical correction, while DE results tests broader overlaps (all DE genes) without correction.

Both approaches are scientifically valid but answer different questions:
- **Signature nomination**: "Do pathway-disrupting genes overlap more than expected?"
- **DE results**: "Do all significant genes overlap more than expected?"

## Status: INVESTIGATION ONGOING
Next step: Numerical validation with actual gene counts and p-value calculations.