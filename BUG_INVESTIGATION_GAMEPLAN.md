# Bug Investigation Game Plan for iSCORE-PDecipher
**Date:** January 16, 2025
**Version:** v0.2.7

## Overview
This document tracks the investigation and resolution plan for 9 identified bugs in the iSCORE-PDecipher Shiny app.

---

## Bug #1: Gene-pair Analysis - Investigate Shared DE Genes/Terms
**Priority:** High  
**Status:** In Progress

### Issue
When the gene-pair analysis table shows shared DE genes or functional enrichment terms between MAST and CRISPRi, users cannot investigate these further.

### Investigation Findings
**GOOD NEWS:** The gene lists are already available in the analysis results!
- `signature_analysis.R` returns `overlap_genes` in the results
- Both intersection and union approaches store the actual gene lists
- The data is available but not exposed in the UI

**Current Table Shows:**
- "Shared DE Genes*" (count only)
- "Enriched Pathways" (count only)
- Jaccard Index
- p-values

**Missing:**
- Actual list of shared genes
- Shared pathway terms
- Gene-specific statistics from both methods

### Proposed Solution
1. **Add "Details" button** to each row in the gene_pair_table
2. **Create modal dialog** that shows:
   - **Shared Genes Tab**: 
     - List of overlapping genes
     - MAST statistics (log2FC, p-value) for each gene
     - CRISPRi statistics (log2FC, p-value) for each gene
     - Sortable/filterable table
   - **Shared Pathways Tab**:
     - List of enriched pathways found in both methods
     - p-values from both methods
     - Gene counts
3. **Export functionality**: Download buttons for both gene and pathway lists

### Implementation Location
- `inst/shiny/modules/mod_signature_nomination.R` - add modal UI and server logic
- Access overlap data from `values$analysis_results`
- Use `overlap_genes` from intersection/union approach results
- May need to modify `calculate_cross_method_signatures` to also return pathway overlaps

---

## Bug #2: Fisher's Exact Test P-value Discrepancy ⚠️ CRITICAL
**Priority:** High  
**Status:** In Progress

### Issue
Fisher's exact test p-values differ drastically between:
- Gene-pair analysis page (signature nomination module)
- DE Genes page (DE results module)

### Investigation Findings
**KEY DISCOVERY:** The two modules test fundamentally different things!

1. **Signature Nomination Module** (Gene-pair analysis):
   - Tests overlap of genes contributing to **significantly enriched pathways**
   - Uses filtered gene sets from enrichment analysis
   - Background: Genes associated with enriched terms
   - More stringent, pathway-focused analysis

2. **DE Results Module** (DE Genes page):
   - Tests overlap of **ALL differentially expressed genes**
   - Uses direct DE gene lists (p < 0.05, |log2FC| > 0.25)
   - Background: All tested genes in the experiment
   - Broader, unfiltered DE analysis

### Explanation
The discrepancy is **expected and correct**! The signature nomination module has:
- Smaller gene sets (only pathway-associated genes)
- More biologically meaningful overlaps
- Lower p-values due to focused comparison

The DE results module has:
- Larger gene sets (all DE genes)
- More noise from non-pathway genes
- Higher p-values due to broader comparison

### Proposed Solution
1. Add clear documentation/tooltips explaining the difference
2. Label the tests distinctly:
   - Signature Nomination: "Pathway-Focused Gene Overlap"
   - DE Results: "All DE Gene Overlap"
3. Consider adding the pathway-focused test to DE Results page for consistency

### Trust Assessment
**Both are correct**, but serve different purposes:
- **Signature Nomination**: Better for finding biologically meaningful convergence
- **DE Results**: Better for overall DE similarity assessment

---

## Bug #3: Cross-cluster DE Gene Reporting
**Priority:** Medium  
**Status:** Pending

### Issue
No way to see if the same DE genes appear across multiple clusters for a given mutation/perturbation.

### Investigation Findings
- Current tables show gene-cluster pairs independently
- No aggregation across clusters

### Proposed Solution
1. Add a new tab or section showing:
   - Genes that are DE in multiple clusters
   - Heatmap of gene x cluster showing DE status
   - Count of clusters where each gene is DE
2. Filter by mutation/perturbation
3. Sort by number of clusters affected

### Implementation Location
- Add new module or extend `mod_de_results.R`
- Query DE results across all clusters for selected gene

---

## Bug #4: DE Genes Heatmap Using All DE Genes
**Priority:** Medium  
**Status:** Pending

### Issue
Request to implement Fisher's test heatmap on DE Genes page using ALL DE genes (not just enriched ones).

### Investigation Findings
- Current signature heatmap uses enrichment-filtered genes
- No equivalent visualization for all DE genes

### Proposed Solution
1. Adapt the heatmap code from `mod_signature_nomination.R`
2. Calculate Fisher's tests for all gene pairs using ALL DE genes
3. Create matrix of p-values or overlap counts
4. Display as interactive heatmap similar to signature heatmap

### Implementation Location
- `mod_de_results.R` - add new heatmap visualization
- Reuse heatmap logic from signature nomination module

---

## Bug #5: Dotplot Hover "No Genes Available" Issue
**Priority:** High  
**Status:** In Progress

### Issue
Functional enrichment dotplots sometimes show "no genes available" on hover for PerturbSeq results.

### Investigation Findings
**ROOT CAUSE IDENTIFIED:** The gene lookup function uses a strict composite key matching system:
```
composite_key = analysis_type|gene|cluster|enrichment_type|direction|default|term_id
```

For PerturbSeq data:
1. The exact composite key match often fails due to naming differences
2. The fallback (term_id only match) also fails if gene associations weren't properly stored
3. The matching is case-sensitive and format-sensitive

**Code Analysis:**
- Lines 197-232 in `mod_visualization_enhanced.R`
- Creates composite key with 7 fields separated by "|"
- Falls back to term_id match if exact match fails
- Returns "No genes found" if both fail

### Proposed Solution
1. **Immediate fix**: Add debug logging to identify which terms/keys are failing
2. **Better matching**: Make the composite key matching more flexible:
   - Try partial key matches (e.g., term_id + enrichment_type)
   - Handle case variations
   - Account for PerturbSeq-specific naming conventions
3. **Data validation**: Check if gene associations are properly loaded for PerturbSeq
4. **Fallback improvement**: Use a more intelligent fallback strategy

### Implementation Location
- `inst/shiny/modules/mod_visualization_enhanced.R` - lines 197-232
- Focus on `get_genes_for_term_enhanced` function
- May need to check data loading in app initialization

---

## Bug #6: Heatmap Color Scale for P-values
**Priority:** Medium  
**Status:** Pending

### Issue
Fisher's p-value color scale (-log10(p)) is not intuitive. Hard to identify significant values.

### Investigation Findings
- Current scale is continuous
- No clear threshold indication for p < 0.05

### Proposed Solution
1. Add significance threshold line or annotation at -log10(0.05) = 1.3
2. Consider discrete color bins:
   - Not significant (p > 0.05): gray/white
   - Significant (p < 0.05): color gradient
3. Add legend explaining the scale
4. Option to show raw p-values instead of -log10

### Implementation Location
- `mod_signature_nomination.R` - signature heatmap
- Any other heatmaps showing p-values

---

## Bug #7: DE Gene Table Search Restriction
**Priority:** Medium  
**Status:** Pending

### Issue
Search function searches all columns, not just gene names. Searching "PINK1" returns rows where PINK1 is the mutation, not the DE gene.

### Investigation Findings
- DT::datatable search is global by default
- Need column-specific search

### Proposed Solution
1. Disable global search
2. Add column-specific search filters
3. Or modify search to only target the gene name column
4. Use DT options: `searchCols` parameter

### Implementation Location
- `mod_de_results.R` - gene overlap table configuration
- Modify DT::datatable options

---

## Bug #8: DE Gene Table - Add Experiment Metadata
**Priority:** Low  
**Status:** Pending

### Issue
CRISPRi DE genes appear multiple times (from different experiments) without indication of which experiment.

### Investigation Findings
- Multiple CRISPRi experiments can identify the same gene as DE
- No experiment ID shown in current table

### Proposed Solution
1. Add "Experiment" column to the gene table
2. For MAST: show batch ID
3. For CRISPRi: show experiment ID (e.g., "C12_FPD-23")
4. Extract from existing data structure

### Implementation Location
- `mod_de_results.R` - modify table generation
- Add experiment metadata extraction

---

## Bug #9: Gene Selection Dropdown Consolidation
**Priority:** Medium  
**Status:** Pending

### Issue
When "all datasets" selected, dropdown shows variants separately (VPS13C, VPS13C_A444P, VPS13C_W395C) instead of consolidated gene names.

### Investigation Findings
- Current implementation shows all unique mutation/perturbation names
- Need intelligent grouping by base gene name

### Proposed Solution
1. When "all datasets" selected:
   - Show base gene names only (VPS13C, SNCA, PRKN)
   - Query should include all variants automatically
2. When specific dataset selected:
   - Show dataset-specific names (mutations for iSCORE-PD, genes for CRISPRi)
3. Handle special cases:
   - PARK2/PRKN consolidation
   - VPS13C variants
   - SNCA variants

### Implementation Location
- `mod_de_results.R` - gene selection dropdown logic
- May need helper function to extract base gene names

---

## Next Steps
1. Complete investigation of Bug #2 (DONE - documented findings)
2. Create test cases for each bug
3. Implement fixes in priority order
4. Update documentation for each fix
5. Test cross-module interactions

## Testing Requirements
- Verify fixes don't break existing functionality
- Test with all three dataset options
- Confirm performance isn't degraded
- Check cross-platform compatibility