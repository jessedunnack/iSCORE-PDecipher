# Fisher's Exact Test Improvements - Technical Documentation

## Overview

This document describes the critical improvements made to Fisher's exact test implementation for gene overlap significance testing between MAST (genetic mutations) and CRISPRi (gene knockdowns) experiments.

## Problem Identified

### Original Issue
The original Fisher's exact test implementation had a **fundamental statistical flaw**:
- **Background Gene Universe**: Used the union of significant genes from both methods as the background
- **Incorrect Logic**: This meant asking "Are these overlapping genes more than expected, given that they're all significant?"
- **Result**: Circular logic that inflated significance and used inappropriate background sizes

### Example of the Problem
```r
# WRONG (Original Implementation)
mast_significant_genes <- c("GENE1", "GENE2", "GENE3")
crispri_significant_genes <- c("GENE2", "GENE3", "GENE4") 
background_genes <- union(mast_significant_genes, crispri_significant_genes)  # Only 4 genes!
# Fisher's test asks: Is overlap of 2 genes significant in universe of 4 genes?
```

## Solution Implemented

### Proper Background Gene Handling
The improved implementation uses **genes actually tested in the differential expression analysis** as the background:

```r
# CORRECT (New Implementation) 
mast_background_genes <- all_genes_tested_in_mast_de_analysis    # ~20,000 genes
crispri_background_genes <- all_genes_tested_in_crispri_de_analysis  # ~18,000 genes

# Two approaches for combining backgrounds:
intersection_background <- intersect(mast_background_genes, crispri_background_genes)  # Conservative
union_background <- union(mast_background_genes, crispri_background_genes)  # Liberal
```

### Dual Approach Implementation

The new system provides **both intersection and union approaches** to give users comprehensive statistical information:

#### 1. **Intersection Approach (Conservative)**
- **Background**: Genes tested in BOTH methods
- **Logic**: More stringent test, focuses on genes comparable across methods
- **Use Case**: When you want to be confident both methods can detect the genes
- **Typical Background Size**: ~15,000-18,000 genes

#### 2. **Union Approach (Liberal)** 
- **Background**: Genes tested in EITHER method
- **Logic**: More inclusive test, captures broader gene universe
- **Use Case**: When you want to include method-specific capabilities
- **Typical Background Size**: ~20,000-22,000 genes

## Technical Implementation

### Backend Changes

#### 1. Core Function: `calculate_gene_overlap_significance_proper()`
**File**: `R/signature_analysis.R`

```r
calculate_gene_overlap_significance_proper <- function(mast_genes, crispri_genes, 
                                                      mast_background_genes = NULL,
                                                      crispri_background_genes = NULL,
                                                      alternative = "greater") {
  
  # Returns both approaches
  results <- list(
    intersection_approach = list(
      overlap_count = ...,
      fisher_p = ...,
      background_size = ...,
      background_type = "intersection"
    ),
    union_approach = list(
      overlap_count = ...,
      fisher_p = ...,
      background_size = ...,
      background_type = "union"
    )
  )
}
```

#### 2. Results Integration: `compute_signature_rankings()`
**File**: `R/manuscript_signature_discovery.R`

The function now includes both approaches in results:
```r
signature_entry <- data.frame(
  # ... existing columns ...
  # NEW: Intersection approach
  intersection_overlap_count = ...,
  intersection_fisher_p = ...,
  intersection_background_size = ...,
  # NEW: Union approach  
  union_overlap_count = ...,
  union_fisher_p = ...,
  union_background_size = ...
)
```

#### 3. DE Data Integration: `analyze_gene_pair_signatures()`
**File**: `R/signature_analysis.R`

Enhanced to extract proper background genes from DE data:
```r
# Extract background genes from DE analysis results
mast_background_genes <- de_data$iSCORE_PD_MAST[[gene]][[cluster]]$background_genes
crispri_background_genes <- de_data$CRISPRi_Mixscale[[gene]][[cluster]]$background_genes

# Use proper function for Fisher's exact test
overlap_stats <- calculate_gene_overlap_significance_proper(
  mast_genes, crispri_genes, 
  mast_background_genes, crispri_background_genes
)
```

### Frontend Changes

#### 1. UI Display Enhancement
**File**: `inst/shiny/modules/mod_signature_nomination.R`

The Gene Pair Analysis table now shows both approaches:

| Column | Description |
|--------|-------------|
| Fisher p (Intersection) | P-value using conservative background |
| Background (Intersection) | Number of genes tested in BOTH methods |
| Fisher p (Union) | P-value using liberal background |
| Background (Union) | Number of genes tested in EITHER method |
| Intersection Sig | Significance level for intersection approach |
| Union Sig | Significance level for union approach |

#### 2. Backwards Compatibility
The UI automatically detects whether new or legacy results are available:
```r
has_intersection_union <- all(c("intersection_fisher_p", "union_fisher_p") %in% colnames(data))

if (has_intersection_union) {
  # Display new dual approach results
} else {
  # Display legacy single approach results
}
```

## Expected User Experience

### Before (Legacy)
```
Gene Pair: LRRK2_vs_LRRK2
Cluster: cluster_0
Shared DE Genes: 12
DE Overlap p-value: 0.023
Background Genes: 20,207
Test Type: legacy
```

### After (Improved)
```
Gene Pair: LRRK2_vs_LRRK2  
Cluster: cluster_0
Shared DE Genes: 12
Fisher p (Intersection): 0.045    Background (Intersection): 15,823
Fisher p (Union): 0.012          Background (Union): 20,891
Intersection Sig: *              Union Sig: *
```

## Statistical Interpretation Guide

### When to Use Each Approach

#### Intersection Approach (Conservative)
- **Use when**: You want high confidence in cross-method comparability
- **Interpretation**: "Among genes that both methods can reliably detect, is this overlap significant?"
- **Best for**: Identifying robust, reproducible signatures
- **Limitation**: May miss method-specific effects

#### Union Approach (Liberal)  
- **Use when**: You want to capture the full scope of both methods
- **Interpretation**: "Considering all genes either method can detect, is this overlap significant?"
- **Best for**: Exploratory analysis and comprehensive screening
- **Limitation**: May include method-specific artifacts

### Combined Interpretation
- **Both significant**: Strong evidence of genuine biological overlap
- **Only union significant**: Possible method-specific effects, investigate further
- **Only intersection significant**: Robust core signature, methods may have different sensitivities
- **Neither significant**: No evidence of meaningful overlap

## Validation and Testing

### Files Modified
1. `R/signature_analysis.R` - Core Fisher's exact test functions
2. `R/manuscript_signature_discovery.R` - Results processing and ranking
3. `inst/shiny/modules/mod_signature_nomination.R` - UI display
4. `inst/shiny/modules/mod_signature_nomination.R` - DE data loading

### Testing Requirements
1. **Backend Testing**: Verify intersection/union calculations
2. **UI Testing**: Confirm both approaches display correctly
3. **Integration Testing**: Ensure DE data loads and flows through pipeline
4. **Edge Case Testing**: Handle missing DE data gracefully

## Benefits of This Implementation

### 1. **Statistical Rigor**
- Proper background gene universe from actual DE analysis
- Eliminates circular logic of original implementation
- Provides conservative and liberal statistical approaches

### 2. **User Transparency**
- Shows both background sizes so users understand the test
- Clear labeling of intersection vs union approaches
- Color-coded significance levels for easy interpretation

### 3. **Flexibility**
- Users can choose appropriate approach based on analysis goals
- Backwards compatibility with legacy results
- Supports both exploratory and confirmatory analysis

### 4. **Reproducibility**
- Background genes are extracted consistently from DE analysis
- Methodology is clearly documented and traceable
- Results include all necessary information for replication

## Future Enhancements

### Potential Improvements
1. **Method-Specific Background Genes**: Consider experiment-specific backgrounds
2. **Multiple Testing Correction**: Implement FDR correction across gene pairs
3. **Effect Size Reporting**: Add odds ratios and confidence intervals
4. **Visualization**: Create plots showing background gene overlaps
5. **Export Functionality**: Allow users to download detailed Fisher's test results

### User Education
1. **Tutorial**: Interactive guide explaining intersection vs union
2. **Examples**: Worked examples with real data
3. **Best Practices**: Guidelines for choosing appropriate approach
4. **Troubleshooting**: Common issues and solutions

## Conclusion

This implementation addresses the fundamental statistical flaw in the original Fisher's exact test while providing users with comprehensive, transparent, and flexible tools for cross-method signature analysis. The dual approach (intersection/union) gives researchers the ability to choose the appropriate statistical rigor for their specific research questions.

---

**Document Version**: 1.0  
**Date**: January 11, 2025  
**Author**: Claude Code Assistant  
**Status**: Implemented and Ready for Testing