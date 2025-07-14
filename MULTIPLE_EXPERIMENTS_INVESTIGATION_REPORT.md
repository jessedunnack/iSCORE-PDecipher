# Multiple CRISPRi Experiments Handling Investigation Report

**Date**: January 14, 2025  
**Investigator**: Claude Code Analysis  
**Scope**: Complete analysis of how multiple CRISPRi experiments are handled across all modules

## 🎯 Executive Summary

This investigation reveals **critical inconsistencies** in how multiple CRISPRi experiments (C12_FPD-23, C12_FPD-24, C18_FPD-23) are handled across different analysis modules. While some components properly preserve experiment information, others arbitrarily select only the first experiment, leading to incomplete and potentially misleading results.

## 🚨 Critical Issues Identified

### 1. DE Results Page - MAJOR FLAW ❌
**Location**: `inst/shiny/modules/mod_de_results.R`, line 78  
**Function**: `process_mixscale_for_volcano()`

**Problem**: 
```r
# Use the first log2FC column and corresponding p-value
log2fc_col <- log2fc_cols[1]  # ← ONLY USES FIRST EXPERIMENT!
```

**Impact**:
- Summary statistics show incomplete data (only one experiment)
- Overlap analysis misses genes that are DE in other experiments
- Arbitrary selection based on file order, not scientific criteria
- Users completely unaware that only one experiment is used
- **Scientific validity compromised**

### 2. User Transparency - INSUFFICIENT ⚠️
**Problem**: No clear indication to users about which experiments contribute to analyses

**Impact**:
- Users cannot make informed decisions about result interpretation
- Manuscript reporting may be inaccurate
- Reproducibility compromised

## ✅ Components Working Correctly

### 1. Enrichment Data Processing ✅
**Location**: `R/process_enrichment_results.R`

**Correct Behavior**:
- Lines 272-277: Properly extracts and preserves experiment metadata
- Each CRISPRi experiment creates separate enrichment entries
- Consolidated data includes: `method`, `mutation_perturbation`, `cluster`, `experiment`

### 2. Effect Size Correlation (NEW) ✅
**Location**: `inst/shiny/modules/mod_signature_nomination.R`, lines 1971-2273

**Correct Behavior**:
- Handles ALL experiments properly with user selection
- Transparent about which experiments are compared
- Provides both single-experiment and multi-experiment visualization
- Clear statistical reporting for each comparison

## 🔍 Areas Requiring Further Investigation

### 1. Signature Nomination Analysis ⚠️
**Status**: UNCLEAR

**Questions**:
- How does `analyze_gene_pair_signatures()` handle multiple experiments?
- Are pathway overlaps calculated using all experiments or subset?
- Do Fisher's exact tests account for experiment multiplicity?

**Investigation Needed**:
```r
# Check how this filtering works with multiple experiments:
crispri_data <- enrichment_data[enrichment_data$method == "MixScale" & 
                               enrichment_data$mutation_perturbation == gene_pair$crispri_gene, ]
```

## 🔧 **APPROVED SOLUTION PLAN** (Updated January 14, 2025)

### **Cell-Based Experiment Weighting Approach** ✅

**Primary Experiment**: **C12_FPD-24** (strongest results, highest cell count)

**Weighting Strategy**: Cell-number based weighting
- **Weight Calculation**: `weight = (perturbed_cells + control_cells_in_cluster) / total_cells_across_experiments`
- **Data Source**: Cell count CSV provided + Seurat object analysis if needed
- **Application**: Weighted meta-analysis across experiments

**Implementation**:
```r
# Cell-based experiment weighting
calculate_experiment_weights <- function(cell_counts_data) {
  
  experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
  weights <- list()
  
  for (exp in experiments) {
    for (cluster in paste0("cluster_", 0:11)) {
      perturbed_cells <- get_perturbed_cell_count(exp, cluster)
      control_cells <- cell_counts_data[[paste0("NT_", cluster)]]
      total_cells <- perturbed_cells + control_cells
      
      weights[[paste0(exp, "_", cluster)]] <- total_cells
    }
  }
  
  # Normalize weights within each cluster
  normalized_weights <- normalize_weights_by_cluster(weights)
  return(normalized_weights)
}

# Weighted meta-analysis implementation  
weighted_meta_analysis <- function(experiment_results, weights) {
  
  # For each gene and cluster, combine evidence across experiments
  combined_results <- list()
  
  for (gene in names(experiment_results)) {
    for (cluster in names(experiment_results[[gene]])) {
      
      exp_effects <- c()
      exp_weights <- c()
      
      for (exp in c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")) {
        if (!is.null(experiment_results[[gene]][[cluster]][[exp]])) {
          exp_effects <- c(exp_effects, experiment_results[[gene]][[cluster]][[exp]]$log2FC)
          exp_weights <- c(exp_weights, weights[[paste0(exp, "_", cluster)]])
        }
      }
      
      # Weighted average effect size
      if (length(exp_effects) > 0) {
        weighted_effect <- sum(exp_effects * exp_weights) / sum(exp_weights)
        combined_results[[gene]][[cluster]] <- list(
          weighted_log2FC = weighted_effect,
          experiments_included = length(exp_effects),
          primary_experiment = "C12_FPD-24",  # Highest weight expected
          weighting_method = "cell_count_based"
        )
      }
    }
  }
  
  return(combined_results)
}
```

#### Add Transparency Indicators
```r
# In UI, clearly show which experiments are included
div(class = "alert alert-info",
    "Analysis includes: ", paste(included_experiments, collapse = ", "),
    " (", length(included_experiments), " of 3 experiments)")
```

### Medium-Term Improvements (Priority 2)

#### Comprehensive Experiment Handling Framework
1. **Standardize**: Create consistent experiment handling across all modules
2. **Document**: Clear documentation of how each analysis handles experiments
3. **Validate**: Ensure statistical methods account for multiple testing
4. **Test**: Verify results consistency across different experiment selections

#### User Control Features
1. **Global Settings**: Allow users to specify experiment inclusion preferences
2. **Comparison Tools**: Enable side-by-side comparison of different experiment combinations
3. **Export Options**: Allow export of experiment-specific results

## 📊 Data Structure Analysis

### Enrichment Data Structure (CORRECT)
```
Columns: ID, Description, GeneRatio, BgRatio, pvalue, p.adjust, qvalue, geneID, Count, 
         method, mutation_perturbation, cluster, experiment, direction, enrichment_type, source_file
         
Example rows:
- method="MixScale", mutation_perturbation="LRRK2", cluster="cluster_0", experiment="C12_FPD-23"
- method="MixScale", mutation_perturbation="LRRK2", cluster="cluster_0", experiment="C12_FPD-24"  
- method="MixScale", mutation_perturbation="LRRK2", cluster="cluster_0", experiment="C18_FPD-23"
```

### DE Results Structure (PROBLEMATIC)
```
Current: Only includes data from first experiment found
Needed: Either all experiments as separate rows, or user-selectable experiment
```

## 🧪 Testing Protocol

### Validation Steps
1. **Verify Current Behavior**: 
   - Load DE results and check which experiment is actually used
   - Compare summary statistics with manual calculation including all experiments

2. **Test Proposed Fixes**:
   - Implement experiment selector
   - Verify overlap calculations include correct data
   - Test with different experiment selections

3. **Cross-Module Consistency**:
   - Ensure signature nomination and DE results use compatible approaches
   - Verify effect size correlation matches other analyses

## 📈 Impact Assessment

### Scientific Validity
- **High Risk**: Current DE results may miss important gene overlaps
- **Medium Risk**: Signature nomination may have similar issues
- **Low Risk**: Enrichment processing and effect size correlation handle experiments correctly

### User Experience
- **Confusion**: Users don't know which experiments are included
- **Inconsistency**: Different pages may use different experiment handling
- **Trust**: Discovery of arbitrary experiment selection may undermine confidence

### Manuscript Implications
- **Methods Section**: Need clear description of experiment handling
- **Results**: May need to re-analyze with all experiments included
- **Supplementary**: Should include experiment-specific analyses

## 📝 Action Items

### Immediate (This Week)
- [ ] Fix `process_mixscale_for_volcano()` to handle experiment selection
- [ ] Add transparency indicators to DE Results page
- [ ] Investigate signature nomination experiment handling

### Short-Term (Next Sprint)  
- [ ] Implement user-controllable experiment selection
- [ ] Add experiment comparison features
- [ ] Update documentation and help text

### Long-Term (Future Releases)
- [ ] Standardize experiment handling framework
- [ ] Add meta-analysis capabilities
- [ ] Implement experiment-specific statistical corrections

## 🔗 Related Files

### Files Modified in This Investigation
- `inst/shiny/modules/mod_signature_nomination.R` - Added effect size correlation with multi-experiment support

### Files Requiring Attention
- `inst/shiny/modules/mod_de_results.R` - **CRITICAL**: Fix experiment handling in lines 64-101
- `R/signature_analysis.R` - **INVESTIGATE**: Check experiment handling in lines 608-700

### Files Working Correctly
- `R/process_enrichment_results.R` - Proper experiment metadata preservation
- `R/calculate_effect_size_correlation` - Handles experiments correctly when called properly

---

**Last Updated**: January 14, 2025  
**Status**: Active Investigation - Critical Issues Identified  
**Next Review**: After implementing immediate fixes