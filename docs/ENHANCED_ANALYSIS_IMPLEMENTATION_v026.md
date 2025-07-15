# Enhanced Cross-Method Analysis Implementation Plan v0.2.6

**Date**: January 14, 2025  
**Target Release**: v0.2.6  
**Previous Version**: v0.2.5 (backup)  
**Scope**: Complete implementation of experiment weighting + direction-aware analysis

## 🎯 Executive Summary

This document outlines the comprehensive implementation plan for version 0.2.6, which addresses two critical statistical gaps in cross-method comparison:

1. **Multiple CRISPRi Experiments**: Implement cell-based weighting with C12_FPD-24 as primary
2. **Direction-Aware Analysis**: Unified statistical framework with biological context (LRRK2 vs SNCA)

## 📊 User-Approved Design Decisions

### Experiment Handling
- **Primary Experiment**: C12_FPD-24 (highest cell count, strongest results)
- **Weighting Method**: Cell-count based (perturbed + control cells per cluster)  
- **Meta-Analysis**: Weighted combination across experiments
- **Default Behavior**: New enhanced analysis (with legacy compatibility)

### Direction Analysis  
- **LRRK2**: Expect opposing effects (gain-of-function mutation vs loss-of-function CRISPRi)
- **SNCA Variants**: Expect same-direction effects (aggregation-related in both cases)
- **Other Mutations**: Generally same-direction (loss-of-function mutations)
- **Statistical Approach**: Single unified test with biological weighting
- **Pathway Analysis**: Direction-aware maintained

### Multiple Testing Correction
- **Research-Based Approach**: Investigate factorial design FDR methods
- **Implementation**: Enhanced hierarchical correction for experiments × directions × gene pairs
- **Transparency**: Clear documentation of correction methodology

### User Experience
- **Single Analysis Approach**: Unified enhanced method (not separate tabs)
- **Default Settings**: New method as default for all new analyses
- **Backward Compatibility**: Preserve option to revert to legacy approach
- **Transparency**: Clear indicators of analysis parameters

## 🧪 Implementation Phases

### Phase 1: Research & Core Functions (Days 1-2)

#### 1.1 FDR Correction Research ✅
**Objective**: Research appropriate FDR correction for factorial design
**Tasks**:
- [ ] Literature review: Multiple testing in factorial designs
- [ ] Evaluate Benjamini-Hochberg vs other methods for nested structure
- [ ] Document chosen approach with scientific rationale
- [ ] Implement enhanced FDR correction function

**Research Questions**:
- Should experiments be treated as separate hypothesis families?
- How to handle Direction × Experiment × Gene Pair factorial structure?
- What are best practices for weighted meta-analysis p-value correction?

#### 1.2 Cell Count Data Integration ✅
**Objective**: Extract and process cell count data for experiment weighting
**Tasks**:
- [ ] Load and analyze provided cell count CSV
- [ ] Extract perturbed cell counts from Seurat object if needed
- [ ] Calculate cluster-specific experiment weights
- [ ] Validate that C12_FPD-24 receives highest weights

**Technical Implementation**:
```r
# Load cell count data and calculate weights
load_experiment_weights <- function(csv_path = NULL, seurat_path = NULL) {
  
  # Load non-targeting cell counts from CSV
  nt_counts <- read.csv(csv_path)  # Columns: cluster_0 through cluster_11
  
  # Extract perturbed cell counts from Seurat if needed
  if (!is.null(seurat_path)) {
    seurat_obj <- readRDS(seurat_path)
    perturbed_counts <- extract_perturbed_cell_counts(seurat_obj)
  }
  
  # Calculate weights by cluster and experiment
  weights <- calculate_cluster_experiment_weights(nt_counts, perturbed_counts)
  return(weights)
}
```

#### 1.3 Direction-Aware Statistical Functions ✅
**Objective**: Implement core direction-aware Fisher's tests
**Tasks**:
- [ ] Create `calculate_same_direction_overlap()` function
- [ ] Create `calculate_opposite_direction_overlap()` function  
- [ ] Implement `combine_direction_pvalues()` with biological weighting
- [ ] Add biological expectation function for each gene

### Phase 2: Core Statistical Implementation (Days 2-3)

#### 2.1 Enhanced Fisher's Exact Tests ✅
**File**: `R/signature_analysis.R`
**Tasks**:
- [ ] Modify `calculate_gene_overlap_significance_proper()` to include direction analysis
- [ ] Add experiment weighting to Fisher's test calculations
- [ ] Implement biological context weighting (LRRK2 vs SNCA)
- [ ] Add comprehensive error handling and validation

**Key Function**:
```r
calculate_enhanced_overlap_significance <- function(mast_genes, crispri_experiments_data, 
                                                   gene_name, cluster, experiment_weights,
                                                   background_genes = NULL) {
  
  # Determine biological direction expectation
  direction_exp <- get_biological_direction_expectation(gene_name)
  
  # Calculate direction-specific overlaps for each experiment
  experiment_results <- list()
  
  for (exp in names(crispri_experiments_data)) {
    exp_data <- crispri_experiments_data[[exp]]
    
    same_dir <- calculate_same_direction_overlap(mast_genes, exp_data, background_genes)
    opposite_dir <- calculate_opposite_direction_overlap(mast_genes, exp_data, background_genes)
    
    experiment_results[[exp]] <- list(
      same_direction = same_dir,
      opposite_direction = opposite_dir,
      weight = experiment_weights[[paste0(exp, "_", cluster)]]
    )
  }
  
  # Weighted meta-analysis across experiments
  meta_analysis <- weighted_direction_meta_analysis(experiment_results, direction_exp)
  
  return(meta_analysis)
}
```

#### 2.2 Enhanced FDR Correction ✅
**File**: `R/manuscript_signature_discovery.R`
**Tasks**:
- [ ] Extend `apply_hierarchical_fdr_correction()` for experiments and directions
- [ ] Implement three-level correction: within gene pairs → across gene pairs → across experiments
- [ ] Add clear documentation of correction rationale
- [ ] Validate that FDR remains controlled

**Implementation Structure**:
```r
apply_enhanced_fdr_correction <- function(signature_rankings) {
  
  # Level 1: Within gene pair correction (experiments × directions)
  signature_rankings <- apply_within_gene_pair_correction(signature_rankings)
  
  # Level 2: Across gene pairs correction  
  signature_rankings <- apply_across_gene_pairs_correction(signature_rankings)
  
  # Level 3: Global experiment family correction
  signature_rankings <- apply_experiment_family_correction(signature_rankings)
  
  # Document correction pathway
  attr(signature_rankings, "fdr_method") <- "enhanced_hierarchical_three_level"
  attr(signature_rankings, "correction_levels") <- c("gene_pair", "across_pairs", "experiment_family")
  
  return(signature_rankings)
}
```

### Phase 3: Integration & UI Updates (Days 4-5)

#### 3.1 Signature Analysis Integration ✅
**File**: `R/manuscript_signature_discovery.R`
**Tasks**:
- [ ] Enable direction statistics in line 332 (currently NULL)
- [ ] Integrate experiment weighting into signature scoring
- [ ] Update composite score calculation with enhanced statistics
- [ ] Add biological context to signature ranking

**Critical Update**:
```r
# Line 332 modification - ENABLE DIRECTION STATISTICS
composite_score <- calculate_composite_signature_score(
  overlap_stats = overlap_stats_display,
  correlation_stats = NULL,  # Still not available from enrichment data directly
  direction_stats = enhanced_direction_statistics,  # ← NEW: Enable direction analysis
  pathway_overlap_stats = pathway_stats,
  experiment_weights = experiment_weights  # ← NEW: Add experiment weighting
)
```

#### 3.2 DE Results Module Fix ✅
**File**: `inst/shiny/modules/mod_de_results.R`
**Tasks**:
- [ ] Fix `process_mixscale_for_volcano()` line 78 (currently only uses first experiment)
- [ ] Implement experiment selection or weighted meta-analysis
- [ ] Add transparency indicators for experiment inclusion
- [ ] Update summary statistics to reflect all experiments

**Before (Broken)**:
```r
log2fc_col <- log2fc_cols[1]  # Only first experiment
```

**After (Fixed)**:
```r
# Option 1: Weighted meta-analysis approach
weighted_results <- calculate_weighted_de_across_experiments(de_data, experiment_weights)
log2fc_values <- weighted_results$weighted_log2fc
p_values <- weighted_results$weighted_p_values

# Option 2: User-selectable experiment with default to C12_FPD-24
selected_experiment <- input$selected_experiment %||% "C12_FPD-24"
log2fc_col <- paste0("log2FC_", selected_experiment)
```

#### 3.3 UI Enhancement ✅
**File**: `inst/shiny/modules/mod_signature_nomination.R`
**Tasks**:
- [ ] Add experiment handling transparency indicators
- [ ] Include direction pattern information in results tables
- [ ] Add biological expectation vs observed pattern comparisons
- [ ] Update help text and documentation

**UI Elements**:
```r
# Experiment transparency
div(class = "alert alert-info",
    icon("info-circle"), " Analysis Method: ",
    "Enhanced cross-method comparison with cell-weighted experiment meta-analysis. ",
    "Primary experiment: C12_FPD-24 (highest cell count). ",
    "Direction analysis: Biologically-informed (LRRK2 opposing, SNCA same-direction)."
)

# Direction pattern display
DT::datatable(display_data) %>%
  DT::formatStyle("Direction_Pattern",
    backgroundColor = DT::styleEqual(
      c("Convergent ↑↑/↓↓", "Divergent ↑↓", "Both patterns **"),
      c("#d4edda", "#fff3cd", "#f8d7da")
    )
  ) %>%
  DT::formatStyle("Biological_Relevance", 
    color = DT::styleEqual(
      c("Expected", "Unexpected", "Mixed"),
      c("#155724", "#856404", "#721c24")
    )
  )
```

### Phase 4: Testing & Validation (Day 6)

#### 4.1 LRRK2 Validation ✅
**Objective**: Verify LRRK2 opposing effects detection
**Tasks**:
- [ ] Test LRRK2 MAST vs CRISPRi analysis with new method
- [ ] Compare results before/after direction-aware analysis
- [ ] Verify that opposite-direction overlaps receive higher statistical weight
- [ ] Document specific biological signatures detected

**Test Case**:
```r
# LRRK2 validation test
test_lrrk2_opposing_effects <- function() {
  
  # Load LRRK2 data
  lrrk2_mast <- load_mast_data("LRRK2", cluster = "cluster_0")
  lrrk2_crispri <- load_crispri_data("LRRK2", cluster = "cluster_0")
  
  # Run enhanced analysis
  enhanced_results <- calculate_enhanced_overlap_significance(
    lrrk2_mast, lrrk2_crispri, gene_name = "LRRK2", cluster = "cluster_0"
  )
  
  # Validate expectations
  expect_true(enhanced_results$biological_expectation == "opposing")
  expect_true(enhanced_results$opposite_direction$fisher_p < enhanced_results$same_direction$fisher_p)
  expect_true(enhanced_results$primary_pattern == "opposite")
  
  return(enhanced_results)
}
```

#### 4.2 SNCA Validation ✅  
**Objective**: Verify SNCA same-direction analysis
**Tasks**:
- [ ] Test SNCA variants analysis with new method
- [ ] Verify that same-direction patterns receive priority
- [ ] Confirm that opposite-direction is not artificially forced
- [ ] Compare results across SNCA_A30P and SNCA_A53T variants

#### 4.3 Experiment Weighting Validation ✅
**Objective**: Verify experiment weighting functions correctly
**Tasks**:
- [ ] Confirm C12_FPD-24 receives highest weights across clusters
- [ ] Test that meta-analysis results change appropriately with weighting
- [ ] Validate that cell count data integration works correctly
- [ ] Compare weighted vs unweighted results

#### 4.4 FDR Correction Validation ✅
**Objective**: Ensure FDR correction maintains statistical validity
**Tasks**:
- [ ] Test that false discovery rate remains controlled
- [ ] Compare p-value distributions before/after correction
- [ ] Validate correction factors are appropriate
- [ ] Document any inflation or deflation in significance levels

## 📚 Files Modified

### Core Statistical Functions
- `R/signature_analysis.R` - Enhanced Fisher's tests with direction and experiment weighting
- `R/manuscript_signature_discovery.R` - Enable direction stats, enhanced FDR correction
- `R/gene_harmonization.R` - Add biological direction expectations

### UI Modules  
- `inst/shiny/modules/mod_signature_nomination.R` - UI transparency and direction display
- `inst/shiny/modules/mod_de_results.R` - Fix experiment handling in DE analysis

### New Utility Functions
- `R/experiment_weighting.R` - Cell count analysis and experiment weighting
- `R/enhanced_direction_analysis.R` - Direction-aware statistical functions
- `R/enhanced_fdr_correction.R` - Three-level FDR correction methodology

## 🧪 Success Metrics

### Statistical Improvements
- [ ] LRRK2 analysis detects previously missed opposing-direction signatures
- [ ] SNCA analysis appropriately handles same-direction expectations
- [ ] Experiment weighting increases statistical power through proper meta-analysis
- [ ] FDR correction maintains Type I error control while increasing power

### User Experience  
- [ ] Complete transparency about which experiments and directions are analyzed
- [ ] Clear biological context provided for each gene's direction expectations
- [ ] Enhanced results tables with direction pattern information
- [ ] Backward compatibility preserved for legacy analysis approach

### Scientific Validity
- [ ] All statistical methods properly documented for manuscript reporting
- [ ] Biological expectations match known mutation mechanisms
- [ ] Cross-method validation strengthened through enhanced analysis
- [ ] Results interpretable and actionable for downstream analysis

## 🔄 Backward Compatibility

### Legacy Analysis Preservation
- [ ] Add `use_enhanced_analysis = TRUE` parameter to main functions
- [ ] Preserve old functions with `_legacy` suffix
- [ ] Create wrapper functions for switching between approaches
- [ ] Document migration path from legacy to enhanced analysis

### Version Control Strategy
- **v0.2.5**: Current stable backup before changes
- **v0.2.6**: Enhanced analysis implementation
- **Fallback**: Ability to revert to v0.2.5 if needed
- **Migration**: Clear documentation for users upgrading

## 📋 Implementation Checklist

### Phase 1: Research & Core Functions
- [ ] Complete FDR correction methodology research
- [ ] Implement cell count data loading and processing
- [ ] Create direction-aware statistical functions
- [ ] Validate biological expectation assignments

### Phase 2: Core Statistical Implementation  
- [ ] Enhance Fisher's exact tests with direction and weighting
- [ ] Implement three-level FDR correction
- [ ] Integrate experiment weighting throughout analysis pipeline
- [ ] Add comprehensive error handling and validation

### Phase 3: Integration & UI Updates
- [ ] Enable direction statistics in signature discovery
- [ ] Fix DE results experiment handling
- [ ] Update UI with transparency indicators
- [ ] Add direction pattern visualization

### Phase 4: Testing & Validation
- [ ] Validate LRRK2 opposing effects detection
- [ ] Validate SNCA same-direction analysis
- [ ] Test experiment weighting functionality
- [ ] Verify FDR correction performance

### Phase 5: Documentation & Release
- [ ] Update all documentation with new methodology
- [ ] Create comprehensive testing report
- [ ] Generate migration guide for users
- [ ] Tag v0.2.6 release with full feature documentation

---

**Status**: Implementation Plan Approved  
**Next Step**: Begin Phase 1 Research & Core Functions  
**Target Completion**: 6 days from approval  
**Quality Standard**: Scientific rigor, user transparency, backward compatibility