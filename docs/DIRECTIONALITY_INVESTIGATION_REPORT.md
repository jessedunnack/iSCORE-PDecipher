# Directionality Analysis Investigation Report

**Date**: January 14, 2025  
**Investigator**: Claude Code Analysis  
**Scope**: Investigation of directionality consideration in cross-method overlap analyses (MAST vs CRISPRi)

## 🎯 Executive Summary

This investigation reveals a **critical gap** in the current cross-method comparison analysis: **Fisher's exact tests and overlap analyses are completely direction-agnostic**, ignoring whether gene expression changes occur in the same or opposite directions between mutation and perturbation experiments. This is particularly problematic for gain-of-function mutations like LRRK2 G2019S that should show **opposing effects** to loss-of-function CRISPRi perturbations.

## 🧬 Mutation Classification Analysis

### Gene Categories Discovered ✅
**Location**: `R/gene_harmonization.R`, `get_mutation_categories()` function (lines 126-192)

#### 1. Point Mutations (Missense Variants)
**Genes**: SNCA_A30P, SNCA_A53T, **LRRK2**, VPS13C_W395C, VPS13C_A444P, SYNJ1
- **Expected Expression Effect**: "Minimal_Change" 
- **Biological Effect**: Functional changes without major expression level changes
- **Key Example**: **LRRK2 G2019S** = hyperactive kinase (gain-of-function)

#### 2. Loss-of-Function Mutations
**Nonsense/Truncating**: PINK1, FBXO7 → "Moderate_Severe_Reduction"  
**Large Deletions**: PRKN, PARK7 → "Strong_Reduction"  
**Frameshift**: ATP13A2, DNAJC6 → "Moderate_Severe_Reduction"  
**Splice Site**: GBA → "Reduced_Expression"

### **UPDATED Biological Context** 🧬 (January 14, 2025)

#### **LRRK2 - Opposing Effects Expected** ✅
- **LRRK2 G2019S mutation**: Creates **hyperactive kinase** → **gain-of-function**
- **LRRK2 CRISPRi**: Reduces **LRRK2 expression** → **loss-of-function**
- **Expected Result**: **OPPOSING effects** on downstream gene expression
- **Analysis Approach**: **Prioritize opposite-direction overlap tests**

#### **SNCA Variants - NON-Opposing Effects** ⚠️ 
- **SNCA A30P/A53T mutations**: **Aggregation-prone variants** → **Lewy body formation**
- **SNCA CRISPRi**: Reduces **SNCA expression** → **Less protein aggregation**
- **Expected Result**: **SAME-DIRECTION effects** (both reduce aggregation pathways)
- **Analysis Approach**: **Do NOT force opposite-direction analysis**

#### **Other Mutations - Generally Same Direction** 📊
- **Loss-of-function mutations** (PINK1, PRKN, PARK7, etc.): Expected to show **same-direction** effects as CRISPRi
- **Analysis Approach**: **Default to same-direction analysis, but test both patterns**

## 🔍 Current Directionality Handling Analysis

### Direction Analysis Function EXISTS ✅
**Location**: `R/signature_analysis.R`, `calculate_direction_consistency()` (lines 372-414)

**Capabilities**:
```r
# Can detect three patterns:
same_direction_up     # Both methods increase gene expression
same_direction_down   # Both methods decrease gene expression  
opposite_direction    # One increases, other decreases (LRRK2 case!)
```

**Statistical Output**:
- Counts for each direction pattern
- Consistency percentage
- Detailed gene-by-gene analysis

### BUT: Function is NOT USED ❌
**Location**: `R/manuscript_signature_discovery.R`, line 332

**Problem**:
```r
composite_score <- calculate_composite_signature_score(
  overlap_stats = overlap_stats_display,
  correlation_stats = NULL,  # Not available from enrichment data directly
  direction_stats = NULL,    # Would need DE data ← CRITICAL ISSUE!
  pathway_overlap_stats = pathway_stats
)
```

**Impact**: Direction analysis is completely bypassed in main signature discovery

## 🚨 Critical Issues Identified

### 1. Direction-Agnostic Fisher's Exact Tests ❌
**Current Behavior**:
- Tests only whether genes are significantly DE in **both** methods
- **Ignores direction** of expression changes completely
- Treats same-direction and opposite-direction as equivalent

**Problem for LRRK2**:
- LRRK2 mutation upregulates Gene X → CRISPRi downregulates Gene X
- **Current analysis**: NO overlap detected (different directions)
- **Biological reality**: **STRONG** inverse correlation (opposing mechanisms)

### 2. Missing Biological Context ❌
**Point Mutations vs Loss-of-Function**:
- Point mutations (LRRK2, SNCA variants): May have **opposing** effects to CRISPRi
- Loss-of-function mutations (PINK1, PRKN): May have **similar** effects to CRISPRi
- **Current analysis**: Treats all mutations identically

### 3. Statistical Power Loss ❌
**Consequence**: Missing significant biological relationships
- True biological convergence through opposing mechanisms goes undetected
- False negatives for important cross-method validations
- Reduced power for detecting meaningful signatures

## 📊 Detailed Technical Analysis

### Current Fisher's Test Implementation
**Location**: Multiple locations in `R/signature_analysis.R`

**Current Logic**:
```r
# Simplified representation
mast_significant_genes <- genes[p < 0.05 & abs(log2FC) > 0.25]
crispri_significant_genes <- genes[p < 0.05 & abs(log2FC) > 0.25]
overlap_genes <- intersect(mast_significant_genes, crispri_significant_genes)
fisher_test(overlap_genes, background_genes)
```

**Missing Logic**:
```r
# What SHOULD be implemented
mast_up <- genes[p < 0.05 & log2FC > 0.25]
mast_down <- genes[p < 0.05 & log2FC < -0.25]
crispri_up <- genes[p < 0.05 & log2FC > 0.25]  
crispri_down <- genes[p < 0.05 & log2FC < -0.25]

# Test multiple patterns:
same_direction_overlap <- c(intersect(mast_up, crispri_up), 
                           intersect(mast_down, crispri_down))
opposite_direction_overlap <- c(intersect(mast_up, crispri_down),
                               intersect(mast_down, crispri_up))
```

### Available Data Sources for Direction Analysis
**DE Results Structure**: Available in `full_DE_results.rds`
- **MAST**: `avg_log2FC` column with direction information
- **CRISPRi**: `log2FC_[experiment]` columns with direction information
- **Accessibility**: Already loaded in signature analysis

**Implementation Readiness**: ✅ **HIGH** - all required data and functions exist

## 🔧 **APPROVED IMPLEMENTATION STRATEGY** (Updated January 14, 2025)

### **Unified Direction-Aware Analysis** ✅

**Approach**: Single test that accounts for direction with biological context

**Implementation Strategy**:
```r
# Unified direction-aware Fisher's test with biological expectations
enhanced_direction_analysis <- function(mast_data, crispri_data, gene_name) {
  
  # Determine biological expectation
  direction_expectation <- get_biological_direction_expectation(gene_name)
  
  # Calculate direction-aware overlap statistics
  same_direction_overlap <- calculate_same_direction_overlap(mast_data, crispri_data)
  opposite_direction_overlap <- calculate_opposite_direction_overlap(mast_data, crispri_data)
  
  # Apply biological weighting
  if (direction_expectation == "opposing") {
    # LRRK2 case: prioritize opposite direction test
    primary_result <- opposite_direction_overlap
    secondary_result <- same_direction_overlap
    primary_weight <- 0.8
    secondary_weight <- 0.2
  } else if (direction_expectation == "same") {
    # SNCA variants, loss-of-function mutations: prioritize same direction
    primary_result <- same_direction_overlap
    secondary_result <- opposite_direction_overlap  
    primary_weight <- 0.8
    secondary_weight <- 0.2
  } else {
    # Unknown/mixed: equal weighting
    primary_result <- same_direction_overlap
    secondary_result <- opposite_direction_overlap
    primary_weight <- 0.5
    secondary_weight <- 0.5
  }
  
  # Combined direction-aware p-value using Fisher's method
  combined_p <- combine_direction_pvalues(
    primary_result$fisher_p, secondary_result$fisher_p,
    weights = c(primary_weight, secondary_weight)
  )
  
  return(list(
    combined_p = combined_p,
    same_direction = same_direction_overlap,
    opposite_direction = opposite_direction_overlap,
    biological_expectation = direction_expectation,
    primary_pattern = ifelse(direction_expectation == "opposing", "opposite", "same")
  ))
}

# Biological expectation function
get_biological_direction_expectation <- function(gene_name) {
  case_when(
    gene_name == "LRRK2" ~ "opposing",  # Hyperactive kinase vs knockdown
    gene_name %in% c("SNCA_A30P", "SNCA_A53T", "SNCA") ~ "same",  # Aggregation in both
    gene_name %in% c("PINK1", "PRKN", "PARK7", "ATP13A2", "DNAJC6", "FBXO7") ~ "same",  # Loss-of-function
    TRUE ~ "unknown"  # Default: test both patterns equally
  )
}
```

**Key Features**:
- **Biologically-informed weighting**: LRRK2 prioritizes opposite-direction tests
- **SNCA-appropriate analysis**: Does not force opposite-direction for aggregation variants  
- **Unified statistical framework**: Single test incorporating direction information
- **Pathway analysis**: Direction-aware approach maintained for pathway overlaps

### Priority 2: Enhanced UI Transparency (HIGH)

#### Direction Indicators in Results Tables
```r
# Add columns showing direction patterns
display_data$Direction_Pattern <- c("Same ↑↑", "Same ↓↓", "Opposite ↑↓", "Mixed")
display_data$Expected_Pattern <- get_expected_pattern(mutation_type)
display_data$Biological_Relevance <- assess_biological_relevance(observed, expected)
```

#### Visual Direction Analysis
```r
# Enhanced correlation plots showing direction
create_direction_correlation_plot(
  data = correlation_data,
  highlight_opposing = TRUE,
  mutation_type = mutation_type,
  expected_pattern = expected_pattern
)
```

### Priority 3: Mutation-Aware Analysis Framework (MEDIUM)

#### Biological Context Integration
```r
# Load mutation classifications in analysis
mutation_info <- get_mutation_categories()

# Enhance signature scoring
enhanced_signature_score <- calculate_signature_score(
  overlap_stats = overlap_stats,
  direction_stats = direction_stats,
  mutation_context = mutation_info,
  biological_expectation = derive_expectation(mutation_type)
)
```

## 📈 Expected Improvements

### Scientific Validity
- **Detect Previously Missed Signatures**: LRRK2-type opposing effects
- **Improved Biological Interpretation**: Context-appropriate direction analysis
- **Enhanced Statistical Power**: More sensitive to true biological relationships

### User Experience  
- **Clear Direction Information**: Users understand same vs opposite effects
- **Biological Context**: Mutation type influences interpretation
- **Enhanced Visualizations**: Direction-aware correlation plots

### Manuscript Impact
- **Stronger Cross-Method Validation**: Both convergent and divergent mechanisms
- **Biologically Informed Results**: Mutation-type-specific analysis
- **Novel Biological Insights**: Opposing mechanism discovery

## 🧪 Implementation Strategy

### Phase 1: Core Direction Analysis (Week 1)
1. **Implement Direction-Aware Fisher's Tests**
   - Modify `calculate_gene_overlap_significance_proper()` to include direction options
   - Add `test_same_direction` and `test_opposite_direction` parameters
   - Update signature ranking to include direction results

2. **Enable Direction Stats in Main Analysis**
   - Modify line 332 in `manuscript_signature_discovery.R`:
   ```r
   direction_stats = calculate_direction_consistency(mast_data, crispri_data)
   ```
   - Update composite scoring to use direction information

### Phase 2: UI Integration (Week 2)
3. **Add Direction Information to Tables**
   - Include direction pattern columns in gene pair analysis
   - Add expected vs observed direction comparison
   - Color-code by biological relevance

4. **Enhance Correlation Plots**
   - Highlight opposing vs convergent patterns
   - Show mutation type context
   - Add biological interpretation tooltips

### Phase 3: Biological Context (Week 3)
5. **Mutation-Type-Specific Analysis**
   - Load mutation categories in signature analysis
   - Adjust statistical expectations based on mutation type
   - Provide mutation-specific interpretation

6. **Comprehensive Direction Framework**
   - Standardize direction analysis across all modules
   - Add user controls for direction analysis preferences
   - Create direction-focused visualization suite

## 🔬 Validation Protocol

### Test Cases
1. **LRRK2 Analysis**: Verify detection of opposing effects
2. **PINK1 Analysis**: Verify detection of convergent effects (loss-of-function)
3. **Cross-Mutation Comparison**: Ensure mutation-type differences are captured

### Statistical Validation
1. **Power Analysis**: Compare detection rates with/without direction analysis
2. **False Discovery Rate**: Ensure direction awareness doesn't inflate false positives
3. **Biological Validation**: Confirm results match known biological mechanisms

## 📋 Implementation Files

### Files Requiring Modification
1. **`R/signature_analysis.R`** - Add direction-aware Fisher's tests
2. **`R/manuscript_signature_discovery.R`** - Enable direction stats (line 332)
3. **`inst/shiny/modules/mod_signature_nomination.R`** - Add direction UI elements
4. **`inst/shiny/modules/mod_de_results.R`** - Add direction info to summary stats

### New Functions Needed
1. **`test_same_direction_overlap()`** - Fisher's test for same-direction genes
2. **`test_opposite_direction_overlap()`** - Fisher's test for opposite-direction genes  
3. **`get_expected_direction_pattern()`** - Mutation-type-based expectations
4. **`assess_biological_relevance()`** - Context-appropriate interpretation

### Files Working Correctly
1. **`R/gene_harmonization.R`** - ✅ Mutation classification available
2. **`R/signature_analysis.R`** - ✅ Direction analysis function exists
3. **Effect Size Correlation** - ✅ Already shows direction relationships

## 🎯 Success Metrics

### Immediate (Phase 1)
- [ ] Direction-aware Fisher's tests implemented and functional
- [ ] LRRK2 analysis detects opposing effects (should show significance in opposite direction test)
- [ ] Main signature analysis includes direction statistics

### Short-term (Phase 2)
- [ ] UI displays direction patterns clearly
- [ ] Users can distinguish same vs opposite direction overlaps
- [ ] Correlation plots highlight direction relationships

### Long-term (Phase 3) 
- [ ] Mutation-type-specific analysis framework operational
- [ ] Comprehensive direction analysis across all modules
- [ ] Enhanced biological interpretations for all cross-method comparisons

---

**Last Updated**: January 14, 2025  
**Status**: Critical Direction Analysis Gap Identified  
**Next Steps**: Implement direction-aware Fisher's exact tests and enable direction statistics in main analysis pipeline