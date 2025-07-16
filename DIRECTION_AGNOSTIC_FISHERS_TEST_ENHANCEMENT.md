# 🔄 ENHANCEMENT PROPOSAL: Direction-Agnostic Fisher's Tests

**Date:** January 16, 2025  
**Status:** 📋 PROPOSED - For Investigation After Current Bug Fixes  
**Priority:** HIGH - Addresses Core User Concern About Missing Signatures

## 🎯 **Core Motivation**

After fixing the directionality inflation issue in Fisher's exact tests, the user expressed concern:

> **"I am worried that we will no longer pick up significant signatures... we MUST find significant signatures that are elicited by both experimental approaches!"**

This enhancement would address that concern by offering **both direction-aware and direction-agnostic Fisher's tests**.

## 🔬 **Scientific Rationale**

### Current Approach (Direction-Aware)
- **UP only**: Tests genes upregulated in both conditions
- **DOWN only**: Tests genes downregulated in both conditions  
- **ALL (deduplicated)**: Tests all DE genes but prevents double-counting

### Proposed Addition (Direction-Agnostic)
- **Any direction**: Tests overlap of DE genes regardless of their direction
- **Biological value**: Captures "antagonistic" effects where both approaches affect the same genes but in opposite directions

## 🔍 **Key Biological Scenarios**

### Example: LRRK2 Analysis
- **LRRK2 mutation** (gain-of-function): Upregulates Gene X
- **LRRK2 CRISPRi** (loss-of-function): Downregulates Gene X
- **Current direction-aware test**: Would NOT detect this overlap
- **Proposed direction-agnostic test**: Would detect this as significant overlap

### Why This Matters
Both perturbations affect the same biological pathway/gene but in opposite directions - this is **biologically meaningful convergence** that we currently miss.

## 🛠️ **Implementation Strategy**

### 1. **UI Enhancement**
Add a **Direction Mode** selector with two options:

```
Direction Mode: 
○ Direction-aware (default) - Conservative, prevents inflation
○ Direction-agnostic - Broader detection, includes opposite effects

[If Direction-aware selected:]
Direction Filter:
○ All directions (deduplicated)  
○ Up-regulated only
○ Down-regulated only
```

### 2. **Statistical Implementation**

#### Direction-Aware (Current)
```r
# UP only
if (selected_direction == "UP") {
  genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC > 0, ]$gene_name
}
```

#### Direction-Agnostic (New)
```r
# All DE genes regardless of direction (but deduplicated)
if (direction_mode == "agnostic") {
  genes <- unique(gene_subset[gene_subset$p_val_adj < 0.05, ]$gene_name)
}
```

### 3. **Implementation Locations**

#### A. Signature Nomination Module
- **File**: `inst/shiny/modules/mod_signature_nomination.R`
- **Function**: Update Fisher's test calculation
- **Impact**: Could detect more significant gene pairs

#### B. DE Gene Overlap Heatmap  
- **File**: `inst/shiny/modules/mod_de_results.R` 
- **Function**: `output$de_overlap_heatmap`
- **Impact**: Could show more significant overlaps in heatmap

#### C. Signature Analysis Core
- **File**: `R/signature_analysis.R`
- **Function**: `analyze_gene_pair_signatures()`
- **Impact**: Core logic for cross-method comparison

## 📊 **Expected Impact**

### Advantages
1. **Addresses user concern** about missing signatures after directionality fix
2. **Captures biological antagonism** (gain-of-function vs loss-of-function effects)
3. **Maintains statistical rigor** by offering both conservative and exploratory options
4. **Default to conservative** approach (direction-aware) while enabling broader exploration

### Potential Results
- **More significant overlaps** detected in direction-agnostic mode
- **Better capture of LRRK2-type scenarios** where mutation and knockdown have opposite effects
- **Enhanced biological interpretation** of cross-method convergence

## ⚠️ **Important Considerations**

### 1. **Not the Same as Original Inflation Issue**
- **Original problem**: Same gene counted multiple times in UP+DOWN+ALL enrichment categories
- **This approach**: Same gene counted once, but direction doesn't matter

### 2. **Statistical Validity**
- Direction-agnostic tests are **scientifically valid**
- They answer a different biological question: "Do both approaches affect the same genes?"
- vs direction-aware: "Do both approaches affect the same genes in the same direction?"

### 3. **Default Behavior**
- **Direction-aware should remain the default** (conservative)
- **Direction-agnostic as optional** for exploratory analysis

## 🎯 **Implementation Plan**

### Phase 1: Core Logic (High Priority)
1. **Update signature analysis functions** to support direction mode parameter
2. **Modify Fisher's test calculations** for both modes
3. **Test with known gene pairs** to validate approach

### Phase 2: UI Integration (Medium Priority)  
1. **Add direction mode selectors** to signature nomination module
2. **Add direction mode selectors** to DE overlap heatmap
3. **Update help text** to explain the difference

### Phase 3: Validation (High Priority)
1. **Test with established positive controls** (known convergent pathways)
2. **Compare results** between direction-aware and direction-agnostic approaches
3. **Validate that direction-agnostic doesn't reintroduce inflation**

## 📋 **Success Metrics**

1. **✅ User Satisfaction**: Address concern about missing signatures
2. **✅ Scientific Rigor**: Maintain statistical accuracy in both modes
3. **✅ Biological Insight**: Capture meaningful antagonistic effects
4. **✅ Usability**: Clear UI distinguishing between the two approaches

## 🔗 **Related Issues**

- **Bug #2**: Original directionality inflation fix
- **Bug #4**: DE gene overlap heatmap (just implemented)
- **User Concern**: [CRITICAL_CONCERN_SIGNATURE_DETECTION.md](CRITICAL_CONCERN_SIGNATURE_DETECTION.md)

---

## 🚨 **IMMEDIATE NEXT STEPS**

1. **Complete current bug fixes** (Bugs #6-#9)
2. **Implement this enhancement** as priority follow-up
3. **Test thoroughly** with known positive/negative controls
4. **Validate user satisfaction** with signature detection capability

**This enhancement could be the key to maintaining both statistical rigor AND biological discovery power.**