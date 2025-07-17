# Comprehensive Bug Fix Analysis for iSCORE-PDecipher
**Date**: 2025-07-17  
**Purpose**: Document all bug fix attempts, failures, and minimal fixes needed  
**Baseline Commit**: d6c08c0 (2025-07-15 13:03) - "fix: Resolve 'subscript out of bounds' error"  
**First Bug Fix Attempt**: 4b10a6b (2025-07-16 13:05) - "feat: Comprehensive bug fixes...9 bugs resolved"
**Current Progress**: v0.2.9-fixes branch - Implementing fixes from baseline

## Executive Summary

This document provides a comprehensive analysis of bug fix attempts in the iSCORE-PDecipher Shiny app. The original attempt to fix 9 bugs on 2025-07-16 introduced numerous tangential issues due to incomplete implementations, wrong assumptions about data structures, and fixes applied in incorrect locations. This analysis identifies the core fixes needed and recommends a minimal intervention approach.

## Critical Distinction: Original vs Tangential Bugs

**Original Bugs**: The 9 bugs listed in commit 4b10a6b that were the initial targets for fixing.

**Tangential Bugs**: New issues that arose from failed attempts to fix the original bugs. These are what the user reported as "NONE OF THE THINGS ARE FIXED!!" and are the bugs referenced in subsequent debugging sessions.

### Mapping of Original → Tangential Issues:
- Original Bug #1 (Details modal) → Tangential Bug #1 (modal shows count only)
- Original Bug #3 (Cross-cluster reporting) → Tangential Bug #3 ([object Object] errors)
- Original Bug #5 (Dotplot hover) → Tangential Bug #5 (CRISPRi dotplots not rendering)
- Original Bug #6 (P-value viz) → Tangential Bug #6 (heatmap not interactive)
- Original Bug #8 (Experiment metadata) → Tangential Bug #8 (CRISPRa not excluded)
- Original Bug #9 (Variant consolidation) → Tangential Bug #9 (still shows "X variants")

## Original 9 Bugs (from commit 4b10a6b)

### Bug #1: Gene-pair Analysis Details Modal
**Original Intent**: Add a modal dialog showing shared genes and pathways when users click for details in gene pair analysis

**What Was Attempted**:
- Added `extract_shared_genes()` and `extract_shared_pathways()` functions in `mod_signature_nomination.R`
- These functions were supposed to be called by a modal dialog

**What Went Wrong**:
- **CRITICAL**: No modal dialog UI was actually implemented
- The extraction functions were added but never connected to any UI element
- User feedback: "Details modal showed generic message instead of actual genes"

**Tangential Issues Created**:
- Later attempts to fix this introduced early returns that prevented gene extraction
- Multiple revisions trying to fix "enrichment_data not found" errors
- Functions were looking for wrong data structures

**Root Cause**: Incomplete implementation - backend functions added without frontend UI

---

### Bug #2: Fisher's Test Directionality Inflation (14.158x)
**Original Intent**: Fix critical statistical bug where genes were counted multiple times across UP/DOWN/ALL categories

**What Was Attempted**:
- Added direction filtering logic in Fisher's exact test calculations
- Implemented proper background gene universe calculations

**What Went Wrong**:
- This fix appears to have worked correctly based on code analysis
- However, user was concerned about losing significant signatures

**Tangential Issues Created**:
- Uncertainty about whether to use direction-aware vs direction-agnostic tests
- Questions about proper background gene selection (intersection vs union)

**Root Cause**: Fix was correct but lacked clear documentation of statistical choices

---

### Bug #3: Cross-cluster DE Gene Reporting with Heatmap
**Original Intent**: Add comprehensive cross-cluster DE gene analysis with heatmap visualization

**What Was Attempted**:
- Added new section in `mod_de_results.R` with cross-cluster analysis
- Implemented heatmap visualization
- Added summary statistics table

**What Went Wrong**:
- **[object Object] errors**: JavaScript serialization issues in data tables
- Column name mismatches when 'experiment' column was included
- Data type conversion issues between R and JavaScript

**Tangential Issues Created**:
- Multiple attempts to fix column names
- Confusion about where the error was occurring (user: "not under signature nomination tab!!")
- Wrong fixes applied in wrong modules

**Root Cause**: Improper handling of data serialization for JavaScript/DT tables

---

### Bug #4: Fisher's Test Heatmap for ALL DE Genes
**Original Intent**: Create heatmap showing Fisher's test results for all DE genes with direction filtering

**What Was Attempted**:
- Added heatmap generation logic
- Implemented direction filtering options

**What Went Wrong**:
- Wrong column names used (`p_val_adj` instead of `pvalue`)
- Data structure assumptions didn't match actual data

**Tangential Issues Created**:
- Heatmap generation failures
- Missing data due to column name mismatches

**Root Cause**: Incorrect assumptions about data structure and column names

---

### Bug #5: Dotplot Hover Functionality Enhancement
**Original Intent**: Improve hover information on dotplots with progressive gene lookup fallbacks

**What Was Attempted**:
- Added 4-level progressive gene lookup strategy
- Enhanced hover text generation

**What Went Wrong**:
- CRISPRi dotplots stopped rendering entirely
- Missing 'modality' column in enrichment data
- Filtering logic too restrictive

**Tangential Issues Created**:
- Multiple debug logging additions
- Experiment filter issues with 'default' values
- Complete loss of CRISPRi visualization functionality

**Root Cause**: Over-engineered solution broke existing functionality

---

### Bug #6: Fisher's P-value Visualization
**Original Intent**: Improve p-value visualization with significance thresholds

**What Was Attempted**:
- Enhanced color scales for p-value visualization
- Added discrete significance bins

**What Went Wrong**:
- Initially appears to have worked
- Later user wanted heatmap to be interactive (not part of original bug)

**Tangential Issues Created**:
- Confusion between static and interactive heatmap requirements
- Multiple attempts to convert from renderPlot to renderPlotly

**Root Cause**: Scope creep - original fix worked but requirements expanded

---

### Bug #7: Gene Table Search Restriction
**Original Intent**: Restrict DataTable search to gene names column only

**What Was Attempted**:
- Set `searchable = TRUE` only for gene column
- Set `searchable = FALSE` for other columns

**What Went Wrong**:
- This fix appears to have worked correctly

**Tangential Issues Created**:
- None identified

**Root Cause**: N/A - Fix was successful

---

### Bug #8: Experiment ID Metadata for CRISPRi
**Original Intent**: Add experiment metadata to CRISPRi DE results

**What Was Attempted**:
- Added experiment column extraction
- Displayed experiment info in results

**What Went Wrong**:
- Initially worked but got confused with CRISPRa exclusion
- User wanted complete exclusion of CRISPRa data

**Tangential Issues Created**:
- CRISPRa vs CRISPRi filtering confusion
- Multiple attempts to exclude CRISPRa data

**Root Cause**: Conflation of two different requirements

---

### Bug #9: Consolidated Mutation Variants in Dropdown
**Original Intent**: Show "VPS13C" instead of "VPS13C_W395C, VPS13C_A444P" in dropdowns

**What Was Attempted**:
- Added variant consolidation logic
- Grouped variants by base gene name

**What Went Wrong**:
- Still showed "PRKN (2 variants)" in dropdown
- Fix was applied in wrong location initially

**Tangential Issues Created**:
- "All Datasets" option added that broke filtering
- Multiple attempts to fix in different modules
- User frustration: "you are completely wrong about the bug 9 fix!!"

**Root Cause**: Fix applied in wrong module (global settings vs DE analysis)

---

## Tangential Issues Summary

### Issues Created by Failed Fixes:
1. **"All Datasets" Option**: Added to global settings, broke data filtering
2. **Missing Modal Implementation**: Bug #1 backend without frontend
3. **Column Name Mismatches**: Assumptions about data structure
4. **Over-filtering**: CRISPRi data completely filtered out
5. **Wrong Module Fixes**: Fixes applied in incorrect locations

### User-Confirmed Status of Tangential Issues:

**CONFIRMED WORKING**:
- Bug #8 (Tangential): CRISPRa exclusion - "the only bug you fixed properly"
- Bug #6 (Tangential): DE Gene Overlap Heatmap made interactive (worked after conversion to plotly)
- Bug #9 (Tangential): Gene dropdown cleaned up (after fixing in correct location)

**STILL BROKEN / UNVERIFIED**:
- Bug #1 (Tangential): Gene pair details modal - "still only shows the count! not the actual shared DE genes!"
- Bug #3 (Tangential): [object Object] error - NOT CONFIRMED FIXED (user hasn't verified)
- Bug #5 (Tangential): CRISPRi dotplots - Not rendering due to missing modality column

**Important Note**: These tangential bugs arose from attempting to fix the original 9 bugs. The original bugs themselves may have had different outcomes:
- Original Bug #2 (Directionality): Appears correctly implemented but user concerned about losing signatures
- Original Bug #7 (Search restriction): Likely working but not explicitly confirmed

---

## Minimal Fix Recommendations

Based on this analysis and user feedback, here are the minimal fixes needed for the legacy code (commit d6c08c0):

### CRITICAL - Still Broken (Must Fix):

1. **Bug #1 (Tangential) - Gene Details Modal**: ✅ FIXED IN v0.2.9
   - Problem: Modal only shows count, not actual genes
   - Fix: Remove early returns in `extract_shared_genes()` that prevent gene extraction
   - Ensure modal actually displays the extracted gene lists
   - Location: `mod_signature_nomination.R`
   - **UPDATE (2025-07-17)**: Successfully fixed by:
     - Correcting method names from "iSCORE_PD_MAST"/"CRISPRi_Mixscale" to "MAST"/"MixScale"
     - Using enrichment data with geneID column instead of DE results
     - Adding Entrez ID to gene symbol conversion with clusterProfiler
     - Implementing data manager fallback when app_data not available
     - Changing default to "union" approach and 100 genes display

2. **Bug #3 (Tangential) - [object Object] Error**:
   - Problem: JavaScript serialization showing [object Object] in tables
   - Status: NOT VERIFIED - User hasn't confirmed if fixed
   - Fix: Ensure proper column name handling and data serialization
   - Location: `mod_de_results.R`

3. **Bug #5 (Tangential) - CRISPRi Dotplot Rendering**:
   - Problem: Dotplots not rendering due to missing modality column
   - Fix: Check if modality column exists before filtering, or skip modality filter entirely
   - Do NOT over-engineer the solution
   - Location: `global.R` filter function and `mod_visualization_enhanced.R`

### ALREADY WORKING (Keep These Fixes):

4. **Bug #8 (Tangential) - CRISPRa Exclusion**:
   - Status: WORKING - User confirmed this is properly fixed
   - Keep the existing implementation

### STILL NEED TO VERIFY:

5. **Bug #6 (Tangential) - Interactive Heatmap**: ⚠️ NOT YET CONFIRMED
   - Status: NEEDS VERIFICATION - User hasn't confirmed if working
   - Need to convert from static to interactive heatmap
   - Keep the renderPlotly implementation attempt

6. **Bug #9 (Tangential) - Gene Dropdown**: ⚠️ NOT YET CONFIRMED
   - Status: NEEDS VERIFICATION - User hasn't confirmed if working
   - Should remove "X variants" text from dropdown
   - Keep the fix in `mod_de_analysis.R`

### ADDITIONAL FIXES IMPLEMENTED IN v0.2.9:

7. **DE Genes Table Search Restriction**: ✅ FIXED
   - Problem: Search was finding matches in all columns (e.g., searching "PINK1" found matches in Target column)
   - Fix: Used DataTables' built-in `searchable = FALSE` for columns 1-7
   - Location: `mod_de_analysis.R`
   - **UPDATE (2025-07-17)**: Successfully fixed by:
     - Setting `list(searchable = FALSE, targets = 1:7)` in columnDefs
     - This restricts search to column 0 (Gene) only
     - Confirmed working by user

### ORIGINAL BUGS - Evaluate Carefully:

8. **Original Bug #2 - Fisher's Test Directionality**:
   - The fix appears correct but user is concerned about losing signatures
   - Consider offering both direction-aware (default) and direction-agnostic options
   - Add clear documentation about the statistical implications

9. **Original Bugs #4, #7**:
   - These may be working but need verification
   - Lower priority unless user reports issues

### Implementation Strategy:

1. Start with the three CRITICAL fixes (Bugs #1, #3, and #5)
2. Preserve all WORKING fixes (#6, #8, #9)
3. Test thoroughly after each change
4. Document the statistical choices for Original Bug #2

---

## Lessons Learned

1. **Test Incrementally**: Each bug should be fixed and tested independently
2. **Understand Data Structures**: Verify column names and data types before implementing
3. **Fix in Correct Location**: Ensure fixes are applied in the right module
4. **Complete Implementations**: Don't add backend without frontend
5. **Avoid Over-Engineering**: Simple fixes are often better than complex solutions
6. **Document Statistical Choices**: Especially for analysis methods

---

## Next Steps

1. Revert to commit d6c08c0 (baseline)
2. Apply minimal fixes one at a time
3. Test each fix thoroughly before proceeding
4. Avoid scope creep and tangential changes
5. Document each fix clearly in commit messages

---

## Appendix: Specific Code Fixes Needed

### Fix for Bug #1 (Gene Details Modal)

In `mod_signature_nomination.R`, the `extract_shared_genes()` function has early returns that prevent gene extraction:

```r
# WRONG - Current code with early returns:
if (!is.na(overlap_count) && overlap_count > 0) {
  cat("[DETAILS DEBUG] Found overlap count:", overlap_count, "\n")
  return(list(
    mast_genes = character(0),
    mixscale_genes = character(0),
    shared_genes = character(0)
  ))
}

# CORRECT - Remove early returns:
if (!is.na(overlap_count) && overlap_count > 0) {
  cat("[DETAILS DEBUG] Found overlap count:", overlap_count, "but continuing to extract actual genes\n")
  # Don't return early - continue to extract actual gene lists
}
```

### Fix for Bug #3 ([object Object] Error)

In `mod_de_results.R`, ensure proper column handling for the Cross-Cluster DE Gene Analysis:

```r
# WRONG - Mismatch between expected and actual columns:
col_names <- c("Gene", "# Clusters", "Mean log2FC", "Min p-value", "Clusters")
# But data might have 'experiment' column

# CORRECT - Check columns and set names accordingly:
col_names <- c("Gene", "# Clusters", "Mean log2FC", "Min p-value", "Clusters")
if ("experiment" %in% colnames(display_data)) {
  col_names <- c(col_names, "Experiment")
}
colnames(display_data) <- col_names
```

### Fix for Bug #5 (CRISPRi Dotplots)

In `global.R` filter function, check for modality column existence:

```r
# WRONG - Assumes modality column exists:
if (modality == "CRISPRi") {
  filtered_data <- filtered_data[filtered_data$modality == "CRISPRi", ]
}

# CORRECT - Check column existence first:
if (modality == "CRISPRi" && "modality" %in% colnames(filtered_data)) {
  filtered_data <- filtered_data[filtered_data$modality == "CRISPRi", ]
} else if (modality == "CRISPRi") {
  # If modality column missing, use analysis_type to infer
  cat("[FILTER DEBUG] Modality column missing, filtering by analysis type\n")
  filtered_data <- filtered_data[filtered_data$analysis_type == "MixScale", ]
}
```