# Recent Work Documentation - January 16, 2025

## 🔬 Comprehensive Analysis Completed

### Fisher's Exact Test Analysis Framework
- **Created comprehensive_4_test_fisher_analysis.R** - Implements all 4 requested Fisher's tests:
  1. ALL DE genes with intersection background
  2. ALL DE genes with union background  
  3. Enrichment-associated DE genes with intersection background
  4. Enrichment-associated DE genes with union background

### Key Findings
1. **35 statistically significant cross-method signatures** identified via enrichment term overlap
2. **46 significant results** from comprehensive 4-test DE gene analysis
3. **PARK7/DJ-1** emerges as top candidate across multiple frameworks (p=1.97e-19)
4. **Union background consistently outperforms intersection** for statistical power

### Directionality Investigation
- **Created investigate_directionality_experiments.R**
- **Confirmed user skepticism is justified**: 76,181 pathways still show multiple directions
- **137,103 cases** where same direction appears in multiple experiments
- **"ALL" category continues to include "UP"/"DOWN" from same experiments**

### Scripts Created/Run
1. `investigate_directionality_experiments.R` - Confirms directionality inflation persists
2. `simple_signature_analysis.R` - Found 35 significant enrichment term overlaps  
3. `comprehensive_4_test_fisher_analysis.R` - Complete 4-test Fisher framework

## 📦 Package Installation Issue

### User Reports All Bug Fixes Failed After GitHub Reinstall
Despite successful `remotes::install_github("jessedunnack/iSCORE-PDecipher", force = TRUE)`, all previously "fixed" bugs persist:

- Bug #1: Details modal still shows "no DE genes found"
- Bug #3: "[object Object]" still in cross-cluster analysis
- Bug #4: "found: 0 gene set(s)" in Fisher's test heatmap
- Bug #5: CRISPRi dotplots completely missing
- Bug #6: Heatmap color fix not visible
- Bug #8: Experiment column missing from DE results table  
- Bug #9: Gene dropdown consolidation not working

### Critical Investigation Needed
**Question**: Do GitHub commits to `/inst/shiny/` files actually get installed by `remotes::install_github()`?

**Hypothesis**: Shiny app files may not be properly packaged/installed, requiring different installation approach.

## 🚨 Immediate Actions Required
1. Investigate GitHub package installation mechanism for Shiny apps
2. Verify if changes to `/inst/shiny/` are included in package installation
3. Determine if alternative installation method needed
4. Re-implement bug fixes if installation issue identified

## 📋 Status
- **Analysis work**: Complete and documented
- **Bug fixes**: Implemented but apparently not taking effect
- **User concerns**: Valid - installation may be the root issue
- **Next steps**: Investigate package installation mechanism