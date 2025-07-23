# Version Status Summary
**Date**: July 22, 2025
**Current Version**: v0.3.5
**Branch**: main

## Current State

We are on the `main` branch with the latest stable version. Current features include:

### Major Features Added (v0.3.0 - v0.3.5):
1. **v0.3.0**: Comprehensive PD signature analysis suite
   - Full analysis of all 16 PD genes and 15 clusters
   - Created visualization pipeline for thesis committee

2. **v0.3.1**: Gene harmonization and visualization fixes
   - Fixed PRKN→PARK2 mapping
   - Added hierarchical clustering to heatmaps
   - Natural sorting for clusters

3. **v0.3.2**: Label readability improvements
   - Word wrapping for all plot labels
   - Horizontal labels with rounded backgrounds
   - Enhanced readability across all visualizations

4. **v0.3.3**: Convergence plot and pathway count fixes
   - Fixed pathway count bug (showing actual totals: 826/531/934)
   - Improved label visibility in convergence strength plot
   - Added metric recommendations documentation

5. **v0.3.4**: Heatmap UI/UX enhancements (July 20, 2025)
   - All collapsible panels start expanded by default
   - Biological pathway categories enabled by default
   - Auto color scale selection based on gene direction (UP→red, DOWN→blue)
   - Improved workflow with fewer clicks

6. **v0.3.5**: Enrichment type selection moved to visualization module (July 20, 2025)
   - Moved enrichment database dropdown from global settings to visualization module
   - Users can switch between GO_BP, GO_CC, GO_MF, KEGG, Reactome, WikiPathways, STRING, GSEA directly in visualizations
   - Dynamic dropdown updates based on available data
   - Better workflow - no need to leave visualization page

## Current Status

The main branch is now the default with all stable features from the stable-heatmap branch merged. Repository has been reorganized with proper folder structure.

## All Changes Since v0.2.9

### Analysis Scripts Added/Modified:
- `pd_signature_discovery_fast.R` - Core analysis with pathway totals
- `create_comprehensive_visualizations.R` - 12-panel figure suite
- `create_quick_plots.R` - Individual method plots
- `create_thesis_committee_summary.R` - Summary figure
- `pd_signature_comprehensive_viz.R` - Comprehensive analysis
- `pd_signature_by_gene_analysis.R` - Gene-specific analysis
- `pd_signature_cluster_analysis.R` - Cluster-specific analysis

### Key Fixes Applied:
1. Gene harmonization (PRKN→PARK2, VPS13C variants, SNCA variants)
2. Hierarchical clustering on heatmaps
3. Natural cluster sorting
4. Word wrapping for labels
5. Pathway count correction (actual totals vs top 30)
6. Label visibility improvements

### Documentation Created:
- `COMPREHENSIVE_FINDINGS_SUMMARY.md` - Key findings
- `PATHWAY_COUNT_FIX_SUMMARY.md` - Count bug documentation
- `KEY_INSIGHTS_REVISED.md` - Updated insights
- `METRIC_RECOMMENDATIONS.md` - Enrichment metric analysis
- `FIX_SUMMARY.md` - All fixes applied

## Recommendation

The main branch is now the default and contains all features. Users should install directly from GitHub:

```r
remotes::install_github("jessedunnack/iSCORE-PDecipher")
```

Current stable version v0.3.5 includes:
1. Complete PD signature analysis pipeline
2. All visualization improvements and fixes
3. Enhanced UI/UX for better workflow
4. Comprehensive documentation
5. Repository reorganization with proper folder structure