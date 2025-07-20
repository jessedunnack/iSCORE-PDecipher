# Version Status Summary
**Date**: January 20, 2025
**Current Version**: v0.3.3
**Branch**: v0.2.9-fixes

## Current State

We are on the `v0.2.9-fixes` branch, which has diverged from `main`. Our branch contains:

### Major Features Added (v0.3.0 - v0.3.3):
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

## Branch Divergence

- **Our branch (v0.2.9-fixes)**: Contains PD signature analysis features (v0.3.x)
- **Main branch**: At v0.2.8 with different bug fixes

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

The v0.2.9-fixes branch contains all the PD signature analysis work and should be considered the latest version. The main branch has different bug fixes that may need to be evaluated separately.

To consolidate:
1. This branch (v0.2.9-fixes) has the complete analysis pipeline
2. All visualizations are updated and corrected
3. Documentation is comprehensive
4. Version v0.3.3 represents the latest stable state