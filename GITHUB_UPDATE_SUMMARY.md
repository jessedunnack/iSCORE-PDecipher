# GitHub Update Summary - iSCORE-PDecipher v0.1.3

## Ready to Push!

All changes have been staged and committed. You can now push to GitHub with:

```bash
git push origin main
```

## What's Been Added in This Session:

### 1. UMAP Visualization (14MB of data files)
- **Landing Page Redesign**: UMAP on left, summary stats on right
- **Automatic Dataset Selection**: Based on loaded enrichment data
- **Lightweight Data**: Extracted from 60-90GB Seurat objects → 14MB total
- **Clean Interface**: No unnecessary controls, automatic configuration

### 2. GSEA Visualization Support
- **Basic Visualization Enhanced**: Now supports GSEA plots
- **Multiple Plot Types**: Enrichment plots, NES dot plots, ridge plots
- **Dynamic UI**: Automatically switches when GSEA is selected
- **Publication Ready**: Static plots suitable for papers

### 3. Bug Fixes
- **Heatmap Direction Labels**: Now shows "UP/DOWN/ALL" instead of color codes
- **GSEA Support**: Added "RANKED" direction and GSEA colors

### 4. New Files Structure:
```
iSCORE-PDecipher/
├── DESCRIPTION (updated to v0.1.3 with new dependencies)
├── Documentation/
│   ├── GSEA_VISUALIZATION_SUMMARY.md
│   ├── UMAP_INTEGRATION_GUIDE.md
│   ├── UMAP_LANDING_PAGE_V2_SUMMARY.md
│   └── UMAP_SOLUTION_SUMMARY.md
├── inst/
│   ├── extdata/
│   │   └── umap_data/
│   │       ├── Full_Dataset_umap_data.rds (4.7MB)
│   │       ├── iSCORE_PD_CRISPRi_umap_data.rds (4.6MB)
│   │       ├── iSCORE_PD_umap_data.rds (4.1MB)
│   │       ├── all_umap_data_combined.rds (14MB)
│   │       └── umap_data_summary.csv
│   ├── scripts/
│   │   ├── extract_umap_data.R
│   │   └── test_umap_visualization.R
│   └── shiny/
│       ├── app.R (updated)
│       ├── modules/
│       │   ├── mod_heatmap.R (fixed)
│       │   ├── mod_landing_page_with_umap_v2.R (new)
│       │   ├── mod_umap_viewer.R (new)
│       │   ├── mod_umap_viewer_simple.R (new)
│       │   └── mod_visualization_enhanced.R (new)
│       └── www/
│           └── custom.css (updated)
```

## Dependencies Added:
- SingleCellExperiment
- dittoSeq  
- enrichplot
- ggridges (suggested)

## Version Bump:
- 0.1.2 → 0.1.3

## Commits in This Session:
1. `d3c9dda` - Major update: UMAP visualization and GSEA support (v0.1.3)
2. `41149e1` - Major enhancement: Interactive heatmap module v0.1.3

## Testing Before Push:
Consider testing the app locally one more time:
```r
library(iSCORE.PDecipher)
launch_app()
```

Verify:
- [ ] UMAP displays correctly on landing page
- [ ] GSEA plots work in Basic Visualization
- [ ] Heatmap direction labels show correctly
- [ ] All dependencies install properly

## Push Command:
```bash
git push origin main
```

The repository is ready for the update! All changes are committed and organized properly.