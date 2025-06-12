# UMAP Landing Page V2 - Summary of Changes

## Overview

I've created a redesigned landing page (`mod_landing_page_with_umap_v2.R`) that matches your mockup with:
- UMAP on the left (60% width)
- Summary statistics boxes on the right (40% width)
- Automatic dataset selection based on loaded data
- No user controls for UMAP (cleaner interface)
- Removed title text as requested

## Key Features

### 1. Automatic Dataset Detection
The UMAP viewer automatically determines which dataset to display based on the loaded enrichment data:
- If both MAST and MixScale data → Shows "Full_Dataset"
- If only MixScale data → Shows "iSCORE_PD_CRISPRi" 
- If only MAST data → Shows "iSCORE_PD"

### 2. Compact Layout
- **Left Column (60%)**: UMAP visualization
  - No settings/controls
  - Clean, maximized plot area
  - Automatic cluster labeling

- **Right Column (40%)**: Summary statistics in colored boxes
  - Total Cells (aqua)
  - Cell Clusters (yellow)
  - Enrichment Results (blue)
  - Genes Analyzed (green)
  - Analysis Methods (purple)
  - Enrichment Types (orange)
  - Dataset Information panel below

### 3. Visual Improvements
- Removed "PD Enrichment Explorer - Interactive Data Overview" title
- Compact value boxes with icons
- Consistent color scheme
- Better use of screen space
- Professional appearance matching AdminLTE style

### 4. Cell Count Integration
Instead of the text output, cell counts are now displayed in attractive colored boxes:
- Total Cells: Shows actual cell count from UMAP data
- Cell Clusters: Shows cluster count from UMAP or enrichment data
- Clean, visual presentation without text clutter

## Implementation

To use this new version:

1. The app.R has been updated to use `mod_landing_page_with_umap_v2.R`
2. Custom CSS has been added for the compact value boxes
3. The UMAP dataset is automatically selected based on your data

## Benefits

1. **Cleaner Interface**: No unnecessary controls or settings
2. **Automatic Configuration**: Dataset selection based on loaded data
3. **Better Space Usage**: Horizontal layout maximizes screen real estate
4. **Professional Look**: Colored boxes match modern dashboard design
5. **Responsive**: Works well on different screen sizes

## How It Works

When the app loads:
1. It detects which enrichment data is loaded (MAST/MixScale)
2. Automatically loads the corresponding UMAP dataset
3. Displays the UMAP with cluster labels
4. Shows summary statistics in the colored boxes
5. No user intervention required

The UMAP will always show the most appropriate dataset for the enrichment results being analyzed, creating a seamless user experience.