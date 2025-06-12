# GSEA Visualization Implementation Summary

## Overview

I've enhanced the Basic Visualization module to include GSEA-specific visualization capabilities. The module now automatically detects when GSEA data is selected and provides specialized plot types.

## Key Features

### 1. Automatic GSEA Detection
- The module automatically detects when the selected enrichment type is "GSEA"
- UI dynamically switches between standard enrichment controls and GSEA-specific options

### 2. GSEA Plot Types

#### A. Enrichment Plot (gseaplot2-style)
- Shows the classic GSEA enrichment plot
- Displays:
  - Gene set name
  - NES (Normalized Enrichment Score)
  - Adjusted p-value
  - Note: Full gseaplot2 functionality requires the original GSEA object, so this provides a summary view

#### B. GSEA Dot Plot
- Similar to standard dot plots but optimized for GSEA
- X-axis: Normalized Enrichment Score (NES)
- Y-axis: Gene sets ordered by NES
- Size: Set size (number of genes)
- Color: Adjusted p-value
- Vertical line at NES = 0 to show up/down regulation

#### C. Ridge Plot
- Shows distribution of enrichment scores
- Useful for comparing multiple gene sets
- Colored by adjusted p-value
- Requires optional ggridges package

#### D. GSEA Table View
- Detailed tabular view of GSEA results
- Sortable by NES, p-value, set size

### 3. GSEA-Specific Controls
- **Gene Set Selection**: Dropdown to select specific gene set for enrichment plot
- **Show p-value**: Toggle p-value display on plots
- **Show running score**: Toggle running score display
- **Base font size**: Adjust text size for readability
- **Top terms**: Control number of gene sets shown in multi-set plots

### 4. Standard Enrichment Features (Still Available)
- Dot plots, bar plots, and lollipop plots
- Interactive (plotly) and static versions
- Customizable axes and colors
- Download functionality

## How It Works

### For Standard Enrichment Data:
1. Select any non-GSEA enrichment type (GO, KEGG, etc.)
2. Use standard plot types (dot, bar, lollipop)
3. Interactive plotly version is the default

### For GSEA Data:
1. Select "GSEA" as the enrichment type in global settings
2. UI automatically switches to GSEA-specific options
3. Choose from:
   - **Enrichment Plot**: Individual gene set visualization
   - **Dot Plot**: Compare multiple gene sets by NES
   - **Ridge Plot**: Distribution visualization
   - **Table**: Detailed data view

## Technical Implementation

### Dependencies Added:
- `enrichplot`: For GSEA visualization functions
- `ggridges`: (Suggested) For ridge plots

### Module Structure:
- Conditional UI panels based on enrichment type
- Separate rendering functions for GSEA vs standard plots
- Automatic data detection and UI updates

## Usage Example

```r
# In the Shiny app:
1. Go to "Basic Visualization" tab
2. In global settings, select:
   - Analysis Type: [Your choice]
   - Enrichment Database: GSEA
3. The visualization panel will automatically show GSEA options
4. Select "Enrichment Plot" to see individual gene set plots
5. Or select "Dot Plot" to compare top gene sets by NES
```

## Limitations

1. **Full gseaplot2 functionality**: The complete gseaplot2 visualization requires the original GSEA result object with gene-level statistics. Since we only have summarized data, the enrichment plot shows a simplified version.

2. **Interactive GSEA plots**: Currently, GSEA plots are static (using ggplot2). While dot plots can be made interactive with plotly, the enrichment plots remain static for accuracy.

## Benefits

1. **Seamless Integration**: GSEA plots appear alongside other enrichment visualizations
2. **Automatic Detection**: No need to navigate to a separate module
3. **Consistent Interface**: Uses the same global settings as other visualizations
4. **Multiple Views**: Different plot types for different analytical needs
5. **Publication Ready**: High-quality static plots suitable for publications

The enhanced visualization module now provides comprehensive support for both standard enrichment results and GSEA-specific visualizations, making it a one-stop solution for all enrichment visualization needs.