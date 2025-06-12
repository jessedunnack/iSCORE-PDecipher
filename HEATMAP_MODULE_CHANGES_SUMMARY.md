# iSCORE-PDecipher Heatmap Module Enhancement Summary

## Version Update: 0.1.2 → 0.1.3

### Overview
The heatmap module (`mod_heatmap.R`) has been completely overhauled to provide interactive, publication-quality visualizations with advanced filtering and export capabilities.

## Major Enhancements

### 1. Interactive Heatmaps with heatmaply
- **Before**: Static heatmaps with limited interactivity
- **After**: Fully interactive heatmaps with:
  - Hierarchical clustering with visible dendrograms
  - Zoom and pan functionality
  - Hover tooltips showing full term descriptions
  - Row and column reordering based on clustering

### 2. Data Type Support
Added support for multiple data visualization types:
- **P-values**: -log10 transformation for better visualization
- **Fold Enrichment**: Direct visualization of enrichment ratios
- **Z-scores**: For normalized comparisons
- **Gene Counts**: Raw count data
- **GSEA NES**: Normalized Enrichment Scores with diverging color scales

### 3. Enhanced Filtering Options

#### Direction Filtering
- **ALL_DIR**: All gene directions
- **UP**: Up-regulated genes only
- **DOWN**: Down-regulated genes only  
- **BOTH**: Both UP and DOWN (excluding ALL)

#### Method Filtering
- All methods
- MAST only
- MixScale only
- Intersection (terms found in both methods)
- Union (terms from either method)

### 4. Color Scale Options
Multiple color scales to suit different data types:
- **Red (0 to max)**: Sequential scale for positive values
- **Blue (0 to max)**: Alternative sequential scale
- **Red-Blue (diverging)**: For data with positive/negative values
- **Viridis**: Perceptually uniform scale
- **Yellow-Orange-Red**: Traditional heatmap colors

### 5. Color Scaling Methods
Three approaches to handle different data distributions:
- **Linear**: Standard linear mapping
- **Quantile (adaptive)**: Uses quantile breaks for skewed data
- **Fixed breaks**: User-defined breakpoints

### 6. GSEA-Specific Features
Conditional UI panel for GSEA analysis:
- Toggle to show GSEA results only
- NES threshold slider (0-3)
- Automatic detection of NES columns
- Diverging color scales for positive/negative enrichment

### 7. PDF Export Functionality
New "PDF Export" tab with ComplexHeatmap integration:
- Customizable dimensions (width/height in inches)
- Adjustable font sizes
- DPI settings for publication quality
- Maintains same clustering and annotations as interactive view

### 8. Row Annotations
Visual indicators for:
- **Enrichment Type**: GO_BP, GO_CC, GO_MF, KEGG, Reactome, etc.
- **Direction**: UP (red), DOWN (blue), ALL (purple)

### 9. Data Integration Fixes
- Fixed hardcoded file path issue
- Now uses `app_data$consolidated_data` for dynamic data loading
- Handles both "method" and "analysis_type" column names
- Graceful handling of missing data

### 10. Download Options
Two download formats:
- **HTML**: Self-contained interactive heatmap
- **PDF**: Static publication-quality figure

## Technical Implementation Details

### Dependencies Added
```r
- heatmaply      # Interactive heatmaps
- ComplexHeatmap # Static PDF generation
- colourpicker   # Color selection UI
```

### Key Functions Modified
- `mod_heatmap_ui()`: Completely redesigned UI with new controls
- `mod_heatmap_server()`: Rewritten server logic for data processing
- Added PDF generation functionality
- Implemented adaptive data pivoting for matrix creation

### Error Handling Improvements
- Robust error messages with helpful feedback
- Graceful degradation when packages unavailable
- Validation of data before visualization
- Prevention of empty heatmap generation

## Usage Examples

### Basic Usage
1. Click "Generate Heatmap" to create visualization
2. Use sidebar controls to filter and customize
3. Interact with heatmap (zoom, hover, reorder)
4. Export as HTML or PDF

### Advanced Features
- Combine enrichment types for comprehensive view
- Use intersection/union methods for cross-method analysis
- Apply quantile scaling for datasets with outliers
- Generate PDFs with custom dimensions for publications

## Benefits
- **Improved Data Exploration**: Interactive features enable deeper insights
- **Publication Ready**: PDF export with customizable settings
- **Flexible Visualization**: Multiple data types and scaling options
- **Better Performance**: Efficient data handling and rendering
- **User Friendly**: Intuitive controls and clear feedback

## Future Enhancements (Suggested)
- Save/load heatmap configurations
- Batch export functionality
- Additional clustering algorithms
- Integration with other visualization modules
- Custom color palette creation

---
Last updated: January 6, 2025