# Enhanced Volcano Plot Implementation

## Overview
Implemented dual-mode volcano plots in the "DE Genes > UMAP & Volcano Plots" panel with radio button selection between Static (Publication) and Interactive modes.

## User Interface Changes

### Plot Type Selection
- **New Radio Button**: "Plot Type" with options:
  - "Static (Publication)" - **DEFAULT** selection
  - "Interactive" - Original plotly-based plots

### Conditional Controls
- **Color Options**: Only visible when "Interactive" mode is selected
- **Settings Hidden**: Color controls automatically hidden for static plots (not applicable)

## Static Plot Features (EnhancedVolcano Package)

### Publication Quality
- **Gene Labeling**: Top 15-20 most significant genes with white background boxes
- **Font Styling**: Bold, readable gene labels with connectors
- **Professional Layout**: Clean borders, grids, and publication-ready aesthetics
- **High Resolution**: Vector-based plots suitable for manuscripts

### Descriptive Titles
- **MAST Data**: `"[Gene] vs batch-specific eWT (MAST) - Cluster [X]"`
- **MixScale Data**: `"[Gene] CRISPRi knockdown (MixScale - [Experiment]) - Cluster [X]"`
- **Automatic Context**: Titles update based on selected gene, cluster, and experiment

### Intelligent Gene Selection
- **Primary**: Top 15-20 significant genes (p < 0.05 & |log2FC| > 1) ranked by combined significance and fold change
- **Fallback**: If no significant genes, shows top genes by absolute fold change
- **Smart Ranking**: Uses rank score = |log2FC| × (-log10(p-value)) for optimal label selection

## Technical Implementation

### Package Management
- **Automatic Installation**: Attempts to install EnhancedVolcano via BiocManager if missing
- **Graceful Fallback**: Shows installation instructions if package unavailable
- **Conditional Loading**: Only loads EnhancedVolcano when available

### Rendering Architecture
```
UI Container (uiOutput) 
    ↓
renderUI() switches between:
    ├── Static: plotOutput() → renderPlot() → generate_static_volcano_plot()
    └── Interactive: plotlyOutput() → renderPlotly() → generate_volcano_plot()
```

### Data Processing
- **Gene Name Handling**: Robust extraction of gene names from various data formats
- **Experiment Detection**: Automatic experiment information extraction for titles
- **Error Handling**: Comprehensive error catching with informative fallback plots

## Enhanced Features

### Static Plot Customization
```r
EnhancedVolcano(
  # Gene selection and labeling
  selectLab = genes_to_label,           # Top 15-20 genes
  boxedLabels = TRUE,                   # White background boxes
  drawConnectors = TRUE,                # Clean connectors
  
  # Publication styling
  labSize = 3.5,                       # Readable font size
  labCol = 'black',                    # High contrast
  labFace = 'bold',                    # Bold text
  
  # Color scheme
  col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
  
  # Layout
  gridlines.major = TRUE,              # Professional grid
  border = 'partial',                  # Clean borders
  axisLabSize = 12,                    # Readable axes
  titleLabSize = 14                    # Prominent title
)
```

### Title Generation Logic
- **Context-Aware**: Combines gene, analysis type, experiment, and cluster information
- **Flexible**: Handles missing information gracefully
- **Descriptive**: Clear indication of comparison type and methodology

## Files Modified

### `/inst/shiny/modules/mod_de_results.R`
1. **Package Loading**: Added EnhancedVolcano conditional import
2. **UI Enhancement**: Added plot type radio button with conditional controls
3. **New Function**: `generate_static_volcano_plot()` with full EnhancedVolcano integration
4. **Container Rendering**: Added `renderUI()` functions for dynamic plot type switching
5. **Plot Rendering**: Added separate static and interactive rendering functions
6. **Error Handling**: Comprehensive error management for both plot types

## User Experience

### Default Behavior
1. **App Launch**: Volcano plots start in "Static (Publication)" mode
2. **Gene Labeling**: Automatically shows top significant genes with clear labels
3. **Professional Appearance**: Publication-ready plots with proper styling

### Interactive Mode
1. **Manual Selection**: Users can switch to "Interactive" mode via radio button
2. **Enhanced Controls**: Color options become available for interactive plots
3. **Original Functionality**: All existing plotly features preserved

### Dynamic Switching
- **Real-time**: Plot type changes immediately when radio button is toggled
- **State Preservation**: Settings maintained when switching between modes
- **Responsive**: Both MAST and MixScale plots update simultaneously

## Benefits

### For Publication
- **Manuscript-Ready**: High-quality static plots suitable for journals
- **Professional Labels**: Clear gene identification with proper formatting
- **Consistent Styling**: Standardized appearance across all volcano plots
- **Vector Graphics**: Scalable plots that maintain quality at any size

### For Analysis
- **Intelligent Labeling**: Automatically highlights most important genes
- **Context Information**: Descriptive titles provide analysis context
- **Flexible Display**: Choice between static and interactive based on needs
- **Robust Performance**: Graceful handling of missing data or packages

## Installation Requirements

### Required Package
```r
# EnhancedVolcano installation
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("EnhancedVolcano")
```

### Automatic Handling
- App attempts automatic installation if package missing
- Shows clear instructions if installation fails
- Continues to function with interactive plots if EnhancedVolcano unavailable

## Testing Checklist

- [ ] Radio button switches between static and interactive modes
- [ ] Static plots show 15-20 labeled genes with white boxes
- [ ] Titles correctly indicate analysis type and experiment
- [ ] Color controls only visible in interactive mode
- [ ] Both MAST and MixScale plots work in both modes
- [ ] Error handling works when data or packages missing
- [ ] Gene labels are readable and well-positioned
- [ ] Professional styling applied to static plots

## Future Enhancements

### Potential Improvements
- **Custom Gene Lists**: Allow users to specify which genes to label
- **Export Options**: Direct export of static plots in various formats
- **Styling Options**: User-configurable color schemes and themes
- **Batch Export**: Generate static plots for all clusters automatically

### Performance Optimizations
- **Caching**: Cache static plots to improve switching speed
- **Lazy Loading**: Only generate plots when tab is viewed
- **Progressive Enhancement**: Load EnhancedVolcano features progressively