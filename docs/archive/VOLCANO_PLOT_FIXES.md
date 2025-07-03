# Volcano Plot Fixes - User Feedback Implementation

## Issues Addressed

### 1. ✅ Terminology Fix
**Problem**: Radio button showed "Static (Publication)" 
**Solution**: Changed to just "Static" as requested

### 2. ✅ Viewport Zero Dimensions Error
**Problem**: Static plots showing "Viewport has zero dimension(s)" error
**Solutions Applied**:
- Added explicit `width = "100%"` and `height = "400px"` to all `plotOutput()` containers
- Added `width = 600, height = 400` parameters to `renderPlot()` functions
- Added comprehensive debugging with console output (`cat()` statements)
- Added explicit `plot.background = element_rect(fill = "white", color = NA)` to all ggplot objects
- Enhanced error handling with proper plot object returns

### 3. ✅ Descriptive Titles Implementation
**Problem**: Titles not descriptive enough about specific mutations/perturbations
**Solution**: Implemented highly descriptive titles for both static and interactive plots:

#### MAST Titles:
- **With Gene**: `"LRRK2 mutation vs isogenic eWT controls (MAST) - Cluster 0"`
- **Without Gene**: `"MAST mutation analysis vs isogenic eWT controls - Cluster 0"`

#### MixScale Titles:
- **With Gene & Experiment**: `"LRRK2 CRISPRi knockdown vs Non-Targeting (C12_FPD-23) - Cluster 0"`
- **With Gene Only**: `"LRRK2 CRISPRi knockdown vs Non-Targeting controls - Cluster 0"`
- **Experiment Only**: `"CRISPRi knockdown vs Non-Targeting (C12_FPD-23) - Cluster 0"`
- **Generic**: `"CRISPRi knockdown vs Non-Targeting controls - Cluster 0"`

### 4. ✅ Legend Optimization
**Problem**: Interactive plots had too many legends taking up space

#### Static Plot Legends (EnhancedVolcano):
- **Position**: Changed from right side to bottom (`legendPosition = 'bottom'`)
- **Size**: Reduced legend label size (`legendLabSize = 8`)
- **Icons**: Smaller legend icons (`legendIconSize = 3.0`)

#### Interactive Plot Legends (Plotly):
- **Orientation**: Changed from vertical to horizontal (`orientation = "h"`)
- **Position**: Centered at bottom (`x = 0.5, y = -0.15, xanchor = "center"`)
- **Styling**: Added background, border, and smaller font
- **Layout**: Added bottom margin (`margin = list(b = 80)`) to accommodate legend

## Technical Implementation

### Enhanced Function Signatures
Updated `generate_volcano_plot()` to accept additional parameters:
```r
generate_volcano_plot(de_data, analysis_type, selected_cluster, color_by, current_gene = NULL, experiment_info = NULL)
```

### Debugging Features
Added comprehensive console logging:
```r
cat("[STATIC PLOT] Starting MAST volcano plot generation\n")
cat("[STATIC PLOT] Current gene:", current_gene, "| Cluster:", values$selected_cluster, "\n")
cat("[STATIC PLOT] Filtered data rows:", nrow(filtered_data), "\n")
cat("[STATIC PLOT] Generated static plot, class:", class(result_plot), "\n")
```

### Error Handling
- Enhanced `tryCatch()` blocks with detailed error messages
- Fallback plots for all error conditions
- Package availability checks for EnhancedVolcano

### Container Improvements
- Dynamic UI containers with proper dimensions
- Conditional rendering based on plot type selection
- Explicit width/height specifications to prevent viewport errors

## Expected Results

### Static Plots
- ✅ No more "Viewport has zero dimension(s)" errors
- ✅ Professional publication-quality plots with EnhancedVolcano
- ✅ 15-20 gene labels with white background boxes
- ✅ Descriptive titles showing exact mutation/perturbation
- ✅ Compact bottom legends

### Interactive Plots  
- ✅ Descriptive titles matching static plot format
- ✅ Horizontal legends at bottom for space efficiency
- ✅ Enhanced hover information with mutation/perturbation details
- ✅ Preserved all existing interactive functionality

### User Experience
- ✅ Default "Static" option (no "Publication" terminology)
- ✅ Seamless switching between static and interactive modes
- ✅ Clear indication of what comparison is being shown
- ✅ More efficient use of screen space with optimized legends

## Files Modified
- `inst/shiny/modules/mod_de_results.R`: Complete volcano plot enhancement with all fixes
- `VOLCANO_PLOT_FIXES.md`: This documentation file

## Testing Checklist
- [ ] Static plots render without viewport errors
- [ ] Titles show specific mutation/perturbation information
- [ ] Legends positioned efficiently (bottom for both plot types)
- [ ] Radio button shows "Static" and "Interactive" options
- [ ] Gene labels visible with white background boxes on static plots
- [ ] Interactive plots maintain all original functionality
- [ ] Debugging console output helps troubleshoot issues

## Debug Instructions
If static plots still don't render:
1. Check console for `[STATIC PLOT]` debug messages
2. Verify EnhancedVolcano package installation
3. Check that plotOutput containers have proper dimensions
4. Ensure ggplot objects have explicit background settings