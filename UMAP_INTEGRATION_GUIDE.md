# UMAP Integration Guide for iSCORE-PDecipher

## Overview

This guide explains how to extract UMAP data from large Seurat objects and integrate interactive UMAP visualizations into the iSCORE-PDecipher Shiny app.

## Step 1: Extract UMAP Data

### Prerequisites
- Access to the Seurat RDS files (20-30 GB each)
- R with Seurat and SingleCellExperiment packages installed
- At least 32 GB RAM recommended

### Running the Extraction Script

1. **Adjust paths in the extraction script** (`inst/scripts/extract_umap_data.R`):
   ```r
   SEURAT_DATA_DIR <- "/path/to/your/seurat/files"
   OUTPUT_DIR <- "/path/to/iSCORE-PDecipher/inst/extdata/umap_data"
   ```

2. **Run the extraction script**:
   ```bash
   cd /path/to/iSCORE-PDecipher
   Rscript inst/scripts/extract_umap_data.R
   ```

   Or run interactively in R:
   ```r
   source("inst/scripts/extract_umap_data.R")
   main()
   ```

3. **Expected output**:
   - Individual files: `GWAS_umap_data.rds`, `FPD_umap_data.rds`, `CRISPRi_umap_data.rds`
   - Combined file: `all_umap_data_combined.rds`
   - Summary CSV: `umap_data_summary.csv`
   - Total size: ~100-500 MB (vs 20-30 GB per original file)

### What the Script Extracts
- UMAP coordinates (umap.cca reduction)
- Cell metadata:
  - seurat_clusters (coarse)
  - seurat_clusters_fine (if available)
  - scMAGeCK_gene_assignment (perturbations)
  - experiments
  - orig.ident
- Minimal SingleCellExperiment objects for dittoSeq compatibility

## Step 2: Integrate UMAP Viewer into Shiny App

### Option A: Use the Enhanced Landing Page (Recommended)

1. **Update app.R** to use the new landing page with UMAP:
   ```r
   # Replace this line:
   source("modules/mod_landing_page.R")
   
   # With:
   source("modules/mod_landing_page_with_umap.R")
   
   # And update the UI call:
   # Replace:
   mod_landing_page_ui("landing")
   
   # With:
   landingPageWithUmapUI("landing")
   
   # And update the server call:
   # Replace:
   mod_landing_page_server("landing", app_data)
   
   # With:
   landingPageWithUmapServer("landing", app_data)
   ```

### Option B: Add UMAP as a Separate Tab

1. **Add to sidebar menu** in app.R:
   ```r
   menuItem("UMAP Explorer", 
            tabName = "umap", 
            icon = icon("chart-scatter"))
   ```

2. **Add to tab items**:
   ```r
   tabItem(tabName = "umap",
           mod_umap_viewer_ui("umap_viewer"))
   ```

3. **Add to server**:
   ```r
   mod_umap_viewer_server("umap_viewer", app_data)
   ```

## Step 3: Features of the UMAP Viewer

### Interactive Controls
- **Dataset Selection**: Switch between GWAS, FPD, and CRISPRi datasets
- **Color By**: 
  - Coarse clusters
  - Fine clusters
  - Perturbations (gene knockdowns)
- **Visual Options**:
  - Toggle cluster labels
  - Adjust point size (0.1-2)
  - Adjust label size (2-8)
  - Show/hide legend

### Performance Optimizations
- Lightweight data (100-500 MB vs 20-30 GB)
- Fast loading and rendering
- No need to load full Seurat objects
- Cached data for quick dataset switching

## Step 4: Customization Options

### Modify Color Schemes
Edit in `mod_umap_viewer.R`:
```r
# dittoSeq uses automatic color assignment
# To use custom colors:
p <- dittoDimPlot(...) + 
  scale_color_manual(values = your_color_vector)
```

### Add Additional Metadata
1. Update extraction script to include more metadata columns
2. Add to `metadata_cols` list in `extract_minimal_umap_data()`
3. Update UI choices in `mod_umap_viewer.R`

### Change Reduction Method
To use a different reduction (e.g., standard UMAP instead of umap.cca):
```r
# In extraction script:
umap_coords <- Embeddings(seurat_obj, reduction = "umap")

# In viewer module:
reduction.use = "UMAP"  # Must match reducedDims name
```

## Troubleshooting

### Common Issues

1. **"UMAP data not found" error**:
   - Ensure extraction script has been run
   - Check file paths in `mod_umap_viewer_server`
   - Verify `inst/extdata/umap_data/` directory exists

2. **"dittoSeq package is required" error**:
   ```r
   BiocManager::install("dittoSeq")
   ```

3. **Memory issues during extraction**:
   - Process one dataset at a time
   - Clear memory between datasets: `rm(seurat_obj); gc()`
   - Consider using a high-memory compute node

4. **Missing metadata columns**:
   - Check available columns in Seurat object
   - Update `metadata_cols` list in extraction script
   - Some datasets may not have all metadata

## Benefits of This Approach

1. **Performance**: 50-100x smaller data files
2. **Speed**: Instant dataset switching
3. **Flexibility**: Easy to add/remove datasets
4. **Compatibility**: Works with dittoSeq visualization
5. **Integration**: Seamless with existing Shiny app

## Next Steps

1. Consider adding:
   - Dataset comparison views (side-by-side UMAPs)
   - Cell type annotations
   - Gene expression overlay
   - Export functionality for publication figures

2. For production deployment:
   - Pre-generate UMAP data during package build
   - Include sample data for testing
   - Add progress indicators for initial load

---

Last updated: January 6, 2025