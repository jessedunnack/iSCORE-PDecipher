# UMAP Integration Solution Summary

## Overview

I've created a complete solution for integrating interactive UMAP visualizations into the iSCORE-PDecipher Shiny app without loading the massive (20-30 GB) Seurat objects.

## Solution Components

### 1. Data Extraction Script (`inst/scripts/extract_umap_data.R`)

**Purpose**: Extracts minimal data from large Seurat objects

**Key Features**:
- Creates lightweight SingleCellExperiment objects (~100-500 MB vs 20-30 GB)
- Preserves UMAP coordinates and essential metadata
- Maintains compatibility with dittoSeq visualization
- Processes all three datasets (GWAS, FPD, CRISPRi)
- Generates summary statistics

**What it extracts**:
- UMAP coordinates (umap.cca reduction)
- Cell cluster assignments (coarse and fine)
- Perturbation information (scMAGeCK_gene_assignment)
- Experiment identifiers
- Dataset metadata

### 2. UMAP Viewer Module (`inst/shiny/modules/mod_umap_viewer.R`)

**Interactive Features**:
- Dataset switching (GWAS/FPD/CRISPRi)
- Multiple coloring options:
  - Coarse clusters (with labels)
  - Fine clusters
  - Perturbations (gene knockdowns)
- Visual customization:
  - Adjustable point size
  - Toggle cluster labels
  - Configurable label size
  - Show/hide legend
- Real-time plot updates

**Smart Features**:
- Automatic detection of available metadata
- Dynamic UI updates based on dataset
- Handles high-dimensional perturbation data gracefully
- Error handling for missing data

### 3. Enhanced Landing Page (`mod_landing_page_with_umap.R`)

**Integration Approach**:
- UMAP viewer placed prominently at top of homepage
- Collapsible design to save space
- Maintains all existing landing page functionality
- Seamless visual integration

**Benefits**:
- Immediate visual context for users
- No navigation required to see data overview
- Optional - can be toggled on/off

### 4. Test Scripts

**`test_umap_visualization.R`**:
- Validates the approach with synthetic data
- Tests actual extracted data
- Generates example plots
- Verifies dittoSeq compatibility

## Implementation Steps

### Step 1: Extract UMAP Data
```bash
# Run from command line
Rscript inst/scripts/extract_umap_data.R

# Or in R
source("inst/scripts/extract_umap_data.R")
main()
```

### Step 2: Update app.R
Apply the changes from `app_umap_integration.patch`:
- Replace landing page module source
- Update UI and server calls

### Step 3: Install Dependencies
```r
BiocManager::install("dittoSeq")
install.packages("SingleCellExperiment")
```

## Technical Advantages

### Memory Efficiency
- **Original Seurat objects**: 20-30 GB each
- **Extracted SCE objects**: 100-500 MB each
- **Reduction**: ~50-100x smaller

### Performance
- Instant dataset switching
- No lag when changing visualization parameters
- Efficient rendering with dittoSeq

### Flexibility
- Easy to add new metadata columns
- Simple to update with new datasets
- Modular design for easy maintenance

## Visualization Examples

### 1. Cluster Analysis
```r
dittoDimPlot(sce, "seurat_clusters", 
             reduction.use = "UMAP", 
             do.label = TRUE)
```

### 2. Perturbation Effects
```r
dittoDimPlot(sce, "scMAGeCK_gene_assignment",
             reduction.use = "UMAP",
             size = 0.5)
```

### 3. Fine Clustering
```r
dittoDimPlot(sce, "seurat_clusters_fine",
             reduction.use = "UMAP",
             legend.show = FALSE)  # Too many clusters
```

## Customization Options

### Add New Metadata
1. Update extraction script `metadata_cols`
2. Re-run extraction
3. Update UI choices in viewer module

### Change Visual Style
- Modify dittoSeq parameters in `mod_umap_viewer.R`
- Add color palettes
- Adjust default sizes

### Add Export Features
- PNG/PDF export buttons
- High-resolution output for publications
- Batch export for all datasets

## Future Enhancements

1. **Side-by-side comparisons**: Show multiple datasets simultaneously
2. **Gene expression overlay**: Color by specific gene expression
3. **Cell type annotations**: Add predicted cell type labels
4. **Trajectory analysis**: Show developmental trajectories
5. **Selection tools**: Interactive cell selection for downstream analysis

## Files Created

```
inst/
├── scripts/
│   ├── extract_umap_data.R          # Data extraction script
│   └── test_umap_visualization.R    # Test and validation
├── shiny/
│   └── modules/
│       ├── mod_umap_viewer.R              # UMAP viewer module
│       └── mod_landing_page_with_umap.R   # Enhanced landing page
└── extdata/
    └── umap_data/                    # Output directory for extracted data
```

## Package Updates

- **Version**: Bumped to 0.1.3
- **New Dependencies**:
  - SingleCellExperiment
  - dittoSeq

## Summary

This solution provides a performant, user-friendly way to visualize UMAP embeddings in the iSCORE-PDecipher app without the computational burden of loading full Seurat objects. The modular design allows for easy maintenance and future enhancements while maintaining the existing app functionality.

---

Last updated: January 6, 2025