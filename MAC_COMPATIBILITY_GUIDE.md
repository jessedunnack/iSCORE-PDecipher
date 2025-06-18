# Mac Compatibility Guide for iSCORE-PDecipher

## Overview

The iSCORE-PDecipher R package has been updated to be fully Mac compatible with cross-platform path management. This guide explains how to set up the package on a Mac and what files are needed for full functionality.

## New Features for Cross-Platform Compatibility

### 1. Automatic Configuration Management
- **First Launch Setup**: On first launch, the app prompts for the parent directory containing your dataset folders
- **Persistent Configuration**: Settings are saved in a platform-appropriate location:
  - **Mac**: `~/Library/Application Support/iSCORE.PDecipher/`
  - **Windows**: `%APPDATA%/iSCORE.PDecipher/`
- **No Hard-coded Paths**: All dataset paths are now configurable and stored in user config

### 2. Required Directory Structure

Your Mac should have a parent directory containing these subdirectories:
```
your_data_folder/
├── iSCORE-PD/
├── iSCORE-PD_plus_CRISPRi/
└── iSCORE-PD_plus_CRISPRi_plus_CRISPRa/
```

The app will automatically detect which datasets are available and only show options for existing directories.

## Installation Instructions for Mac

### 1. Install the R Package
```r
# Install from GitHub
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

remotes::install_github("jessedunnack/iSCORE-PDecipher")
```

### 2. Install Required Dependencies
Some Bioconductor packages may need manual installation:
```r
# Install BiocManager if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Install Bioconductor dependencies
BiocManager::install(c(
  "clusterProfiler",
  "ReactomePA", 
  "DOSE",
  "org.Hs.eg.db",
  "pathview",
  "ComplexHeatmap",
  "SingleCellExperiment",
  "dittoSeq",
  "enrichplot"
))

# Install CRAN dependencies
install.packages(c(
  "rappdirs",
  "jsonlite",
  "heatmaply"
))
```

### 3. First Launch
```r
library(iSCORE.PDecipher)
launch_iscore_app()
```

On first launch, you'll see:
```
=== iSCORE-PDecipher First-Time Setup ===

Please specify the parent directory that contains your dataset folders.
This directory should contain the following subdirectories:
  - iSCORE-PD/
  - iSCORE-PD_plus_CRISPRi/
  - iSCORE-PD_plus_CRISPRi_plus_CRISPRa/

Enter the full path to your parent data directory: 
```

## Required Files for Flash Drive Transfer

### Essential Files (Minimum for App Launch)

#### **Option 1: Minimal Setup (~280 MB)**
For basic functionality, copy these files to each dataset folder:

```
dataset_folder/
├── all_enrichment_padj005_complete_with_direction.rds  # 266 MB
└── inst/extdata/umap_data/
    └── all_umap_data_combined.rds                      # 14 MB
```

#### **Option 2: Complete Setup (~470 MB per dataset)**
For full functionality including volcano plots:

```
dataset_folder/
├── all_enrichment_padj005_complete_with_direction.rds  # 266 MB - Main enrichment data
├── full_DE_results.rds                                 # 190 MB - For volcano plots
└── inst/extdata/umap_data/                            # UMAP visualization data
    ├── {dataset_name}_umap_data_30pc.rds              # ~5 MB each
    ├── {dataset_name}_umap_data_50pc.rds              
    ├── {dataset_name}_umap_data_100pc.rds             
    ├── {dataset_name}_cluster_markers.rds             # ~1 MB
    └── all_umap_data_combined.rds                      # 14 MB (fallback)
```

### Dataset-Specific UMAP Files

Replace `{dataset_name}` with:
- **iSCORE-PD**: `iSCORE_PD_umap_data_*`
- **iSCORE-PD + CRISPRi**: `iSCORE_PD_CRISPRi_umap_data_*`  
- **Full Dataset**: `Full_Dataset_umap_data_*`

### Optional Source Data (For Auto-Generation)

If you want the ability to regenerate files on Mac:

```
dataset_folder/
├── iSCORE-PD_MAST_analysis/           # MAST source data
│   └── mutation_*.rds
├── PerturbSeq_MixScale_analysis*/     # MixScale source data
│   ├── CRISPRi/
│   │   └── *_DEGs.rds
│   └── CRISPRa/
│       └── *_DEGs.rds
└── enrichment_results/                # Individual enrichment files (14K+ files)
    └── (individual result files)
```

## File Size Summary

| Component | Size | Required | Purpose |
|-----------|------|----------|---------|
| Main enrichment data | 266 MB | ✅ Essential | All visualizations |
| DE results | 190 MB | 🔄 Auto-generated | Volcano plots |
| UMAP data (combined) | 14 MB | ✅ Essential | Cell visualization |
| UMAP data (specific) | ~20 MB | 🔄 Optional | Better performance |
| Source MAST data | Variable | 🔄 Optional | Auto-generation |
| Source MixScale data | Variable | 🔄 Optional | Auto-generation |
| **Total minimum** | **~280 MB** | per dataset | Basic functionality |
| **Total complete** | **~470 MB** | per dataset | Full functionality |

## Flash Drive Preparation Checklist

### For Your Boss's Mac Setup

1. **Create the directory structure** on flash drive:
   ```
   iSCORE_datasets/
   ├── iSCORE-PD/
   │   ├── all_enrichment_padj005_complete_with_direction.rds
   │   ├── full_DE_results.rds
   │   └── inst/extdata/umap_data/
   │       ├── iSCORE_PD_umap_data_30pc.rds
   │       ├── iSCORE_PD_umap_data_50pc.rds
   │       ├── iSCORE_PD_umap_data_100pc.rds
   │       └── iSCORE_PD_cluster_markers.rds
   ├── iSCORE-PD_plus_CRISPRi/
   │   └── (same structure with CRISPRi files)
   └── iSCORE-PD_plus_CRISPRi_plus_CRISPRa/
       └── (same structure with full dataset files)
   ```

2. **Include this guide** as `MAC_SETUP_INSTRUCTIONS.md`

3. **Test on your system** before transferring:
   ```r
   # Test the setup
   library(iSCORE.PDecipher)
   launch_iscore_app(data_dir = "/path/to/test/dataset")
   ```

## Mac-Specific Setup Steps for Your Boss

1. **Copy data from flash drive** to desired location on Mac (e.g., `~/Documents/iSCORE_datasets/`)

2. **Install R and RStudio** (if not already installed)

3. **Install the package and dependencies** (see Installation Instructions above)

4. **Launch the app**:
   ```r
   library(iSCORE.PDecipher)
   launch_iscore_app()
   ```

5. **When prompted for parent directory**, enter the path where the data was copied (e.g., `/Users/username/Documents/iSCORE_datasets`)

6. **The app will validate** the directory structure and show available datasets

## Troubleshooting

### Path Issues
- **Use full paths**: `/Users/username/Documents/iSCORE_datasets` not `~/Documents/iSCORE_datasets`
- **Check permissions**: Ensure read access to all data files
- **Verify structure**: Use `ls -la` in Terminal to check directory contents

### Missing Dependencies
```r
# Check for missing packages
required_packages <- c("rappdirs", "jsonlite", "heatmaply", "clusterProfiler")
missing <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(missing)) install.packages(missing)
```

### Configuration Reset
If you need to reconfigure the parent directory:
```r
# Reset configuration
iSCORE.PDecipher:::set_parent_data_dir("/new/path/to/data")

# Or delete config to force re-setup
unlink(iSCORE.PDecipher:::get_config_path())
```

## New Configuration Functions

The package now includes these functions for path management:

- `is_first_launch()`: Check if this is the first time running
- `get_parent_data_dir()`: Get configured parent directory
- `set_parent_data_dir(path)`: Set parent directory manually
- `setup_parent_dir()`: Interactive setup with validation

These functions handle cross-platform compatibility automatically.

## Version Compatibility

- **Package Version**: 0.1.3+
- **R Version**: >= 4.0.0
- **Tested on**: macOS Big Sur, macOS Monterey, Windows 10/11
- **Dependencies**: Cross-platform compatible (rappdirs, jsonlite)

This update ensures the iSCORE-PDecipher package works seamlessly on both Mac and Windows systems without any hard-coded path dependencies.