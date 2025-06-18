# Mac Setup Checklist for iSCORE-PDecipher

## Quick Summary

✅ **Mac Compatibility**: Complete - no hard-coded Windows paths  
✅ **Automatic Setup**: First-launch prompts for data location  
✅ **File Transfer**: Helper scripts identify required files  
✅ **Cross-Platform**: Uses standard R config directories  

## For You (Data Preparation)

### Step 1: Generate Transfer Report
```r
# Load the updated package
library(iSCORE.PDecipher)

# Check what files exist and their sizes
generate_transfer_report()
```

This will show you:
- Which files exist in each dataset
- Total size for transfer planning
- Missing essential files (if any)

### Step 2: Prepare Files for Transfer
```r
# Option A: Copy to flash drive (recommended files only)
prepare_mac_transfer_copy(
  source_parent_dir = "E:/ASAP/scRNASeq/PerturbSeq/final",
  dest_parent_dir = "F:/iSCORE_for_mac",  # Flash drive path
  copy_mode = "recommended"  # ~470 MB per dataset
)

# Option B: Minimal copy (essential files only) 
prepare_mac_transfer_copy(
  source_parent_dir = "E:/ASAP/scRNASeq/PerturbSeq/final",
  dest_parent_dir = "F:/iSCORE_for_mac", 
  copy_mode = "minimal"  # ~280 MB per dataset
)
```

### Step 3: Create Flash Drive
1. Copy the prepared directories to flash drive
2. Include `MAC_COMPATIBILITY_GUIDE.md` 
3. Include this checklist as `MAC_SETUP_CHECKLIST.md`

## For Your Boss (Mac Installation)

### Step 1: Install R Dependencies
```r
# Install package manager
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "clusterProfiler", "ReactomePA", "DOSE", "org.Hs.eg.db",
  "pathview", "ComplexHeatmap", "SingleCellExperiment", 
  "dittoSeq", "enrichplot"
))

# Install additional dependencies
install.packages(c("rappdirs", "jsonlite", "heatmaply"))
```

### Step 2: Install iSCORE-PDecipher
```r
remotes::install_github("jessedunnack/iSCORE-PDecipher")
```

### Step 3: Copy Data from Flash Drive
Copy the flash drive contents to Mac, e.g.:
```bash
cp -r /Volumes/FlashDrive/iSCORE_datasets ~/Documents/
```

### Step 4: First Launch
```r
library(iSCORE.PDecipher)
launch_iscore_app()
```

When prompted, enter the full Mac path:
```
/Users/username/Documents/iSCORE_datasets
```

The app will:
- ✅ Validate the directory structure
- ✅ Show available datasets  
- ✅ Save the configuration for future launches

## Required Files Summary

| File Type | Size | Required | Purpose |
|-----------|------|----------|---------|
| **Essential** | | | |
| `all_enrichment_padj005_complete_with_direction.rds` | 266 MB | ✅ | Main enrichment data |
| **UMAP Data** | | | |
| `all_umap_data_combined.rds` | 14 MB | ✅ | Cell visualization |
| **Recommended** | | | |
| `full_DE_results.rds` | 190 MB | 🔄 | Volcano plots |
| Dataset-specific UMAP files | ~20 MB | 🔄 | Better performance |

**Total per dataset:**
- **Minimal**: ~280 MB (basic functionality)
- **Recommended**: ~470 MB (full functionality)

## Auto-Generation Capabilities

If files are missing, the app can generate them:

1. **`full_DE_results.rds`** (5-10 min) - if source MAST/MixScale data available
2. **Enrichment analysis** (hours) - if source data available  
3. **UMAP data** - requires original Seurat objects

## Troubleshooting

### If Setup Fails
```r
# Reset configuration
iSCORE.PDecipher:::set_parent_data_dir("/new/path")

# Or force re-setup
unlink(iSCORE.PDecipher:::get_config_path())
launch_iscore_app()
```

### Check Configuration
```r
# See current config
iSCORE.PDecipher:::get_parent_data_dir()

# Check if first launch
iSCORE.PDecipher:::is_first_launch()
```

### Validate Data
```r
# Check specific dataset
iSCORE.PDecipher:::check_dataset_files("/path/to/dataset")

# Generate new report
generate_transfer_report("/path/to/parent")
```

## Key Changes Made

1. **Added `rappdirs` + `jsonlite` dependencies** for cross-platform config
2. **Created `config_manager.R`** with persistent settings
3. **Updated `dataset_validator.R`** to use stored paths
4. **Modified `launch_app.R`** for first-launch setup
5. **Added transfer helper functions** for file preparation

## Testing Checklist

Before giving to your boss:

- [ ] Test `generate_transfer_report()` shows correct files
- [ ] Test `prepare_mac_transfer_copy()` creates proper structure  
- [ ] Verify flash drive has complete directory structure
- [ ] Test config reset: `unlink(get_config_path())` then `launch_iscore_app()`
- [ ] Confirm app works with new path on your system

## File Locations on Mac

**Config file**: `~/Library/Application Support/iSCORE.PDecipher/config.json`
**Data location**: User-specified (e.g., `~/Documents/iSCORE_datasets/`)

The config file persists the data location so your boss only needs to set it up once.

---

**Result**: Your boss can install the R package, copy the flash drive data, and run the app without any path modifications or code changes. The app will automatically detect available datasets and work seamlessly on Mac.