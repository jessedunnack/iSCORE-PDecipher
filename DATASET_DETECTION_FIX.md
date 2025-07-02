# DATASET DETECTION LOGIC FIX

## Root Cause Discovered ✅

The 35-cluster issue was **NOT** due to wrong UMAP files or cluster data. The problem was **incorrect dataset detection logic** in the landing page module.

## What Was Wrong

### Console Evidence
```
[DE Results] Determined dataset to load: Full_Dataset 
[DE Results] Successfully loaded UMAP (30 PCs) from: Full_Dataset_umap_data_30pc.rds
```

**The app was loading the wrong dataset entirely!**

- ❌ **User selected**: Dataset 2 (iSCORE-PD + CRISPRi) → should load `iSCORE_PD_CRISPRi_umap_data_30pc.rds` (15 clusters)
- ❌ **App actually loaded**: `Full_Dataset_umap_data_30pc.rds` (36 clusters)

### Faulty Logic (Before Fix)
```r
if (has_crispri && has_mutations) {
  dataset_to_load <- "Full_Dataset"  // WRONG!
}
```

This logic treated ANY dataset with both MAST and CRISPRi as "Full_Dataset", but:
- **Dataset 2**: MAST + CRISPRi (no CRISPRa) → should load "iSCORE_PD_CRISPRi"
- **Dataset 3**: MAST + CRISPRi + CRISPRa → should load "Full_Dataset"

## The Fix Applied ✅

### Updated Logic (After Fix)
```r
# Check for specific markers in the data
has_crispri <- any(grepl("MixScale", data$consolidated_data$method))
has_mutations <- any(grepl("MAST", data$consolidated_data$method))
has_crispa <- any(grepl("CRISPRa", data$consolidated_data$method))  // NEW!

# Determine which dataset to load based on actual content
if (has_crispri && has_mutations && has_crispa) {
  dataset_to_load <- "Full_Dataset"                    // Dataset 3
} else if (has_crispri && has_mutations) {
  dataset_to_load <- "iSCORE_PD_CRISPRi"              // Dataset 2 ✓
} else if (has_crispri) {
  dataset_to_load <- "iSCORE_PD_CRISPRi"              // CRISPRi only
} else {
  dataset_to_load <- "iSCORE_PD"                      // Dataset 1
}
```

### Key Changes
1. **Added CRISPRa detection**: `has_crispa <- any(grepl("CRISPRa", data$consolidated_data$method))`
2. **Fixed dataset 2 logic**: Only load "Full_Dataset" if ALL three modalities present
3. **Added debug logging**: Console will show which dataset is selected and why

## Expected Result After Fix

When you select Dataset 2 and relaunch the app:

```
[DATASET DETECTION] MAST: TRUE CRISPRi: TRUE CRISPRa: FALSE
[DATASET DETECTION] Loading iSCORE_PD_CRISPRi (MAST + CRISPRi only)
[DE Results] Successfully loaded UMAP (30 PCs) from: iSCORE_PD_CRISPRi_umap_data_30pc.rds
```

And the overview page should show **15 clusters** instead of 36.

## Files Modified
- `inst/shiny/modules/mod_landing_page_with_umap_v2.R` - Fixed dataset detection logic

## Why This Took So Long to Find
1. **Misleading symptoms**: The issue appeared to be cluster count/data problems
2. **Correct file names**: All UMAP files existed with correct timestamps
3. **Working cluster data**: The Seurat object had correct 15 clusters
4. **Deep caching**: Package installation was working, but loading wrong dataset

The debug output finally revealed the app was loading an entirely different dataset file than expected.

## Next Steps
1. **Reinstall package**: `devtools::install()`
2. **Test app**: `launch_iscore_app()` → select Dataset 2
3. **Verify debug output**: Should show "Loading iSCORE_PD_CRISPRi" 
4. **Confirm UMAP**: Should display 15 clusters in overview