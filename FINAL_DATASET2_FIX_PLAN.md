# FINAL DATASET 2 UMAP FIX PLAN

## Executive Summary
After systematic investigation, we've created a streamlined solution to fix the persistent 35-cluster issue in Dataset 2. The problem appears to be related to file structure or caching, not cluster naming (which is correct).

## Root Cause Analysis Completed

### ✅ Confirmed Correct
1. **Seurat Object**: Has 15 clusters with correct "cluster_0" format ✓
2. **Cluster Naming**: Matches DE results and enrichment data expectations ✓
3. **App Access Pattern**: Uses `var = "seurat_clusters"` from SCE objects ✓
4. **Data Expectations**: All downstream data expects "cluster_0" format ✓

### 🔍 Potential Issues Identified
1. **File Structure**: Our SCE objects may not have correct internal structure
2. **Caching**: App may be loading cached or wrong file paths
3. **Package vs Local**: Mismatch between installed package and local files

## Complete Solution Workflow

### Step 1: Archive Existing Files
```r
source("archive_dataset2_files.R")
```
This will backup and remove all existing iSCORE_PD_CRISPRi UMAP files.

### Step 2: Generate Fresh UMAP Data (Streamlined - No Seurat Saving)
```r
source("generate_dataset2_umap_only.R")
```
This streamlined script:
- ✅ Skips reclustering (saves 1 hour)
- ✅ Uses existing 15 clusters from Seurat object
- ✅ Creates proper SingleCellExperiment objects
- ✅ Generates all PC variants (30pc, 50pc, 100pc)
- ✅ Includes comprehensive validation

### Step 3: Reinstall Package
```r
remove.packages("iSCORE.PDecipher")
devtools::install()
library(iSCORE.PDecipher)
```

### Step 4: Test App
```r
launch_iscore_app()
```
Check that Dataset 2 overview page shows 15 clusters.

### Step 5: Push to GitHub (If Successful)
```bash
git add inst/extdata/umap_data/iSCORE_PD_CRISPRi_*
git commit -m "Fix Dataset 2 UMAP clustering - 15 clusters with correct format

- Archive old 35-cluster UMAP files
- Generate fresh SingleCellExperiment objects with proper structure
- Maintain cluster_0 format consistent with DE/enrichment data
- Update all PC variants (30pc, 50pc, 100pc)
- Verify seurat_clusters column in SCE colData"

git push origin main
```

## Files Created

### Core Scripts
1. **`generate_dataset2_umap_only.R`** - Streamlined UMAP generation (no Seurat saving)
2. **`archive_dataset2_files.R`** - Archive existing files before regeneration
3. **`DATASET2_UMAP_INVESTIGATION_PLAN.md`** - Detailed investigation documentation
4. **`FINAL_DATASET2_FIX_PLAN.md`** - This action plan

### Debug Scripts (For Investigation)
1. **`check_seurat_clusters.R`** - Verify Seurat object cluster status
2. **`inspect_umap_structure.R`** - Examine SCE object structure
3. **`check_umap_clusters_simple.R`** - Quick cluster count verification

## Key Improvements in Streamlined Script

### Performance Optimizations
- ⚡ **No Seurat object saving** (saves ~1 hour)
- ⚡ **Read-only Seurat access** (faster startup)
- ⚡ **Efficient marker calculation** with downsampling

### Data Structure Compliance
- 🔧 **Proper SCE object creation** matching extract_umap_data.R
- 🔧 **Correct metadata columns** (seurat_clusters, mutation_tidy, etc.)
- 🔧 **All PC variants generated** (30pc, 50pc, 100pc)
- 🔧 **Comprehensive validation** of cluster counts and labels

### Debug & Verification
- 🔍 **Detailed logging** of cluster counts and labels
- 🔍 **SCE structure verification** including colData columns
- 🔍 **File path confirmation** for all generated files
- 🔍 **Expected vs actual cluster comparison**

## Expected Outcome

After running this workflow:
1. **App Launch**: Dataset 2 overview should show 15 clusters instead of 35
2. **Cluster Labels**: Should display cluster_0, cluster_1, ..., cluster_14
3. **PC Selection**: All variants (30pc, 50pc, 100pc) should work correctly
4. **Performance**: No regression - should load as fast as before

## Troubleshooting

### If Still Shows 35 Clusters
1. **Check debug output** from the landing page module
2. **Verify file locations** using system.file() commands
3. **Clear R session** completely and restart
4. **Check browser cache** if using RStudio viewer

### If Script Fails
1. **Check SingleCellExperiment package** installation
2. **Verify file paths** are correct for your platform
3. **Check disk space** for large file operations
4. **Verify Seurat object** loads correctly

## Questions Before Proceeding

1. **Should we archive existing files first** or try the streamlined script directly?
2. **Do you see any debug output** when launching the current app?
3. **Would you prefer to test locally first** before pushing to GitHub?

## Next Action
Ready to execute Step 1 (archiving) when you give the go-ahead.