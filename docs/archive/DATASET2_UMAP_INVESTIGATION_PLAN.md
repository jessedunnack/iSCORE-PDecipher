# Dataset 2 UMAP Issue Investigation Plan

## Problem Summary
User reports that after running the clustering fix script and reinstalling the package, the app still shows 35 clusters instead of the expected 15 clusters when loading Dataset 2 (iSCORE-PD + CRISPRi) with 30 PCs.

## Key Findings So Far

### ✅ Confirmed Working
1. **Seurat Object**: Has correct 15 clusters with "cluster_0" format labels ✓
2. **Cluster Naming**: "cluster_0" format matches DE results and enrichment data expectations ✓  
3. **File Updates**: UMAP files were updated July 1, 23:34 with correct timestamps ✓
4. **App Access**: App uses `var = "seurat_clusters"` to access cluster column from SCE object ✓

### ❓ Still Unknown
1. **File Loading**: Which specific file path is the app actually loading from?
2. **Data Structure**: Does our SCE object have the correct "seurat_clusters" column?
3. **Cluster Content**: What cluster values are actually in our updated UMAP files?
4. **Cache Issues**: Is there any caching preventing the new data from loading?

## Investigation Plan

### Phase 1: Verify Current State (Immediate)
1. **Check App Debug Output**: Look for the debug messages we added to see which file is loaded
2. **Verify File Contents**: Confirm our UMAP files contain 15 clusters in correct format
3. **Check Package Installation**: Ensure installed package has updated files

### Phase 2: Root Cause Analysis (If files are correct)
1. **Cache Investigation**: Check for any app caching mechanisms
2. **Path Analysis**: Verify all possible file loading paths in the app
3. **Data Structure Validation**: Ensure SCE object structure matches expectations

### Phase 3: Solution Implementation
1. **Archive Old Files**: Move existing UMAP files to backup location
2. **Generate Fresh Files**: Create new UMAP files from scratch using extract_umap_data.R
3. **Streamlined Script**: Create optimized version without Seurat object saving
4. **Push to GitHub**: Update repository with corrected files

## Detailed Action Items

### A. Create Streamlined Clustering Script
```r
# Skip Seurat object saving (saves 1 hour)
# Focus only on UMAP data generation
# Use extract_umap_data.R approach for consistency
```

### B. Debug File Loading
```r
# Add detailed debug messages to landing page
# Log exact file path loaded
# Log cluster count and names from loaded data
# Check system.file() vs local file paths
```

### C. Fresh File Generation Workflow
1. **Backup existing files**
   ```bash
   mkdir inst/extdata/umap_data/backup_$(date +%Y%m%d)
   mv inst/extdata/umap_data/iSCORE_PD_CRISPRi* inst/extdata/umap_data/backup_$(date +%Y%m%d)/
   ```

2. **Generate new files using established script**
   ```r
   # Use extract_umap_data.R with dataset 2
   # Force recalculation of all PC variants
   # Verify cluster counts match DE expectations
   ```

3. **Install and test**
   ```r
   devtools::install()
   launch_iscore_app()  # Test with dataset 2
   ```

### D. Verification Checklist
- [ ] App loads correct number of clusters (15)
- [ ] UMAP visualization shows proper clustering
- [ ] Cluster labels match DE data expectations
- [ ] All PC variants (30pc, 50pc, 100pc) work correctly
- [ ] No performance regression from file changes

## Questions for User
1. **Debug Output**: When you launch the app, do you see the debug messages we added showing file paths and cluster counts?
2. **Cache Clearing**: Have you tried clearing any browser cache or restarting R session completely?
3. **File Verification**: Would you prefer we focus on archiving old files and generating completely fresh ones?

## Next Steps
Based on user feedback, proceed with:
1. Immediate debug output analysis
2. Fresh file generation if needed
3. Streamlined script creation (no Seurat saving)
4. Comprehensive testing and GitHub push