# Checkpoint v0.1.4-stable

**Created:** $(date)  
**Status:** ✅ **STABLE & TESTED**  
**Tag:** `v0.1.4-stable`

## 🎯 What This Checkpoint Includes

### **Core Functionality Verified:**
- ✅ **App Launch:** `launch_iscore_app()` works without errors
- ✅ **UMAP Visualization:** Overview page with cluster selection
- ✅ **Markers Table:** Top 25 genes with MAST methodology info
- ✅ **Multi-PC Support:** 30/50/100 PC UMAP selector (when data available)
- ✅ **Clean Repository:** Professional R package structure

### **Major Features Complete:**
1. **Repository Organization**
   - Essential R package files only in root
   - Documentation organized in `docs/`
   - Development files organized in `dev/`

2. **Enhanced UMAP Interface**
   - Multi-PC selector (30/50/100 principal components)
   - Dynamic data loading based on selection
   - Perfect panel alignment (720px height)

3. **Improved Markers Display**
   - Clear column names: "P-adj", "% this clust", "% in others"
   - MAST methodology info integrated in controls area
   - Fixed table height utilization (385px)

4. **Robust App Launch**
   - Enhanced path resolution for all scenarios
   - Proper sourcing of helper functions
   - Better error messages

### **Repository Structure:**
```
iSCORE-PDecipher/
├── DESCRIPTION, NAMESPACE, README.md  # R package essentials
├── R/                                 # Core functions
├── inst/                             # Data & Shiny app
├── man/                              # Documentation
├── docs/                             # User guides & sharing info
└── dev/                              # Test scripts & utilities
```

### **Recent Improvements (since v0.1.3):**
- **Repository cleanup:** Removed 9 outdated files, organized 14 files
- **Multi-PC UMAP:** Added PC selector and data generation script
- **Layout fixes:** Perfect panel alignment and table space utilization
- **Methodology clarity:** Clear MAST parameters display
- **Launch reliability:** Fixed path resolution issues

## 🚀 How to Use This Checkpoint

### **To Revert to This Version:**
```bash
git checkout v0.1.4-stable
```

### **To Compare with Current State:**
```bash
git diff v0.1.4-stable HEAD
```

### **To Create a Branch from This Checkpoint:**
```bash
git checkout -b feature-branch v0.1.4-stable
```

## 🧪 Verification Commands

Test that everything works:
```r
# Install and test
remotes::install_github("jessedunnack/iSCORE-PDecipher", ref = "v0.1.4-stable")
library(iSCORE.PDecipher)

# Should launch without errors
launch_iscore_app()
```

## ⚠️ Important Notes

1. **Multi-PC Data:** PC selector requires running `generate_multipc_umaps.R` first
2. **Path Requirements:** Launch app from package directory or with proper installation
3. **Data Files:** Requires dataset with `all_enrichment_padj005_complete_with_direction.rds`

## 🎉 Safe Development

This checkpoint preserves a **fully functional state** before any experimental features. Use this as a safe fallback for:
- Testing new features
- Rolling back breaking changes  
- Sharing stable version with collaborators
- Creating feature branches

**Last verified:** App launches and displays UMAP visualization correctly