# Google Drive Sharing Checklist for Labmates

## 📁 Files to Upload to Google Drive

### Essential Files (Required)
Upload these **2 files** to a shared Google Drive folder:

1. **`all_enrichment_padj005_complete_with_direction.rds`** (266 MB)
   - **Location**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/data_files/`
   - **Contains**: 767,337 enrichment terms (MAST + CRISPRi + CRISPRa)
   - **Purpose**: Main data for all visualizations and analysis

2. **`all_umap_data_combined.rds`** (14 MB)
   - **Location**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/inst/extdata/umap_data/`
   - **Contains**: UMAP coordinates for all three datasets
   - **Purpose**: Interactive cell cluster visualization

**Total size:** ~280 MB

## 📋 Sharing Instructions

### Step 1: Create Shared Folder
1. Create a new Google Drive folder: `iSCORE-PDecipher-Data`
2. Share with lab members (view/download permissions)
3. Add folder description: "Pre-computed data for iSCORE-PDecipher Shiny app"

### Step 2: Upload Files
1. Upload both RDS files to the shared folder
2. Verify file sizes match expected values
3. Test download to ensure files aren't corrupted

### Step 3: Share Instructions
1. Send labmates the `LABMATE_QUICKSTART.md` guide
2. Include Google Drive folder link
3. Recommend they create a local folder like `~/PD_Analysis/` for the data

## 🎯 What Labmates Get

### Immediate Access To:
- **767,337 enrichment terms** across all PD experiments
- **Interactive UMAP visualizations** of 201,679 cells
- **13 genetic mutations** (LRRK2, GBA, PINK1, PRKN, etc.)
- **10 CRISPR perturbations** (knockdown + activation)
- **6 enrichment databases** (GO, KEGG, Reactome, etc.)

### Analysis Capabilities:
- Cross-method comparisons (mutations vs perturbations)
- Cell cluster-specific responses
- Pathway enrichment visualization
- Interactive heatmaps and plots
- Publication-quality figure export

## ⚡ Alternative: Smaller Datasets

If 280 MB is too large, you can create focused subsets:

### Option A: Single Method (MAST only)
- Filter `all_enrichment_padj005_complete_with_direction.rds` to MAST only
- ~65 MB instead of 266 MB
- Still includes all mutations and clusters

### Option B: Top Pathways Only  
- Filter to top 100K most significant terms
- ~100 MB instead of 266 MB
- Covers most important biology

### Code to Create Subsets:
```r
# Load full data
data <- readRDS("all_enrichment_padj005_complete_with_direction.rds")

# Option A: MAST only
mast_data <- data[data$method == "MAST", ]
saveRDS(mast_data, "mast_only_enrichment.rds")

# Option B: Top 100K terms
top_data <- data[order(data$p.adjust), ][1:100000, ]
saveRDS(top_data, "top_100k_enrichment.rds")
```

## 📞 Support Plan

### For Labmates:
1. **Installation issues**: Point to `LABMATE_QUICKSTART.md`
2. **Data download**: Verify Google Drive permissions
3. **App crashes**: Check R/package versions
4. **Analysis questions**: Direct to Jesse

### For Jesse:
1. Monitor Google Drive usage/downloads
2. Be available for first-time setup help
3. Consider lab meeting demo if multiple people interested
4. Update quickstart guide based on feedback

---

**Ready to share!** Your labmates will have a powerful PD analysis tool with just a 5-minute setup.