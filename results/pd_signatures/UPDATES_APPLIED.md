# Updates Applied to Original Files - July 20, 2025

## ✅ All Fixes Applied In-Place

### 1. Gene Harmonization
**Files Updated in `/by_gene/`:**
- `all_genes_summary.csv` - Now shows 11 harmonized genes instead of 16 variants
- Individual CSV files renamed to harmonized names:
  - `PARK2_*_pathways.csv` (includes data from both PRKN and PARK2)
  - `VPS13C_*_pathways.csv` (includes A444P and W395C variants)
  - `SNCA_*_pathways.csv` (includes A30P and A53T variants)
- All plots in `/by_gene/plots/` updated with harmonized names

**Key Results:**
- PARK2 now shows 106 convergent pathways (was 0 for PRKN alone)
- VPS13C shows 160 convergent pathways (variants combined)
- SNCA shows 223 convergent pathways (variants combined)

### 2. Clustering Added to Heatmaps
**Files Updated:**
- `/by_gene/all_genes_pathway_heatmap.pdf` - Now hierarchically clustered by gene similarity
- `/by_cluster/cluster_pathway_heatmap.pdf` - Pathways clustered, clusters in natural order

### 3. Natural Cluster Sorting
**Files Updated in `/by_cluster/`:**
- `cluster_method_breakdown.csv` - Clusters sorted 0,1,2...9,10,11,12,13,14
- `cluster_method_distribution.pdf` - X-axis properly ordered
- `cluster_pathway_heatmap.pdf` - Columns in natural order
- All individual `cluster_*_top_pathways.csv` files maintained

## 📊 Summary Statistics After Updates

### Harmonized Genes (11 total):
1. **SYNJ1**: 301 convergent pathways (highest)
2. **ATP13A2**: 224 convergent pathways
3. **SNCA**: 223 convergent (combined A30P, A53T)
4. **LRRK2**: 166 convergent pathways
5. **VPS13C**: 160 convergent (combined A444P, W395C)
6. **FBXO7**: 130 convergent pathways
7. **PARK7**: 121 convergent pathways
8. **PARK2**: 106 convergent (combined PRKN, PARK2)
9. **PINK1**: 63 convergent pathways
10. **DNAJC6**: 19 convergent pathways
11. **GBA**: 0 convergent (MAST only)

### Clusters (15 total, naturally sorted):
- cluster_0 through cluster_14 in proper numerical order
- cluster_4 has highest enrichment (1,152 pathways)

## 🗂️ No New Directories Created
All updates were applied directly to original files in:
- `/results/pd_signatures/by_gene/`
- `/results/pd_signatures/by_cluster/`
- `/results/pd_signatures/visualizations/comprehensive/`

Original file structure preserved, content improved with fixes.