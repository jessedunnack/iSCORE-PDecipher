# iSCORE-PDecipher Quick Start Guide for Lab Members

## 📋 What You Need

### Google Drive Files to Download
Download these **2 essential files** from Jesse's Google Drive:

1. **`all_enrichment_padj005_complete_with_direction.rds`** (266 MB)
   - Main enrichment results: 767,337 terms across all experiments
   - Contains MAST mutations + CRISPRi + CRISPRa data

2. **`all_umap_data_combined.rds`** (14 MB)  
   - UMAP visualizations for all three datasets
   - Lightweight extracts from 20-30GB Seurat objects

**Total download:** ~280 MB

## 🚀 Quick Setup (5 minutes)

### Step 1: Install the Package from GitHub

Open RStudio and run:

```r
# Install devtools if needed
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install iSCORE-PDecipher
devtools::install_github("jessedunnack/iSCORE-PDecipher")
```

### Step 2: Install Dependencies

```r
# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required packages (this may take 5-10 minutes)
BiocManager::install(c(
    "clusterProfiler", "ReactomePA", "DOSE", "org.Hs.eg.db", 
    "pathview", "SingleCellExperiment", "dittoSeq", 
    "enrichplot", "ComplexHeatmap"
))

install.packages(c(
    "heatmaply", "colourpicker", "shinyWidgets", 
    "shinycssloaders", "shinyjs", "ggridges"
))
```

### Step 3: Set Up Your Data Directory

1. Create a folder on your computer (e.g., `~/PD_Analysis/`)
2. Place the 2 downloaded RDS files in this folder
3. Note the full path to this folder

### Step 4: Launch the App

```r
library(iSCORE.PDecipher)

# Launch with your data (replace with your actual path)
launch_iscore_app("~/PD_Analysis/all_enrichment_padj005_complete_with_direction.rds")
```

**That's it!** The app will open in your web browser.

## 🎯 What You Can Explore

### Overview Tab
- **Interactive UMAP**: Cell clusters colored by type
- **Dataset metrics**: 767K terms, 201K cells, 13 mutations, 10 CRISPR genes
- **Cluster markers**: Click clusters to see marker genes

### Basic Visualization Tab
- **Dot plots**: Size = gene count, color = significance
- **Bar plots**: Top enriched pathways
- **GSEA plots**: For gene set enrichment analysis

### Heatmap Visualization Tab  
- **Interactive heatmaps**: Compare across conditions
- **Multiple data types**: P-values, fold enrichment, z-scores
- **Filtering**: UP/DOWN genes, specific databases
- **Export**: HTML (interactive) and PDF (publication)

### Method Comparison Tab
- **MAST vs CRISPRi**: Compare mutations with knockdowns
- **Convergent pathways**: Find overlapping results

## 🔍 Key Things to Try

### 1. Compare Different Genes
- **Global Settings sidebar**: Change gene from LRRK2 to GBA, PINK1, etc.
- See how different PD genes affect different pathways

### 2. Explore Cell Clusters  
- **Overview tab**: Different clusters respond differently
- **Cluster markers**: See what makes each cluster unique

### 3. Find Convergent Biology
- **Method Comparison**: Do mutations and knockdowns hit the same pathways?
- **Different directions**: UP vs DOWN regulated genes

### 4. Export Results
- **Heatmaps**: Download interactive HTML or publication PDF
- **Data tables**: Export filtered results to Excel
- **Plots**: Save high-resolution figures

## 🛠️ Troubleshooting

### App Won't Launch
```r
# Check if packages installed correctly
library(iSCORE.PDecipher)
library(shiny)
library(heatmaply)
```

### Memory Issues
```r
# Increase memory if needed
options(java.parameters = "-Xmx8g")
gc()  # Garbage collection
```

### File Path Issues
```r
# Check your file path
file.exists("~/PD_Analysis/all_enrichment_padj005_complete_with_direction.rds")
# Should return TRUE

# If FALSE, find the correct path:
file.choose()  # Opens file browser
```

### Missing Dependencies
```r
# Check for missing packages
required_pkgs <- c("clusterProfiler", "heatmaply", "dittoSeq")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) {
    BiocManager::install(missing)
}
```

## 📊 Understanding the Data

### Three Experimental Approaches
1. **iSCORE-PD Mutations** (211K terms)
   - 13 mutations: LRRK2, GBA, PINK1, PRKN, etc.
   - Mutant vs wild-type comparisons
   - 14 cell clusters

2. **CRISPRi Knockdowns** (549K terms)  
   - 10 genes knocked down with CRISPR
   - Perturbed vs non-targeting controls
   - 10 cell clusters, 3 experiments

3. **CRISPRa Activations** (7K terms)
   - 10 genes activated with CRISPR
   - 1 specialized cluster

### Enrichment Databases
- **GO**: Biological Process, Cellular Component, Molecular Function
- **KEGG**: Metabolic and signaling pathways  
- **Reactome**: Curated biological pathways
- **WikiPathways**: Community-curated pathways
- **STRING**: Protein interaction networks
- **GSEA**: Gene set enrichment analysis

## 💡 Research Tips

### Start Simple
1. Pick one gene (e.g., LRRK2)
2. Look at one cluster (e.g., cluster_0)  
3. Try GO_BP enrichment with UP-regulated genes

### Ask Biological Questions
- Which pathways are consistently affected across methods?
- Do mutations and CRISPR perturbations show similar effects?
- Are there cluster-specific responses?

### Use Filters Strategically
- **P-value threshold**: Start with 0.05, try 0.01 for stringency
- **Direction**: UP/DOWN for specific hypotheses, ALL for overview
- **Enrichment type**: GO_BP for processes, KEGG for pathways

## 📞 Getting Help

**Jesse Dunnack**  
📧 jessedunnack@berkeley.edu | jessedunnack@gmail.com  
🔗 GitHub Issues: [Report bugs](https://github.com/jessedunnack/iSCORE-PDecipher/issues)

### Quick Questions
- Slack Jesse for immediate help
- Email for detailed analysis questions

### Bug Reports  
- Use GitHub Issues for reproducible problems
- Include error messages and what you were trying to do

---

**Happy exploring!** 🧬 This tool contains 767,337 enrichment terms from comprehensive PD research - there's a lot to discover!