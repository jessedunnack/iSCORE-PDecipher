# iSCORE-PDecipher

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/Version-0.1.3-brightgreen.svg)](https://github.com/jessedunnack/iSCORE-PDecipher)

**Integrated Analysis of Parkinson's Disease Mutations and Perturbations**

iSCORE-PDecipher is an R package for comprehensive analysis of Parkinson's disease research data, integrating three distinct experimental approaches: genetic mutations (iSCORE-PD), gene knockdowns (CRISPRi), and gene activations (CRISPRa). It provides tools for differential expression analysis, functional enrichment analysis, and professional interactive visualizations.

## 🧬 Three Distinct Dataset Collections

### 1. **iSCORE-PD Genetic Mutations**
- **13 PD-associated mutations** across 14 cell clusters
- **Genes**: ATP13A2, DNAJC6, FBXO7, GBA, LRRK2, PARK7, PINK1, PRKN, SNCA variants, SYNJ1, VPS13C variants
- **Analysis**: Mutant vs. isogenic wild-type comparisons using MAST
- **Results**: 211,470 enrichment terms

### 2. **CRISPRi Gene Knockdowns (PerturbSeq)**
- **10 PD genes** across 10 cell clusters, 3 experiments
- **Target genes**: ATP13A2, DNAJC6, FBXO7, LRRK2, PARK2, PARK7, PINK1, SNCA, SYNJ1, VPS13C
- **Analysis**: Perturbed vs. non-targeting controls using MixScale
- **Results**: 548,921 enrichment terms

### 3. **CRISPRa Gene Activations (PerturbSeq)**
- **10 PD genes** in 1 specialized cluster
- **Target genes**: Same as CRISPRi experiments
- **Analysis**: Activated vs. non-targeting controls using MixScale
- **Results**: 6,946 enrichment terms

## Installation

### From GitHub (Recommended)

```r
# Install devtools if needed
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install iSCORE-PDecipher
devtools::install_github("jessedunnack/iSCORE-PDecipher")
```

### Prerequisites

Install required Bioconductor and CRAN packages:

```r
# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install(c(
    "clusterProfiler", "ReactomePA", "DOSE", "org.Hs.eg.db", 
    "pathview", "SingleCellExperiment", "dittoSeq", 
    "enrichplot", "ComplexHeatmap"
))

# Install CRAN packages
install.packages(c(
    "heatmaply", "colourpicker", "shinyWidgets", 
    "shinycssloaders", "shinyjs", "ggridges"
))
```

## Quick Start

### 1. Launch Interactive Shiny App

```r
library(iSCORE.PDecipher)

# Launch with your enrichment data
launch_iscore_app("path/to/enrichment_results.rds")

# Or launch with file upload interface
launch_iscore_app()
```

### 2. Run Complete Analysis Pipeline

```r
# Process from raw differential expression to visualization
results <- run_complete_pipeline(
    mast_directory = "./iSCORE-PD_MAST_analysis/",
    mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
    output_directory = "./analysis_output/",
    launch_app = TRUE
)
```

### 3. Extract UMAP Visualizations

```r
# Extract lightweight UMAP data from large Seurat objects
source("inst/scripts/extract_umap_data.R")
main()  # Creates interactive UMAP visualizations
```

## 🖥️ Shiny App Features

### **Collapsible Global Settings** ⭐ NEW
- Streamlined leftward-collapsing sidebar for all settings
- Clean interface with icon representations when collapsed
- Synchronizes with all visualization modules

### **DE Results Page** ⭐ NEW
- **Interactive UMAP**: Click clusters to update volcano plots
- **Dual Volcano Plots**: MAST and MixScale results side-by-side
- **Dynamic Coloring**: By significance, experiment, or gene
- **Real-time Updates**: Instant synchronization between panels

### **Overview Dashboard**
- Interactive UMAP plots with cell cluster visualization
- Dataset metrics and cluster marker tables
- Method comparison charts

### **Visualizations**
- Interactive heatmaps with hierarchical clustering
- Dot plots and bar plots for pathway enrichment
- GSEA plots with NES visualizations

### **Data Exploration**
- Multi-level filtering by gene, cluster, experiment, direction
- Real-time updates and customizable statistical thresholds
- Export options for figures and data tables

## Data Format Requirements

### MAST Results (Genetic Mutations)
```
Required columns:
- mutation_tidy: Mutation identifier (e.g., "LRRK2", "GBA")
- cluster: Cell cluster (e.g., "cluster_0")
- gene: Gene symbol
- log2FC: Log2 fold change
- p_val_adj: Adjusted p-value
```

### MixScale Results (CRISPR Perturbations)
```
Required columns:
- scMAGeCK_gene_assignment: Target gene
- cluster: Cell cluster identifier
- experiment: Experiment ID (e.g., "C12_FPD-23")
- Dynamic columns: log2FC_*, p_cell_type*:weight
```

## Analysis Workflow

```
1. Data Import
   ┌─ MAST Results ──┐
   │                 ├─→ Combined DE Results
   └─ MixScale Results┘     (full_DE_results.rds)
                            
2. Enrichment Analysis                            
   DE Results ──→ GO/KEGG/Reactome/WikiPathways/STRING/GSEA
               ├─→ 14,052 individual result files
               └─→ Consolidated dataset (767K+ terms)

3. UMAP Processing
   Large Seurat Objects ──→ Lightweight SCE Objects
   (20-30GB each)           (100-500MB each)

4. Interactive Visualization
   Consolidated Data + UMAP ──→ Shiny App
                              ├─→ Interactive plots
                              ├─→ Statistical analysis
                              └─→ Export capabilities
```

## Troubleshooting

### Memory Issues
```r
# Increase memory limits for large datasets
options(java.parameters = "-Xmx16g")
gc()  # Garbage collection
```

### Missing Dependencies
```r
# Check for required packages
required_pkgs <- c("clusterProfiler", "heatmaply", "dittoSeq")
missing <- setdiff(required_pkgs, rownames(installed.packages()))
if (length(missing) > 0) BiocManager::install(missing)
```

## Citation

**Manuscript in preparation**

Please cite this repository until the manuscript is published:

```
Dunnack J. (2025). iSCORE-PDecipher: Integrated Analysis of 
Parkinson's Disease Mutations and Perturbations. 
GitHub: https://github.com/jessedunnack/iSCORE-PDecipher
```

## Contact

**Jesse Dunnack**  
PhD Student, Molecular and Cell Biology Department  
University of California, Berkeley  
Hockemeyer Lab + Bateup Lab  

📧 **Email:** jessedunnack@berkeley.edu | jessedunnack@gmail.com  
🔗 **ORCID:** [0000-0002-0387-0090](https://orcid.org/0000-0002-0387-0090)  
📍 **GitHub Issues:** [Report bugs or request features](https://github.com/jessedunnack/iSCORE-PDecipher/issues)

---

## ✨ Key Features (v0.1.3)

### 🎯 **Interactive UMAP Visualization**
- **Cell cluster exploration** with automatic dataset detection
- **Cluster marker genes** with interactive tables and statistical analysis
- **Lightweight data** (14MB vs 20-30GB Seurat objects)
- **Publication-quality plots** using dittoSeq integration

### 📊 **Advanced Interactive Heatmaps**
- **heatmaply integration** with hierarchical clustering and dendrograms
- **Multiple data types**: P-values, fold enrichment, z-scores, GSEA NES
- **Direction filtering**: ALL/UP/DOWN/BOTH regulated genes
- **Color customization**: 5 color scales with 3 scaling methods
- **Export options**: Interactive HTML and publication PDF formats

### 🧪 **GSEA Visualization Support**
- **Normalized Enrichment Score (NES)** heatmaps and plots
- **enrichplot integration** for static GSEA visualizations
- **Interactive filtering** with NES threshold controls
- **Ridge plots and dot plots** for gene set analysis

### ⚡ **Performance Optimizations**
- **50x faster startup** with centralized data management
- **Eliminated UI flickering** through reactive optimization
- **Memory efficient** data processing pipeline
- **Professional error handling** with informative feedback

### 🎨 **Enhanced User Interface**
- **Responsive design** with optimized space utilization
- **Professional styling** with consistent visual themes
- **Intuitive navigation** between analysis modules
- **Real-time feedback** and progress indicators