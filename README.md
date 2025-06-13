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

### **Overview Dashboard**
- **Interactive UMAP plots** with cell cluster visualization
- **Dataset metrics** showing analysis scope and coverage
- **Cluster marker tables** with differential expression statistics
- **Method comparison** charts (MAST vs MixScale results)

### **Advanced Visualizations**
- **Interactive heatmaps**: Cross-condition enrichment comparisons with clustering
- **Dot plots**: Gene set enrichment with size and color mapping
- **Bar plots**: Top enriched pathways with statistical significance
- **GSEA plots**: Enrichment plots, ridge plots, and NES visualizations

### **Data Exploration**
- **Multi-level filtering**: Gene, cluster, experiment, direction, enrichment type
- **Real-time updates**: Instant visualization refresh on parameter changes
- **Statistical thresholds**: Customizable p-value and effect size cutoffs
- **Cross-modal analysis**: Compare mutations vs perturbations

### **Export & Reports**
- **Publication figures**: High-resolution PNG/PDF export
- **Interactive downloads**: HTML heatmaps for presentations
- **Data tables**: Filtered results in CSV/Excel formats
- **Analysis reports**: Comprehensive summaries with statistics

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

## Example Analysis Results

**Current Dataset (v0.1.3):**
- **767,337 total enrichment terms** across all experiments
- **MAST mutations**: 13 mutations × 14 clusters = 211,470 terms
- **CRISPRi knockdowns**: 10 genes × 10 clusters × 3 experiments = 548,921 terms
- **CRISPRa activations**: 10 genes × 1 cluster × 1 experiment = 6,946 terms
- **Enrichment databases**: GO (BP/CC/MF), KEGG, Reactome, WikiPathways, STRING, GSEA
- **Cell populations**: 201,679 total cells analyzed across all datasets

## Research Applications

### Hypothesis Testing
- **Convergent pathways**: Identify pathways affected by both mutations and perturbations
- **Dosage effects**: Compare knockdown vs activation of the same genes
- **Cell-type specificity**: Analyze cluster-specific responses to perturbations
- **Disease mechanisms**: Map PD gene network interactions and dependencies

### Comparative Analysis
- **Method validation**: Cross-validate findings between genetic and CRISPR approaches
- **Effect magnitude**: Quantify differential responses across experimental conditions
- **Pathway enrichment**: Identify consistently dysregulated biological processes
- **Therapeutic targets**: Prioritize genes based on multi-modal evidence

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

### Performance Optimization
```r
# Use subset of data for testing
data_subset <- data[sample(nrow(data), 10000), ]
launch_iscore_app(data_subset)
```

## Citation

If you use iSCORE-PDecipher in your research, please cite:

```
Dunnack J, et al. (2025). iSCORE-PDecipher: Integrated Analysis of 
Parkinson's Disease Mutations and Perturbations. 
GitHub: https://github.com/jessedunnack/iSCORE-PDecipher
```

## Contributing

We welcome contributions! Please:
1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description
4. Ensure all tests pass

## License

GPL-3 License - see [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/jessedunnack/iSCORE-PDecipher/issues)
- **Documentation**: Package vignettes and help files
- **Updates**: Follow repository for latest features and fixes

## Acknowledgments

- **iSCORE-PD Consortium**: Genetic mutation data and experimental design
- **PerturbSeq Community**: CRISPR perturbation methodology and protocols  
- **Bioconductor Project**: Enrichment analysis infrastructure and tools
- **Shiny Community**: Interactive visualization framework and best practices

---

**Version 0.1.3** | Last updated: January 2025 | 🧬 Powered by R and Bioconductor