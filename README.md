# iSCORE-PDecipher

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![Version](https://img.shields.io/badge/Version-0.2.8-brightgreen.svg)](https://github.com/jessedunnack/iSCORE-PDecipher)

**Integrated Analysis of Parkinson's Disease Mutations and Perturbations**

iSCORE-PDecipher is an R package for comprehensive analysis of Parkinson's disease research data, integrating genetic mutations (iSCORE-PD) and gene knockdown perturbations (CRISPRi). The package provides tools for differential expression analysis, functional enrichment analysis, cross-method signature discovery, and interactive visualizations.

## Latest Updates (v0.2.8)

**Critical bug fixes for improved Shiny app stability:**
- ✅ Gene pair details modal now displays actual DE genes and enrichment terms
- ✅ Fixed [object Object] errors in Cross-Cluster DE Gene Analysis
- ✅ Enhanced debugging for CRISPRi dotplot rendering issues
- ✅ DE Gene Overlap Heatmap now interactive with heatmaply
- ✅ CRISPRa data properly excluded from all analyses
- ✅ Gene dropdown formatting cleaned (removed "X variants" text)

## Enhanced Interactive Correlation Plots (v0.2.7)

**Major visualization improvements for correlation analysis:**
- **Vertical stacking layout**: Replaced compressed grid with scrollable vertical plots for better readability
- **Correlation statistics display**: r-values, p-values, and sample sizes shown in plot titles
- **Bold reference lines**: x=0 and y=0 lines for easy identification of effect directions
- **Gene filtering methodology**: Top N most changed genes (default: 200) from both methods

**Performance**: Gene filtering improves correlations by 11x (mean |r| = 0.593 vs 0.053)
**Strong correlations**: 61 gene pairs with |r| ≥ 0.5, including DNAJC6 (r=0.99)

[See correlation analysis documentation](docs/correlation_analysis_guide.md) for methodology.

## Cross-Platform Compatibility

The package supports Windows, Mac, and Linux systems with automatic setup and data transfer utilities for seamless cross-platform deployment.

## Primary Dataset Collections

### 1. iSCORE-PD Genetic Mutations
- 13 PD-associated mutations across 14 cell clusters
- Target genes: ATP13A2, DNAJC6, FBXO7, GBA, LRRK2, PARK7, PINK1, PRKN, SNCA variants, SYNJ1, VPS13C variants
- Analysis methodology: Mutant vs. isogenic wild-type comparisons using MAST
- Results: ~211,470 enrichment terms

### 2. CRISPRi Gene Knockdowns
- 10 PD genes across 10 cell clusters, 3 experiments
- Target genes: ATP13A2, DNAJC6, FBXO7, LRRK2, PARK2, PARK7, PINK1, SNCA, SYNJ1, VPS13C
- Analysis methodology: Perturbed vs. non-targeting controls using MixScale
- Results: ~450,000 enrichment terms

### 3. Cross-Method Signature Analysis (v0.2.0)
- Signature discovery between MAST mutations and CRISPRi knockdowns
- PD-focused biological interpretation with pathway categorization
- Pan-cluster and cluster-specific signatures for comprehensive analysis
- Interactive signature nomination interface for rapid analysis

**Total Dataset**: ~663,000+ significant enrichment terms (p.adjust < 0.05)

Note: CRISPRa activation data included in some analyses but de-prioritized for primary workflows.

## Installation

### From GitHub (Recommended)

```r
# Install remotes if needed (recommended over devtools)
if (!require("remotes", quietly = TRUE))
    install.packages("remotes")

# Install iSCORE-PDecipher (works on Windows, Mac, Linux)
remotes::install_github("jessedunnack/iSCORE-PDecipher")
```

**Platform Requirements:**
- Windows: R 4.0+ and RTools
- Mac: R 4.0+ and Xcode command line tools  
- Linux: R 4.0+ and build-essential
- First Launch: App will prompt for data directory location and save configuration

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

## Shiny App Features

### Collapsible Global Settings (NEW)
- Streamlined leftward-collapsing sidebar for all settings
- Clean interface with icon representations when collapsed
- Synchronizes with all visualization modules

### Signature Nomination Module (v0.2.0)
- Cross-method signature discovery between MAST mutations and CRISPRi knockdowns
- PD Biology Focus tab with biological interpretation and pathway categorization
- Pan-cluster signatures showing effects across multiple cell types
- Cluster-specific signatures for cell type-specific analysis
- Interactive heatmaps for signature visualization and comparison
- Signature Trends Analysis - Data-driven discovery of most frequent and impactful signatures without manual curation

### DE Results Page (Enhanced)
- Interactive UMAP: Click clusters to update volcano plots
- Dual Volcano Plots: MAST and MixScale results side-by-side
- Dynamic Coloring: By significance, experiment, or gene
- Real-time Updates: Instant synchronization between panels

### Overview Dashboard
- Interactive UMAP plots with cell cluster visualization
- Dataset metrics and cluster marker tables
- Method comparison charts

### Visualizations
- Interactive heatmaps with hierarchical clustering
- Dot plots and bar plots for pathway enrichment
- GSEA plots with NES visualizations

### Data Exploration
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

### Cross-Platform Issues
```r
# Reset configuration if setup fails
unlink(iSCORE.PDecipher:::get_config_path())
launch_iscore_app()  # Will prompt for setup again

# Check configuration status
iSCORE.PDecipher:::is_first_launch()
iSCORE.PDecipher:::get_parent_data_dir()
```

## Documentation

### Quick Start Resources
- [Complete Documentation Index](docs/DOCUMENTATION_INDEX.md) - Organized guide to all documentation
- [Labmate Quickstart Guide](docs/LABMATE_QUICKSTART.md) - Step-by-step setup for new users

### Cross-Platform Setup
- [Mac Compatibility Guide](docs/MAC_COMPATIBILITY_GUIDE.md) - Complete setup instructions for Mac users
- [Mac Setup Checklist](docs/MAC_SETUP_CHECKLIST.md) - Quick reference for data transfer and setup
- [Cross-Platform Development Guidelines](docs/CROSS_PLATFORM_DEVELOPMENT_GUIDELINES.md) - For developers contributing to the package

### Analysis & Features
- [PD Signature Analysis Guide](docs/PD_SIGNATURE_ANALYSIS_GUIDE.md) - Guide to signature nomination and biological interpretation
- [Future Enhancements](docs/FUTURE_ENHANCEMENTS.md) - Planned data-driven discovery features

### Advanced Documentation  
- [Project Documentation (CLAUDE.md)](CLAUDE.md) - Complete project overview and implementation details
- [All Documentation](docs/) - Complete documentation directory

## Citation

Please cite this repository:

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

Email: jessedunnack@berkeley.edu | jessedunnack@gmail.com  
ORCID: [0000-0002-0387-0090](https://orcid.org/0000-0002-0387-0090)  
GitHub Issues: [Report bugs or request features](https://github.com/jessedunnack/iSCORE-PDecipher/issues)

---

## Key Features (v0.2.7)

### Enhanced Interactive Correlation Plots (v0.2.7)
- **Revolutionary vertical layout**: Replaced grid with scrollable stacked plots for optimal readability
- **Statistical transparency**: r-values, p-values, and sample sizes displayed in plot titles
- **Visual reference system**: Bold x=0, y=0 lines for effect direction identification
- **Gene filtering methodology**: Top N most changed genes (11x correlation improvement)
- **Independent scaling**: Each cluster plot has optimal axis ranges
- **Dynamic sizing**: Automatic height adjustment based on cluster count
- **Enhanced user experience**: Natural scrolling through clusters 0, 1, 2...

### Cross-Method Signature Analysis (v0.2.0)
- Shared signature discovery between MAST mutations and CRISPRi knockdowns
- PD-focused biological interpretation with automated pathway categorization
- Prioritization of strongest cross-method signatures
- Interactive signature nomination interface for rapid hypothesis generation
- Data-driven signature trends analysis with frequency rankings and impact scoring for unbiased discovery

### Interactive UMAP Visualization
- Cell cluster exploration with automatic dataset detection
- Cluster marker genes with interactive tables and statistical analysis
- Lightweight data (14MB vs 20-30GB Seurat objects)
- Publication-quality plots using dittoSeq integration

### Advanced Interactive Heatmaps
- heatmaply integration with hierarchical clustering and dendrograms
- Multiple data types: P-values, fold enrichment, z-scores, GSEA NES
- Direction filtering: ALL/UP/DOWN/BOTH regulated genes
- Color customization: 5 color scales with 3 scaling methods
- Export options: Interactive HTML and publication PDF formats

### GSEA Visualization Support
- Normalized Enrichment Score (NES) heatmaps and plots
- enrichplot integration for static GSEA visualizations
- Interactive filtering with NES threshold controls
- Ridge plots and dot plots for gene set analysis

### Performance Optimizations
- 50x faster startup with centralized data management
- Eliminated UI flickering through reactive optimization
- Memory efficient data processing pipeline
- Professional error handling with informative feedback

### Enhanced User Interface
- Responsive design with optimized space utilization
- Professional styling with consistent visual themes
- Intuitive navigation between analysis modules
- Real-time feedback and progress indicators