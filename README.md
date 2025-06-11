# iSCORE-PDecipher

[![R](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

**Integrated Analysis of Parkinson's Disease Mutations and Perturbations**

iSCORE-PDecipher is an R package for comprehensive analysis of Parkinson's disease research data, combining iSCORE-PD mutation studies with PerturbSeq (CRISPRi/CRISPRa) experiments. It provides tools for differential expression analysis, functional enrichment analysis, and interactive visualization.

## Features

### 🧬 Multi-Modal Data Integration
- **iSCORE-PD mutations**: 13 mutations across 14 cell clusters
- **CRISPRi perturbations**: Gene knockdown experiments across 10 clusters  
- **CRISPRa perturbations**: Gene activation experiments

### 📊 Comprehensive Analysis Pipeline
- **Differential Expression**: MAST and MixScale analysis
- **Functional Enrichment**: GO, KEGG, Reactome, WikiPathways, STRING, GSEA
- **Result Consolidation**: 767,337+ enrichment terms processed
- **Quality Control**: Automated validation and error handling

### 🎯 Interactive Visualization
- **Shiny Application**: Professional web interface
- **Heatmap Visualizations**: Cross-condition comparisons
- **Export Capabilities**: Publication-ready figures
- **Modality Comparison**: Compare mutations vs. perturbations

## Installation

### From GitHub (Recommended)

```r
# Install devtools if you haven't already
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")

# Install iSCORE-PDecipher
devtools::install_github("jessedunnack/iSCORE-PDecipher")
```

### Prerequisites

The package requires several Bioconductor packages:

```r
# Install BiocManager if needed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install required Bioconductor packages
BiocManager::install(c(
    "clusterProfiler",
    "ReactomePA", 
    "DOSE",
    "org.Hs.eg.db",
    "pathview"
))
```

## Quick Start

### 1. Launch the Shiny App

```r
library(iSCORE.PDecipher)

# Launch with your data file
launch_iscore_app("path/to/your/enrichment_results.rds")

# Or launch with file upload interface
launch_iscore_app()
```

### 2. Run Complete Analysis Pipeline

```r
# Run the full pipeline from DE results to visualization
results <- run_complete_pipeline(
    mast_directory = "./iSCORE-PD_MAST_analysis/",
    mixscale_directory = "./PerturbSeq_MixScale_analysis_full_dataset/",
    output_directory = "./analysis_output/",
    launch_app = TRUE
)
```

### 3. Import and Process Data

```r
# Import MAST results
mast_data <- import_mast_data("./iSCORE-PD_MAST_analysis/")

# Import CRISPRi results
crispi_data <- import_mixscale_data("./PerturbSeq_data/", modality = "CRISPRi")

# Run enrichment analysis
enrichment_results <- run_enrichment_analysis(
    input_file = "full_DE_results.rds",
    output_dir = "./enrichment_results/",
    run_methods = c("GO", "KEGG", "Reactome")
)
```

## Data Format

### Input Requirements

The package expects differential expression results in RDS format with specific column names:

**MAST Results:**
- `mutation_tidy`: Mutation identifier
- `cluster`: Cell cluster identifier
- `gene`: Gene symbol
- `log2FC`: Log2 fold change
- `p_val_adj`: Adjusted p-value

**MixScale Results:**
- `scMAGeCK_gene_assignment`: Perturbed gene
- `cluster`: Cell cluster identifier  
- `experiment`: Experiment identifier
- Log2FC columns (e.g., `log2FC_C12_FPD-23`)
- P-value columns (e.g., `p_cell_typeC12_FPD-23:weight`)

### Output Structure

The package generates:
- `full_DE_results.rds`: Combined differential expression results
- `enrichment_results/`: Directory with enrichment analysis outputs
- `all_enrichment_complete.rds`: Consolidated enrichment results (767K+ terms)

## Shiny App Features

### Data Exploration
- **Gene Selection**: Choose genes/mutations of interest
- **Cluster Filtering**: Select specific cell clusters
- **Enrichment Types**: Filter by GO, KEGG, Reactome, etc.
- **Direction Analysis**: Up/down-regulated gene sets

### Visualizations
- **Heatmaps**: Cross-condition enrichment comparisons
- **Bar Plots**: Top enriched pathways
- **Network Plots**: Pathway relationship networks
- **KEGG Pathview**: Pathway visualization with gene overlays

### Export Options
- **Data Export**: CSV/Excel format results
- **Figure Export**: PNG/PDF publication figures
- **Report Generation**: Comprehensive analysis reports

## Analysis Workflow

### 1. Data Import
```
MAST Results → import_mast_data()
MixScale Results → import_mixscale_data()
                ↓
Combined DE Results (full_DE_results.rds)
```

### 2. Enrichment Analysis
```
DE Results → run_enrichment_analysis()
           ↓
GO/KEGG/Reactome/WikiPathways/STRING/GSEA
           ↓
Enrichment Results (14,052 files)
```

### 3. Result Processing
```
Individual Results → process_enrichment_results()
                   ↓
Consolidated Dataset (767,337 terms)
```

### 4. Visualization
```
Consolidated Data → Shiny App
                  ↓
Interactive Visualizations & Exports
```

## Example Dataset

The package has been tested with:
- **767,337 total enrichment terms**
- **MAST**: 13 mutations × 14 clusters (211,470 terms)
- **CRISPRi**: 10 genes × 10 clusters (548,921 terms)  
- **CRISPRa**: 10 genes × 1 cluster (6,946 terms)

## Troubleshooting

### Common Issues

**Memory Errors:**
```r
# Increase memory limit if needed
options(java.parameters = "-Xmx8g")
```

**Missing Dependencies:**
```r
# Check for missing packages
missing_pkgs <- setdiff(
    c("clusterProfiler", "ReactomePA", "org.Hs.eg.db"),
    rownames(installed.packages())
)
if (length(missing_pkgs) > 0) {
    BiocManager::install(missing_pkgs)
}
```

**File Path Issues:**
- Use absolute paths for data files
- Ensure proper directory structure for input data
- Check file permissions

## Citation

If you use iSCORE-PDecipher in your research, please cite:

```
[Your Publication Citation Here]
```

## Contributing

We welcome contributions! Please see our [Contributing Guidelines](CONTRIBUTING.md) for details.

## License

This project is licensed under the GPL-3 License - see the [LICENSE](LICENSE) file for details.

## Support

- **Issues**: [GitHub Issues](https://github.com/jessedunnack/iSCORE-PDecipher/issues)
- **Documentation**: See package vignettes
- **Contact**: jesse.dunnack@example.com

## Acknowledgments

- iSCORE-PD consortium for mutation data
- PerturbSeq methodology for CRISPR perturbation data
- Bioconductor community for enrichment analysis tools