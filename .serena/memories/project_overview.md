# iSCORE-PDecipher Project Overview

## Purpose
iSCORE-PDecipher is an R package (v0.3.5) for integrated analysis of Parkinson's disease mutations and perturbations. It combines iSCORE-PD mutation data with PerturbSeq (CRISPRi/CRISPRa) experiments for comprehensive differential expression and functional enrichment analysis.

## Tech Stack
- **Language**: R (>= 4.0.0)
- **Framework**: Standard R package with roxygen2 documentation
- **UI**: Shiny dashboard applications
- **Dependencies**: Seurat, clusterProfiler, ReactomePA, DOSE, plotly, DT, and 30+ other packages
- **Testing**: testthat (>= 3.0.0) listed in Suggests but NOT implemented

## Core Functionality
1. **MAST Analysis**: Differential expression for mutations vs controls
2. **MixScale Analysis**: PerturbSeq experiment analysis
3. **Enrichment Analysis**: GO, KEGG, Reactome, WikiPathways
4. **Interactive Visualization**: Shiny apps for exploration
5. **Signature Discovery**: PD-relevant pathway identification

## Repository Structure
- `R/`: 25+ source files with core functions
- `inst/shiny/`: Interactive Shiny applications  
- `inst/extdata/`: Data files and UMAP data
- `man/`: Generated documentation
- `dev/`: Development scripts and informal tests
- `docs/`: Extensive project documentation
- `results/`: Analysis outputs and visualizations

## Development Status
- **CRITICAL ISSUE**: No formal tests/ directory - major stability risk
- **NAMESPACE Issue**: Core functions tagged @export but missing from NAMESPACE
- **Version**: 0.3.5 with active development
- **Documentation**: Extensive but informal testing only