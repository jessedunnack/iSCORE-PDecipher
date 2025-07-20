# PD Signature Analysis Scripts

This directory contains analysis scripts for identifying compelling Parkinson's Disease signatures across mutation and CRISPRi perturbation datasets.

## Scripts Overview

### Primary Analysis Scripts

1. **pd_signature_discovery_fast.R**
   - Fast discovery of PD-relevant pathways from enrichment data
   - Filters 663K+ enrichments to ~45K PD-relevant pathways
   - Identifies method-specific and convergent signatures

2. **pd_signature_by_gene_analysis.R**
   - Analyzes each PD gene individually
   - Generates gene-specific reports and visualizations
   - Creates summary heatmap across all genes

3. **pd_signature_cluster_analysis.R**
   - Analyzes enrichment patterns across all clusters
   - Identifies cluster-specific vs ubiquitous pathways
   - Creates cluster × pathway matrices

4. **pd_signature_comprehensive_viz.R**
   - Creates multi-panel visualizations
   - Generates manuscript-ready figures
   - Produces summary statistics

### Visualization Scripts

5. **create_quick_plots.R**
   - Generates presentation-ready plots
   - Uses updated labeling: "Mutation - iSCORE-PD" and "CRISPRi Perturbation"
   - Creates top pathway bar charts and comparison figures

6. **pd_signature_visualization.R**
   - Original visualization script with detailed plots
   - Creates comprehensive visual summaries

### Reporting Scripts

7. **pd_signature_comprehensive_report.Rmd**
   - Complete R Markdown report with all findings
   - Interactive tables and figures
   - Executive summary and recommendations

8. **render_comprehensive_report.R**
   - Helper script to render the Rmd report
   - Handles pandoc requirements

## Usage

All scripts expect data to be in:
```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/all_enrichment_padj005_complete_with_direction.rds
```

Results are saved to:
```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/update_analysis_scripts/iSCORE-PDecipher/results/pd_signatures/
```

## Key Findings

- **16 PD genes analyzed** (including SNCA variants)
- **15 clusters analyzed** (cluster_0 through cluster_14)
- **30 mutation-only pathways** identified
- **30 CRISPRi-only pathways** identified
- **30 convergent pathways** validated across methods

## Version

Created for iSCORE-PDecipher v0.3.0