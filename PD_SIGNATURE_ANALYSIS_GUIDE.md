# Parkinson's Disease Signature Analysis Guide

## Overview

This guide explains how to use the enhanced PD-focused signature analysis tools that provide clearer, more biologically meaningful interpretations of your cross-method comparison results.

## What's New

### 1. PD-Focused Biological Interpretation
- **Filters results for PD-relevant pathways** (mitochondrial, autophagy, dopamine, etc.)
- **Categorizes biological processes** by their relevance to Parkinson's disease
- **Generates manuscript-ready summaries** with biological insights

### 2. Enhanced Shiny App Interface
- **New "PD Biology Focus" tab** in the Signature Nomination module
- **Interactive visualizations** showing pathway frequency and biological categories
- **Downloadable reports** for manuscript preparation

### 3. Standalone Analysis Scripts
- **Command-line analysis** for batch processing
- **Comprehensive output files** (CSV, TXT) for further analysis

## How to Use

### Option 1: Run Standalone Analysis (Recommended for Initial Exploration)

```r
# Navigate to your package directory
setwd("/path/to/iSCORE-PDecipher")

# Run the comprehensive PD analysis
source("analyze_current_results.R")

# This will automatically:
# 1. Load your enrichment data
# 2. Run signature discovery
# 3. Perform PD biological interpretation
# 4. Generate summary reports
# 5. Save results to 'pd_signature_analysis_results/'
```

**Expected Output Files:**
- `detailed_signature_analysis.csv` - Complete signature table with PD relevance scores
- `pd_pathway_summary.csv` - Most frequently disrupted PD pathways
- `manuscript_ready_summary.txt` - Publication-ready biological interpretation
- `individual_reports/` - Detailed reports for each signature

### Option 2: Use Enhanced Shiny App

```r
# Install latest package version
remotes::install_github("jessedunnack/iSCORE-PDecipher", force = TRUE)

# Launch app
library(iSCORE.PDecipher)
launch_iscore_app()

# Navigate to "Signature Nomination" tab
# Run analysis with "PD-relevant focus" checked
# View results in new "PD Biology Focus" tab
```

## Key Features Explained

### 1. PD Relevance Scoring

Each signature receives a **PD Relevance Score** based on:
- **Mitochondrial pathways** (weight: 3.0) - Core PD mechanism
- **Dopamine pathways** (weight: 3.0) - Central to PD pathophysiology  
- **Protein quality control** (weight: 2.5) - Major PD feature
- **Autophagy/lysosomal** (weight: 2.0) - Important PD mechanism
- **Other pathways** (weights: 1.5-2.0)

### 2. Biological Category Analysis

Results are automatically categorized into:
- **Mitochondrial dysfunction** - Energy production defects
- **Protein quality control** - Misfolding and aggregation
- **Autophagy/lysosomal** - Cellular cleanup mechanisms
- **Dopamine metabolism** - Direct PD relevance
- **Synaptic function** - Neurotransmission defects
- **Oxidative stress** - Cellular damage mechanisms
- **Neuronal function** - Neurodegenerative processes

### 3. Cross-Method Validation

The analysis identifies:
- **Shared pathways** between mutation and knockdown approaches
- **Biological convergence** across different perturbation types
- **Core disease mechanisms** vs method-specific artifacts

## Interpreting Results

### Key Metrics to Focus On:

1. **Signature Strength** (2.5-2.7 range is meaningful)
   - Measures overall similarity between MAST and CRISPRi
   - Higher values indicate stronger cross-method agreement

2. **PD Relevance Score** (>2.0 is biologically significant)
   - Weighted by importance to PD biology
   - Helps prioritize signatures for follow-up

3. **Biological Categories** (count > 0 indicates pathway involvement)
   - Shows which PD mechanisms are most affected
   - Guides mechanistic interpretation

### Example Interpretation:

```
Gene Pair: LRRK2_vs_LRRK2 (Cluster 3)
- Signature Strength: 2.7 (Strong cross-method agreement)
- PD Relevance: 4.2 (High PD significance)
- Mitochondrial: 3 pathways (Energy defects)
- Protein Quality: 2 pathways (Misfolding issues)
- Autophagy: 1 pathway (Clearance defects)

Biological Interpretation:
LRRK2 mutations and knockdowns show convergent disruption of mitochondrial 
function coupled with protein quality control defects, suggesting unified 
mechanisms where energy deficits and protein misfolding reinforce each other.
```

## Manuscript Integration

### For Methods Section:
"Cross-method signature analysis was performed using iSCORE-PDecipher to identify shared biological signatures between genetic mutations (MAST) and CRISPR knockdowns (CRISPRi). Signatures were scored for PD relevance based on pathway involvement in mitochondrial function, protein quality control, autophagy, and dopaminergic signaling."

### For Results Section:
Use the generated `manuscript_ready_summary.txt` which includes:
- **Cross-method validation results**
- **Biological insights** with pathway categories
- **Manuscript implications** 
- **Recommended follow-up experiments**

### For Figures:
The analysis generates publication-ready visualizations:
- **Pathway frequency plots** - Most disrupted PD pathways
- **Biological category breakdowns** - Process-level analysis
- **Signature strength heatmaps** - Cross-method comparisons

## Troubleshooting

### If Analysis Returns No Results:
1. **Check data availability** - Ensure both MAST and CRISPRi data present
2. **Lower thresholds** - Reduce min_cluster_breadth to 6 or lower
3. **Check gene pairs** - Verify comparable genes exist in both methods

### If PD Analysis Fails:
1. **Ensure "PD-relevant focus" is checked** in analysis scope
2. **Check package installation** - Reinstall if functions missing
3. **Verify data format** - Ensure enrichment data has required columns

### Low PD Relevance Scores:
- **Normal for some genes** - Not all PD genes affect core pathways equally
- **Focus on relative ranking** - Compare signatures within your dataset
- **Consider cluster-specific effects** - Some signatures may be cell-type specific

## Advanced Usage

### Custom PD Pathway Lists:
```r
# Modify PD-relevant pathways
custom_pathways <- c("your_custom_pathway_terms")
pd_analysis <- analyze_pd_signatures(
  signature_results, 
  enrichment_data,
  custom_pd_pathways = custom_pathways
)
```

### Focused Gene Pair Analysis:
```r
# Analyze specific gene pairs only
target_pairs <- c("LRRK2_vs_LRRK2", "PINK1_vs_PINK1")
filtered_results <- signature_results[signature_results$gene_pair %in% target_pairs, ]
```

### Export for External Tools:
```r
# Export for pathway analysis tools
write.csv(detailed_table, "for_ingenuity_pathway_analysis.csv")
write.csv(pathway_summary, "for_david_functional_annotation.csv")
```

## Summary

The enhanced PD signature analysis provides:

✅ **Clear biological interpretation** of technical results  
✅ **PD-specific pathway focus** for relevant insights  
✅ **Manuscript-ready summaries** for publication  
✅ **Interactive visualizations** for exploration  
✅ **Cross-method validation** for robust conclusions  
✅ **Downloadable reports** for further analysis  

This transforms your technical signature strengths (2.5-2.7 range) into actionable biological insights about mitochondrial dysfunction, protein quality control defects, and other core PD mechanisms.