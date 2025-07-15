# Cross-Method Correlation Analysis Guide

## Overview

This guide describes the gene filtering methodology for improving correlations between MAST (mutations) and CRISPRi (knockdown) results.

## Methodology

### Traditional Approach
- **Method**: Correlate log2FC values for all genes
- **Results**: mean |r| = 0.053
- **Issue**: Thousands of unchanged genes dilute the signal

### Improved Approach
- **Method**: Filter to top N most-changed genes per method
- **Results**: mean |r| = 0.593 (11x improvement with top 200 genes)
- **Advantage**: Focuses on genes most affected by each perturbation

## Gene Filtering Options

### Top 100 genes
- **Correlation**: mean |r| = 0.601
- **Advantage**: Highest correlation strength
- **Limitation**: Fewer overlapping genes

### Top 200 genes (recommended)
- **Correlation**: mean |r| = 0.593
- **Advantage**: Optimal balance of correlation strength and data points
- **Usage**: Default setting in interactive plots

### Top 500 genes
- **Correlation**: mean |r| = 0.487
- **Advantage**: More genes included in analysis
- **Trade-off**: Slightly weaker correlations

## Implementation

### Interactive Plots
Gene filtering is available in the Signature Nomination module:
1. Navigate to "Gene Pair Analysis" 
2. Select filtering approach from dropdown
3. View updated correlation plots with trend lines

### Analysis Scripts
Correlation analysis scripts are located in:
- `inst/scripts/correlation_analysis/comprehensive_correlation_analysis.R`
- `inst/scripts/correlation_analysis/test_correlation_approaches.R`

## Results

### Validation
- **Combinations tested**: 180 (12 gene pairs × 3 experiments × 5 clusters)
- **Approaches compared**: 7 different filtering strategies
- **Strong correlations**: 61 gene pairs with |r| ≥ 0.5

### Notable Findings
- **DNAJC6**: r = 0.99 (strongest correlation observed)
- **VPS13C variants**: Consistent strong correlations across experiments
- **SNCA variants**: Strong correlations with directional consistency

## Interpretation

Strong correlations indicate biological pathway convergence between:
- Genetic mutations (MAST analysis)
- Gene knockdowns (CRISPRi analysis)

This validates that both perturbation methods affect similar downstream pathways, supporting the biological relevance of observed effects.

## References

Complete analysis results are available in:
- `inst/results/comprehensive_correlation_results.csv`
- `inst/results/correlation_quality_analysis.csv`