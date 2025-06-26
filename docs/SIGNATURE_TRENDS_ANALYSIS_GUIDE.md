# Signature Trends Analysis Guide

## Overview

The Signature Trends Analysis module provides **data-driven discovery** of signature patterns without manual curation bias. Instead of pre-defining biological pathways to focus on, this module analyzes all discovered signatures to identify the most frequently occurring patterns and highest impact signatures across all conditions.

## Key Features

### 🔍 **Frequency Analysis**
- **Most frequent signatures**: Identifies gene pairs that consistently appear across multiple conditions
- **Pattern breadth**: Shows how many clusters each signature appears in
- **Frequency scoring**: Normalized scores from 0-1 indicating relative frequency

### ⚡ **Impact Analysis**
- **Statistical strength ranking**: Orders signatures by composite strength scores
- **Impact scoring**: Normalized metrics showing relative statistical significance
- **Cross-method prioritization**: Highlights signatures with strongest cross-method agreement

### 📊 **Term Pattern Discovery**
- **Most common enrichment terms**: Identifies biological themes appearing most frequently
- **Pattern categorization**: Automatic classification into biological categories:
  - **Mitochondrial**: Respiratory chain, oxidative phosphorylation
  - **Autophagy**: Lysosomal function, cellular cleanup
  - **Protein Quality**: Protein folding, proteasome, chaperones
  - **Neurotransmission**: Dopamine, neurotransmitter metabolism
  - **Synaptic**: Synaptic function, vesicle transport

## User Interface

### Analysis Parameters
- **Minimum Frequency**: Sets threshold for signature inclusion (recommended: 2-5)
- **Top Results**: Number of top-ranked signatures to display (recommended: 50-100)

### Interactive Visualizations
- **Frequency Distribution**: Bar plots showing signature occurrence patterns
- **Impact Distribution**: Distribution of signature strength scores
- **Pattern Categories**: Breakdown of enrichment term categories
- **Frequency vs Impact Scatter**: Interactive comparison of frequency and impact metrics

### Metric Explanations

#### Signature Strength
**Definition**: Composite score combining gene overlap count, pathway overlap count, and statistical significance (Fisher's p-value)
**Range**: 0-10+
**Interpretation**: 
- >2.0 = Strong cross-method agreement
- >1.0 = Moderate agreement  
- <1.0 = Weak agreement

#### Frequency Score  
**Definition**: How often a signature pattern appears across different cluster/condition combinations, normalized to 0-1 scale
**Range**: 0-1
**Interpretation**:
- >0.8 = Very common patterns
- >0.5 = Common patterns
- <0.5 = Rare patterns

#### Impact Score
**Definition**: Statistical strength normalized to the strongest signature in the dataset, based on significance and effect size
**Range**: 0-1  
**Interpretation**:
- >0.9 = Highest impact signatures
- >0.7 = High impact signatures
- <0.5 = Moderate impact signatures

## Analysis Workflow

### 1. **Input Data Requirements**
- **Signature discovery results** from the Signature Nomination module
- **Consolidated enrichment data** from the main app dataset
- Minimum of 10+ discovered signatures for meaningful trends analysis

### 2. **Running the Analysis**
1. Navigate to **Signature Nomination** → **Signature Trends Analysis**
2. Set analysis parameters (frequency threshold, top N results)
3. Click **"Run Trends Analysis"** button
4. Wait for processing completion (typically 30-60 seconds)

### 3. **Interpreting Results**

#### Frequency Tab
- Shows signatures that appear most consistently across conditions
- **Use for**: Identifying robust, reproducible patterns for manuscript focus
- **Key metric**: Frequency Count and Cluster Breadth

#### Impact Tab  
- Shows signatures with highest statistical significance and effect size
- **Use for**: Finding the most statistically compelling signatures
- **Key metric**: Signature Strength and Impact Score

#### Term Patterns Tab
- Shows most commonly occurring enrichment terms across all signatures
- **Use for**: Understanding dominant biological themes without bias
- **Key metric**: Term frequency and pattern categories

#### Visualizations Tab
- Interactive plots for comprehensive pattern exploration
- **Use for**: Discovering relationships between frequency and impact
- **Key plots**: Scatter plots and heatmaps

### 4. **Export and Documentation**
- All tables can be downloaded as CSV files
- Plots can be exported as high-resolution images
- Results include analysis metadata and timestamp

## Best Practices

### Parameter Selection
- **Start with frequency = 2**: Captures most patterns while filtering noise
- **Use top N = 50**: Balances comprehensiveness with focus
- **Increase frequency for robust patterns**: Use 3-5 for only very consistent signatures

### Result Interpretation
- **Focus on frequency + impact intersection**: Signatures high in both metrics are most valuable
- **Consider biological relevance**: Validate data-driven findings with literature
- **Use for hypothesis generation**: Let data guide rather than confirm existing hypotheses

### Manuscript Development
- **Report top 10-20 signatures**: From combined frequency and impact rankings
- **Document methodology**: Emphasize unbiased, data-driven approach
- **Compare with literature**: Validate novel findings with existing PD research

## Technical Details

### Algorithm Overview
1. **Data Validation**: Comprehensive input checking and error handling
2. **Frequency Calculation**: Gene pair occurrence counting across conditions
3. **Impact Scoring**: Statistical strength normalization and ranking
4. **Pattern Discovery**: Enrichment term frequency analysis and categorization
5. **Visualization**: Interactive plot generation with hover explanations

### Performance Notes
- Analysis typically completes in 30-60 seconds
- Handles datasets with 100-1000+ signatures efficiently
- Memory usage scales linearly with signature count
- Interactive visualizations optimized for datasets up to 10,000 terms

### Error Handling
- Comprehensive validation of input data structure
- Graceful handling of missing or incomplete data
- Informative error messages for troubleshooting
- Fallback options for various edge cases

## Integration with Other Modules

### Signature Nomination
- **Input source**: Uses results from signature discovery analysis
- **Complementary analysis**: Provides quantitative ranking for qualitative PD Biology Focus

### Export Module
- **Results integration**: All trends analysis results can be exported
- **Report generation**: Integrates with comprehensive analysis reports

### Visualization Modules
- **Data sharing**: Trends results can inform other visualization choices
- **Cross-validation**: Compare trends with individual signature details

## Future Enhancements

### Planned Features
- **Network analysis**: Signature relationship mapping
- **Temporal analysis**: Tracking signature trends over analysis iterations
- **Cross-dataset comparison**: Comparing trends across different experimental conditions
- **Machine learning integration**: Predictive modeling of signature importance

### User-Requested Features
- **Custom categorization**: User-defined biological categories
- **Advanced filtering**: Complex multi-parameter signature selection
- **Batch analysis**: Processing multiple datasets simultaneously

---

**Version**: 0.2.0  
**Last Updated**: January 2025  
**Module**: `mod_signature_trends.R`  
**Core Functions**: `analyze_signature_trends()`, `compute_signature_frequency_analysis()`, `compute_signature_impact_analysis()`