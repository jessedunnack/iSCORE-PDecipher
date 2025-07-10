# Signature Metrics Explained

## Overview
The Signature Nomination system uses several key metrics to identify and rank cross-method signatures (MAST mutations vs CRISPRi knockdowns). These metrics help researchers identify the most biologically relevant and statistically robust signatures.

## Core Metrics

### 1. **Signature Strength** (Score: 0-10+)
**What it measures:** The overall strength of agreement between MAST and CRISPRi methods for a given gene/cluster combination.

**How it's calculated:**
```
Signature Strength = (Gene Overlap Count × 0.4) + 
                    (Pathway Overlap Count × 0.3) + 
                    (-log10(Fisher p-value) × 0.3)
```

**Components:**
- **Gene Overlap Count (40% weight):** Number of differentially expressed genes shared between methods
- **Pathway Overlap Count (30% weight):** Number of enriched pathways/terms shared between methods  
- **Statistical Significance (30% weight):** Fisher's exact test p-value for overlap significance

**Interpretation:**
- **> 2.0**: Strong signature - High confidence in cross-method agreement
- **1.0 - 2.0**: Moderate signature - Good agreement, worth investigating
- **< 1.0**: Weak signature - Limited cross-method support

### 2. **Frequency Score** (0-1)
**What it measures:** How often a signature pattern appears across different clusters or experimental conditions.

**How it's calculated:**
```
Frequency Score = (Number of clusters where signature appears) / (Total clusters analyzed)
```

**Interpretation:**
- **> 0.8**: Very Common - Signature appears in most clusters (pan-cluster effect)
- **0.5 - 0.8**: Common - Signature in multiple clusters  
- **< 0.5**: Cluster-specific - Limited to few clusters

### 3. **Impact Score** (0-1)
**What it measures:** The relative statistical strength compared to the strongest signature in the dataset.

**How it's calculated:**
```
Impact Score = (Current signature strength) / (Maximum signature strength in dataset)
```

**Interpretation:**
- **> 0.9**: Highest Impact - Among the strongest signatures
- **0.7 - 0.9**: High Impact - Statistically robust
- **< 0.5**: Moderate Impact - Significant but not top-tier

### 4. **PD Relevance Score** (0-1)
**What it measures:** Enrichment for Parkinson's disease-related pathways and biological processes.

**How it's calculated:**
The PD relevance score considers enrichment in key PD-associated categories:
- Mitochondrial dysfunction pathways
- Protein aggregation/degradation (proteasome, autophagy)
- Synaptic/neurotransmitter pathways
- Lysosomal function
- Oxidative stress response

**Scoring:**
- Each PD-relevant pathway adds to the score
- Weighted by statistical significance
- Normalized to 0-1 scale

**Interpretation:**
- **> 0.7**: Highly PD-relevant - Strong enrichment for known PD pathways
- **0.4 - 0.7**: Moderately relevant - Some PD pathway involvement
- **< 0.4**: Low PD relevance - May represent general cellular processes

## Additional Metrics

### Gene Overlap Metrics
- **Gene Overlap Count**: Absolute number of shared DE genes
- **Gene Jaccard Index**: Overlap size relative to union size (0-1)
- **Gene Fisher p-value**: Statistical significance of overlap

### Pathway Overlap Metrics  
- **Pathway Overlap Count**: Number of shared enriched pathways
- **Pathway Categories**: Distribution across GO, KEGG, Reactome, etc.
- **Direction Concordance**: Agreement in up/down regulation

## Using These Metrics

### For Manuscript Priorities:
1. **Focus on High Signature Strength** (>2.0) for main findings
2. **Check Frequency Score** to distinguish pan-cluster vs cluster-specific effects
3. **Use PD Relevance** to prioritize disease-relevant signatures
4. **Validate with Impact Score** for statistical robustness

### For Biological Interpretation:
- **High strength + High frequency** = Core cellular response
- **High strength + Low frequency** = Context-specific response
- **High PD relevance** = Disease-relevant mechanism

### Quality Control:
- Signatures with **very low p-values** (<1e-10) but **low overlap counts** (<5) may be statistical artifacts
- Always check the actual gene lists and pathways, not just scores
- Consider biological plausibility alongside statistical metrics