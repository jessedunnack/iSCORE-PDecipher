# Comprehensive Fisher's Exact Test Analysis Summary
**Date:** January 17, 2025
**Analysis:** DE Gene Overlap Between CRISPRi and iSCORE-PD Mutations

## Executive Summary

The most statistically significant overlaps between CRISPRi perturbations and iSCORE-PD mutations reveal strong convergent signatures in Parkinson's disease genes, with proper handling of directionality to avoid double-counting.

## Top 5 Most Significant Overlaps

### 1. **SYNJ1** (Synaptojanin 1)
- **Best p-value:** 6.47 × 10⁻⁴² (cluster_7, C12_FPD-24)
- **Overlap:** 571 genes (280 same direction, 291 opposite)
- **PD Relevance:** 1,541 PD-relevant pathways
  - Mitophagy: 491 pathways
  - Mitochondrial: 481 pathways
  - Synaptic: 367 pathways
- **Known Function:** Synaptic vesicle uncoating, phosphoinositide metabolism
- **Key Finding:** Strongest statistical evidence of convergent dysfunction

### 2. **ATP13A2** (PARK9)
- **Best p-value:** 1.62 × 10⁻³⁵ (cluster_6, A15_FPD-24)
- **Overlap:** 118 genes (96 same direction, 22 opposite)
- **PD Relevance:** 1,540 PD-relevant pathways
  - Mitophagy: 506 pathways
  - Mitochondrial: 510 pathways
- **Known Function:** Lysosomal ATPase, α-synuclein clearance
- **Key Finding:** Strongest mitochondrial/lysosomal convergence

### 3. **FBXO7** (PARK15)
- **Best p-value:** 4.60 × 10⁻³⁵ (cluster_3, C12_FPD-24)
- **Overlap:** 942 genes (417 same direction, 525 opposite)
- **PD Relevance:** 1,085 PD-relevant pathways
- **Known Function:** Mitophagy regulation, PINK1-Parkin pathway
- **Key Finding:** Broadest impact across 15 significant comparisons

### 4. **LRRK2** (PARK8)
- **Best p-value:** 1.37 × 10⁻²⁶ (cluster_7, C12_FPD-23)
- **Overlap:** 148 genes (66 same direction, 82 opposite)
- **PD Relevance:** 981 PD-relevant pathways
- **Known Function:** Kinase activity, autophagy regulation
- **Key Finding:** Most clinically relevant PD mutation

### 5. **PARK7** (DJ-1)
- **Best p-value:** 4.90 × 10⁻²⁶ (cluster_8, C12_FPD-23)
- **Overlap:** 672 genes (255 same direction, 417 opposite)
- **PD Relevance:** 814 PD-relevant pathways
- **Known Function:** Oxidative stress response, mitochondrial protection
- **Key Finding:** Strong oxidative stress signature

## Directionality Analysis

The analysis properly accounts for directionality to avoid double-counting:

1. **Same Direction Effects:** Genes up/down-regulated in both methods
2. **Opposite Direction Effects:** Genes with opposing regulation
3. **Union Approach:** All overlapping genes counted once, regardless of direction

### Key Observations:
- Most genes show mixed directionality (not purely same/opposite)
- PRKN/PARK2 shows strongest directional consistency (48 same vs 8 opposite)
- Synaptic genes (SYNJ1) show balanced directionality

## Convergent PD Pathways

### Mitochondrial Dysfunction Cluster
All top genes show strong mitochondrial signatures:
- **ATP13A2:** 510 mitochondrial pathways
- **SYNJ1:** 481 mitochondrial pathways
- **PRKN:** 272 mitochondrial pathways
- **FBXO7:** 241 mitochondrial pathways

### Protein Quality Control Network
Strong enrichment for protein homeostasis:
- **LRRK2:** 102 protein aggregation pathways
- **ATP13A2:** 96 protein aggregation pathways
- **FBXO7:** 95 protein aggregation pathways
- **PRKN:** 95 protein aggregation pathways

### Synaptic Dysfunction
Major synaptic involvement:
- **SYNJ1:** 367 synaptic pathways (expected)
- **FBXO7:** 360 synaptic pathways
- **LRRK2:** 346 synaptic pathways
- **PRKN:** 316 synaptic pathways

## Statistical Robustness

### Multiple Testing Correction
- Hierarchical FDR correction applied
- Results remain significant after correction
- Multiple experiments analyzed independently

### Avoiding Double-Counting
- Each gene counted once per comparison
- Experiments analyzed separately (C12_FPD-23, C12_FPD-24, etc.)
- Union of overlapping genes used for Fisher's test

## Biological Interpretation

### Convergent Mechanisms
1. **Mitochondrial-Lysosomal Axis:** ATP13A2, FBXO7, PRKN
2. **Synaptic-Mitochondrial Crosstalk:** SYNJ1, LRRK2
3. **Oxidative Stress Response:** PARK7, PINK1

### Clinical Relevance
- **LRRK2:** Most common genetic cause of PD
- **ATP13A2:** Early-onset parkinsonism
- **PRKN:** Second most common recessive PD gene
- **SYNJ1:** Early-onset atypical parkinsonism

## Recommendations for Manuscript

### Strongest Candidates:
1. **SYNJ1** - Unprecedented statistical significance (p < 10⁻⁴⁰)
2. **ATP13A2** - Perfect mitochondrial/lysosomal convergence
3. **FBXO7** - Broadest cross-cluster impact
4. **LRRK2** - Clinical relevance + strong statistics

### Key Messages:
1. Orthogonal methods (mutation vs knockdown) converge on same pathways
2. Mitochondrial dysfunction is central across all top genes
3. Direction-aware analysis strengthens biological interpretation
4. Statistical rigor with proper FDR correction validates findings

## Data Files Generated
- `fisher_test_results_detailed.rds` - Complete results object
- `fisher_test_results_flat.csv` - All comparisons table
- `fisher_test_results_summary.csv` - Gene pair summary
- `pd_biological_relevance.csv` - PD pathway enrichment counts