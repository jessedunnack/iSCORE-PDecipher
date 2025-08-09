# Cell Type Annotation Criteria and Expansion Opportunities

## Current Confidence Level Criteria

### High Confidence (15 clusters)
- **Well-established marker genes**: Multiple canonical markers that clearly define the cell type
- **Literature support**: Strong evidence from published studies
- **Functional coherence**: Markers belong to related pathways/functions
- **Examples**:
  - Cluster 5 (DA neurons): TH, DDC, SLC6A3, SLC18A2 - complete dopamine synthesis/transport pathway
  - Cluster 4 (Choroid plexus): TTR, FOLR1, HTR2C - classic choroid plexus markers
  - Cluster 22/35 (Proliferating): TOP2A, MKI67, CDK1 - canonical cell cycle markers

### Medium Confidence (7 clusters)
- **Some canonical markers present** but incomplete signature
- **Mixed or transitional phenotypes**
- **Regional identity markers** without clear cell type
- **Examples**:
  - Cluster 3 (Floor plate): EN2+ but missing some classic markers
  - Cluster 31 (Caudal): HOX genes indicate position but not specific cell type
  - Cluster 30 (Hypothalamic): POU3F2+ but limited other markers

### Low Confidence (14 clusters)
- **Novel or uncharacterized markers**
- **Mixed transcriptional signatures**
- **Technical concerns** (stress, artifacts)
- **Insufficient markers** for definitive classification
- **Examples**:
  - Cluster 11: PRSS56+ (function unknown in neural context)
  - Cluster 14: Melanocyte markers (unexpected in neural cultures)
  - Cluster 2: Mixed markers without clear identity

## Marker Gene Information We Can Expand

### 1. Functional Annotation of Novel Markers
```r
# For each "Unknown" cluster, we could:
- PRSS56 (Cluster 11): Serine protease, role in eye development
- IGFBPL1 (Cluster 12): IGF binding protein-like, neural development
- RGCC (Cluster 16): Cell cycle regulator, stress response
- RTL1 (Multiple clusters): Retrotransposon-like, imprinted gene
```

### 2. Pathway Enrichment Analysis
- Run GO/KEGG enrichment for each cluster's top 50-100 markers
- Identify functional themes beyond individual markers
- Could reveal cell state/function even without clear cell type

### 3. Cross-Reference with Atlases
```r
# Potential reference datasets:
- Allen Brain Atlas (developing human brain)
- La Manno et al. 2016 (human midbrain development)
- Tiklová et al. 2019 (mouse dopaminergic subtypes)
- Kamath et al. 2022 (human midbrain scRNA-seq)
```

### 4. Developmental Stage Markers
Many "unknown" clusters might represent specific developmental stages:
- Cluster 0: Immature neurons (DCX+, TUBB3+)
- Cluster 8: Early progenitors (PTN+, PTPRZ1+)
- Cluster 26: Neural stem cells (NES+, VIM+)

### 5. Subtype Resolution for Known Types

#### Dopaminergic Neurons (Current: 3 clusters)
Could expand to identify:
- **A8 (RRF)**: High ALDH1A1, low GIRK2
- **A9 (SNc)**: High GIRK2/KCNJ6, vulnerable in PD
- **A10 (VTA)**: High CALB1, OTX2+
- **Pioneer DA**: High EN1/2, early born

#### Floor Plate Progenitors (Current: 3 clusters)
Could distinguish:
- **Early FP**: FOXA2+/SHH+/CORIN+
- **Late FP**: LMX1A+/MSX1+
- **FP-derived glia**: SOX9+/NKX6.1+

## Reasoning Behind Classifications

### Strong Classifications (High Confidence)
1. **Dopaminergic neurons**: Complete dopamine pathway (TH→DDC→VMAT2/DAT)
2. **Choroid plexus**: TTR (transthyretin) is nearly exclusive to CP
3. **Oligodendrocytes**: Complete myelin program (MBP, PLP1, MOG)
4. **Proliferating**: Universal cell cycle markers

### Uncertain Classifications (Low Confidence)
1. **Mixed markers**: e.g., Cluster 2 has stress + neural markers
2. **Novel combinations**: e.g., Cluster 14 has melanocyte + neural markers
3. **Sparse markers**: Some clusters have <5 strong markers
4. **Position without identity**: HOX+ clusters indicate where, not what

## Proposed Annotation Improvements

### 1. Hierarchical Classification
```
Level 1: Neural vs Non-neural
Level 2: Neurons vs Progenitors vs Glia
Level 3: Neuron subtypes (DA, GABA, etc.)
Level 4: Molecular subtypes (A9-like, A10-like)
```

### 2. Confidence Scoring (0-100)
Instead of Low/Medium/High, use numeric scores based on:
- Number of canonical markers (0-40 points)
- Literature support (0-30 points)
- Functional coherence (0-20 points)
- Absence of contradictory markers (0-10 points)

### 3. Alternative Annotations
For each cluster, provide:
- Primary annotation (most likely)
- Alternative annotations (other possibilities)
- Key evidence for/against each

### 4. Expanded Marker Lists
Currently showing top 3-5 markers, could expand to:
- Top 10 positive markers
- Top 5 negative markers (not expressed)
- Marker specificity scores

### 5. Cell State Annotations
Beyond cell type, annotate:
- Maturation state (progenitor/immature/mature)
- Activation state (resting/activated)
- Stress level (healthy/stressed/dying)
- Cell cycle phase

## Data Available for Expansion

1. **Full marker gene lists**: Currently using top 25, have access to all DE genes
2. **Expression levels**: Could use expression magnitude, not just presence/absence
3. **Co-expression patterns**: Identify gene modules/programs
4. **Negative markers**: Genes specifically NOT expressed
5. **Cluster relationships**: Use UMAP proximity and hierarchical clustering

## Recommended Next Steps

1. **Rerun marker analysis with top 50-100 genes per cluster**
2. **Perform pathway enrichment to identify functional themes**
3. **Score marker specificity** (how unique to each cluster)
4. **Compare with published atlases** for validation
5. **Create detailed annotation reports** for each cluster

Would you like me to implement any of these expansions?