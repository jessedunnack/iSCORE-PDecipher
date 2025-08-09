# iSCORE-PDecipher Comprehensive Data Structure Documentation

## Overview
The iSCORE-PDecipher project contains a complex multi-dimensional data structure for analyzing Parkinson's Disease single-cell RNA-seq data with both mutations (MAST) and CRISPR perturbations (MixScale).

## Data Dimensions

### 1. Resolution Levels
- **COARSE**: 16 clusters (resolution 0.2)
  - Supports all MixScale approaches (Pooled, Split, Combined)
  - Complete FDR-corrected enrichment available
- **FINE**: 40 clusters (higher resolution)
  - Currently supports Pooled approach only
  - FDR-corrected enrichment incomplete (~1200 jobs needed)

### 2. FDR Status
- **WITH FDR**: Benjamini-Hochberg corrected p-values
  - Located in `enrichment_FDR_corrected/` directory
  - MAST data already corrected (p_val_adj column)
  - MixScale requires manual correction
- **WITHOUT FDR**: Raw p-values
  - Used for initial exploration
  - Higher sensitivity but more false positives

### 3. MixScale Approaches
- **Pooled**: All experiments combined before analysis
  - Available for both COARSE and FINE resolutions
  - Most comprehensive approach
- **Split**: Individual experiments analyzed separately
  - C12_FPD-23, C12_FPD-24, C18_FPD-23
  - COARSE resolution only
- **Combined**: UNION of Split results with consistent direction
  - Requires genes significant in same direction across experiments
  - COARSE resolution only

### 4. Direction
- **ALL**: All differentially expressed genes
- **UP**: Only upregulated genes
- **DOWN**: Only downregulated genes

### 5. Genes/Mutations
- **MAST**: 13 mutations
  - SNCA_A30P, SNCA_A53T
  - LRRK2_G2019S, LRRK2_R1441C
  - VPS13C_A444P, VPS13C_W395C
  - VPS35_D620N
  - GBA_IVS2+1
  - PARK7_M26I
  - ATP13A2_compound
  - RAB29_enlarged
  - PINK1_Q456X
  - PRKN_various
- **MixScale**: 10 genes
  - PARK7, ATP13A2, SNCA, LRRK2
  - DNAJC6, FBXO7, PINK1, SYNJ1
  - VPS13C, PARK2 (same as PRKN)

## File Locations and Formats

### Primary Data Files
```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PD_plus_CRISPRi/
├── iSCORE-PD_plus_CRISPRi_final.rds          # Seurat object (~200K cells)
├── full_DE_results.rds                       # Original DE results
├── full_DE_results_with_pooled.rds           # DE results with pooled approach
├── all_enrichment_padj005_complete_with_direction.rds  # 627,959 enrichment terms
└── pooled_mixscale_results.rds               # Pooled MixScale-specific results
```

### DE Results Structure

#### MAST Format
```r
$iSCORE_PD_MAST$<mutation>$<cluster>$results
# Columns:
# - avg_log2FC: Average log2 fold change
# - p_val: Raw p-value
# - p_val_adj: FDR-corrected p-value
# - pct.1: Percent expression in mutation
# - pct.2: Percent expression in control
# - gene: Gene names (in rownames)
```

#### MixScale Format
```r
$CRISPRi_Mixscale$<gene>$<cluster>$results
# Columns:
# - gene_ID: Gene identifier
# - log2FC_<experiment>: Log2 fold change per experiment
# - p_cell_type<experiment>:weight: Weighted p-value per experiment
# Experiments: C12_FPD-24, C12_FPD-23, C18_FPD-23
```

### Enrichment Data Organization

#### FDR-Corrected Enrichment Hierarchy
```
enrichment_FDR_corrected/
├── COARSE/
│   ├── Pooled/
│   │   └── <Gene>/
│   │       └── cluster_<n>/
│   │           └── <Direction>/
│   │               ├── GO_BP_ALL.rds
│   │               ├── GO_CC_ALL.rds
│   │               ├── GO_MF_ALL.rds
│   │               ├── KEGG_ALL.rds
│   │               ├── Reactome_ALL.rds
│   │               ├── WikiPathways_ALL.rds
│   │               ├── STRING_ALL.rds
│   │               └── job_metadata.rds
│   ├── Split/
│   │   └── <Experiment>/
│   │       └── <Gene>/...
│   └── Combined/
│       └── <Gene>/...
└── FINE/
    └── Pooled/
        └── <Gene>/...
```

#### Gene Lists
```
gene_lists/
├── COARSE_Pooled_<Gene>_cluster_<n>_<Direction>.txt
├── COARSE_Split_<Experiment>_<Gene>_cluster_<n>_<Direction>.txt
└── FINE_Pooled_<Gene>_cluster_<n>_<Direction>.txt
```

### Convergence Analysis Results

#### Gene Convergence
```
coarse_cluster_convergence_results.rds
# Columns:
# - mast_gene: MAST mutation name
# - mixscale_gene: MixScale gene name
# - cluster: Cluster identifier
# - approach: MixScale approach used
# - n_mast: Number of MAST DE genes
# - n_mixscale: Number of MixScale DE genes
# - n_overlap: Number of overlapping genes
# - jaccard: Jaccard similarity index
```

#### Pathway Convergence
```
coarse_pathway_convergence_results.rds
fine_cluster_pathway_convergence_results.rds
# Similar structure but for pathway-level comparisons
# Shows much higher convergence (37.2% COARSE, 55.6% FINE)
```

## Data Import Recommendations

### 1. Create Unified Import Function
```r
import_all_data <- function(base_dir, resolution = "COARSE", 
                           fdr_status = "WITH", 
                           approach = "Pooled") {
  data <- list(
    de_results = import_de_results(base_dir, resolution),
    enrichment = import_enrichment(base_dir, resolution, fdr_status, approach),
    convergence = import_convergence(base_dir, resolution, fdr_status),
    metadata = import_metadata(base_dir)
  )
  return(data)
}
```

### 2. Handle Missing Data Gracefully
```r
safe_import <- function(file_path, default = NULL) {
  if (!file.exists(file_path)) {
    message(sprintf("File not found: %s", file_path))
    return(default)
  }
  tryCatch(
    readRDS(file_path),
    error = function(e) {
      warning(sprintf("Error reading %s: %s", file_path, e$message))
      return(default)
    }
  )
}
```

### 3. Standardize Data Access Patterns
```r
get_enrichment_path <- function(base_dir, resolution, approach, 
                               gene, cluster, direction, database) {
  file.path(base_dir, "enrichment_FDR_corrected", 
           resolution, approach, gene, 
           paste0("cluster_", cluster), 
           direction, 
           paste0(database, "_ALL.rds"))
}
```

### 4. Create Data Validation Functions
```r
validate_data_dimensions <- function(data) {
  required_dims <- list(
    resolution = c("COARSE", "FINE"),
    fdr_status = c("WITH", "WITHOUT"),
    approach = c("Pooled", "Split", "Combined"),
    direction = c("ALL", "UP", "DOWN")
  )
  
  # Check each dimension
  for (dim in names(required_dims)) {
    if (!all(data[[dim]] %in% required_dims[[dim]])) {
      warning(sprintf("Invalid %s values found", dim))
    }
  }
}
```

### 5. Implement Lazy Loading
```r
create_data_accessor <- function(base_dir) {
  # Return a function that loads data on demand
  function(data_type, ...) {
    switch(data_type,
      "de" = import_de_results(base_dir, ...),
      "enrichment" = import_enrichment(base_dir, ...),
      "convergence" = import_convergence(base_dir, ...),
      stop("Unknown data type")
    )
  }
}
```

## Data Quality Considerations

### Known Issues
1. **FINE resolution gaps**: FDR-corrected enrichment incomplete
2. **Approach limitations**: Split/Combined only for COARSE
3. **CRISPRa coverage**: Limited to A15_FPD-24 experiment
4. **Timeout failures**: Some enrichment analyses exceed 600s limit
5. **Missing combinations**: Not all gene-cluster pairs have results

### Validation Checklist
- [ ] Check file existence before loading
- [ ] Validate data dimensions match expected values
- [ ] Handle empty/NULL results gracefully
- [ ] Verify background genes are present
- [ ] Check for consistent gene naming (PRKN vs PARK2)
- [ ] Validate cluster numbering consistency
- [ ] Ensure FDR correction status is clear

## Usage Examples

### Load COARSE Pooled FDR-corrected data
```r
data <- import_all_data(
  base_dir = "/mnt/e/ASAP/scRNASeq/PerturbSeq/final",
  resolution = "COARSE",
  fdr_status = "WITH",
  approach = "Pooled"
)
```

### Access specific enrichment result
```r
enrichment <- get_enrichment_result(
  gene = "LRRK2",
  cluster = 6,
  direction = "UP",
  database = "GO_BP",
  resolution = "COARSE",
  approach = "Pooled"
)
```

### Compare convergence across approaches
```r
convergence_comparison <- compare_convergence(
  mast_gene = "SNCA_A53T",
  mixscale_gene = "SNCA",
  cluster = "cluster_6",
  approaches = c("Pooled", "Split_C12_FPD-24", "Combined")
)
```

## Future Improvements

1. **Complete FINE FDR enrichment**: Run remaining ~1200 jobs
2. **Add Split/Combined for FINE**: Extend approaches to fine resolution
3. **Integrate CRISPRa fully**: Standardize A15_FPD-24 handling
4. **Create data catalog**: Auto-generate available data inventory
5. **Implement caching**: Speed up repeated data access
6. **Add data versioning**: Track analysis iterations
7. **Create data validation suite**: Comprehensive QC checks

## Contact
For questions about data structure or access patterns, refer to the iSCORE-PDecipher documentation or check the agent coordinator system.