# Complete File Locations Reference for iSCORE-PDecipher Update

**Purpose**: Comprehensive listing of all files the agent needs to access for integrating FDR-corrected MixScale data

**Last Updated**: October 23, 2025

---

## 📂 Directory Structure Overview

```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/
├── iSCORE-PDecipher/                          # Agent working directory
│   ├── CLAUDE.md                               # Main instructions (THIS IS KEY!)
│   ├── FILE_LOCATIONS.md                       # This file
│   ├── R/                                      # Package R functions to update
│   └── inst/shiny/modules/                     # Shiny modules to update
│
├── final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/
│   ├── CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/     # FPD data
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds  # FDR-corrected DE results
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster1/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_1_mixscale_DEGs.rds
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster2/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_2_mixscale_DEGs.rds
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster3/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_3_mixscale_DEGs.rds
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster4/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_4_mixscale_DEGs.rds
│   │   ├── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster5/
│   │   │   └── all_FPD_no_multiplets_noExptSplit_clust_5_mixscale_DEGs.rds
│   │   └── all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster6/
│   │       └── all_FPD_no_multiplets_noExptSplit_clust_6_mixscale_DEGs.rds
│   │
│   └── CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/  # CRISPRi data
│       ├── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster0/
│       │   └── all_CRISPRi_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds  # FDR-corrected DE results
│       ├── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster1/
│       │   └── all_CRISPRi_no_multiplets_noExptSplit_clust_1_mixscale_DEGs.rds
│       ├── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster2/
│       │   └── all_CRISPRi_no_multiplets_noExptSplit_clust_2_mixscale_DEGs.rds
│       ├── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster3/
│       │   └── all_CRISPRi_no_multiplets_noExptSplit_clust_3_mixscale_DEGs.rds
│       ├── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster4/
│       │   └── all_CRISPRi_no_multiplets_noExptSplit_clust_4_mixscale_DEGs.rds
│       └── all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster5/
│           └── all_CRISPRi_no_multiplets_noExptSplit_clust_5_mixscale_DEGs.rds
│
├── enrichment_analysis_simple.R               # Core enrichment library (reference)
├── run_enrichment_p_weight.R                  # Wrapper script (reference)
├── run_enrichment_p_weight_BH.R               # Wrapper script (reference)
├── run_enrichment_p_weight_bonferroni.R       # Wrapper script (reference)
├── ENRICHMENT_SIMPLE_README.md                # Enrichment pipeline documentation
├── HPC_AGENT_PROMPT.md                        # HPC deployment guide
│
└── [FUTURE] enrichment_results_*/             # Will exist after HPC jobs complete
    ├── enrichment_results_FPD_p_weight/
    ├── enrichment_results_FPD_p_weight_BH/
    ├── enrichment_results_FPD_p_weight_bonferroni/
    ├── enrichment_results_CRISPRi_p_weight/
    ├── enrichment_results_CRISPRi_p_weight_BH/
    └── enrichment_results_CRISPRi_p_weight_bonferroni/
```

---

## 🎯 Key File Groups

### 1. Agent Working Directory
**Location**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/`

**Key Files**:
- `CLAUDE.md` - **MAIN INSTRUCTIONS** (read this first!)
- `FILE_LOCATIONS.md` - This reference document
- `R/data_import_functions.R` - **UPDATE REQUIRED** (add pooled data import)
- `R/data_import_functions_optimized.R` - Reference for optimization patterns
- `R/import_functions_enhanced.R` - Additional import logic
- `inst/shiny/modules/mod_data_loader.R` - **UPDATE REQUIRED** (add p-value selector)
- `inst/shiny/modules/mod_de_results.R` - **UPDATE REQUIRED** (handle pooled structure)
- `inst/shiny/modules/mod_enrichment_gene_display.R` - **UPDATE REQUIRED** (comparison view)
- `inst/shiny/modules/mod_heatmap.R` - **UPDATE REQUIRED** (pooled columns)

### 2. FPD MixScale DE Results (7 files with FDR corrections)

**Base Directory**: `../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/`

**Files** (relative to base):
1. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0/all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds`
2. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster1/all_FPD_no_multiplets_noExptSplit_clust_1_mixscale_DEGs.rds`
3. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster2/all_FPD_no_multiplets_noExptSplit_clust_2_mixscale_DEGs.rds`
4. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster3/all_FPD_no_multiplets_noExptSplit_clust_3_mixscale_DEGs.rds`
5. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster4/all_FPD_no_multiplets_noExptSplit_clust_4_mixscale_DEGs.rds`
6. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster5/all_FPD_no_multiplets_noExptSplit_clust_5_mixscale_DEGs.rds`
7. `all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster6/all_FPD_no_multiplets_noExptSplit_clust_6_mixscale_DEGs.rds`

**Absolute Paths** (for testing):
```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0/all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds
(etc. for Cluster1-6)
```

### 3. CRISPRi MixScale DE Results (6 files with FDR corrections)

**Base Directory**: `../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/`

**Files** (relative to base):
1. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster0/all_CRISPRi_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds`
2. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster1/all_CRISPRi_no_multiplets_noExptSplit_clust_1_mixscale_DEGs.rds`
3. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster2/all_CRISPRi_no_multiplets_noExptSplit_clust_2_mixscale_DEGs.rds`
4. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster3/all_CRISPRi_no_multiplets_noExptSplit_clust_3_mixscale_DEGs.rds`
5. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster4/all_CRISPRi_no_multiplets_noExptSplit_clust_4_mixscale_DEGs.rds`
6. `all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster5/all_CRISPRi_no_multiplets_noExptSplit_clust_5_mixscale_DEGs.rds`

**Absolute Paths** (for testing):
```
/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster0/all_CRISPRi_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds
(etc. for Cluster1-5)
```

### 4. Enrichment Pipeline Scripts (Reference Only)

**Location**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/`

**Files**:
- `enrichment_analysis_simple.R` - Core enrichment library for pooled data
- `run_enrichment_p_weight.R` - Wrapper for uncorrected p-values
- `run_enrichment_p_weight_BH.R` - Wrapper for BH-corrected p-values
- `run_enrichment_p_weight_bonferroni.R` - Wrapper for Bonferroni-corrected p-values
- `submit_enrichment_p_weight.sh` - SLURM script (uncorrected)
- `submit_enrichment_p_weight_BH.sh` - SLURM script (BH)
- `submit_enrichment_p_weight_bonferroni.sh` - SLURM script (Bonferroni)
- `submit_all_enrichment_parallel.sh` - Master launcher

**Purpose**: Reference for understanding how enrichment data was generated

### 5. Documentation Files

**Primary Documentation**:
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/CLAUDE.md` - **MAIN AGENT INSTRUCTIONS**
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/FILE_LOCATIONS.md` - This file

**Supporting Documentation**:
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/ENRICHMENT_SIMPLE_README.md` - Enrichment pipeline guide
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/HPC_AGENT_PROMPT.md` - HPC deployment instructions
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/README.md` - Original package README
- `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/NEWS.md` - Package change log

### 6. Future Enrichment Results (After HPC Jobs Complete)

**Base Directory**: `/mnt/e/ASAP/scRNASeq/PerturbSeq/final/`

**Directories** (will be created by HPC jobs):
1. `enrichment_results_FPD_p_weight/` - FPD with original p-values
2. `enrichment_results_FPD_p_weight_BH/` - FPD with BH correction
3. `enrichment_results_FPD_p_weight_bonferroni/` - FPD with Bonferroni correction
4. `enrichment_results_CRISPRi_p_weight/` - CRISPRi with original p-values
5. `enrichment_results_CRISPRi_p_weight_BH/` - CRISPRi with BH correction
6. `enrichment_results_CRISPRi_p_weight_bonferroni/` - CRISPRi with Bonferroni correction

**Structure** (example):
```
enrichment_results_FPD_p_weight_BH/
├── cluster_0/
│   ├── perturbation1_enrichment.rds
│   ├── perturbation2_enrichment.rds
│   └── ...
├── cluster_1/
└── ...
```

---

## 🔍 File Content Verification

### Verify FDR Columns Exist

**Test on FPD Cluster0**:
```r
# Load data
de_results <- readRDS("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/all_CRISPRi_all_FPD_no_multiplets_noExptSplit_Cluster0/all_FPD_no_multiplets_noExptSplit_clust_0_mixscale_DEGs.rds")

# Check structure
names(de_results)  # Should show perturbation names

# Check first perturbation
first_pert <- de_results[[1]]
colnames(first_pert)  # Should include: p_weight, p_weight_BH, p_weight_bonferroni

# Verify FDR columns exist
required_cols <- c("gene_ID", "log2FC", "p_weight", "p_weight_BH", "p_weight_bonferroni")
all(required_cols %in% colnames(first_pert))  # Should be TRUE
```

**Test on CRISPRi Cluster3**:
```r
# Load data
de_results <- readRDS("/mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_Cluster3/all_CRISPRi_no_multiplets_noExptSplit_clust_3_mixscale_DEGs.rds")

# Check structure
length(de_results)  # Should be 338 perturbations

# Verify FDR columns
required_cols <- c("gene_ID", "log2FC", "p_weight", "p_weight_BH", "p_weight_bonferroni")
all(required_cols %in% colnames(de_results[[1]]))  # Should be TRUE
```

---

## 📦 Data Summary

### Total Files with FDR Corrections
- **FPD**: 7 RDS files (Cluster0-6)
- **CRISPRi**: 6 RDS files (Cluster0-5)
- **Total**: 13 RDS files

### File Sizes (Approximate)
- **FPD files**: 18-19 MB each
- **CRISPRi files**: 213 MB (Cluster0), varies by cluster
- **Total dataset size**: ~1.5 GB

### Perturbations per Cluster
- **FPD Cluster0**: ~500+ perturbations (exact number varies)
- **CRISPRi Cluster3**: 338 perturbations (verified)

### Required Columns in All Files
All 13 RDS files contain these columns:
- `gene_ID` - Gene identifier
- `log2FC` - Log2 fold change (pooled)
- `p_weight` - Original uncorrected p-value
- `p_weight_BH` - Benjamini-Hochberg FDR corrected
- `p_weight_bonferroni` - Bonferroni corrected
- `beta_Intercept`, `beta_weight`, `beta_log_ct` - Regression coefficients
- `p_Intercept`, `p_log_ct` - Additional p-values
- `DE_method` - "wmvReg" for all

---

## 🗺️ Relative Path Reference

### From iSCORE-PDecipher Directory

**To FPD Data**:
```
../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/all_CRISPRi_all_FPD_no_multiplets_noExptSplit_ClusterX/all_FPD_no_multiplets_noExptSplit_clust_X_mixscale_DEGs.rds
```

**To CRISPRi Data**:
```
../final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/all_CRISPRi_all_CRISPRi_no_multiplets_noExptSplit_ClusterX/all_CRISPRi_no_multiplets_noExptSplit_clust_X_mixscale_DEGs.rds
```

**To Enrichment Pipeline**:
```
../enrichment_analysis_simple.R
../run_enrichment_p_weight*.R
../ENRICHMENT_SIMPLE_README.md
```

**To Future Enrichment Results**:
```
../enrichment_results_FPD_p_weight/
../enrichment_results_FPD_p_weight_BH/
../enrichment_results_FPD_p_weight_bonferroni/
../enrichment_results_CRISPRi_p_weight/
../enrichment_results_CRISPRi_p_weight_BH/
../enrichment_results_CRISPRi_p_weight_bonferroni/
```

---

## ✅ Quick Access Commands

### Navigate to Agent Working Directory
```bash
cd /mnt/e/ASAP/scRNASeq/PerturbSeq/final/iSCORE-PDecipher/
```

### List All FPD Cluster Directories
```bash
ls -d /mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit/all_CRISPRi_all_FPD_*/
```

### List All CRISPRi Cluster Directories
```bash
ls -d /mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit/all_CRISPRi_all_CRISPRi_*/
```

### Find All FDR-Corrected RDS Files
```bash
find /mnt/e/ASAP/scRNASeq/PerturbSeq/final/final_hdWGCNA_results/testing_hdWGCNA/mixscale_for_paper/ -name "*_mixscale_DEGs.rds"
```

---

**End of File Locations Reference**
