# iSCORE-PDecipher Refactoring Plan
## Strategic Pivot: MAST Removal & Perturb-seq Focus

**Date:** 2025-11-18
**Goal:** Transform package from dual-method (MAST + MixScale) to **Perturb-seq exploration tool**
**User Workflow:** Download package → Download published dataset → Explore interactively
**Reversibility:** All MAST code archived on separate branch, easily recoverable

---

## 📊 **Current State Audit (PHASE 1)**

### **Codebase Statistics:**
- **Total R files:** 75
- **MAST/mutation references:** 1,315 instances across 49 files
- **Module sizes:**
  - `mod_signature_nomination.R`: 3,213 lines (too large)
  - `mod_de_results.R`: 3,088 lines (too large)
  - `mod_heatmap.R`: 2,154 lines (too large)

### **Files Categorized by MAST Dependency:**

#### **Category A: MAST-ONLY (Archive Entirely)**
Files that exclusively handle MAST analysis and can be archived:

```
R/MAST_analysis.R                      (4.9KB)
R/MAST_analysis_optimized.R            (15KB)
```

#### **Category B: CROSS-METHOD (Heavy MAST Integration)**
Files that compare MAST vs MixScale - need strategic refactoring:

```
R/signature_analysis.R                 (Cross-method signature discovery)
R/manuscript_signature_discovery.R     (MAST vs CRISPRi comparison)
R/gene_harmonization.R                 (Maps MAST gene names to CRISPRi names)
R/comprehensive_fishers_analysis.R     (Fisher's test for gene overlap)
R/enhanced_direction_analysis.R        (Direction concordance MAST vs CRISPRi)
R/pd_signature_interpretation.R        (Interprets MAST+CRISPRi convergence)

inst/shiny/modules/mod_signature_nomination.R  (3,213 lines - MAST vs CRISPRi UI)
inst/shiny/modules/mod_comparison.R            (MAST vs MixScale comparison UI)
```

**Decision:** These need to be **completely reconceptualized** for Perturb-seq only:
- Signature discovery → Perturbation effect signatures (cross-cluster consistency)
- Gene harmonization → Perturbation name standardization
- Direction analysis → Effect direction across clusters

#### **Category C: SHARED (Selective MAST Removal)**
Files that handle both MAST and MixScale with separable code paths:

```
R/data_import_functions.R              (import_de_results handles both)
R/data_import_functions_optimized.R    (detect_mast_format functions)
R/dataset_validator.R                  (validates MAST and MixScale directories)
R/enrichment_analysis.R                (run_enrichment for both methods)
R/launch_app.R                         (dataset selection includes MAST options)

inst/shiny/modules/mod_de_results.R     (Volcano plots for MAST and MixScale)
inst/shiny/modules/mod_heatmap.R        (Heatmaps for both methods)
inst/shiny/modules/mod_landing_page_with_umap_v2.R  (Overview with mutation metadata)
inst/shiny/R/data_manager.R            (Manages MAST and MixScale data)
```

**Decision:** These need **surgical MAST removal** with code paths commented/archived

#### **Category D: MixScale/Perturb-seq ONLY (Keep As-Is)**
Files already focused on Perturb-seq:

```
R/import_pooled_mixscale_functions.R   (Pooled MixScale import - KEEP)
inst/shiny/modules/mod_perturbseq_only.R  (Perturb-seq module - KEEP & EXPAND)
```

---

## 🎯 **STRATEGIC REMOVAL PLAN (PHASE 2-3)**

### **Step 1: Create Archive Branch**
```bash
git checkout -b archive/mast-legacy
git push -u origin archive/mast-legacy

# Tag current state
git tag -a v0.4.0-mast-legacy -m "Last version with full MAST support"
git push origin v0.4.0-mast-legacy
```

### **Step 2: Create Feature Flag System (Temporary)**
Add to `R/config_manager.R`:
```r
#' Check if MAST functionality is enabled
#' @keywords internal
is_mast_enabled <- function() {
  getOption("iscore.enable_mast", default = FALSE)
}
```

This allows code to remain but be disabled:
```r
if (is_mast_enabled()) {
  # MAST code paths
} else {
  # Skip MAST analysis
}
```

### **Step 3: Move MAST-ONLY Files to Archive**
```bash
mkdir -p archive/mast_functions/
mv R/MAST_analysis.R archive/mast_functions/
mv R/MAST_analysis_optimized.R archive/mast_functions/
```

### **Step 4: Refactor Cross-Method Files**

#### **4A: Signature Discovery → Perturbation Signatures**
**OLD CONCEPT:** Find convergent pathways between MAST mutations and CRISPRi knockdowns

**NEW CONCEPT:** Find perturbation effect signatures across clusters
- **Pan-cluster signatures:** Perturbations with consistent effects across all/most clusters
- **Cluster-specific signatures:** Perturbations with strong effects in specific cell types
- **Pathway convergence:** Multiple perturbations converging on same pathway

**Refactored Functions:**
```r
# REMOVE:
discover_cross_method_signatures()
compare_mast_vs_crispri()
calculate_gene_overlap_fisher()

# REPLACE WITH:
discover_perturbation_signatures(
  mixscale_results,
  min_clusters = 3,           # Minimum clusters showing effect
  consistency_threshold = 0.7  # Direction consistency across clusters
)

identify_pathway_convergence(
  perturbations,               # Multiple genes perturbed
  enrichment_data,
  pathway_overlap_threshold = 0.3
)
```

#### **4B: Gene Harmonization → Perturbation Standardization**
**OLD:** Map `SNCA_A30P` (MAST) ↔ `SNCA` (CRISPRi)

**NEW:** Standardize perturbation names across experiments
```r
# NEW FUNCTION:
standardize_perturbation_names(
  perturbation_list,
  naming_scheme = c("HGNC", "Ensembl", "custom")
)

# Handle cases like:
# - PARK2 vs PRKN (same gene, different names)
# - Non-Targeting vs NT-Control vs Scramble (all controls)
```

#### **4C: Data Import - Remove MAST Paths**
In `R/data_import_functions.R`:
```r
# COMMENT OUT:
# import_mast_results()
# detect_mast_format()
# consolidate_mast_mixscale()

# KEEP & ENHANCE:
import_pooled_mixscale_data()
import_enrichment_results()
```

---

## 🔨 **SIMPLIFICATION PLAN (PHASE 4)**

### **Principle 1: Follow Google's R Style Guide**
Reference: https://google.github.io/styleguide/Rguide.html

**Key Rules:**
1. **File names:** `snake_case.R` (already following ✓)
2. **Function names:** `verb_noun()` format
3. **Variable names:** `snake_case` not `camelCase`
4. **Line length:** Max 80 characters
5. **Documentation:** Always use roxygen2
6. **No global state:** Pass parameters explicitly

### **Principle 2: DRY (Don't Repeat Yourself)**
**Problem:** Multiple functions doing same thing with slight variations

**Example - CURRENT:**
```r
# Found in codebase:
import_de_results()
import_mast_results()
import_mixscale_results()
import_pooled_mixscale_data()
import_enrichment_results()
```

**REFACTORED:**
```r
#' Universal Data Import Function
#' @param source Path to data file/directory
#' @param type One of: "de_results", "enrichment", "seurat"
#' @param format Auto-detected or specify: "pooled", "experiment_split"
import_data <- function(source, type, format = "auto") {
  switch(type,
    "de_results" = import_de_results_internal(source, format),
    "enrichment" = import_enrichment_internal(source),
    "seurat" = import_seurat_internal(source),
    stop("Unknown type: ", type)
  )
}
```

### **Principle 3: Single Responsibility**
**Problem:** 3,000-line modules doing everything

**Solution:** Break into focused sub-modules

**CURRENT:**
```
inst/shiny/modules/mod_signature_nomination.R  (3,213 lines)
  ├─ UI definition                              (500 lines)
  ├─ Server logic                               (800 lines)
  ├─ Data processing functions                  (600 lines)
  ├─ Fisher's exact test calculations           (400 lines)
  ├─ Visualization helpers                      (500 lines)
  └─ PD pathway interpretation                  (413 lines)
```

**REFACTORED:**
```
inst/shiny/modules/signature/
  ├─ mod_sig_ui.R                    (UI only, 200 lines)
  ├─ mod_sig_server.R                (Server logic, 300 lines)
  ├─ sig_processing.R                (Data processing, 400 lines)
  └─ sig_visualization.R             (Plotting, 300 lines)

R/
  ├─ perturbation_signatures.R       (Core algorithm, 400 lines)
  └─ pathway_convergence.R           (Pathway analysis, 300 lines)
```

### **Principle 4: Tidyverse Best Practices**
Reference: https://style.tidyverse.org/

**Use:**
- `%>%` pipe for readability
- `dplyr` verbs: `filter()`, `mutate()`, `summarize()`
- `purrr::map()` instead of `lapply()`
- `rlang` for tidy evaluation

**Example Refactor:**
```r
# OLD (base R):
results <- list()
for (i in seq_along(clusters)) {
  cluster_data <- data[data$cluster == clusters[i], ]
  enriched <- cluster_data[cluster_data$p_val_adj < 0.05, ]
  results[[clusters[i]]] <- enriched
}

# NEW (tidyverse):
results <- data %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  group_split() %>%
  set_names(unique(data$cluster))
```

---

## 📦 **PACKAGE BEST PRACTICES (PHASE 5)**

### **Directory Structure (Simplified)**

**CURRENT (Messy):**
```
iSCORE-PDecipher/
├── R/                           (35 files, many duplicates)
├── inst/shiny/
│   ├── R/                       (9 helper files)
│   ├── modules/                 (23 modules, some 3000+ lines)
│   ├── www/                     (CSS, images)
│   └── app.R
├── man/                         (Auto-generated)
├── vignettes/                   (4 vignettes)
└── DESCRIPTION
```

**REFACTORED (Clean):**
```
iSCORE-PDecipher/
├── R/                           (Core functions, organized)
│   ├── import.R                 (All import functions)
│   ├── enrichment.R             (Enrichment analysis)
│   ├── perturbation_analysis.R  (Perturb-seq specific)
│   ├── visualization.R          (Core plotting functions)
│   ├── utils.R                  (Helper functions)
│   └── zzz.R                    (.onLoad, .onAttach hooks)
├── inst/
│   ├── shiny/
│   │   ├── app.R                (Main app entry)
│   │   ├── ui/                  (UI components, modular)
│   │   ├── server/              (Server logic, modular)
│   │   └── www/                 (Assets: CSS, JS, images)
│   └── extdata/                 (Example datasets)
├── data/                        (R data objects for package)
│   └── example_perturbseq.rda   (Small example for vignettes)
├── man/                         (Documentation)
├── vignettes/                   (2-3 focused vignettes)
├── tests/testthat/              (Unit tests)
├── _pkgdown.yml                 (Website config)
├── DESCRIPTION                  (Clean dependencies)
├── NAMESPACE                    (Auto-generated by roxygen2)
└── README.md                    (Clear getting started)
```

### **DESCRIPTION File (Cleaned)**

**REMOVE:**
- Dependencies only used for MAST
- Redundant packages (pheatmap if using ComplexHeatmap)
- Rarely used packages

**KEEP:**
- Core: dplyr, ggplot2, shiny
- Enrichment: clusterProfiler, org.Hs.eg.db
- Visualization: plotly, ComplexHeatmap, heatmaply
- Import: jsonlite, readr

**Example DESCRIPTION:**
```
Package: iSCORE.PDecipher
Title: Interactive Exploration of CRISPRi Perturb-seq Results
Version: 1.0.0
Authors@R: person("Jesse", "Dunnack", email = "your@email.com", role = c("aut", "cre"))
Description: Comprehensive analysis and visualization of CRISPRi Perturb-seq experiments
    in the context of Parkinson's disease. Import MixScale differential expression results,
    perform pathway enrichment analysis, discover perturbation effect signatures, and
    explore results through an interactive Shiny application.
License: MIT + file LICENSE
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
RoxygenNote: 7.2.3
Depends: R (>= 4.1.0)
Imports:
    shiny (>= 1.7.0),
    dplyr (>= 1.1.0),
    ggplot2 (>= 3.4.0),
    plotly (>= 4.10.0),
    clusterProfiler (>= 4.0.0),
    org.Hs.eg.db,
    ComplexHeatmap,
    DT,
    shinydashboard,
    bslib (>= 0.5.0)
Suggests:
    testthat (>= 3.0.0),
    knitr,
    rmarkdown
VignetteBuilder: knitr
```

---

## 📚 **EXAMPLE DATASET STRATEGY (PHASE 6)**

### **Two-Tier Dataset Approach:**

#### **Tier 1: Embedded Example (Small)**
**Purpose:** Package vignettes and testing
**Size:** ~1MB (1000 genes, 3 clusters, 10 perturbations)
**Location:** `data/example_perturbseq.rda`

```r
# Create example
example_perturbseq <- list(
  de_results = list(
    LRRK2 = tibble(gene_ID, log2FC, p_weight_BH, ...),
    SNCA = tibble(...),
    ...
  ),
  enrichment = tibble(perturbation, cluster, pathway, p.adjust, ...),
  metadata = list(
    n_cells = 5000,
    n_perturbations = 10,
    clusters = c("Dopaminergic", "Astrocytes", "Microglia")
  )
)

usethis::use_data(example_perturbseq, overwrite = TRUE)
```

**Vignette usage:**
```r
library(iSCORE.PDecipher)
data("example_perturbseq")
launch_app(example_perturbseq)
```

#### **Tier 2: Full Published Dataset (External)**
**Purpose:** User analysis of complete data
**Size:** ~500MB-2GB (full 338 perturbations, 6 clusters, 230k cells)
**Location:** Zenodo/Figshare/GEO
**DOI:** 10.xxxx/xxxxx

**User workflow:**
```r
# 1. Install package
install.packages("remotes")
remotes::install_github("jessedunnack/iSCORE-PDecipher")

# 2. Download full dataset (one-time)
library(iSCORE.PDecipher)
download_published_dataset(
  dest_dir = "~/iSCORE_data/",
  doi = "10.xxxx/xxxxx"
)

# 3. Launch app with full data
launch_app(data_dir = "~/iSCORE_data/perturbseq_full/")
```

---

## 🔄 **SIMPLIFIED USER WORKFLOW (PHASE 7)**

### **Goal: 5-Minute Getting Started**

**README.md Example:**
```markdown
# iSCORE-PDecipher

> Interactive exploration of CRISPRi Perturb-seq results in Parkinson's disease

## Quick Start (5 minutes)

### Install
install.packages("remotes")
remotes::install_github("jessedunnack/iSCORE-PDecipher")

### Explore Example Data
library(iSCORE.PDecipher)
data("example_perturbseq")
launch_app(example_perturbseq)

### Explore Full Published Dataset
# Download once (2GB, ~10 min)
download_published_dataset(doi = "10.xxxx/xxxxx")

# Launch with full data
launch_app(data_dir = "~/iSCORE_data/")

## Features
- **338 CRISPRi perturbations** across 10 PD genes
- **6 cell clusters** from 230k cells
- **Pathway enrichment** (GO, KEGG, Reactome, WikiPathways)
- **Interactive visualizations** (volcano plots, heatmaps, UMAP)
- **Perturbation signatures** (cross-cluster effects)

## Citation
[Your paper citation]
```

---

## ✅ **PHASE IMPLEMENTATION CHECKLIST**

### **Phase 1: Audit ✅ (COMPLETE)**
- [x] Identify all MAST-dependent files (49 files, 1315 references)
- [x] Categorize by dependency type (A/B/C/D)
- [x] Create refactoring plan document

### **Phase 2: Safe Removal Strategy**
- [ ] Create `archive/mast-legacy` branch
- [ ] Tag current version as `v0.4.0-mast-legacy`
- [ ] Add feature flag system (`is_mast_enabled()`)
- [ ] Document reversal procedure

### **Phase 3: Remove MAST Functionality**
- [ ] Archive MAST-only files (Category A: 2 files)
- [ ] Refactor cross-method files (Category B: 8 files)
  - [ ] Redesign signature discovery for Perturb-seq only
  - [ ] Convert gene harmonization to perturbation standardization
  - [ ] Remove Fisher's test for MAST vs MixScale
- [ ] Remove MAST code paths from shared files (Category C: 15+ files)
  - [ ] Update `data_import_functions.R`
  - [ ] Update `dataset_validator.R`
  - [ ] Clean `mod_de_results.R` (remove MAST volcano plots)
  - [ ] Clean `mod_heatmap.R` (remove MAST heatmap options)
  - [ ] Clean `mod_landing_page_with_umap_v2.R` (remove mutation metadata)
- [ ] Delete unused modules
  - [ ] `mod_comparison.R` (MAST vs MixScale comparison)

### **Phase 4: Simplify & Consolidate**
- [ ] Break up monolithic modules (3 files >2000 lines)
  - [ ] Split `mod_signature_nomination.R` (3,213 → 4 files ~300-400 lines each)
  - [ ] Split `mod_de_results.R` (3,088 → 3 files ~400-500 lines each)
  - [ ] Split `mod_heatmap.R` (2,154 → 3 files ~300-400 lines each)
- [ ] Consolidate duplicate functions
  - [ ] Merge `data_import_functions.R` + `data_import_functions_optimized.R`
  - [ ] Remove redundant visualization functions
- [ ] Apply Google R Style Guide
  - [ ] Consistent naming conventions
  - [ ] Line length <80 characters
  - [ ] Proper documentation
- [ ] Apply tidyverse patterns
  - [ ] Replace loops with `purrr::map()`
  - [ ] Use pipes for readability
  - [ ] Consistent dplyr verbs

### **Phase 5: Package Best Practices**
- [ ] Clean DESCRIPTION file
  - [ ] Remove MAST-only dependencies
  - [ ] Remove redundant packages
  - [ ] Specify minimum versions
- [ ] Reorganize R/ directory
  - [ ] Group related functions in single files
  - [ ] Consistent file naming
- [ ] Add example dataset
  - [ ] Create `data/example_perturbseq.rda`
  - [ ] Document with roxygen2
- [ ] Update NAMESPACE
  - [ ] Export only user-facing functions
  - [ ] Mark internal functions with `@keywords internal`
- [ ] Improve README.md
  - [ ] Clear 5-minute quick start
  - [ ] Screenshot of app
  - [ ] Citation information

### **Phase 6: Create Example Dataset & Workflow**
- [ ] Create small embedded dataset (1MB)
  - [ ] 10 perturbations, 3 clusters, 1000 genes
  - [ ] Full enrichment results
  - [ ] Metadata
- [ ] Create `download_published_dataset()` function
  - [ ] Downloads from Zenodo/Figshare
  - [ ] Validates download integrity
  - [ ] Sets up directory structure
- [ ] Test complete user workflow
  - [ ] Install → Load example → Explore
  - [ ] Download full → Load full → Explore

### **Phase 7: Update Documentation**
- [ ] Update vignettes (remove MAST references)
  - [ ] `getting-started.Rmd` → Perturb-seq focus
  - [ ] `perturbseq-workflow.Rmd` → Keep & expand
  - [ ] `signature-discovery.Rmd` → Redesign for Perturb-seq signatures
  - [ ] `advanced-customization.Rmd` → Keep & simplify
- [ ] Update comprehensive documentation
  - [ ] Remove MAST sections
  - [ ] Update workflows
- [ ] Update `_pkgdown.yml`
  - [ ] Remove MAST function references
  - [ ] Reorganize reference sections

### **Phase 8: Testing & Quality**
- [ ] Create unit tests (50+ tests targeting 80% coverage)
  - [ ] Test import functions
  - [ ] Test enrichment analysis
  - [ ] Test perturbation signatures
  - [ ] Test visualization functions
- [ ] Integration tests
  - [ ] Test Shiny modules with testServer()
  - [ ] End-to-end workflow tests
- [ ] Performance benchmarks
  - [ ] Load time for full dataset
  - [ ] Memory usage
  - [ ] Visualization rendering speed
- [ ] User testing
  - [ ] Fresh install on clean system
  - [ ] Complete getting-started workflow
  - [ ] Test all major app features

---

## 🔄 **REVERSIBILITY STRATEGY**

### **If MAST Needs to be Restored:**

1. **Checkout archive branch:**
   ```bash
   git checkout archive/mast-legacy
   ```

2. **Cherry-pick improvements:**
   ```bash
   git cherry-pick <commit-hash>  # Pick specific refactoring improvements
   ```

3. **Merge strategy:**
   ```bash
   git checkout main
   git merge --no-ff archive/mast-legacy
   # Resolve conflicts, keeping both MAST and MixScale
   ```

4. **Feature flag activation:**
   ```r
   options(iscore.enable_mast = TRUE)
   ```

---

## 📊 **SUCCESS METRICS**

### **Code Quality:**
- [ ] All modules <500 lines
- [ ] No duplicate functions (DRY principle)
- [ ] 80%+ test coverage
- [ ] Zero linter warnings (lintr package)

### **Package Quality:**
- [ ] Passes `R CMD check` with zero errors, zero warnings
- [ ] All examples run successfully
- [ ] All vignettes render without errors
- [ ] pkgdown website builds successfully

### **User Experience:**
- [ ] Install to first visualization: <5 minutes
- [ ] README understandable by R beginners
- [ ] Example dataset loads instantly
- [ ] Full dataset loads in <30 seconds

### **Performance:**
- [ ] App launch: <5 seconds
- [ ] Heatmap rendering: <3 seconds
- [ ] Memory usage: <2GB for full dataset
- [ ] No blocking UI operations

---

## 🎯 **TIMELINE ESTIMATE**

| Phase | Duration | Cumulative |
|-------|----------|------------|
| Phase 1: Audit | ✅ 0.5 days | 0.5 days |
| Phase 2: Archive Strategy | 0.5 days | 1 day |
| Phase 3: Remove MAST | 3 days | 4 days |
| Phase 4: Simplify Code | 4 days | 8 days |
| Phase 5: Package Practices | 2 days | 10 days |
| Phase 6: Example Dataset | 1 day | 11 days |
| Phase 7: Update Docs | 2 days | 13 days |
| Phase 8: Testing | 2 days | 15 days |
| **TOTAL** | **~3 weeks** | **15 work days** |

---

**Next Action:** Proceed to Phase 2 (Create archive branch and feature flags)
