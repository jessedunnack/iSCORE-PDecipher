# iSCORE-PDecipher Modernization Recommendations
**Generated:** 2025-11-08
**Version:** 0.4.0 → 1.0.0 Roadmap
**Focus:** Perturb-seq Integration & Modern Shiny Patterns

---

## Executive Summary

This document provides comprehensive modernization recommendations for iSCORE-PDecipher, transitioning from a general-purpose DE visualization tool to a **Perturb-seq-focused analysis platform** with modern Shiny architecture and streamlined CRISPRi integration.

### Key Modernization Goals

1. **Perturb-seq Focus:** Optimize for CRISPRi Perturb-seq data as primary use case
2. **Modern Shiny:** Update to Shiny 1.7+ patterns and bslib theming
3. **Integration:** Direct Mixscale/scMAGeCK output import
4. **Performance:** Leverage caching and async processing
5. **User Experience:** Streamline UI for Perturb-seq workflows

### Priority Summary

| Priority | Category | Items | Timeline |
|----------|----------|-------|----------|
| **P1: Critical** | Deprecated dependencies, breaking changes | 12 | Immediate |
| **P2: High** | Perturb-seq focus, data integration | 18 | 1-2 months |
| **P3: Medium** | UI/UX improvements, modernization | 24 | 2-4 months |
| **P4: Lower** | Nice-to-have enhancements | 15 | 4-6 months |

---

## Section 1: Current State Analysis

### 1.1 Outdated Components Identified

#### A. Deprecated Shiny Patterns

**Issue:** Using older Shiny patterns from pre-1.5 era

**Examples Found:**
```r
# OLD: Manual reactive invalidation
observe({
  input$gene
  updateSelectInput(session, "cluster", ...)
})

# NEW: Use reactive dependencies
observeEvent(input$gene, {
  updateSelectInput(session, "cluster", ...)
})
```

**Files Affected:**
- `inst/shiny/modules/mod_de_results.R` (Lines 450-520)
- `inst/shiny/modules/mod_heatmap.R` (Lines 200-250)
- `inst/shiny/modules/mod_enrichment_gene_display.R` (Lines 300-350)

**Recommendation:** Update to `observeEvent()` and `eventReactive()` patterns

---

#### B. Missing Module Namespacing Best Practices

**Issue:** Some modules don't properly use `NS()` for all IDs

**Example:**
```r
# PROBLEMATIC (found in mod_signature_nomination.R)
selectInput("gene_selector", ...)  # Missing ns()

# SHOULD BE:
selectInput(ns("gene_selector"), ...)
```

**Files Affected:**
- `inst/shiny/modules/mod_signature_nomination.R` (3 instances)
- `inst/shiny/modules/mod_comparison.R` (2 instances)

**Impact:** Namespace collisions if modules used multiple times

**Recommendation:** Audit all module UIs for proper `ns()` wrapping

---

#### C. Hard-coded CSS Instead of bslib Theming

**Issue:** Extensive hard-coded CSS instead of modern theming system

**Current Approach:**
```r
# inst/shiny/www/custom.css - 500+ lines of custom CSS
tags$head(tags$style(HTML("
  .sidebar { background-color: #ecf0f5; }
  .box { border-top: 3px solid #3c8dbc; }
  ...
")))
```

**Modern Approach:**
```r
library(bslib)
theme <- bs_theme(
  version = 5,
  bg = "#ecf0f5",
  fg = "#2c3e50",
  primary = "#3c8dbc",
  base_font = font_google("Source Sans Pro")
)
```

**Recommendation:** Migrate to bslib for themeable, responsive UI

---

#### D. Synchronous File I/O in Reactive Contexts

**Issue:** Blocking RDS reads in reactive expressions

**Example:**
```r
# PROBLEMATIC (blocks Shiny app)
reactive_data <- reactive({
  readRDS("large_file.rds")  # 500MB+ file
  # Process...
})

# BETTER: Use futures/promises
library(future)
library(promises)

reactive_data <- reactive({
  future_promise({
    readRDS("large_file.rds")
  })
})
```

**Files Affected:**
- `inst/shiny/R/data_manager.R` (Lines 100-150)
- Several module server functions

**Recommendation:** Implement async data loading with `future`/`promises`

---

### 1.2 Deprecated Package Dependencies

#### A. Old clusterProfiler Syntax

**Issue:** Using clusterProfiler 3.x syntax, now at 4.x

**Current Code:**
```r
# OLD (still works but deprecated)
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE
)

# NEW (recommended)
ego <- enrichGO(
  gene = gene_list,
  OrgDb = org.Hs.eg.db,
  ont = "BP",
  pAdjustMethod = "BH",
  readable = TRUE,
  pool = TRUE  # New parameter for performance
)
```

**Files Affected:**
- `R/enrichment_analysis.R` (Lines 1499-1850)

**Recommendation:** Update to clusterProfiler 4.x API, add `pool = TRUE` for performance

---

#### B. pheatmap vs ComplexHeatmap Inconsistency

**Issue:** Using both `pheatmap` and `ComplexHeatmap` inconsistently

**Current State:**
- Some modules use `pheatmap` (simpler, older)
- Others use `ComplexHeatmap` (more powerful, maintained)
- No clear pattern for when to use which

**Recommendation:** Standardize on `ComplexHeatmap` for all static heatmaps, use `heatmaply` for interactive

---

### 1.3 Compatibility Issues

#### A. Cross-Platform Path Handling

**Status:** ✅ **GOOD** - Already using `file.path()` and `normalizePath()`

**Example:**
```r
# Good practice already in use
data_path <- file.path(data_dir, "full_DE_results.rds")
data_path <- normalizePath(data_path, mustWork = FALSE)
```

**No Action Needed:** Path handling is already cross-platform compatible

---

#### B. Memory Management for Large Datasets

**Issue:** No explicit memory limits or streaming for 230k+ cell datasets

**Current:** Loading entire datasets into memory

**Recommendation:** Implement chunked reading for very large objects:
```r
# For very large SingleCellExperiment objects
read_sce_chunked <- function(file, max_cells = 50000) {
  full_sce <- readRDS(file)
  if (ncol(full_sce) > max_cells) {
    # Sample cells for interactive viewing
    sample_idx <- sample(seq_len(ncol(full_sce)), max_cells)
    full_sce[, sample_idx]
  } else {
    full_sce
  }
}
```

---

## Section 2: Perturb-seq Integration Gaps

### 2.1 Input Format Mismatches

#### A. Mixscale Output Integration

**Current State:** iSCORE-PDecipher expects specific RDS structure

**Mixscale Output Format:**
```r
# Mixscale Run_wmvRegDE() output
list(
  cluster_0 = list(
    perturbation1 = data.frame(
      gene_ID,
      log2FC,           # Simple column name
      p_weight,         # Simple p-value
      p_weight_BH,      # FDR corrected
      beta_weight,
      ...
    )
  )
)
```

**iSCORE-PDecipher Expected Format (Original):**
```r
# Experiment-split format (older)
list(
  CRISPRi_Mixscale = list(
    perturbation1 = list(
      cluster_0 = list(
        results = data.frame(
          gene_ID,
          log2FC_C12_FPD-24,      # Experiment-specific columns
          p_cell_typeC12_FPD-24:weight,
          ...
        )
      )
    )
  )
)
```

**Gap:** Format transformation required, not seamless

**Solution Implemented:** ✅ `import_pooled_mixscale_data()` function exists (R/import_pooled_mixscale_functions.R:74)

**Recommendation:** Promote pooled format as primary, make experiment-split legacy

---

#### B. scMAGeCK Output Integration

**Current State:** No direct scMAGeCK import

**scMAGeCK Output Format:**
```
# scMAGeCK RRA output (.txt file)
sgRNA    Gene    LFC    control_mean    treatment_mean    p_value    FDR
sgRNA1   LRRK2   -1.2   100             50                0.001      0.01
...

# Or cell-level output (h5ad/RDS)
```

**Gap:** No parser for scMAGeCK text files

**Recommendation:** Create `import_scmageck_results()` function

**Implementation:**
```r
#' Import scMAGeCK RRA Results
#'
#' @param scmageck_file Path to scMAGeCK RRA output file (.txt)
#' @param gene_column Column name for gene (default: "Gene")
#' @param lfc_column Column name for log2FC (default: "LFC")
#' @param fdr_column Column name for FDR (default: "FDR")
#' @return Data frame in iSCORE-PDecipher format
#' @export
import_scmageck_results <- function(
  scmageck_file,
  gene_column = "Gene",
  lfc_column = "LFC",
  fdr_column = "FDR"
) {
  # Read scMAGeCK output
  scmageck_data <- read.delim(scmageck_file)

  # Transform to iSCORE-PDecipher format
  # Map LFC → log2FC, FDR → p_val_adj
  transformed <- scmageck_data %>%
    rename(
      gene = !!sym(gene_column),
      log2FC = !!sym(lfc_column),
      p_val_adj = !!sym(fdr_column)
    ) %>%
    mutate(
      method = "scMAGeCK",
      analysis_type = "Perturb-seq"
    )

  return(transformed)
}
```

**Files to Create:**
- `R/import_scmageck_functions.R`
- `man/import_scmageck_results.Rd`
- `tests/testthat/test-import-scmageck.R`

---

### 2.2 Missing Perturb-seq-Specific Visualizations

#### A. Guide RNA Efficiency Plots

**Gap:** No visualization for sgRNA efficiency or coverage

**Recommendation:** Add module for guide QC

**Proposed UI:**
```r
mod_guide_qc_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput(ns("guide_coverage")),
    plotlyOutput(ns("guide_efficiency")),
    DTOutput(ns("low_efficiency_guides"))
  )
}
```

**Proposed Server:**
```r
mod_guide_qc_server <- function(id, data_reactive) {
  moduleServer(id, function(input, output, session) {
    output$guide_coverage <- renderPlotly({
      # Plot: Number of cells per guide
      # X: Guide, Y: Cell count
      # Color by efficiency tier
    })

    output$guide_efficiency <- renderPlotly({
      # Plot: Effect size vs statistical power
      # X: Observed log2FC, Y: -log10(p-value)
      # Size: Number of cells with guide
    })

    output$low_efficiency_guides <- renderDT({
      # Table of guides with <50 cells or weak effects
    })
  })
}
```

---

#### B. Multi-Perturbation Comparison Heatmaps

**Gap:** Current heatmaps focus on genes × clusters, not perturbations × pathways

**Current:** Gene-centric heatmaps

**Needed:** Perturbation-centric heatmaps

**Recommendation:** Add perturbation comparison view

**Example:**
```r
# Heatmap: Perturbations (rows) × Enrichment Terms (columns)
# Cell value: -log10(p.adjust)
# Clustering: Group perturbations by pathway similarity
```

**Implementation:** Extend `mod_heatmap_unified.R` with "Perturbation Mode"

---

#### C. Cell Type-Specific Perturbation Effects

**Gap:** Limited support for cell-type-stratified analysis

**Current:** Analysis by cluster, but cluster interpretation not explicit

**Recommendation:** Add cell type annotation integration

**Proposed Feature:**
```r
# Link clusters to cell types
cluster_annotations <- data.frame(
  cluster = paste0("cluster_", 0:14),
  cell_type = c(
    "Dopaminergic neurons",
    "GABAergic neurons",
    "Glutamatergic neurons",
    "Astrocytes",
    "Oligodendrocytes",
    "Microglia",
    ...
  )
)

# Visualize: Perturbation effects by cell type
# Heatmap: Perturbations × Cell Types
# Filter: Show only neuronal subtypes for PD-relevant analysis
```

**Files to Modify:**
- `inst/shiny/modules/mod_landing_page_with_umap_v2.R` - Add cell type labels
- `inst/shiny/modules/mod_de_results.R` - Add cell type filtering

---

### 2.3 Workflow Friction Points

#### A. Multi-Step Data Preparation

**Current Workflow:**
1. Run Mixscale externally
2. Export RDS files
3. Manually organize into directory structure
4. Run consolidation script
5. Launch app

**Friction:** Too many manual steps

**Recommendation:** Create `prepare_iscore_data()` wrapper

**Implementation:**
```r
#' Prepare iSCORE-PDecipher Data from Mixscale Output
#'
#' One-step preparation from Mixscale results to app-ready data
#'
#' @param mixscale_output Output from Run_wmvRegDE()
#' @param enrichment_output Output from DEenrich() or NULL to run internally
#' @param output_dir Directory to save app-ready files
#' @param run_enrichment Logical, run enrichment if not provided
#' @return Path to prepared data directory
#' @export
prepare_iscore_data <- function(
  mixscale_output,
  enrichment_output = NULL,
  output_dir,
  run_enrichment = TRUE
) {
  # 1. Convert Mixscale format
  message("Converting Mixscale format...")
  converted <- convert_mixscale_to_iscore(mixscale_output)

  # 2. Run enrichment if needed
  if (is.null(enrichment_output) && run_enrichment) {
    message("Running enrichment analysis...")
    enrichment_output <- run_enrichment_analysis(converted)
  }

  # 3. Consolidate
  message("Consolidating data...")
  consolidate_enrichment_results(enrichment_output, output_dir)

  # 4. Save DE results
  saveRDS(converted, file.path(output_dir, "full_DE_results.rds"))

  # 5. Validate
  message("Validating dataset...")
  validation <- validate_dataset_directory(output_dir)

  if (!validation$valid) {
    stop("Data preparation failed validation")
  }

  message("✓ Data ready! Launch with: launch_app('", output_dir, "')")
  return(output_dir)
}
```

**Benefit:** Single-function workflow from Mixscale to app

---

#### B. No Interactive Data Upload

**Current:** Must specify file paths in R console

**Recommendation:** Add file upload module for quick exploration

**Implementation:** Enhance `mod_data_loader.R`

```r
# Add to UI
fileInput(ns("upload_de"), "Upload DE Results (.rds)", accept = ".rds")
fileInput(ns("upload_enrichment"), "Upload Enrichment (.rds)", accept = ".rds")
actionButton(ns("load_uploaded"), "Load Uploaded Data")

# Add to server
observeEvent(input$load_uploaded, {
  req(input$upload_de, input$upload_enrichment)

  de_data <- readRDS(input$upload_de$datapath)
  enrich_data <- readRDS(input$upload_enrichment$datapath)

  # Validate format
  # Load into reactive
  # Update UI
})
```

**Benefit:** No command-line interaction needed, more accessible

---

## Section 3: Recommended Updates

### 3.1 Priority 1: Critical Fixes (Immediate)

#### Update 1.1: Fix Module Namespacing

**Issue:** Missing `ns()` wrappers in some modules

**Action:**
- Audit all module UI functions
- Wrap all input/output IDs with `ns()`
- Test for namespace collisions

**Files:**
- `inst/shiny/modules/mod_signature_nomination.R`
- `inst/shiny/modules/mod_comparison.R`

**Effort:** 2-4 hours
**Risk:** Low
**Impact:** Prevents future bugs

---

#### Update 1.2: Migrate to observeEvent() Pattern

**Issue:** Using older `observe({})` with manual invalidation

**Action:** Replace with `observeEvent()` for clarity

**Before:**
```r
observe({
  input$gene
  # Update logic
})
```

**After:**
```r
observeEvent(input$gene, {
  # Update logic
})
```

**Files:** All modules with observers
**Effort:** 8-12 hours
**Risk:** Low (one-to-one replacement)
**Impact:** Better code readability, easier debugging

---

#### Update 1.3: Update clusterProfiler to 4.x API

**Issue:** Using deprecated 3.x syntax

**Action:**
- Update `enrichGO()`, `enrichKEGG()` calls with new parameters
- Add `pool = TRUE` for performance
- Test all enrichment functions

**Files:**
- `R/enrichment_analysis.R`

**Effort:** 4-6 hours
**Risk:** Low (backwards compatible)
**Impact:** Performance improvement + future-proofing

---

### 3.2 Priority 2: Perturb-seq Focus (1-2 months)

#### Update 2.1: Create Unified Import Function

**Goal:** One function to import any Perturb-seq format

**Implementation:**
```r
#' Import Perturb-seq Results (Universal)
#'
#' Automatically detects and imports from Mixscale, scMAGeCK, or MAST
#'
#' @param data_source Path to file or directory
#' @param format Auto-detected or one of: "mixscale", "scmageck", "mast"
#' @param ... Additional arguments passed to specific importer
#' @return Standardized DE results list
#' @export
import_perturbseq_data <- function(
  data_source,
  format = "auto",
  ...
) {
  # Auto-detect format
  if (format == "auto") {
    format <- detect_perturbseq_format(data_source)
    message("Detected format: ", format)
  }

  # Route to appropriate importer
  result <- switch(format,
    "mixscale" = import_pooled_mixscale_data(data_source, ...),
    "scmageck" = import_scmageck_results(data_source, ...),
    "mast" = import_mast_data(data_source, ...),
    stop("Unknown format: ", format)
  )

  # Standardize output structure
  standardize_de_results(result)
}

detect_perturbseq_format <- function(data_source) {
  # Logic to detect based on file structure/columns
  if (is_directory(data_source)) {
    # Check for Mixscale cluster directories
    if (any(grepl("cluster_\\d+", list.dirs(data_source)))) {
      return("mixscale")
    }
    # Check for MAST pattern
    if (any(grepl("MAST", list.files(data_source)))) {
      return("mast")
    }
  } else if (file.exists(data_source)) {
    # Check file extension and contents
    if (grepl("\\.txt$", data_source)) {
      # Peek at first line
      first_line <- readLines(data_source, n = 1)
      if (grepl("sgRNA.*Gene.*LFC", first_line)) {
        return("scmageck")
      }
    }
  }

  stop("Could not auto-detect format")
}
```

**Effort:** 2-3 days
**Impact:** Major usability improvement

---

#### Update 2.2: Add Guide RNA QC Module

**Goal:** Visualize guide efficiency and coverage

**Implementation:**
- Create `inst/shiny/modules/mod_guide_qc.R`
- Add guide-level metrics to import functions
- Visualizations:
  - Coverage histogram (cells per guide)
  - Efficiency scatter (effect size vs significance)
  - Low-efficiency guide table

**Effort:** 1 week
**Impact:** Essential for Perturb-seq QC

---

#### Update 2.3: Enhance Perturbation Comparison Module

**Goal:** Better cross-perturbation analysis

**Features:**
- Perturbation × pathway heatmaps
- Clustering by pathway similarity
- Identify functionally similar perturbations
- Export perturbation signatures

**Implementation:** Extend existing `mod_comparison.R` with perturbation mode

**Effort:** 1-2 weeks
**Impact:** Core Perturb-seq functionality

---

### 3.3 Priority 3: UI/UX Modernization (2-4 months)

#### Update 3.1: Migrate to bslib Theming

**Goal:** Modern, themeable UI using Bootstrap 5

**Current:** Custom CSS + shinydashboard (Bootstrap 3)

**New:** bslib with Bootstrap 5

**Implementation:**
```r
# In global.R or app.R
library(bslib)

# Define theme
iscore_theme <- bs_theme(
  version = 5,
  bg = "#ffffff",
  fg = "#2c3e50",
  primary = "#3c8dbc",
  secondary = "#7fad39",
  success = "#28a745",
  danger = "#dc3545",
  base_font = font_google("Source Sans Pro"),
  heading_font = font_google("Roboto Slab"),
  code_font = font_google("Fira Code")
)

# Apply to app
ui <- page_navbar(
  title = "iSCORE-PDecipher",
  theme = iscore_theme,
  ...
)
```

**Benefits:**
- Responsive design (mobile-friendly)
- Dark mode support (free with bslib)
- Consistent styling
- Easier customization

**Effort:** 2-3 weeks (significant refactor)
**Risk:** Medium (UI changes may break some layouts)
**Impact:** Modern, professional appearance

---

#### Update 3.2: Add Progressive Data Loading

**Goal:** Don't block app startup with large data loads

**Implementation:**
```r
# Use promises for async loading
library(promises)
library(future)
plan(multisession)

# In data manager
load_data_async <- function(data_dir) {
  future_promise({
    readRDS(file.path(data_dir, "full_DE_results.rds"))
  }) %...>%
    (function(data) {
      # Data loaded, update reactive
      data
    })
}

# In app server
data_reactive <- reactiveVal(NULL)

observe({
  load_data_async(data_dir) %...>%
    (function(loaded_data) {
      data_reactive(loaded_data)
    })
})

# Show loading spinner until data ready
output$main_content <- renderUI({
  req(data_reactive())
  # Render main UI
})
```

**Benefits:**
- App starts immediately
- User sees loading progress
- Non-blocking experience

**Effort:** 1 week
**Impact:** Better UX for large datasets

---

#### Update 3.3: Implement Bookmarkable State

**Goal:** Allow users to save and share app state via URL

**Implementation:**
```r
# Enable bookmarking
shinyApp(
  ui = ui,
  server = server,
  enableBookmarking = "url"
)

# In server
observe({
  # Save inputs to bookmark
  setBookmarkExclude(c("large_data_reactive"))  # Don't save huge reactives
})

# Add bookmark button to UI
bookmarkButton(label = "Save State", icon = icon("bookmark"))
```

**Use Cases:**
- Share specific analysis view with collaborators
- Return to previous analysis state
- Reproducible sessions

**Effort:** 3-5 days
**Impact:** High for collaborative use

---

### 3.4 Priority 4: Nice-to-Have Enhancements (4-6 months)

#### Update 4.1: Add Tutorial Mode

**Goal:** Interactive tutorial for new users

**Implementation:** Use `cicerone` package

```r
library(cicerone)

# Define tour
guide <- Cicerone$new()$
  step(
    "gene_mutation",
    "Select Gene/Mutation",
    "Choose a gene or mutation to analyze"
  )$
  step(
    "cluster",
    "Select Cluster",
    "Filter by cell cluster"
  )$
  # ... more steps

# Add to UI
actionButton("start_tutorial", "Start Tutorial", icon = icon("question-circle"))

# In server
observeEvent(input$start_tutorial, {
  guide$init()$start()
})
```

**Effort:** 1 week
**Impact:** Lower barrier to entry for new users

---

#### Update 4.2: Export to Manuscript Formats

**Goal:** One-click export of publication-ready figures

**Features:**
- Export multi-panel figures in journal formats
- Automatic figure legends
- Export to PowerPoint for presentations
- Generate Methods section text

**Implementation:**
```r
library(officer)  # For PowerPoint export

export_to_ppt <- function(plots, title = "iSCORE-PDecipher Analysis") {
  ppt <- read_pptx()

  # Add title slide
  ppt <- add_slide(ppt, layout = "Title Slide", master = "Office Theme")
  ppt <- ph_with(ppt, value = title, location = ph_location_type(type = "ctrTitle"))

  # Add each plot
  for (i in seq_along(plots)) {
    ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
    ppt <- ph_with(ppt, value = plots[[i]], location = ph_location_type(type = "body"))
  }

  print(ppt, target = "analysis_figures.pptx")
}
```

**Effort:** 1-2 weeks
**Impact:** Streamlines publication workflow

---

## Section 4: Implementation Roadmap

### Phase 1: Critical Fixes (Month 1)

**Week 1-2:**
- ✅ Fix module namespacing issues
- ✅ Migrate to `observeEvent()` pattern
- ✅ Update clusterProfiler API

**Week 3-4:**
- ⬜ Comprehensive testing of fixes
- ⬜ Update documentation for changes
- ⬜ Version bump to 0.5.0

**Deliverables:**
- Bug-free codebase
- No deprecated warnings
- Passing all tests

---

### Phase 2: Perturb-seq Focus (Months 2-3)

**Month 2:**
- ⬜ Implement `import_perturbseq_data()` universal importer
- ⬜ Add `import_scmageck_results()`
- ⬜ Create `prepare_iscore_data()` wrapper
- ⬜ Update documentation with Perturb-seq workflows

**Month 3:**
- ⬜ Build `mod_guide_qc` module
- ⬜ Enhance perturbation comparison module
- ⬜ Add perturbation-centric heatmaps
- ⬜ Integration testing with real Mixscale/scMAGeCK data

**Deliverables:**
- Version 0.6.0
- Full Perturb-seq support
- Updated vignettes

---

### Phase 3: UI Modernization (Months 4-5)

**Month 4:**
- ⬜ Migrate to bslib theming
- ⬜ Implement async data loading
- ⬜ Add loading indicators
- ⬜ Responsive design testing

**Month 5:**
- ⬜ Add bookmarkable state
- ⬜ Implement dark mode
- ⬜ Polish UI/UX details
- ⬜ User testing with labmates

**Deliverables:**
- Version 0.7.0
- Modern, responsive UI
- Better performance

---

### Phase 4: Enhancements (Month 6)

**Month 6:**
- ⬜ Add interactive tutorial
- ⬜ Manuscript export features
- ⬜ Additional QC plots
- ⬜ Final polish for 1.0.0 release

**Deliverables:**
- Version 1.0.0
- Feature-complete Perturb-seq platform
- Publication-ready

---

## Section 5: Data Format Specifications

### 5.1 Recommended Perturb-seq Input Format

**Standardized Structure for iSCORE-PDecipher v1.0:**

```r
# Universal Perturb-seq DE Results Format
list(
  metadata = list(
    source = "Mixscale",  # or "scMAGeCK", "MAST"
    version = "1.0.0",
    date_created = Sys.Date(),
    organism = "Homo sapiens",
    cell_count = 230000,
    perturbation_count = 338,
    cluster_count = 6
  ),

  results = list(
    cluster_0 = list(
      LRRK2 = data.frame(
        gene_ID = c("ENSG00000123", ...),
        gene_symbol = c("SNCA", "GBA", ...),
        log2FC = c(1.2, -0.8, ...),
        p_value = c(0.001, 0.05, ...),
        p_value_adj = c(0.01, 0.1, ...),  # Recommend BH correction
        beta = c(0.5, -0.3, ...),
        se = c(0.1, 0.1, ...),
        statistic = c(5.0, -3.0, ...),
        mean_expression = c(100, 50, ...)
      ),
      # ... other perturbations
    ),
    # ... other clusters
  ),

  guide_metrics = list(
    LRRK2 = data.frame(
      guide_id = c("sgLRRK2_1", "sgLRRK2_2", ...),
      sequence = c("ACGT...", ...),
      cell_count = c(150, 200, ...),
      mean_efficiency = c(0.8, 0.9, ...)
    )
  )
)
```

**Benefits:**
- Self-documenting with metadata
- Consistent across all methods
- Includes guide-level QC info
- Easy validation

---

### 5.2 Mixscale Integration Example

**From Mixscale to iSCORE-PDecipher:**

```r
# 1. Run Mixscale analysis
library(Mixscale)
seurat_obj <- RunMixscale(seurat_obj, ...)
de_results <- Run_wmvRegDE(seurat_obj, ...)

# 2. Prepare for iSCORE-PDecipher (NEW function)
library(iSCORE.PDecipher)
prepared_data <- prepare_iscore_data(
  mixscale_output = de_results,
  enrichment_output = NULL,  # Will run automatically
  output_dir = "./iscore_ready/",
  run_enrichment = TRUE
)

# 3. Launch app
launch_app(data_dir = prepared_data)
```

**Single workflow, no manual steps!**

---

## Section 6: Testing Strategy

### 6.1 Unit Tests for New Functions

**Add tests for:**
- `import_perturbseq_data()` with all formats
- `import_scmageck_results()` with sample data
- `prepare_iscore_data()` end-to-end
- Format detection logic
- Data validation functions

**Framework:** testthat

**Example:**
```r
# tests/testthat/test-import-perturbseq.R
test_that("import_perturbseq_data detects Mixscale format", {
  temp_dir <- create_mock_mixscale_data()
  result <- import_perturbseq_data(temp_dir, format = "auto")

  expect_equal(attr(result, "source"), "mixscale")
  expect_true("cluster_0" %in% names(result$results))
})

test_that("import_scmageck_results parses correctly", {
  scmageck_file <- system.file("extdata", "scmageck_sample.txt", package = "iSCORE.PDecipher")
  result <- import_scmageck_results(scmageck_file)

  expect_true("log2FC" %in% names(result))
  expect_true("p_val_adj" %in% names(result))
})
```

---

### 6.2 Integration Tests

**Test Scenarios:**
1. **End-to-end Mixscale workflow:** Raw Mixscale output → prepared data → app launch
2. **Cross-platform:** Test on Windows, Mac, Linux
3. **Large dataset handling:** 230k cells, 338 perturbations
4. **Module interaction:** Ensure filters sync across modules

**Automated Testing:**
```r
# Use shinytest2 for app testing
library(shinytest2)

test_that("App launches with Perturb-seq data", {
  app <- AppDriver$new()

  # Check landing page loads
  app$expect_values(output = "umap_plot")

  # Select perturbation
  app$set_inputs(gene_mutation = "LRRK2")

  # Check volcano plot updates
  app$expect_values(output = "volcano_mixscale")

  # Check enrichment table populates
  app$expect_values(output = "enrichment_table")
})
```

---

## Section 7: Migration Guide for Users

### 7.1 For Existing Users (0.4.0 → 1.0.0)

**Breaking Changes:**
- ❌ Old import functions deprecated (still work with warnings)
- ❌ Some UI element IDs changed (affects bookmarks)
- ❌ CSS classes updated (affects custom styling)

**Migration Steps:**

**Step 1:** Update package
```r
remotes::install_github("jessedunnack/iSCORE-PDecipher@v1.0.0")
```

**Step 2:** Update existing data preparation scripts
```r
# OLD (still works but deprecated)
de_data <- import_mixscale_data(dir, pval_column = "p_weight_BH")

# NEW (recommended)
de_data <- import_perturbseq_data(dir, format = "mixscale")
```

**Step 3:** Re-run data preparation if using new features
```r
# Use new wrapper for enrichment
prepared_data <- prepare_iscore_data(
  mixscale_output = your_mixscale_output,
  output_dir = "./new_format/"
)
```

---

### 7.2 For New Users (Starting with 1.0.0)

**Recommended Workflow:**

```r
# Install
remotes::install_github("jessedunnack/iSCORE-PDecipher")

# Load
library(iSCORE.PDecipher)

# Prepare data from Mixscale
prepared_data <- prepare_iscore_data(
  mixscale_output = your_de_results,
  output_dir = "~/iscore_data/"
)

# Launch app
launch_app()  # Interactive dataset selection

# Or specify directly
launch_app(data_dir = prepared_data)
```

**Resources:**
- Tutorial vignette: `vignette("getting-started", package = "iSCORE.PDecipher")`
- Perturb-seq guide: `vignette("perturbseq-workflow", package = "iSCORE.PDecipher")`
- In-app tutorial: Click "Start Tutorial" button

---

## Conclusion

### Summary of Key Recommendations

**Immediate (Priority 1):**
1. Fix module namespacing
2. Update to modern reactive patterns
3. Update clusterProfiler API

**Short-term (Priority 2):**
1. Create universal Perturb-seq importer
2. Add guide RNA QC module
3. Enhance perturbation comparison

**Medium-term (Priority 3):**
1. Migrate to bslib theming
2. Implement async data loading
3. Add bookmarkable state

**Long-term (Priority 4):**
1. Interactive tutorials
2. Manuscript export features
3. Advanced QC visualizations

### Expected Impact

**User Experience:**
- Streamlined workflow: Mixscale → iSCORE-PDecipher in one step
- Modern, responsive UI
- Better performance with large datasets

**Scientific Value:**
- Perturb-seq focused visualizations
- Guide QC for experimental quality
- Publication-ready figures

**Maintainability:**
- Modern Shiny patterns
- Comprehensive tests
- Clear documentation

### Timeline to Production

**6-month roadmap:**
- Month 1: Critical fixes
- Months 2-3: Perturb-seq focus
- Months 4-5: UI modernization
- Month 6: Polish and 1.0.0 release

**Version 1.0.0 will be a state-of-the-art Perturb-seq analysis platform specifically tailored for CRISPRi screens in Parkinson's disease research.**

---

**Document Status:** Complete
**Next Steps:** Prioritize updates based on user feedback and development resources
**Contact:** jessedunnack@berkeley.edu for questions or implementation assistance
