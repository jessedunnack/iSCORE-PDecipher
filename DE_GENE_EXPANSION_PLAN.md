# DE Gene Expansion Plan for iSCORE-PDecipher
## Implementation Plan Document - June 18, 2025

### 📋 **COMPREHENSIVE DE GENE EXPANSION PLAN**

#### **🔍 Current App Structure Analysis (Completed)**

##### **Data Architecture**
1. **Enrichment Data**: `consolidated_data` (pathway-level results, 767K+ terms)
   - Columns: `mutation_perturbation`, `cluster`, `enrichment_type`, `direction`, `p.adjust`, `Description`, `method`
   - Used by: Basic Visualization, Heatmap, Method Comparison, KEGG Pathview

2. **DE Data**: `full_DE_results.rds` (gene-level statistics: log2FC, p-values)
   - Structure: `iSCORE_PD_MAST[gene][cluster]$results` with `avg_log2FC`, `p_val_adj`
   - Structure: `CRISPRi_Mixscale[gene][cluster]$results` with `log2FC_[experiment]`, `p_cell_type[experiment]:weight`
   - Currently used by: DE Results module (volcano plots only)

3. **Metadata**: `final_dataset_metadata.rds` (222K cells, 105 variables)
   - Critical columns: `dataset`, `gene`, `crispr_modality`, `mutation_tidy`, `seurat_clusters`

##### **Reactive Data Flow Pattern**
```
Global Sidebar → global_data_selection() → filtered_data() → Individual Modules
```

##### **Current Tabs Structure**
- Overview (UMAP), DE Results (volcano), Basic Visualization (enrichment), Method Comparison, Heatmap (enrichment), KEGG Pathview, Export

#### **🎯 User Requirements (From Q&A Session)**

1. **Tab Organization**: Option B - Two new dedicated tabs
2. **Visualizations**: Gene overlaps (Venn, UpSet), correlation plots, scatter plots, bar charts  
3. **Controls**: Module-specific settings (independent of global sidebar)
4. **Data Filtering**: Maintain current pre-filter pattern
5. **Performance**: Large scale (thousands of genes, all perturbations) - need sampling approach
6. **Integration**: Interactive hover/click features for cross-referencing
7. **⚠️ CRITICAL**: Distinguish mutation source (iSCORE-PD vs PerturbSeq) in ALL heatmaps

#### **🧬 Metadata Insights (Key Discovery)**

**Critical Columns for Implementation:**
- `dataset`: "PerturbSeq" vs "iSCORE-PD" (PERFECT for source distinction!)
- `gene`/`scMAGeCK_gene_assignment`: Target genes
- `crispr_modality`: "CRISPRa" vs "CRISPRi"  
- `mutation_tidy`: iSCORE-PD mutations
- `seurat_clusters`: Cluster assignments
- `nCount_RNA`, `nFeature_RNA`, `percent.mt`: Quality metrics for summaries

### **🏗️ IMPLEMENTATION PLAN**

#### **Phase 1: Critical Infrastructure Changes**

##### **1.1 Source Distinction System (URGENT - Affects Existing Features)**
**Problem**: Current heatmaps show "ATP13A2" without distinguishing iSCORE-PD mutation vs CRISPRi perturbation

**Solution**: 
```r
# New function to create source-aware labels
create_source_labels <- function(gene, method, dataset) {
  case_when(
    method == "MAST" ~ paste0(gene, " (iSCORE-PD)"),
    method == "MixScale" & dataset == "CRISPRi" ~ paste0(gene, " (CRISPRi)"),
    method == "MixScale_CRISPRa" ~ paste0(gene, " (CRISPRa)"),
    TRUE ~ gene
  )
}
```

**Files to Update:**
- `modules/mod_heatmap.R` - Add source distinction to existing enrichment heatmaps
- All new DE modules - Apply from start

##### **1.2 Enhanced Data Processing Pipeline**
**New Core Functions Needed:**
```r
# Enhanced DE data processing
process_de_data_with_metadata <- function(de_results, metadata) {
  # Combine DE results with cell metadata
  # Add source distinction
  # Calculate summary statistics
}

# Multi-condition DE comparison
prepare_multi_condition_de <- function(conditions, cutoffs) {
  # Apply filters per condition
  # Combine for cross-comparison
  # Handle thousands of genes efficiently
}

# Sampling strategy for large datasets
sample_de_for_preview <- function(de_data, max_genes = 500) {
  # Intelligent sampling (top by significance, random subset)
  # Preserve representative gene sets
}
```

#### **Phase 2: New Tab Implementation**

##### **2.1 Tab 1: "DE Gene Analysis"**
**Purpose**: Cross-condition comparisons, overlaps, correlations

**UI Components:**
```r
mod_de_analysis_ui <- function(id) {
  # Control Panel:
  - Dataset selection (iSCORE-PD only, PerturbSeq only, Both)
  - Gene/mutation multi-select
  - Cluster selection (all, specific, range)
  - Significance cutoffs (log2FC, p-value)
  - Analysis type selector

  # Visualization Panel:
  - tabsetPanel(
      "Gene Overlaps" - Venn diagrams, UpSet plots
      "Correlations" - Cross-condition log2FC correlations  
      "Scatter Plots" - Pairwise condition comparisons
      "Top Genes" - Bar charts of most significant
    )
}
```

**Server Logic:**
```r
mod_de_analysis_server <- function(id, app_data) {
  # Reactive data processing
  de_comparison_data <- reactive({
    # Process selected conditions
    # Apply cutoffs
    # Calculate overlaps and correlations
  })
  
  # Visualization outputs
  output$overlap_plot <- renderPlotly({ ... })
  output$correlation_plot <- renderPlotly({ ... })
  output$scatter_plot <- renderPlotly({ ... })
  output$top_genes_plot <- renderPlotly({ ... })
}
```

##### **2.2 Tab 2: "DE Gene Heatmaps"**
**Purpose**: Interactive heatmap builder with custom controls

**UI Components:**
```r
mod_de_heatmap_ui <- function(id) {
  # Advanced Control Panel:
  - Data scope: Full dataset vs specific clusters
  - Source filter: iSCORE-PD, CRISPRi, CRISPRa, All
  - Gene selection: All, specific genes, top N by significance
  - Cutoff controls: log2FC threshold, p-value threshold
  - Display options: Clustering, color scales, annotations
  - Performance: Preview (sample) vs Full render

  # Heatmap Panel:
  - Dynamic heatmap with source-distinguished column names
  - Interactive features: zoom, hover, click
  - Export options: PNG, PDF, data table
}
```

**Server Logic with Performance Optimization:**
```r
mod_de_heatmap_server <- function(id, app_data) {
  # Smart data processing
  de_heatmap_data <- reactive({
    # Apply user filters
    # Create source-aware labels: "ATP13A2 (iSCORE-PD)" vs "ATP13A2 (CRISPRi)"
    # Implement sampling for large datasets
  })
  
  # Progressive rendering
  preview_data <- reactive({
    sample_de_for_preview(de_heatmap_data(), max_genes = 500)
  })
  
  output$heatmap_preview <- renderHeatmaply({ ... })
  output$heatmap_full <- renderHeatmaply({ ... })  # On-demand
}
```

#### **Phase 3: Interactive Integration Features**

##### **3.1 Cross-Module Communication**
```r
# Shared reactive values for cross-module communication
shared_selections <- reactiveValues(
  selected_genes = NULL,
  selected_pathways = NULL,
  current_comparison = NULL
)

# Click handlers for gene-pathway integration
observeEvent(input$gene_click, {
  # Find pathways associated with clicked gene
  # Highlight in enrichment modules
})

observeEvent(input$pathway_click, {
  # Find genes in clicked pathway
  # Highlight in DE modules
})
```

##### **3.2 Summary Statistics Integration**
```r
# Metadata-driven summary statistics
generate_condition_summary <- function(condition, metadata) {
  # Cell counts per cluster
  # Quality metrics (nCount_RNA, percent.mt)
  # Perturbation efficiency metrics
  # Experimental metadata (age, duration, etc.)
}
```

#### **Phase 4: Enhanced User Experience**

##### **4.1 Progressive Loading Strategy**
1. **Initial Load**: Show preview with sampled data (fast)
2. **User Decision**: Option to render full dataset 
3. **Background Processing**: Full analysis with progress indicators
4. **Smart Caching**: Cache large computations

##### **4.2 Intuitive Interactions (Non-Clunky)**
- **Hover**: Show gene details, pathway membership, statistics
- **Click**: Select genes for cross-module highlighting
- **Brush/Select**: Multi-gene selection for further analysis
- **Context Menus**: Right-click for quick actions (export, analyze)

### **📊 Technical Architecture**

#### **New Module Structure**
```
inst/shiny/modules/
├── mod_de_analysis.R          # NEW - Cross-condition analysis
├── mod_de_heatmap.R           # NEW - Interactive DE heatmaps  
├── mod_heatmap.R              # UPDATE - Add source distinction
└── shared/
    ├── de_processing.R        # NEW - DE data utilities
    ├── source_labeling.R      # NEW - Source distinction functions
    └── interactive_utils.R    # NEW - Cross-module communication
```

#### **Data Flow Enhancement**
```r
# Enhanced app_data structure
app_data <- reactiveValues(
  consolidated_data = NULL,      # Existing enrichment data
  de_data = NULL,                # Enhanced DE data with metadata
  metadata = NULL,               # Cell-level metadata
  source_labels = NULL,          # Source-aware gene labels
  shared_selections = NULL       # Cross-module communication
)
```

### **⚠️ CRITICAL IMMEDIATE CHANGES REQUIRED**

#### **1. Fix Existing Enrichment Heatmaps**
**Current Problem**: Shows "ATP13A2" without distinguishing source
**Required Fix**: Update `mod_heatmap.R` to use source-aware labels

#### **2. Update Tab Structure in app.R**
**Add New Tabs**:
```r
# Add after existing tabs
tabPanel(
  "DE Gene Analysis",
  icon = icon("project-diagram"),
  value = "de_analysis",
  br(),
  mod_de_analysis_ui("de_analysis_module")
),

tabPanel(
  "DE Gene Heatmaps", 
  icon = icon("fire"),
  value = "de_heatmap",
  br(),
  mod_de_heatmap_ui("de_heatmap_module")
)
```

### **🚀 Implementation Priority**

#### **Priority 1 (Immediate)**
1. ✅ Source distinction system implementation
2. ✅ Fix existing enrichment heatmaps  
3. ✅ Enhanced DE data processing pipeline

#### **Priority 2 (Core Features)**
1. ✅ DE Gene Analysis tab (overlaps, correlations, scatter plots)
2. ✅ DE Gene Heatmaps tab (interactive builder)
3. ✅ Performance optimization (sampling strategy)

#### **Priority 3 (Enhanced Experience)**
1. ✅ Cross-module communication
2. ✅ Interactive hover/click features
3. ✅ Metadata-driven summaries

### **📝 Future Considerations**

#### **Global Sidebar Evolution**
- User indicated potential removal of global sidebar
- New modules designed with independent controls
- Easy transition if sidebar architecture changes

#### **Scalability Features**
- Data export for external analysis
- API endpoints for programmatic access
- Advanced filtering and search capabilities

---

## **🎯 IMPLEMENTATION STEPS**

### **Step 1: Source Distinction System**
1. Create `shared/source_labeling.R` with labeling functions
2. Update existing `mod_heatmap.R` to use source-aware labels
3. Test with existing enrichment heatmaps

### **Step 2: Enhanced DE Data Pipeline**
1. Create `shared/de_processing.R` with enhanced processing functions
2. Integrate metadata loading and processing
3. Implement sampling strategy for performance

### **Step 3: DE Gene Analysis Tab**
1. Create `mod_de_analysis.R` with UI and server functions
2. Implement gene overlap analysis (Venn, UpSet)
3. Add correlation and scatter plot functionality
4. Add top genes bar chart visualization

### **Step 4: DE Gene Heatmaps Tab**
1. Create `mod_de_heatmap.R` with advanced controls
2. Implement progressive loading (preview → full)
3. Add interactive features and export options
4. Apply source distinction to heatmap columns

### **Step 5: Integration Features**
1. Create `shared/interactive_utils.R` for cross-module communication
2. Implement hover/click interactions
3. Add metadata-driven summary statistics
4. Test cross-module integration

---

**Created**: June 18, 2025  
**Status**: Ready for Implementation  
**Priority**: Phase 1 (Source Distinction) → Phase 2 (New Tabs) → Phase 3 (Integration)