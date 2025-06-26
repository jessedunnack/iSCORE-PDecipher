# Bidirectional Sync Implementation Plan
## iSCORE-PDecipher App Enhancement

### Date: January 18, 2025

## 🎯 Overview
Implement bidirectional synchronization between global settings panel and local module interactions, ensuring consistent state across all app components while preventing circular updates.

### Key Requirements
1. Local interactions (like clicking UMAP clusters) update global settings
2. When switching pages, the new page shows data based on updated global settings
3. Analysis Type dropdown shows dataset-specific options (MAST, CRISPRi, CRISPRa)
4. Invalid gene/method combinations are prevented in dropdowns

## 📋 Implementation Phases

### Phase 1: Data Infrastructure (30 mins)
**Goal**: Establish reactive data structures for tracking available options and current state

#### 1.1 Enhance app_data structure in app.R
```r
app_data <- reactiveValues(
  consolidated_data = NULL,
  data_loaded = FALSE,
  # NEW: Track available options dynamically
  available_methods = list(),      # e.g., list(MAST=TRUE, CRISPRi=TRUE, CRISPRa=FALSE)
  available_genes_by_method = list(),  # e.g., list(MAST=c("LRRK2",...), CRISPRi=c(...))
  # NEW: Global state tracking
  update_in_progress = FALSE,      # Flag to prevent circular updates
  last_update_source = NULL        # Track where update originated
)
```

#### 1.2 Create helper functions
- `detect_available_methods()` - Scan consolidated data for available analysis types
- `get_valid_genes_for_method()` - Return genes available for selected method
- `is_valid_combination()` - Check if gene/method combo exists in data

### Phase 2: Dynamic Dropdowns (45 mins)
**Goal**: Make Analysis Type and Gene dropdowns respond to available data

#### 2.1 Update Analysis Type dropdown
```r
# Replace static choices with dynamic detection
observe({
  req(app_data$data_loaded)
  
  # Detect what's actually in the data
  methods_available <- detect_available_methods(app_data$consolidated_data)
  
  # Create user-friendly labels
  choices <- c()
  if (methods_available$MAST) choices["iSCORE-PD (MAST)"] <- "MAST"
  if (methods_available$CRISPRi) choices["PerturbSeq (CRISPRi)"] <- "MixScale_CRISPRi"
  if (methods_available$CRISPRa) choices["PerturbSeq (CRISPRa)"] <- "MixScale_CRISPRa"
  
  updateSelectInput(session, "global_analysis_type", choices = choices)
})
```

#### 2.2 Filter gene dropdown based on analysis type
```r
observe({
  req(input$global_analysis_type)
  
  # Get only genes valid for selected method
  valid_genes <- app_data$available_genes_by_method[[input$global_analysis_type]]
  
  updateSelectInput(session, "global_gene", 
                   choices = valid_genes,
                   selected = ifelse(input$global_gene %in% valid_genes, 
                                   input$global_gene, 
                                   valid_genes[1]))
})
```

### Phase 3: Bidirectional Communication (1 hour)
**Goal**: Enable modules to update global settings and vice versa

#### 3.1 Create session-level message handlers in app.R
```r
# Handler for cluster selection from UMAP clicks
observeEvent(input$update_cluster_from_module, {
  if (!app_data$update_in_progress) {
    app_data$update_in_progress <- TRUE
    app_data$last_update_source <- "module"
    
    updateSelectInput(session, "global_cluster", 
                     selected = input$update_cluster_from_module)
    
    app_data$update_in_progress <- FALSE
  }
})
```

#### 3.2 Update mod_de_results.R
```r
# Add to plotly click event handler (around line 300)
observeEvent(event_data("plotly_click", source = "umap_plot"), {
  clicked_data <- event_data("plotly_click", source = "umap_plot")
  if (!is.null(clicked_data)) {
    selected_cluster <- paste0("cluster_", clicked_data$key)
    values$selected_cluster <- selected_cluster
    
    # NEW: Update global settings
    session$sendInputMessage("update_cluster_from_module", 
                           list(value = selected_cluster))
  }
})

# Re-enable initialization from global but with protection
observe({
  # Only update if not currently processing a local change
  if (!isTRUE(values$updating_global)) {
    values$selected_cluster <- global_selection()$cluster
  }
})
```

### Phase 4: Module Integration (45 mins)
**Goal**: Update all modules to respect bidirectional sync

#### 4.1 Standard module pattern
Each module should follow this pattern:
```r
# At module start
local_updating <- reactiveVal(FALSE)

# When receiving global updates
observe({
  if (!local_updating()) {
    # Update local state from global
    updateSelectInput(session, ns("local_gene"), 
                     selected = global_selection()$gene)
  }
})

# When making local changes
observeEvent(input$local_gene, {
  local_updating(TRUE)
  session$sendInputMessage("update_gene_from_module", 
                         list(value = input$local_gene))
  local_updating(FALSE)
})
```

#### 4.2 Priority modules to update
1. **mod_de_results.R** - UMAP cluster clicks
2. **mod_visualization_enhanced.R** - If any local gene selection
3. **mod_heatmap.R** - If any local filtering
4. **mod_comparison.R** - Method selection updates

### Phase 5: Testing & Validation (30 mins)
**Goal**: Ensure robust operation without circular updates

#### 5.1 Test scenarios
1. Click UMAP cluster → Global updates → Other modules update
2. Change global gene → All modules update → No loops
3. Change analysis type → Gene list updates → Invalid selections handled
4. Rapid clicking → No race conditions

#### 5.2 Debug helpers
```r
# Add debug mode
if (isTRUE(getOption("iscore.debug"))) {
  observe({
    cat("Global state changed:", 
        "Method:", input$global_analysis_type,
        "Gene:", input$global_gene,
        "Cluster:", input$global_cluster, "\n")
  })
}
```

## 🔧 Technical Details

### Preventing Circular Updates
1. Use `update_in_progress` flag in app_data
2. Check flag before processing updates
3. Use `isolate()` when reading values that might trigger updates
4. Implement `local_updating` flags in modules

### State Consistency
1. All state changes go through global settings first
2. Modules update their local state from global
3. Use reactive expressions for derived values
4. Validate state transitions

## 📊 Implementation Order & Risk Assessment

**Low Risk - Start Here:**
1. ✅ Data infrastructure (Phase 1)
2. ✅ Dynamic dropdowns (Phase 2)
3. ✅ Debug logging

**Medium Risk:**
4. ⚠️ DE Results module update (Phase 3.2)
5. ⚠️ Global message handlers (Phase 3.1)

**Higher Risk - Do Last:**
6. ⚡ Full module integration (Phase 4)
7. ⚡ Cross-module testing

## 🔄 Rollback Plan
```bash
# Before starting
git tag v0.1.5-pre-bidirectional-sync

# If issues arise
git checkout v0.1.5-pre-bidirectional-sync
```

## ✅ Success Criteria
- ✓ UMAP clicks update global cluster selection
- ✓ Analysis Type shows only available options
- ✓ Gene dropdown filters by analysis type
- ✓ No circular update loops
- ✓ All modules stay synchronized
- ✓ Performance remains responsive

## ⏱️ Estimated Total Time: 3-4 hours

## 📝 Notes
- Start with low-risk phases to validate approach
- Test incrementally after each phase
- Keep old behavior as fallback during development
- Monitor performance impact of additional reactivity