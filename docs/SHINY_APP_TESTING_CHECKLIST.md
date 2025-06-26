# Shiny App Runtime Testing Checklist

## 🚀 **Priority 1: Basic App Launch & Stability**

### **App Startup**
- [ ] **Launch app successfully** without crashes
  ```r
  library(iSCORE.PDecipher)
  launch_iscore_app()
  ```
- [ ] **Dataset selection works** (single prompt, no double prompts)
- [ ] **Data loads completely** (should show 663,280+ enrichment terms)
- [ ] **No console errors** during startup
- [ ] **All tabs accessible** without crashes

### **Core Navigation**
- [ ] **UMAP tab loads** properly with visualization
- [ ] **Enrichment Analysis tabs** all accessible
- [ ] **DE Analysis tabs** accessible (critical - recent focus of fixes)
- [ ] **Settings panel** responsive to changes

## 🔧 **Priority 2: Critical Module Testing**

### **DE Heatmap Module (Recent Problem Area)**
- [ ] **DE Heatmap tab loads** without hanging
- [ ] **Can select different clusters** without freezing
- [ ] **Data processes within reasonable time** (<30 seconds)
- [ ] **Heatmap renders** or shows appropriate error message
- [ ] **Memory usage stable** during DE processing

### **Interactive Heatmap Module (Known Issues)**
- [ ] **Heatmap tab accessible**
- [ ] **Heatmap generates** with row annotations **OFF** (default setting)
- [ ] **Test row annotations** (may fail with color errors - document if so)
- [ ] **Export functionality works** (HTML/PDF download)
- [ ] **Different heatmap types work** (p-value, fold enrichment, etc.)

### **Core Visualization Modules**
- [ ] **Dotplot generates** properly with all plot types (dot, bar, lollipop)
- [ ] **UMAP shows clusters** correctly with proper dimensions (950×700px)
- [ ] **Volcano plots work** in DE Results module
- [ ] **Global settings update** all visualizations reactively

## 📊 **Priority 3: Data Filtering & Reactivity**

### **Global Settings Panel**
- [ ] **Gene selection updates** all modules
- [ ] **Cluster selection updates** all modules  
- [ ] **Enrichment type filtering** works across modules
- [ ] **Direction filtering** (ALL/UP/DOWN) works correctly
- [ ] **Method filtering** (MAST/CRISPRi/CRISPRa) works if data supports it

### **Cross-Module Communication**
- [ ] **Changes in one module** reflect in others
- [ ] **No circular update loops** causing freezing
- [ ] **Reactive updates complete** within reasonable time
- [ ] **No data corruption** between module switches

## 🔍 **Priority 4: Performance & Edge Cases**

### **Large Dataset Handling**
- [ ] **Memory usage stable** with 200K+ cell datasets
- [ ] **Visualizations render** within reasonable time
- [ ] **App remains responsive** during data processing
- [ ] **No browser tab crashes** due to memory

### **Edge Cases**
- [ ] **Empty results handled gracefully** (e.g., no significant terms)
- [ ] **Invalid combinations handled** (e.g., gene not in selected method)
- [ ] **Error messages informative** and user-friendly
- [ ] **Fallback mechanisms work** when expected

## 🐛 **Known Issues to Verify Status**

### **Issues from Recent Commits**
- [ ] **DE heatmap hanging resolved** (commits 44e3bef, 08c6ff0)
- [ ] **Automatic DE loading disabled** properly
- [ ] **No more massive DE processing** causing timeouts
- [ ] **Targeted cluster extraction working**

### **Issues from CLAUDE.md**
- [ ] **Heatmap row annotations** work or are properly disabled
- [ ] **Color palette errors** resolved or documented workaround
- [ ] **Interactive heatmaply** functions without crashes

## 📝 **Testing Protocol**

### **For Each Test:**
1. **Document the result** (✅ Pass / ❌ Fail / ⚠️ Partial)
2. **Note any error messages** or console warnings
3. **Record timing** for slow operations
4. **Screenshot issues** if visual problems occur

### **Test Data Combinations:**
- [ ] **MAST data**: Test with LRRK2, cluster_0, GO_BP
- [ ] **CRISPRi data**: Test if available in your dataset
- [ ] **CRISPRa data**: Test if available in your dataset
- [ ] **Mixed analysis**: Test intersection/union modes

### **Browser Compatibility:**
- [ ] **Chrome/Chromium**: Primary testing browser
- [ ] **Firefox**: Secondary testing if issues found
- [ ] **Safari/Edge**: Optional if cross-browser issues suspected

## 🎯 **Success Criteria**

### **Minimum Viable:**
- App launches without crashes
- Basic visualization modules work
- Data filtering functional
- No major memory leaks

### **Fully Functional:**
- All modules accessible and responsive
- DE heatmap module stable
- Interactive heatmaps work (with or without annotations)
- Export functionality operational
- Performance acceptable for large datasets

---

## 📋 **Quick Testing Commands**

```r
# Basic launch test
library(iSCORE.PDecipher)
launch_iscore_app()

# Check data loading in R console
data <- readRDS("path/to/enrichment/data.rds")
cat("Loaded", nrow(data), "enrichment terms\n")

# Memory monitoring
gc()  # Check memory usage
```

**Priority Testing Order:** Basic Launch → DE Heatmap → Interactive Heatmap → Other Modules → Performance → Edge Cases