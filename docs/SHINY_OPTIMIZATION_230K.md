# Shiny App Optimization for 230K+ Cells

## Current Status
The iSCORE-PDecipher Shiny app needs optimization to handle the expanded dataset (230K+ cells from Batch 7).

## Key Findings

### 1. **UMAP Visualization (Critical Path)**
- **Current**: Using `renderPlot()` with ggplot2 for static UMAP
- **Data Format**: SingleCellExperiment object with UMAP coordinates
- **Performance**: Already optimized with static plots instead of plotly for large datasets
- **Cache**: Not currently caching UMAP plots

### 2. **Existing Optimizations**
- ✅ Cache Manager already implemented (`R/cache_manager.R`)
- ✅ Static plots for UMAP (ggplot2 instead of plotly)
- ✅ Lazy loading of DE data (loads on demand)
- ✅ SingleCellExperiment for lightweight UMAP data storage
- ✅ Multiple PC options (30pc, 100pc) for different analysis needs

### 3. **Current Bottlenecks**
- [ ] UMAP plot regeneration on every cluster selection
- [ ] No caching of rendered UMAP plots
- [ ] Full dataset loaded even when filtering
- [ ] Volcano plots using plotly (interactive) can be slow

## Optimization Recommendations

### Phase 1: Quick Wins (Immediate)

#### 1.1 Cache UMAP Plots
```r
# Add to mod_de_results.R
umap_plot_cache <- CacheManager$new(max_size = 20, ttl = 3600)

output$umap_plot <- renderPlot({
  cache_key <- paste(dataset_name(), input$pc_selection, input$color_by, collapse = "_")
  
  cached_plot <- umap_plot_cache$get(cache_key)
  if (!is.null(cached_plot)) {
    return(cached_plot)
  }
  
  # Generate plot
  p <- generate_umap_plot(values$umap_data)
  umap_plot_cache$set(cache_key, p)
  return(p)
})
```

#### 1.2 Implement Data Sampling for Preview
```r
# For initial display, use subset
preview_data <- function(full_data, n = 50000) {
  if (nrow(full_data) > n) {
    set.seed(123)
    sample_idx <- sample(nrow(full_data), n)
    return(full_data[sample_idx, ])
  }
  return(full_data)
}
```

#### 1.3 Pre-render Common Views
```r
# On app startup, pre-render common UMAP views
prerender_common_plots <- function() {
  common_views <- list(
    list(dataset = "dataset1", pc = "30", color = "cluster"),
    list(dataset = "dataset2", pc = "30", color = "cluster")
  )
  
  for (view in common_views) {
    plot <- generate_umap_plot(view)
    cache_key <- paste(view, collapse = "_")
    umap_plot_cache$set(cache_key, plot)
  }
}
```

### Phase 2: Architecture Improvements

#### 2.1 Implement Progressive Loading
```r
# Load data in chunks
load_data_progressive <- function(path, chunk_size = 10000) {
  total_cells <- get_cell_count(path)
  
  withProgress(message = 'Loading data', value = 0, {
    for (i in seq(1, total_cells, chunk_size)) {
      chunk <- load_chunk(path, i, min(i + chunk_size - 1, total_cells))
      process_chunk(chunk)
      incProgress(chunk_size / total_cells)
    }
  })
}
```

#### 2.2 Use Data Tables for Large Results
```r
# Replace renderTable with DT for large tables
output$de_results_table <- DT::renderDataTable({
  datatable(de_results,
    extensions = 'Scroller',
    options = list(
      deferRender = TRUE,
      scrollY = 400,
      scroller = TRUE,
      scrollX = TRUE
    )
  )
})
```

#### 2.3 Implement Cluster-Specific Loading
```r
# Only load DE results for selected cluster
load_cluster_specific_de <- function(cluster_id) {
  de_file <- sprintf("de_results_cluster_%s.rds", cluster_id)
  if (file.exists(de_file)) {
    return(readRDS(de_file))
  }
  return(NULL)
}
```

### Phase 3: Advanced Optimizations

#### 3.1 WebGL for Large Scatter Plots
```r
# Use scatterD3 or plotly with webgl
output$umap_interactive <- renderPlotly({
  plot_ly(
    data = values$umap_data,
    x = ~UMAP_1, y = ~UMAP_2,
    color = ~cluster,
    type = 'scattergl',  # WebGL renderer
    mode = 'markers',
    marker = list(size = 2)
  )
})
```

#### 3.2 Implement HDF5 Backend
```r
# Use HDF5 for on-disk data access
library(rhdf5)

# Store large matrices in HDF5
h5createFile("umap_data.h5")
h5write(umap_coords, "umap_data.h5", "coordinates")

# Read only needed portions
read_umap_subset <- function(clusters) {
  h5read("umap_data.h5", "coordinates", 
         index = list(which(metadata$cluster %in% clusters), 1:2))
}
```

#### 3.3 Async Data Loading
```r
# Use promises for non-blocking data loads
library(promises)
library(future)
plan(multisession)

load_data_async <- function(path) {
  future_promise({
    readRDS(path)
  })
}
```

## Implementation Priority

### Immediate (Week 1)
1. ✅ NAMESPACE fixes (COMPLETED)
2. ✅ FDR enrichment data verified (COMPLETED)
3. 🔄 Implement UMAP plot caching
4. 🔄 Add data sampling for preview

### Short-term (Week 2)
1. Progressive loading UI
2. Cluster-specific DE loading
3. Pre-render common views
4. Optimize volcano plots

### Long-term (Month 1)
1. HDF5 backend for 230K+ cells
2. WebGL rendering
3. Async data operations
4. Distributed caching

## Testing Strategy

### Performance Benchmarks
- Target: < 2s UMAP render time
- Target: < 5s initial app load
- Target: < 1s cluster switch

### Test Datasets
1. Small (10K cells) - baseline
2. Medium (50K cells) - optimization check
3. Large (230K cells) - production test

### Monitoring
```r
# Add performance tracking
tictoc::tic()
render_umap()
time_taken <- tictoc::toc()
log_performance("umap_render", time_taken)
```

## Resource Requirements
- **Memory**: 16GB RAM minimum, 32GB recommended
- **CPU**: Multi-core for parallel processing
- **Storage**: SSD for fast data access
- **R Packages**: data.table, fst, arrow (for optimized I/O)

## Next Steps
1. Implement Phase 1 optimizations
2. Benchmark with 230K cell dataset
3. Profile bottlenecks with profvis
4. Iterate based on user feedback

## References
- [Shiny Performance](https://shiny.rstudio.com/articles/performance.html)
- [Large Data in R](https://www.r-bloggers.com/2020/09/how-to-work-with-large-datasets-in-r/)
- [SingleCellExperiment Optimization](https://bioconductor.org/packages/release/bioc/vignettes/SingleCellExperiment/inst/doc/intro.html)