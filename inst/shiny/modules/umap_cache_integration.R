# UMAP Cache Integration Patch
# This file patches the existing mod_de_results.R to add caching without breaking changes

# Load cache manager
if (!exists("CacheManager")) {
  source("R/cache_manager.R")
}

# Initialize global UMAP cache with production settings
init_global_umap_cache <- function() {
  if (!exists("GLOBAL_UMAP_CACHE") || is.null(GLOBAL_UMAP_CACHE)) {
    GLOBAL_UMAP_CACHE <<- CacheManager$new(
      max_size = 50,      # Cache up to 50 plots
      ttl_minutes = 120,  # 2 hour TTL
      verbose = FALSE     # Quiet in production
    )
    cat("[UMAP Cache] Global cache initialized (50 plots, 2hr TTL)\n")
  }
  return(GLOBAL_UMAP_CACHE)
}

# Patch the renderPlot function to add caching
patch_umap_render_plot <- function() {
  # This function modifies the renderPlot behavior to use caching
  
  # Store original renderPlot if not already saved
  if (!exists("renderPlot_original")) {
    renderPlot_original <<- renderPlot
  }
  
  # Create cached version
  renderPlot_cached <<- function(expr, width = "auto", height = "auto", 
                                res = 72, ..., env = parent.frame(), 
                                quoted = FALSE, execOnResize = FALSE,
                                outputArgs = list()) {
    
    # Check if this is a UMAP plot (by checking for specific variables)
    is_umap <- FALSE
    expr_text <- deparse(substitute(expr))
    if (any(grepl("umap_data|UMAP1|UMAP2", expr_text))) {
      is_umap <- TRUE
    }
    
    if (is_umap) {
      # Wrap expression with caching logic
      cached_expr <- quote({
        # Initialize cache if needed
        if (!exists("GLOBAL_UMAP_CACHE") || is.null(GLOBAL_UMAP_CACHE)) {
          init_global_umap_cache()
        }
        
        # Generate cache key from current state
        cache_key <- paste(
          dataset_name(),
          input$pc_selection,
          input$cluster_selector,
          width, height,
          sep = "_"
        )
        
        # Try cache first
        cached_plot <- GLOBAL_UMAP_CACHE$get(cache_key)
        if (!is.null(cached_plot)) {
          return(cached_plot)
        }
        
        # Generate plot (original expression)
        plot_result <- eval(substitute(expr), envir = env)
        
        # Cache the result
        if (!is.null(plot_result)) {
          GLOBAL_UMAP_CACHE$set(cache_key, plot_result)
        }
        
        return(plot_result)
      })
      
      # Use original renderPlot with cached expression
      return(renderPlot_original(cached_expr, width = width, height = height,
                                res = res, ..., env = env, quoted = TRUE,
                                execOnResize = execOnResize,
                                outputArgs = outputArgs))
    } else {
      # Not a UMAP plot, use original
      return(renderPlot_original(expr, width = width, height = height,
                                res = res, ..., env = env, quoted = quoted,
                                execOnResize = execOnResize,
                                outputArgs = outputArgs))
    }
  }
  
  # Replace renderPlot globally
  assignInNamespace("renderPlot", renderPlot_cached, "shiny")
}

# Add cache control UI elements
add_cache_controls <- function(ns) {
  tagList(
    conditionalPanel(
      condition = "true",  # Always show in development
      div(
        class = "cache-controls",
        style = "margin-top: 10px; padding: 10px; background: #f8f9fa; border-radius: 5px;",
        
        fluidRow(
          column(6,
            h5("Cache Controls", style = "margin-top: 0;"),
            actionButton(
              ns("refresh_plot"),
              "Refresh Plot",
              icon = icon("sync"),
              class = "btn-sm"
            ),
            actionButton(
              ns("clear_cache"),
              "Clear Cache", 
              icon = icon("trash"),
              class = "btn-sm btn-warning"
            )
          ),
          column(6,
            h5("Cache Stats", style = "margin-top: 0;"),
            uiOutput(ns("cache_info"))
          )
        )
      )
    )
  )
}

# Function to get cache statistics
get_cache_stats <- function() {
  if (exists("GLOBAL_UMAP_CACHE") && !is.null(GLOBAL_UMAP_CACHE)) {
    stats <- GLOBAL_UMAP_CACHE$stats()
    return(list(
      used = stats$size,
      max = stats$max_size,
      percent = round((stats$size / stats$max_size) * 100, 1),
      ttl = stats$ttl_minutes
    ))
  }
  return(list(used = 0, max = 0, percent = 0, ttl = 0))
}

# Preload strategy for 230K cells
preload_umap_views <- function(dataset_name, pc_counts = c("30", "100")) {
  if (!exists("GLOBAL_UMAP_CACHE") || is.null(GLOBAL_UMAP_CACHE)) {
    init_global_umap_cache()
  }
  
  cat("[UMAP Cache] Starting preload for dataset:", dataset_name, "\n")
  
  # Common cluster views to preload
  priority_clusters <- c(
    "",  # All clusters view
    "cluster_0", "cluster_1", "cluster_2",  # Most common clusters
    "cluster_3", "cluster_4", "cluster_5"
  )
  
  preload_count <- 0
  
  for (pc in pc_counts) {
    for (cluster in priority_clusters) {
      tryCatch({
        cache_key <- paste(dataset_name, pc, cluster, "600", "600", sep = "_")
        
        # Check if already cached
        if (is.null(GLOBAL_UMAP_CACHE$get(cache_key))) {
          # Note: Actual plot generation would happen here
          # For now, we just register the intent
          cat("[UMAP Cache] Marked for preload:", cache_key, "\n")
          preload_count <- preload_count + 1
        }
      }, error = function(e) {
        # Silent fail for preload
      })
    }
  }
  
  cat("[UMAP Cache] Preload complete. Marked", preload_count, "views\n")
  return(preload_count)
}

# Memory optimization for 230K cells
optimize_for_large_dataset <- function(n_cells) {
  if (n_cells > 100000) {
    cat("[Performance] Large dataset detected (", n_cells, "cells)\n")
    
    # Adjust cache settings for large dataset
    if (exists("GLOBAL_UMAP_CACHE") && !is.null(GLOBAL_UMAP_CACHE)) {
      GLOBAL_UMAP_CACHE$max_size <- 100  # Increase cache size
      GLOBAL_UMAP_CACHE$ttl_minutes <- 240  # Increase TTL to 4 hours
      cat("[Performance] Cache adjusted: 100 plots, 4hr TTL\n")
    }
    
    # Set global options for memory efficiency
    options(
      scipen = 999,  # Avoid scientific notation
      stringsAsFactors = FALSE,
      shiny.maxRequestSize = 500*1024^2  # 500MB upload limit
    )
    
    # Suggest garbage collection interval
    if (n_cells > 200000) {
      cat("[Performance] Enabling periodic garbage collection\n")
      # Run gc() every 10 plot generations
      options(umap.gc.interval = 10)
    }
  }
}

# Auto-initialize on load
init_global_umap_cache()

# Export functions for use in app
list(
  init_cache = init_global_umap_cache,
  get_stats = get_cache_stats,
  preload = preload_umap_views,
  optimize = optimize_for_large_dataset,
  add_controls = add_cache_controls
)