# Enhanced DE Results Module with UMAP Plot Caching
# Optimized for 230K+ cells performance

# Source the original module to inherit base functionality
source("mod_de_results.R")

# Create module-specific cache for UMAP plots
# Increased size and TTL for production use with 230K cells
umap_plot_cache <- NULL

#' Initialize UMAP plot cache
#' 
#' @param max_plots Maximum number of plots to cache (default: 50)
#' @param ttl_hours Time-to-live in hours (default: 2)
#' @export
initialize_umap_cache <- function(max_plots = 50, ttl_hours = 2) {
  # Load cache manager if not already loaded
  if (!exists("CacheManager")) {
    source("R/cache_manager.R")
  }
  
  # Create global cache instance
  umap_plot_cache <<- CacheManager$new(
    max_size = max_plots,
    ttl_minutes = ttl_hours * 60,
    verbose = TRUE
  )
  
  cat("UMAP plot cache initialized: max_plots =", max_plots, 
      ", TTL =", ttl_hours, "hours\n")
}

#' Generate cache key for UMAP plot
#' 
#' @param dataset_name Dataset identifier
#' @param pc_count PC count (30pc, 100pc, etc.)
#' @param cluster_selector Selected cluster or NULL
#' @param plot_width Plot width
#' @param plot_height Plot height
#' @return String cache key
generate_umap_cache_key <- function(dataset_name, pc_count, cluster_selector, 
                                   plot_width = 600, plot_height = 600) {
  # Create deterministic key from plot parameters
  key_parts <- c(
    dataset_name,
    pc_count,
    ifelse(is.null(cluster_selector) || cluster_selector == "", "all", cluster_selector),
    paste0(plot_width, "x", plot_height)
  )
  
  paste(key_parts, collapse = "_")
}

#' Generate UMAP plot with caching
#' 
#' This function wraps the UMAP plot generation with caching logic
#' 
#' @param umap_data Data frame with UMAP coordinates
#' @param cluster_selector Selected cluster or NULL
#' @param cache_key Cache key for this plot
#' @param force_refresh Force regeneration even if cached
#' @return ggplot object
generate_cached_umap_plot <- function(umap_data, cluster_selector = NULL, 
                                    cache_key = NULL, force_refresh = FALSE) {
  
  # Check cache first (unless force refresh)
  if (!force_refresh && !is.null(cache_key) && !is.null(umap_plot_cache)) {
    cached_plot <- umap_plot_cache$get(cache_key)
    if (!is.null(cached_plot)) {
      cat("[UMAP Cache] HIT - Returning cached plot for key:", cache_key, "\n")
      return(cached_plot)
    }
  }
  
  cat("[UMAP Cache] MISS - Generating new plot for key:", cache_key, "\n")
  
  # Performance tracking
  start_time <- Sys.time()
  
  # Get cluster colors with proper ordering
  clusters <- natural_sort_clusters(unique(umap_data$cluster))
  n_clusters <- length(clusters)
  ditto_colors <- get_ditto_colors(n_clusters)
  names(ditto_colors) <- clusters
  
  # Create display data with highlighting and proper factor ordering
  plot_data <- umap_data
  plot_data$cluster <- factor(plot_data$cluster, levels = clusters)
  
  # Handle cluster selection for highlighting
  if (!is.null(cluster_selector) && cluster_selector != "") {
    # Create display categories
    plot_data$display_group <- ifelse(
      plot_data$cluster == cluster_selector,
      plot_data$cluster,
      "Background"
    )
    
    # Set colors - selected cluster keeps its color, others gray
    color_values <- c(ditto_colors[cluster_selector], "Background" = "#E8E8E8")
    
    # Set alpha and size values for highlighting
    plot_data$alpha_value <- ifelse(
      plot_data$cluster == cluster_selector,
      0.8,
      0.15
    )
    
    plot_data$size_value <- ifelse(
      plot_data$cluster == cluster_selector,
      0.5,
      0.3
    )
    
    # Calculate cluster centers for labels
    cluster_centers <- plot_data %>%
      filter(cluster == cluster_selector) %>%
      summarise(
        x = median(UMAP1),
        y = median(UMAP2),
        label = unique(cluster)
      )
    
  } else {
    # Show all clusters in full color
    plot_data$display_group <- plot_data$cluster
    color_values <- ditto_colors
    plot_data$alpha_value <- 0.6
    plot_data$size_value <- 0.4
    
    # Calculate all cluster centers
    cluster_centers <- plot_data %>%
      group_by(cluster) %>%
      summarise(
        x = median(UMAP1),
        y = median(UMAP2),
        label = unique(cluster),
        .groups = 'drop'
      )
  }
  
  # Create the base plot
  p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
    geom_point(aes(color = display_group, alpha = alpha_value, size = size_value)) +
    scale_color_manual(values = color_values) +
    scale_alpha_identity() +
    scale_size_identity() +
    theme_minimal() +
    theme(
      legend.position = if(is.null(cluster_selector) || cluster_selector == "") "right" else "none",
      legend.title = element_text(size = 12, face = "bold"),
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(size = 12),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    guides(color = guide_legend(title = "Cluster", override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMAP 1", y = "UMAP 2")
  
  # Add cluster labels based on selection
  if (!is.null(cluster_selector) && cluster_selector != "") {
    # Only label the selected cluster
    p <- p + 
      geom_label(
        data = cluster_centers,
        aes(x = x, y = y, label = label),
        size = 5,
        fontface = "bold",
        fill = "white",
        alpha = 0.8
      )
  } else if (n_clusters <= 20) {
    # Show all labels with repel to avoid overlap
    p <- p + 
      geom_label_repel(
        data = cluster_centers,
        aes(x = x, y = y, label = gsub("cluster_", "", label)),
        size = 4,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.5,
        segment.color = "gray50",
        max.overlaps = Inf,
        fill = "white",
        alpha = 0.8
      )
  }
  
  # Track performance
  end_time <- Sys.time()
  generation_time <- as.numeric(difftime(end_time, start_time, units = "secs"))
  cat("[UMAP Cache] Plot generated in", round(generation_time, 2), "seconds\n")
  
  # Cache the plot if cache is available
  if (!is.null(cache_key) && !is.null(umap_plot_cache)) {
    umap_plot_cache$set(cache_key, p)
    cat("[UMAP Cache] Plot cached with key:", cache_key, "\n")
  }
  
  return(p)
}

#' Enhanced DE Results Server with UMAP caching
#' 
#' Drop-in replacement for de_results_server with caching
#' 
#' @param id Module ID
#' @param de_data_reactive Reactive DE data
#' @param enrichment_data_reactive Reactive enrichment data
#' @param dataset_name Reactive dataset name
#' @export
de_results_server_cached <- function(id, de_data_reactive, enrichment_data_reactive, 
                                    dataset_name = reactive({"dataset1"})) {
  
  moduleServer(id, function(input, output, session) {
    
    # Initialize cache if not already done
    if (is.null(umap_plot_cache)) {
      initialize_umap_cache()
    }
    
    # Call original server logic
    original_server <- de_results_server(id, de_data_reactive, enrichment_data_reactive, dataset_name)
    
    # Override the UMAP plot rendering with cached version
    output$umap_plot <- renderPlot({
      req(values$umap_data)
      
      # Generate cache key
      cache_key <- generate_umap_cache_key(
        dataset_name = dataset_name(),
        pc_count = input$pc_selection,
        cluster_selector = input$cluster_selector,
        plot_width = 600,
        plot_height = 600
      )
      
      # Generate or retrieve cached plot
      p <- generate_cached_umap_plot(
        umap_data = values$umap_data,
        cluster_selector = input$cluster_selector,
        cache_key = cache_key,
        force_refresh = input$force_refresh %||% FALSE
      )
      
      return(p)
      
    }, height = 600, width = 600)
    
    # Add cache statistics to UI
    output$cache_stats <- renderUI({
      req(umap_plot_cache)
      
      stats <- umap_plot_cache$stats()
      
      div(
        class = "cache-stats",
        style = "padding: 10px; background: #f0f0f0; border-radius: 5px; margin-top: 10px;",
        h5("UMAP Cache Statistics"),
        p(
          strong("Cache Usage: "), 
          sprintf("%d / %d plots", stats$size, stats$max_size)
        ),
        p(
          strong("Hit Rate: "),
          if (exists("cache_hits") && exists("cache_misses")) {
            hit_rate <- cache_hits / (cache_hits + cache_misses) * 100
            sprintf("%.1f%%", hit_rate)
          } else {
            "N/A"
          }
        ),
        p(
          strong("TTL: "),
          sprintf("%d minutes", stats$ttl_minutes)
        ),
        actionButton(
          NS(id, "clear_cache"),
          "Clear Cache",
          class = "btn-sm btn-warning"
        )
      )
    })
    
    # Handle cache clearing
    observeEvent(input$clear_cache, {
      if (!is.null(umap_plot_cache)) {
        umap_plot_cache$clear()
        showNotification("UMAP cache cleared", type = "info")
        
        # Force redraw
        values$cache_cleared <- runif(1)
      }
    })
    
    return(original_server)
  })
}

#' Preload common UMAP views
#' 
#' Pregenerate and cache commonly accessed UMAP plots
#' 
#' @param umap_data_list List of UMAP datasets
#' @param common_views List of common view configurations
#' @export
preload_common_umap_views <- function(umap_data_list, common_views = NULL) {
  
  # Default common views if not provided
  if (is.null(common_views)) {
    common_views <- list(
      list(dataset = "dataset1", pc = "30", cluster = NULL),
      list(dataset = "dataset1", pc = "100", cluster = NULL),
      list(dataset = "dataset2", pc = "30", cluster = NULL),
      list(dataset = "dataset2", pc = "100", cluster = NULL)
    )
  }
  
  cat("[UMAP Cache] Starting preload of", length(common_views), "common views\n")
  
  preload_count <- 0
  for (view in common_views) {
    tryCatch({
      # Get data for this view
      if (view$dataset %in% names(umap_data_list)) {
        umap_data <- umap_data_list[[view$dataset]][[view$pc]]
        
        if (!is.null(umap_data)) {
          # Generate cache key
          cache_key <- generate_umap_cache_key(
            dataset_name = view$dataset,
            pc_count = view$pc,
            cluster_selector = view$cluster
          )
          
          # Generate and cache plot
          generate_cached_umap_plot(
            umap_data = umap_data,
            cluster_selector = view$cluster,
            cache_key = cache_key
          )
          
          preload_count <- preload_count + 1
        }
      }
    }, error = function(e) {
      cat("[UMAP Cache] Error preloading view:", e$message, "\n")
    })
  }
  
  cat("[UMAP Cache] Preloaded", preload_count, "views successfully\n")
  
  return(preload_count)
}

#' Export cached DE Results module UI
#' @export
de_results_ui_cached <- de_results_ui