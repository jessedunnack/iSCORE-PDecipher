# Cache Manager for Enrichment Data
# R6 class implementation for efficient caching with TTL and size management

library(R6)

#' CacheManager R6 Class
#' 
#' A cache implementation with time-to-live (TTL) and size management
#' for storing enrichment analysis results
#' 
#' @field cache List storing cached objects
#' @field timestamps List storing cache timestamps
#' @field max_size Maximum number of items to cache
#' @field ttl_minutes Time-to-live in minutes
#' 
#' @examples
#' cache <- CacheManager$new()
#' cache$set("key1", data)
#' cached_data <- cache$get("key1")
CacheManager <- R6::R6Class("CacheManager",
  public = list(
    cache = NULL,
    timestamps = NULL,
    max_size = 10,
    ttl_minutes = 30,
    verbose = FALSE,
    
    #' @description
    #' Initialize a new cache manager
    #' @param max_size Maximum cache size (default: 10)
    #' @param ttl_minutes Time-to-live in minutes (default: 30)
    #' @param verbose Whether to print cache operations (default: FALSE)
    initialize = function(max_size = 10, ttl_minutes = 30, verbose = FALSE) {
      self$cache <- list()
      self$timestamps <- list()
      self$max_size <- max_size
      self$ttl_minutes <- ttl_minutes
      self$verbose <- verbose
      
      if (self$verbose) {
        cat("CacheManager initialized with max_size =", max_size, 
            "and TTL =", ttl_minutes, "minutes\n")
      }
    },
    
    #' @description
    #' Get cached value if valid
    #' @param key Cache key
    #' @return Cached object or NULL
    get = function(key) {
      if (self$is_valid(key)) {
        if (self$verbose) {
          cat("Cache HIT for key:", key, "\n")
        }
        return(self$cache[[key]])
      }
      
      if (self$verbose) {
        cat("Cache MISS for key:", key, "\n")
      }
      return(NULL)
    },
    
    #' @description
    #' Set cache value
    #' @param key Cache key
    #' @param value Object to cache
    set = function(key, value) {
      # Check if we need to evict
      if (length(self$cache) >= self$max_size && !key %in% names(self$cache)) {
        self$evict_oldest()
      }
      
      self$cache[[key]] <- value
      self$timestamps[[key]] <- Sys.time()
      
      if (self$verbose) {
        cat("Cached key:", key, "| Total items:", length(self$cache), "\n")
      }
    },
    
    #' @description
    #' Check if cached value is valid (exists and not expired)
    #' @param key Cache key
    #' @return Logical
    is_valid = function(key) {
      if (!key %in% names(self$cache)) {
        return(FALSE)
      }
      
      age_minutes <- as.numeric(difftime(Sys.time(), self$timestamps[[key]], units = "mins"))
      is_fresh <- age_minutes < self$ttl_minutes
      
      if (!is_fresh) {
        # Remove expired entry
        self$remove(key)
      }
      
      return(is_fresh)
    },
    
    #' @description
    #' Remove cache entry
    #' @param key Cache key
    remove = function(key) {
      if (key %in% names(self$cache)) {
        self$cache[[key]] <- NULL
        self$timestamps[[key]] <- NULL
        
        if (self$verbose) {
          cat("Removed key:", key, "\n")
        }
      }
    },
    
    #' @description
    #' Evict oldest cache entry
    evict_oldest = function() {
      if (length(self$cache) == 0) return()
      
      # Find oldest timestamp
      oldest_time <- min(unlist(self$timestamps))
      oldest_key <- names(self$timestamps)[which(unlist(self$timestamps) == oldest_time)[1]]
      
      if (self$verbose) {
        cat("Evicting oldest key:", oldest_key, "\n")
      }
      
      self$remove(oldest_key)
    },
    
    #' @description
    #' Clear all cache entries
    clear = function() {
      self$cache <- list()
      self$timestamps <- list()
      
      if (self$verbose) {
        cat("Cache cleared\n")
      }
    },
    
    #' @description
    #' Get cache statistics
    #' @return List with cache stats
    stats = function() {
      list(
        size = length(self$cache),
        max_size = self$max_size,
        ttl_minutes = self$ttl_minutes,
        keys = names(self$cache),
        ages = if (length(self$timestamps) > 0) {
          sapply(self$timestamps, function(ts) {
            round(as.numeric(difftime(Sys.time(), ts, units = "mins")), 1)
          })
        } else {
          numeric(0)
        }
      )
    },
    
    #' @description
    #' Print cache information
    print = function() {
      stats <- self$stats()
      cat("CacheManager:\n")
      cat("  Size:", stats$size, "/", stats$max_size, "\n")
      cat("  TTL:", stats$ttl_minutes, "minutes\n")
      if (stats$size > 0) {
        cat("  Keys:", paste(stats$keys, collapse = ", "), "\n")
        cat("  Ages (min):", paste(stats$ages, collapse = ", "), "\n")
      }
    }
  )
)

#' Global cache instance for the Shiny app
#' 
#' This creates a singleton cache instance that can be used
#' throughout the application
#' 
#' @export
create_app_cache <- function(max_size = 20, ttl_minutes = 30, verbose = FALSE) {
  CacheManager$new(max_size = max_size, ttl_minutes = ttl_minutes, verbose = verbose)
}

#' Enhanced load function with caching
#' 
#' Wrapper around load_enrichment_result that uses caching
#' 
#' @param cache CacheManager instance
#' @param analysis_type MAST or MixScale
#' @param gene Gene name
#' @param cluster Cluster ID
#' @param experiment Experiment name
#' @param enrichment_type Enrichment type
#' @param direction UP, DOWN, or ALL
#' @return Enrichment result object or NULL
#' @export
load_enrichment_with_app_cache <- function(cache, analysis_type, gene, cluster, 
                                          experiment, enrichment_type, direction) {
  # Create cache key
  key <- paste(analysis_type, gene, cluster, experiment, 
               enrichment_type, direction, sep = "_")
  
  # Try cache first
  cached_result <- cache$get(key)
  if (!is.null(cached_result)) {
    return(cached_result)
  }
  
  # Load from file
  result <- load_enrichment_result(
    analysis_type = analysis_type,
    gene = gene,
    cluster = cluster,
    experiment = experiment,
    enrichment_type = enrichment_type,
    direction = direction
  )
  
  # Cache if successful
  if (!is.null(result)) {
    cache$set(key, result)
  }
  
  return(result)
}

#' Preload likely next files
#' 
#' Intelligently preload files that are likely to be accessed next
#' 
#' @param cache CacheManager instance
#' @param current_selection Current user selection
#' @param app_data App data containing available options
#' @export
preload_next_files <- function(cache, current_selection, app_data) {
  # Preload adjacent clusters
  if (!is.null(current_selection$cluster)) {
    cluster_num <- as.numeric(gsub("cluster_", "", current_selection$cluster))
    
    # Try next and previous cluster
    for (offset in c(-1, 1)) {
      next_cluster <- paste0("cluster_", cluster_num + offset)
      if (next_cluster %in% app_data$available_clusters) {
        # Preload in background (don't wait for result)
        tryCatch({
          load_enrichment_with_app_cache(
            cache,
            current_selection$analysis_type,
            current_selection$gene,
            next_cluster,
            current_selection$experiment,
            current_selection$enrichment_type,
            current_selection$direction
          )
        }, error = function(e) {
          # Silently ignore preload errors
        })
      }
    }
  }
  
  # Preload other directions for same selection
  all_directions <- c("UP", "DOWN", "ALL")
  other_directions <- setdiff(all_directions, current_selection$direction)
  
  for (dir in other_directions) {
    tryCatch({
      load_enrichment_with_app_cache(
        cache,
        current_selection$analysis_type,
        current_selection$gene,
        current_selection$cluster,
        current_selection$experiment,
        current_selection$enrichment_type,
        dir
      )
    }, error = function(e) {
      # Silently ignore preload errors
    })
  }
}