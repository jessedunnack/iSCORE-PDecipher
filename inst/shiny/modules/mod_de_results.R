# Module for DE Results page with interactive UMAP and volcano plots
# Allows clicking on UMAP clusters to update volcano plots

# Load required packages conditionally
if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
  library(SingleCellExperiment)
}
if (requireNamespace("SummarizedExperiment", quietly = TRUE)) {
  library(SummarizedExperiment)
}
library(ggplot2)
library(dplyr)
library(ggrepel)
library(plotly)  # CRITICAL: Ensure plotly is loaded for volcano plots

# Load EnhancedVolcano for publication-quality static plots
enhancedvolcano_available <- FALSE
if (requireNamespace("EnhancedVolcano", quietly = TRUE)) {
  library(EnhancedVolcano)
  enhancedvolcano_available <- TRUE
} else {
  # Try to install if missing
  tryCatch({
    if (requireNamespace("BiocManager", quietly = TRUE)) {
      BiocManager::install("EnhancedVolcano", quiet = TRUE)
      library(EnhancedVolcano)
      enhancedvolcano_available <- TRUE
    }
  }, error = function(e) {
    warning("EnhancedVolcano package not available. Static volcano plots will show installation message.")
  })
}

# Helper functions to process DE data for volcano plots
process_mast_for_volcano <- function(mast_data) {
  # Convert MAST data structure to volcano plot format
  volcano_data <- data.frame()
  
  for (gene in names(mast_data)) {
    for (cluster in names(mast_data[[gene]])) {
      if (!is.null(mast_data[[gene]][[cluster]]$results)) {
        de_results <- mast_data[[gene]][[cluster]]$results
        
        # Extract log2FC and p-values from MAST results
        if ("avg_log2FC" %in% colnames(de_results) && "p_val_adj" %in% colnames(de_results)) {
          cluster_data <- data.frame(
            gene = gene,
            cluster = cluster,
            gene_name = rownames(de_results),
            log2FC = de_results$avg_log2FC,
            pvalue = de_results$p_val_adj,
            experiment = "default",
            stringsAsFactors = FALSE
          )
          volcano_data <- rbind(volcano_data, cluster_data)
        }
      }
    }
  }
  
  return(volcano_data)
}

process_mixscale_for_volcano_enhanced <- function(mixscale_data, experiment_preference = "C12_FPD-24", 
                                                 use_weighted_meta_analysis = FALSE, 
                                                 experiment_weights = NULL) {
  # Enhanced version that properly handles multiple CRISPRi experiments
  # Fixes critical issue: previous version only used first experiment (arbitrary selection)
  
  volcano_data <- data.frame()
  
  for (gene in names(mixscale_data)) {
    for (cluster in names(mixscale_data[[gene]])) {
      if (!is.null(mixscale_data[[gene]][[cluster]]$results)) {
        de_results <- mixscale_data[[gene]][[cluster]]$results
        
        # Find all available log2FC and p-value columns
        log2fc_cols <- grep("^log2FC_", names(de_results), value = TRUE)
        
        if (length(log2fc_cols) > 0) {
          if (use_weighted_meta_analysis && !is.null(experiment_weights)) {
            # Option 1: Weighted meta-analysis across experiments
            meta_result <- calculate_weighted_volcano_data(
              de_results, log2fc_cols, gene, cluster, experiment_weights
            )
            if (!is.null(meta_result)) {
              volcano_data <- rbind(volcano_data, meta_result)
            }
            
          } else {
            # Option 2: Use preferred experiment or fallback strategy
            selected_experiment <- select_experiment_for_volcano(
              log2fc_cols, experiment_preference
            )
            
            if (!is.null(selected_experiment)) {
              log2fc_col <- paste0("log2FC_", selected_experiment)
              pval_col <- paste0("p_cell_type", selected_experiment, ":weight")
              
              if (log2fc_col %in% names(de_results) && pval_col %in% names(de_results)) {
                cluster_data <- data.frame(
                  gene = gene,
                  cluster = cluster,
                  gene_name = rownames(de_results),
                  log2FC = de_results[[log2fc_col]],
                  pvalue = de_results[[pval_col]],
                  experiment = selected_experiment,
                  experiment_selection_method = ifelse(selected_experiment == experiment_preference, 
                                                     "preferred", "fallback"),
                  total_experiments_available = length(log2fc_cols),
                  stringsAsFactors = FALSE
                )
                volcano_data <- rbind(volcano_data, cluster_data)
              }
            }
          }
        }
      }
    }
  }
  
  # Add metadata about experiment handling
  attr(volcano_data, "experiment_handling") <- list(
    method = ifelse(use_weighted_meta_analysis, "weighted_meta_analysis", "experiment_preference"),
    preference = experiment_preference,
    weighted_analysis = use_weighted_meta_analysis,
    total_experiments_processed = length(unique(volcano_data$experiment))
  )
  
  return(volcano_data)
}

#' Select best available experiment for volcano plot
#'
#' @param log2fc_cols Available log2FC columns
#' @param experiment_preference Preferred experiment (default "C12_FPD-24")
#' @return Selected experiment name or NULL
select_experiment_for_volcano <- function(log2fc_cols, experiment_preference = "C12_FPD-24") {
  
  # Extract experiment names from column names
  available_experiments <- gsub("^log2FC_", "", log2fc_cols)
  
  # Priority order based on user guidance: C12_FPD-24 is strongest
  experiment_priority <- c("C12_FPD-24", "C12_FPD-23", "C18_FPD-23")
  
  # First, try preferred experiment
  if (experiment_preference %in% available_experiments) {
    return(experiment_preference)
  }
  
  # Fall back to priority order
  for (exp in experiment_priority) {
    if (exp %in% available_experiments) {
      return(exp)
    }
  }
  
  # Last resort: use first available
  if (length(available_experiments) > 0) {
    return(available_experiments[1])
  }
  
  return(NULL)
}

#' Calculate weighted meta-analysis volcano data
#'
#' @param de_results DE results data frame
#' @param log2fc_cols Available log2FC columns
#' @param gene Gene name
#' @param cluster Cluster name
#' @param experiment_weights Experiment weights
#' @return Data frame with weighted meta-analysis results
calculate_weighted_volcano_data <- function(de_results, log2fc_cols, gene, cluster, experiment_weights) {
  
  # Extract experiments and their data
  experiments <- gsub("^log2FC_", "", log2fc_cols)
  experiment_data <- list()
  experiment_weights_for_cluster <- list()
  
  for (exp in experiments) {
    log2fc_col <- paste0("log2FC_", exp)
    pval_col <- paste0("p_cell_type", exp, ":weight")
    
    if (log2fc_col %in% names(de_results) && pval_col %in% names(de_results)) {
      # Get weight for this experiment and cluster
      weight_key <- paste0(exp, "_", cluster)
      
      if (!is.null(experiment_weights) && weight_key %in% names(experiment_weights$weights)) {
        exp_weight <- experiment_weights$weights[[weight_key]]
        
        if (exp_weight > 0) {
          experiment_data[[exp]] <- list(
            log2fc = de_results[[log2fc_col]],
            pvalue = de_results[[pval_col]],
            weight = exp_weight
          )
          experiment_weights_for_cluster[[exp]] <- exp_weight
        }
      }
    }
  }
  
  if (length(experiment_data) == 0) {
    return(NULL)
  }
  
  # Calculate weighted averages for each gene
  gene_names <- rownames(de_results)
  weighted_log2fc <- rep(NA, length(gene_names))
  weighted_pvalue <- rep(NA, length(gene_names))
  
  for (i in seq_along(gene_names)) {
    exp_log2fc <- c()
    exp_pvalues <- c()
    exp_weights <- c()
    
    for (exp in names(experiment_data)) {
      if (!is.na(experiment_data[[exp]]$log2fc[i]) && !is.na(experiment_data[[exp]]$pvalue[i])) {
        exp_log2fc <- c(exp_log2fc, experiment_data[[exp]]$log2fc[i])
        exp_pvalues <- c(exp_pvalues, experiment_data[[exp]]$pvalue[i])
        exp_weights <- c(exp_weights, experiment_data[[exp]]$weight)
      }
    }
    
    if (length(exp_weights) > 0) {
      # Weighted average log2FC
      weighted_log2fc[i] <- sum(exp_log2fc * exp_weights) / sum(exp_weights)
      
      # Weighted p-value using Fisher's method
      chi_square_stats <- -2 * log(exp_pvalues)
      weighted_chi_square <- sum(chi_square_stats * exp_weights) / sum(exp_weights)
      weighted_pvalue[i] <- pchisq(weighted_chi_square * 2 * length(exp_pvalues), 
                                  df = 2 * length(exp_pvalues), lower.tail = FALSE)
    }
  }
  
  # Create result data frame
  cluster_data <- data.frame(
    gene = gene,
    cluster = cluster,
    gene_name = gene_names,
    log2FC = weighted_log2fc,
    pvalue = weighted_pvalue,
    experiment = paste(names(experiment_data), collapse = ","),
    experiment_selection_method = "weighted_meta_analysis",
    total_experiments_available = length(experiment_data),
    experiments_included = length(experiment_data),
    stringsAsFactors = FALSE
  )
  
  return(cluster_data)
}

# Legacy function for backward compatibility
process_mixscale_for_volcano <- function(mixscale_data) {
  # Wrapper that uses enhanced version with C12_FPD-24 preference
  # This maintains backward compatibility while fixing the arbitrary selection issue
  
  return(process_mixscale_for_volcano_enhanced(
    mixscale_data = mixscale_data,
    experiment_preference = "C12_FPD-24",  # Use strongest experiment by default
    use_weighted_meta_analysis = FALSE,
    experiment_weights = NULL
  ))
}

# UI function
mod_de_results_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      # Left panel: UMAP with cluster selection
      column(6,
        wellPanel(
          # Header with PC selector
          fluidRow(
            column(9,
              h3("Cell Cluster Analysis", icon("object-group"), style = "margin: 0;")
            ),
            column(3, style = "padding-top: 8px; text-align: right;",
              selectInput(ns("pc_selection"), 
                         label = NULL,
                         choices = c("30 PCs" = "30", "50 PCs" = "50", "100 PCs" = "100"),
                         selected = "30",  # Default to 30 PCs as requested
                         width = "100px")
            )
          ),
          
          # Visual invitation to select cluster
          div(class = "alert alert-info", style = "margin-bottom: 20px;",
            icon("hand-pointer"),
            strong(" Select a cluster to explore its differential expression"),
            br(),
            span("Choose from ", textOutput(ns("n_clusters_text"), inline = TRUE), 
                 " clusters containing ", textOutput(ns("n_cells_text"), inline = TRUE), 
                 " cells", style = "font-size: 0.9em;")
          ),
          
          # Cluster selector dropdown
          div(style = "margin-bottom: 15px;",
            selectInput(ns("cluster_selector"),
                       label = NULL,
                       choices = c("Choose a cluster to analyze..." = ""),
                       selected = "",
                       width = "100%")
          ),
          
          # UMAP plot output
          shinycssloaders::withSpinner(
            plotOutput(ns("umap_plot"), height = "600px"),
            type = 6,
            color = "#3c8dbc"
          ),
          
          # Selected cluster info box
          conditionalPanel(
            condition = "input.cluster_selector != ''",
            ns = ns,
            div(class = "well well-sm", style = "margin-top: 15px; background-color: #f8f9fa;",
              h5("Cluster Information", style = "margin-top: 0;"),
              div(id = ns("cluster_stats"),
                uiOutput(ns("cluster_info"))
              )
            )
          )
        )
      ),
      
      # Right panel: Volcano plots
      column(6,
        wellPanel(
          h3("Volcano Plots", icon("chart-area")),
          
          # Plot type selection
          div(style = "margin-bottom: 15px;",
            radioButtons(ns("plot_type"),
                        "Plot Type:",
                        choices = c("Static" = "static",
                                  "Interactive" = "interactive"),
                        selected = "static",
                        inline = TRUE)
          ),
          
          # Color options (only for interactive plots)
          conditionalPanel(
            condition = "input.plot_type == 'interactive'",
            ns = ns,
            div(style = "margin-bottom: 15px;",
              radioButtons(ns("color_by"),
                          "Color Points By:",
                          choices = c("Significance" = "significance",
                                    "Experiment" = "experiment",
                                    "Gene/Mutation" = "gene"),
                          selected = "significance",
                          inline = TRUE)
            )
          ),
          
          # MAST volcano plot - DYNAMIC: Height adapts based on MixScale availability
          div(style = "margin-bottom: 20px;",
            h4("MAST Results"),
            shinycssloaders::withSpinner(
              uiOutput(ns("mast_volcano_container")),
              type = 6,
              color = "#3c8dbc"
            )
          ),
          
          # MixScale volcano plot - CONDITIONAL: Server-side rendering for reliable hiding
          uiOutput(ns("mixscale_volcano_container"))
        )
      )
    ),
    
    # Bottom panel: Summary statistics, overlapping genes, and cross-cluster analysis
    fluidRow(
      column(12,
        tabsetPanel(
          tabPanel("Summary & Overlaps",
            wellPanel(
              fluidRow(
                # Left half: Compact summary statistics
                column(6,
                  h4("Summary Statistics"),
                  div(id = ns("summary_stats"),
                    uiOutput(ns("stats_content"))
                  )
                ),
                # Right half: Overlapping genes list
                column(6,
                  h4("Overlapping DE Genes"),
                  div(id = ns("overlap_genes"),
                    uiOutput(ns("overlap_content"))
                  )
                )
              )
            )
          ),
          tabPanel("Cross-Cluster Analysis",
            wellPanel(
              fluidRow(
                column(12,
                  h4("Cross-Cluster DE Gene Analysis", icon("layer-group")),
                  p("Identify genes that are differentially expressed across multiple clusters for the selected gene/mutation.", 
                    style = "color: #666; margin-bottom: 20px;"),
                  
                  # Controls for cross-cluster analysis
                  fluidRow(
                    column(4,
                      selectInput(ns("cross_cluster_gene"),
                                 "Select Gene/Mutation:",
                                 choices = c("Choose a gene..." = ""),
                                 selected = "")
                    ),
                    column(4,
                      numericInput(ns("min_clusters"),
                                  "Min. Clusters:",
                                  value = 2,
                                  min = 1,
                                  max = 20)
                    ),
                    column(4,
                      numericInput(ns("pvalue_threshold"),
                                  "P-value Threshold:",
                                  value = 0.05,
                                  min = 0.001,
                                  max = 0.1,
                                  step = 0.01)
                    )
                  ),
                  
                  # Cross-cluster results table
                  div(style = "margin-top: 20px;",
                    h5("Genes DE in Multiple Clusters"),
                    DT::dataTableOutput(ns("cross_cluster_table"))
                  ),
                  
                  # Cross-cluster heatmap
                  div(style = "margin-top: 20px;",
                    h5("Cross-Cluster DE Status Heatmap"),
                    shinycssloaders::withSpinner(
                      plotOutput(ns("cross_cluster_heatmap"), height = "500px"),
                      type = 6,
                      color = "#3c8dbc"
                    )
                  )
                )
              )
            )
          ),
          tabPanel("DE Gene Overlap Heatmap",
            wellPanel(
              fluidRow(
                column(12,
                  h4("Fisher's Test Heatmap - All DE Genes", icon("th")),
                  p("Compare overlapping DE genes between different mutations/perturbations using Fisher's exact test.", 
                    style = "color: #666; margin-bottom: 20px;"),
                  
                  # Controls for DE gene heatmap
                  fluidRow(
                    column(3,
                      selectInput(ns("heatmap_cluster"),
                                 "Select Cluster:",
                                 choices = c("Choose a cluster..." = ""),
                                 selected = "")
                    ),
                    column(3,
                      selectInput(ns("heatmap_method"),
                                 "Analysis Method:",
                                 choices = c("Both MAST & MixScale" = "both",
                                           "MAST only" = "mast",
                                           "MixScale only" = "mixscale"),
                                 selected = "both")
                    ),
                    column(3,
                      selectInput(ns("heatmap_direction"),
                                 "Direction Filter:",
                                 choices = c(
                                   "All directions (deduplicated)" = "ALL",
                                   "Up-regulated only" = "UP",
                                   "Down-regulated only" = "DOWN"
                                 ),
                                 selected = "ALL")
                    ),
                    column(3,
                      numericInput(ns("min_genes_overlap"),
                                  "Min. Genes for Test:",
                                  value = 5,
                                  min = 1,
                                  max = 50)
                    )
                  ),
                  
                  # Direction explanation
                  div(class = "alert alert-info", style = "margin-top: 10px; font-size: 0.85em;",
                    icon("info-circle"),
                    tags$strong(" Direction Filtering: "), 
                    "Controls which gene expression direction to analyze. 'All directions' prevents counting genes multiple times across UP/DOWN tests, ",
                    "ensuring accurate Fisher's exact test calculations (similar to the fix applied to gene-pair analysis)."
                  ),
                  
                  # DE gene overlap heatmap
                  div(style = "margin-top: 20px;",
                    h5("Gene Overlap Significance Matrix"),
                    shinycssloaders::withSpinner(
                      plotOutput(ns("de_overlap_heatmap"), height = "600px"),
                      type = 6,
                      color = "#3c8dbc"
                    )
                  ),
                  
                  # Download button for heatmap data
                  div(style = "margin-top: 15px; text-align: center;",
                    downloadButton(ns("download_overlap_matrix"), 
                                  "Download Overlap Matrix", 
                                  class = "btn btn-primary")
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}

# Server function
mod_de_results_server <- function(id, global_selection, app_data) {
  moduleServer(id, function(input, output, session) {
    
    cat("[DE Results] Module server starting...\n")
    
    # Reactive values
    values <- reactiveValues(
      selected_cluster = NULL,
      de_data_mast = NULL,
      de_data_mixscale = NULL,
      umap_data = NULL,
      sce_list = NULL
    )
    
    # Store dataset name for PC switching
    dataset_name <- reactiveVal(NULL)
    
    # Function to load UMAP data for specific PC count (from Overview page)
    load_umap_data <- function(dataset_name, pc_count) {
      # Construct filename with PC suffix
      filename <- paste0(dataset_name, "_umap_data_", pc_count, "pc.rds")
      
      # Try multiple possible paths
      possible_paths <- c(
        system.file("extdata", "umap_data", filename, package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", filename),
        paste0("../../inst/extdata/umap_data/", filename)
      )
      
      # Also check for legacy filename (without pc suffix) for any PC count
      # This ensures compatibility with older UMAP files that don't have PC suffixes
      legacy_filename <- paste0(dataset_name, "_umap_data.rds")
      possible_paths <- c(possible_paths,
        system.file("extdata", "umap_data", legacy_filename, package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", legacy_filename),
        paste0("../../inst/extdata/umap_data/", legacy_filename)
      )
      
      for (path in possible_paths) {
        if (file.exists(path)) {
          tryCatch({
            sce <- readRDS(path)
            cat("[DE Results] Successfully loaded UMAP (", pc_count, " PCs) from:", path, "\n")
            
            # Extract UMAP coordinates
            if (!is.null(sce)) {
              # Check if SingleCellExperiment functions are available
              if (exists("reducedDim") && exists("colData")) {
                umap_coords <- reducedDim(sce, "UMAP")
                cluster_data <- colData(sce)$seurat_clusters
              } else {
                # Try direct access as a fallback
                cat("[DE Results] SingleCellExperiment functions not available, trying direct access\n")
                # For SCE objects, try accessing slots directly
                if (!is.null(sce@int_colData@listData$reducedDims$UMAP)) {
                  umap_coords <- sce@int_colData@listData$reducedDims$UMAP
                  cluster_data <- sce@colData$seurat_clusters
                } else {
                  stop("Cannot access UMAP data without SingleCellExperiment package")
                }
              }
              
              # Create data frame for plotting
              values$umap_data <- data.frame(
                UMAP1 = umap_coords[, 1],
                UMAP2 = umap_coords[, 2],
                cluster = as.character(cluster_data),
                stringsAsFactors = FALSE
              )
              
              cat("[DE Results] UMAP data extracted:", nrow(values$umap_data), "cells\n")
              return(TRUE)
            }
          }, error = function(e) {
            cat("[DE Results] Failed to load from", path, ":", e$message, "\n")
          })
        }
      }
      
      # If we couldn't load the requested PC count
      showNotification(
        paste0("UMAP with ", pc_count, " PCs not available. Please run generate_multipc_umaps.R"),
        type = "warning",
        duration = 5
      )
      return(FALSE)
    }
    
    # Load initial UMAP data
    observe({
      cat("[DE Results] UMAP data observe block triggered\n")
      cat("[DE Results] app_data$data_loaded =", app_data$data_loaded, "\n")
      
      req(app_data$data_loaded)
      cat("[DE Results] app_data is loaded, attempting to load UMAP data...\n")
      
      # Determine which dataset to load based on app data
      has_crispri <- any(grepl("MixScale", app_data$consolidated_data$method))
      has_mutations <- any(grepl("MAST", app_data$consolidated_data$method))
      has_crispa <- any(grepl("CRISPRa", app_data$consolidated_data$method))
      
      cat("[DE Results] MAST:", has_mutations, "CRISPRi:", has_crispri, "CRISPRa:", has_crispa, "\n")
      
      # Fixed logic: distinguish dataset 2 from dataset 3 using CRISPRa
      if (has_crispri && has_mutations && has_crispa) {
        dataset_to_load <- "Full_Dataset"
      } else if (has_crispri && has_mutations) {
        dataset_to_load <- "iSCORE_PD_CRISPRi"
      } else if (has_crispri) {
        dataset_to_load <- "iSCORE_PD_CRISPRi"
      } else {
        dataset_to_load <- "iSCORE_PD"
      }
      
      cat("[DE Results] Determined dataset to load:", dataset_to_load, "\n")
      dataset_name(dataset_to_load)
      
      # Load default 30 PC UMAP first (as requested by user)
      if (!load_umap_data(dataset_to_load, "30")) {
        # If 30 PC not available, try 100 PC
        if (!load_umap_data(dataset_to_load, "100")) {
          showNotification("UMAP data not found. Please run extract_umap_data.R first.", 
                         type = "error",
                         duration = NULL)
        }
      }
      
      cat("[DE Results] UMAP data population complete. values$umap_data is", 
          ifelse(is.null(values$umap_data), "NULL", "populated"), "\n")
    })
    
    # React to PC selection changes - CLEANED: Removed cat() that could interfere
    observeEvent(input$pc_selection, {
      req(dataset_name())
      
      # Load new UMAP data
      if (load_umap_data(dataset_name(), input$pc_selection)) {
        # Force plot redraw by updating cluster choices
        if (!is.null(values$umap_data)) {
          clusters <- natural_sort_clusters(unique(values$umap_data$cluster))
          cluster_choices <- setNames(clusters, paste("Cluster", gsub("cluster_", "", clusters)))
          
          updateSelectInput(session, "cluster_selector",
                           choices = c("Choose a cluster to analyze..." = "", cluster_choices))
        }
      }
    })
    
    # Track DE loading status
    de_data_loaded <- reactiveVal(FALSE)
    
    # Function to load DE data on demand
    load_de_data <- function() {
      if (!de_data_loaded()) {
        cat("[DE Results] Loading DE data on demand...\n")
        
        # Show loading notification
        showNotification("Loading differential expression data...", type = "message", duration = 3)
        
        # Get dataset directory
        data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
        
        # Look for full_DE_results.rds in the dataset directory
        possible_de_paths <- c(
          file.path(data_dir, "full_DE_results.rds"),
          Sys.getenv("ISCORE_DE_FILE", ""),
          file.path(dirname(Sys.getenv("ISCORE_ENRICHMENT_DIR", "")), "full_DE_results.rds")
        )
        
        # Remove empty paths
        possible_de_paths <- possible_de_paths[possible_de_paths != ""]
        
        de_loaded <- FALSE
        for (path in possible_de_paths) {
          cat("[DE Results] Checking DE path:", path, "\n")
          if (file.exists(path)) {
            tryCatch({
              de_results <- readRDS(path)
              cat("[DE Results] Successfully loaded DE results from:", path, "\n")
              cat("[DE Results] DE results structure:", paste(names(de_results), collapse=", "), "\n")
              
              # Extract MAST and MixScale data
              if ("iSCORE_PD_MAST" %in% names(de_results)) {
                # Convert MAST data to volcano plot format
                mast_data <- de_results$iSCORE_PD_MAST
                values$de_data_mast <- process_mast_for_volcano(mast_data)
                cat("[DE Results] Processed MAST data:", nrow(values$de_data_mast), "rows\n")
                cat("[DE Results] Available MAST genes:", paste(unique(values$de_data_mast$gene), collapse=", "), "\n")
              }
              
              if ("CRISPRi_Mixscale" %in% names(de_results)) {
                # Convert MixScale data to volcano plot format  
                mixscale_data <- de_results$CRISPRi_Mixscale
                values$de_data_mixscale <- process_mixscale_for_volcano(mixscale_data)
                cat("[DE Results] Processed MixScale data:", nrow(values$de_data_mixscale), "rows\n")
                cat("[DE Results] Available MixScale genes:", paste(unique(values$de_data_mixscale$gene), collapse=", "), "\n")
              }
              
              de_loaded <- TRUE
              de_data_loaded(TRUE)
              showNotification("DE data loaded successfully!", type = "message", duration = 2)
              break
            }, error = function(e) {
              cat("[DE Results] Failed to load DE results from", path, ":", e$message, "\n")
              showNotification(paste("Error loading DE data:", e$message), type = "error", duration = 5)
            })
          }
        }
        
        if (!de_loaded) {
          cat("[DE Results] No DE results loaded\n")
          showNotification("No DE results found. Please ensure full_DE_results.rds is available.", type = "warning", duration = 5)
        }
      }
    }
    
    # Track if we're updating to prevent circular updates
    local_updating <- reactiveVal(FALSE)
    
    # Initialize from global settings
    observe({
      if (!local_updating() && !is.null(global_selection()$cluster)) {
        # Only update if different from current selection
        if (is.null(values$selected_cluster) || values$selected_cluster != global_selection()$cluster) {
          values$selected_cluster <- global_selection()$cluster
          # Update the cluster selector dropdown
          updateSelectInput(session, "cluster_selector", selected = values$selected_cluster)
        }
      }
    })
    
    # DYNAMIC LAYOUT: Detect MixScale availability and adjust MAST plot height
    output$has_mixscale_data <- reactive({
      !is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0
    })
    outputOptions(output, "has_mixscale_data", suspendWhenHidden = FALSE)
    
    # Dynamic MAST volcano plot container with adaptive height and plot type support
    output$mast_volcano_container <- renderUI({
      ns <- session$ns
      has_mixscale <- !is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0
      
      # If no MixScale data, expand MAST plot to full height (700px)
      # If MixScale data available, use standard height (350px)
      plot_height <- if (has_mixscale) "350px" else "700px"
      
      div(style = "width: 100%;",  # Consistent container styling
        div(style = paste0("width: 100%; height: ", plot_height, ";"),  # Explicit sizing container
          if (input$plot_type == "static") {
            plotOutput(ns("mast_volcano_static"), height = plot_height, width = "100%")
          } else {
            plotlyOutput(ns("mast_volcano_interactive"), height = plot_height, width = "100%")
          }
        )
      )
    })
    
    # Dynamic MixScale volcano plot container with plot type support - ONLY renders if data available
    output$mixscale_volcano_container <- renderUI({
      ns <- session$ns
      has_mixscale <- !is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0
      
      # Only render MixScale plot if data is available
      if (has_mixscale) {
        div(style = "width: 100%; margin-top: 10px;",
          h4("MixScale Results"),
          div(style = "width: 100%; height: 350px;",  # Explicit sizing container
            shinycssloaders::withSpinner(
              if (input$plot_type == "static") {
                plotOutput(ns("mixscale_volcano_static"), height = "350px", width = "100%")
              } else {
                plotlyOutput(ns("mixscale_volcano_interactive"), height = "350px", width = "100%")
              },
              type = 6,
              color = "#3c8dbc"
            )
          )
        )
      } else {
        # Return NULL to hide the entire section
        NULL
      }
    })

    # Get dittoSeq colors for consistency
    get_ditto_colors <- function(n_colors) {
      if (requireNamespace("dittoSeq", quietly = TRUE)) {
        # Use dittoSeq's color palette
        colors <- dittoSeq::dittoColors()[1:n_colors]
      } else {
        # Fallback to a similar colorblind-friendly palette
        colors <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
                   "#D55E00", "#CC79A7", "#999999", "#000000", "#E7298A",
                   "#66A61E", "#E6AB02", "#A6761D", "#666666", "#7570B3")
        colors <- colors[1:n_colors]
      }
      return(colors)
    }
    
    # Update cluster choices when data is loaded
    observe({
      req(values$umap_data)
      
      clusters <- natural_sort_clusters(unique(values$umap_data$cluster))
      cluster_choices <- setNames(clusters, paste("Cluster", gsub("cluster_", "", clusters)))
      
      updateSelectInput(session, "cluster_selector",
                       choices = c("Choose a cluster to analyze..." = "", cluster_choices))
      
      # Update info text
      output$n_clusters_text <- renderText({ length(clusters) })
      output$n_cells_text <- renderText({ format(nrow(values$umap_data), big.mark = ",") })
    })
    
    # Update selected cluster when dropdown changes
    observeEvent(input$cluster_selector, {
      if (input$cluster_selector != "") {
        # Set flag to prevent circular updates
        local_updating(TRUE)
        
        # Update local state
        values$selected_cluster <- input$cluster_selector
        
        # Send update to global settings
        session$sendInputMessage("update_cluster_from_module", 
                               list(value = input$cluster_selector))
        
        cat("[DE Results] Sent cluster update to global:", input$cluster_selector, "\n")
        
        # Reset flag
        local_updating(FALSE)
      } else {
        values$selected_cluster <- NULL
      }
    })
    
    # Render UMAP plot with ggplot2
    output$umap_plot <- renderPlot({
      req(values$umap_data)
      
      # Get cluster colors with proper ordering
      clusters <- natural_sort_clusters(unique(values$umap_data$cluster))
      n_clusters <- length(clusters)
      ditto_colors <- get_ditto_colors(n_clusters)
      names(ditto_colors) <- clusters
      
      # Create display data with highlighting and proper factor ordering
      plot_data <- values$umap_data
      plot_data$cluster <- factor(plot_data$cluster, levels = clusters)
      
      if (!is.null(input$cluster_selector) && input$cluster_selector != "") {
        # Create display categories
        plot_data$display_group <- ifelse(
          plot_data$cluster == input$cluster_selector,
          plot_data$cluster,
          "Background"
        )
        
        # Set colors - selected cluster keeps its color, others gray
        color_values <- c(ditto_colors[input$cluster_selector], "Background" = "#E8E8E8")
        
        # Set alpha values
        plot_data$alpha_value <- ifelse(
          plot_data$cluster == input$cluster_selector,
          0.8,
          0.15
        )
        
        # Set point sizes
        plot_data$size_value <- ifelse(
          plot_data$cluster == input$cluster_selector,
          0.5,
          0.3
        )
        
        # Calculate cluster centers for labels
        cluster_centers <- plot_data %>%
          filter(cluster == input$cluster_selector) %>%
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
      
      # Create the plot
      p <- ggplot(plot_data, aes(x = UMAP1, y = UMAP2)) +
        geom_point(aes(color = display_group, alpha = alpha_value, size = size_value)) +
        scale_color_manual(values = color_values) +
        scale_alpha_identity() +
        scale_size_identity() +
        theme_minimal() +
        theme(
          legend.position = if(is.null(input$cluster_selector) || input$cluster_selector == "") "right" else "none",
          legend.title = element_text(size = 12, face = "bold"),
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          plot.background = element_rect(fill = "white", color = NA)
        ) +
        guides(color = guide_legend(title = "Cluster", override.aes = list(size = 3, alpha = 1))) +
        labs(x = "UMAP 1", y = "UMAP 2")
      
      # Add cluster labels
      if (!is.null(input$cluster_selector) && input$cluster_selector != "") {
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
      
      p
    }, height = 600, width = 600)
    
    # Render cluster information with cell type breakdown and top markers
    output$cluster_info <- renderUI({
      # Capture namespace function at the start
      ns <- session$ns
      
      req(input$cluster_selector)
      req(values$umap_data)
      
      cluster_cells <- values$umap_data %>%
        filter(cluster == input$cluster_selector)
      
      n_cells <- nrow(cluster_cells)
      pct_cells <- round(100 * n_cells / nrow(values$umap_data), 1)
      
      # Enhanced cell type breakdown
      cell_type_info <- tagList()
      
      # Try to detect cell types based on metadata columns
      if ("mutation_tidy" %in% colnames(cluster_cells)) {
        # MAST data - mutation vs eWT
        mut_cells <- cluster_cells %>% filter(mutation_tidy != "eWT")
        ewt_cells <- cluster_cells %>% filter(mutation_tidy == "eWT")
        
        cell_type_info <- tagList(
          tags$div(style = "margin-bottom: 8px;",
            tags$div(style = "display: flex; justify-content: space-between;",
              tags$div(
                tags$strong("Mutant iSCORE-PD: "),
                tags$span(format(nrow(mut_cells), big.mark = ","), 
                         paste0(" (", round(100 * nrow(mut_cells) / n_cells, 1), "%)"))
              ),
              tags$div(
                tags$strong("Control (eWT): "),
                tags$span(format(nrow(ewt_cells), big.mark = ","),
                         paste0(" (", round(100 * nrow(ewt_cells) / n_cells, 1), "%)"))
              )
            )
          )
        )
      } else if ("scMAGeCK_gene_assignment" %in% colnames(cluster_cells)) {
        # CRISPRi data - perturbation vs Non-Targeting
        pert_cells <- cluster_cells %>% filter(scMAGeCK_gene_assignment != "Non-Targeting")
        ctrl_cells <- cluster_cells %>% filter(scMAGeCK_gene_assignment == "Non-Targeting")
        
        cell_type_info <- tagList(
          tags$div(style = "margin-bottom: 8px;",
            tags$div(style = "display: flex; justify-content: space-between;",
              tags$div(
                tags$strong("CRISPRi PerturbSeq: "),
                tags$span(format(nrow(pert_cells), big.mark = ","), 
                         paste0(" (", round(100 * nrow(pert_cells) / n_cells, 1), "%)"))
              ),
              tags$div(
                tags$strong("Control (ntCTRL): "),
                tags$span(format(nrow(ctrl_cells), big.mark = ","),
                         paste0(" (", round(100 * nrow(ctrl_cells) / n_cells, 1), "%)"))
              )
            )
          )
        )
      }
      
      tagList(
        # Total cell count
        tags$div(style = "margin-bottom: 8px;",
          tags$strong("Total Cells in ", input$cluster_selector, ": "),
          tags$span(format(n_cells, big.mark = ","), 
                   paste0(" (", pct_cells, "% of all cells)"))
        ),
        
        # Cell type breakdown
        cell_type_info,
        
        # Top Cluster Markers section (replicated from Overview page)
        tags$div(style = "margin-top: 15px;",
          tags$h5(style = "margin-bottom: 10px; color: #3c8dbc;",
            icon("dna"), " Top Marker Genes for ", input$cluster_selector
          ),
          tags$div(style = "font-size: 11px; color: #666; margin-bottom: 8px;",
            tags$strong("Method: "), "MAST test, LFC≥0.25, min.pct=0.1, padj<0.05"
          ),
          div(style = "height: 300px; overflow-y: auto; border: 1px solid #ddd; border-radius: 4px;",
            withSpinner(
              DT::dataTableOutput(ns("cluster_markers_table"), height = "285px"),
              type = 1,
              color = "#3c8dbc",
              size = 0.5
            )
          )
        )
      )
    })
    
    # Render cluster markers table (replicated from Overview page)
    output$cluster_markers_table <- DT::renderDataTable({
      req(input$cluster_selector)
      req(values$umap_data_obj$markers)
      
      # Filter markers for selected cluster
      cluster_markers <- values$umap_data_obj$markers %>%
        filter(cluster == input$cluster_selector) %>%
        arrange(desc(avg_log2FC)) %>%
        head(20) %>%  # Show top 20 markers to fit in smaller space
        select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
        mutate(
          avg_log2FC = round(avg_log2FC, 3),
          p_val_adj = formatC(p_val_adj, format = "e", digits = 2),
          pct.1 = round(pct.1, 3),
          pct.2 = round(pct.2, 3)
        )
      
      # Create DataTable with compact settings
      DT::datatable(
        cluster_markers,
        options = list(
          pageLength = 10,  # Show fewer rows for compact space
          scrollY = "260px",
          scrollCollapse = TRUE,
          dom = 't',  # Only show table (no search/pagination)
          autoWidth = FALSE,
          columnDefs = list(
            list(width = '80px', targets = 0),  # Gene column
            list(width = '60px', targets = 1),  # Log2FC
            list(width = '70px', targets = 2),  # P-val
            list(width = '50px', targets = 3),  # % in cluster
            list(width = '50px', targets = 4),  # % in other
            list(className = 'dt-center', targets = 1:4)
          )
        ),
        rownames = FALSE,
        colnames = c('Gene', 'Log2FC', 'P-adj', '% this', '% other')
      ) %>%
        DT::formatStyle(
          'avg_log2FC',
          background = DT::styleColorBar(cluster_markers$avg_log2FC, 'lightblue'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })
    
    # Helper function for user-friendly p-value formatting
    format_pvalue <- function(p_val) {
      if (is.na(p_val)) return("NA")
      if (p_val >= 0.001) {
        # Use decimal format for values >= 0.001
        return(sprintf("%.3f", p_val))
      } else {
        # Use scientific notation for very small values < 0.001
        return(format(p_val, digits = 2, scientific = TRUE))
      }
    }
    
    # Generate volcano plot function - BULLETPROOF: Always returns valid plotly object
    generate_volcano_plot <- function(de_data, analysis_type, selected_cluster, color_by, current_gene = NULL, experiment_info = NULL) {
      # Wrap entire function in tryCatch to ensure we always return a plotly object
      tryCatch({
        
        if (is.null(de_data) || nrow(de_data) == 0) {
          # Create empty plot with message
          p <- plot_ly() %>%
            layout(
              title = paste(analysis_type, "- No DE data available"),
              xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
              yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
              annotations = list(
                x = 0,
                y = 5,
                text = "No differential expression data loaded",
                showarrow = FALSE,
                font = list(size = 16, color = "gray")
              )
            )
          return(p)
        }
      
      # Filter by cluster if selected - CLEANED: Removed interfering cat() statements
      if (!is.null(selected_cluster) && selected_cluster != "All") {
        plot_data <- de_data[de_data$cluster == selected_cluster, ]
        cluster_label <- gsub("cluster_", "", selected_cluster)
      } else {
        plot_data <- de_data
        cluster_label <- "All Clusters"
      }
      
      # Generate descriptive title with multi-line format for clarity
      if (analysis_type == "MAST") {
        if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
          plot_title <- paste0(current_gene, " mutation vs isogenic eWT controls\n", 
                              "MAST Analysis\n",
                              "Cluster ", cluster_label)
        } else {
          plot_title <- paste0("MAST mutation analysis vs isogenic eWT controls\n",
                              "All Genes\n", 
                              "Cluster ", cluster_label)
        }
      } else if (analysis_type == "MixScale") {
        if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
          if (!is.null(experiment_info) && experiment_info != "" && experiment_info != "default") {
            plot_title <- paste0(current_gene, " CRISPRi knockdown vs Non-Targeting\n",
                                "MixScale Analysis (", experiment_info, ")\n",
                                "Cluster ", cluster_label)
          } else {
            plot_title <- paste0(current_gene, " CRISPRi knockdown vs Non-Targeting controls\n",
                                "MixScale Analysis\n",
                                "Cluster ", cluster_label)
          }
        } else {
          if (!is.null(experiment_info) && experiment_info != "" && experiment_info != "default") {
            plot_title <- paste0("CRISPRi knockdown vs Non-Targeting\n",
                                "MixScale Analysis (", experiment_info, ")\n", 
                                "Cluster ", cluster_label)
          } else {
            plot_title <- paste0("CRISPRi knockdown vs Non-Targeting controls\n",
                                "MixScale Analysis\n",
                                "Cluster ", cluster_label)
          }
        }
      } else {
        plot_title <- paste(analysis_type, "\nCluster", cluster_label)
      }
      
      if (nrow(plot_data) == 0) {
        # No data for this cluster - create informative plot
        
        p <- plot_ly() %>%
          layout(
            title = paste(analysis_type, title_suffix, "- No data"),
            xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
            yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
            annotations = list(
              x = 0,
              y = 5,
              text = paste("No data for", selected_cluster, "\nAvailable clusters:", 
                          paste(unique(de_data$cluster), collapse=", ")),
              showarrow = FALSE,
              font = list(size = 14, color = "gray")
            )
          )
        return(p)
      }
      
      # Calculate -log10 p-value
      plot_data$negLog10p <- -log10(plot_data$pvalue + 1e-300)  # Add small value to avoid log(0)
      
      # CRITICAL FIX: Ensure gene_name contains actual gene names, not row numbers
      if ("gene_name" %in% names(plot_data)) {
        # Check if gene_name column contains actual names or just numbers
        if (all(grepl("^[0-9]+$", plot_data$gene_name[1:min(10, nrow(plot_data))]))) {
          # If gene_name contains only numbers, try to use rownames of original data
          cat("[VOLCANO DEBUG] gene_name appears to be row numbers, attempting to fix\n")
          # If we have proper rownames, use them; otherwise create meaningful names
          if (!is.null(rownames(plot_data)) && !all(rownames(plot_data) == seq_len(nrow(plot_data)))) {
            plot_data$gene_name <- rownames(plot_data)
          } else {
            # Create gene names from available information
            plot_data$gene_name <- paste0("Gene_", seq_len(nrow(plot_data)))
          }
        }
      } else {
        plot_data$gene_name <- paste0("Gene_", seq_len(nrow(plot_data)))
      }
      
      # Determine significance
      plot_data$significant <- plot_data$pvalue < 0.05 & abs(plot_data$log2FC) > 0.25
      
      # Color based on selection
      if (color_by == "significance") {
        plot_data$color_group <- ifelse(
          !plot_data$significant, "Not significant",
          ifelse(plot_data$log2FC > 0, "Up-regulated", "Down-regulated")
        )
        color_scale <- c(
          "Not significant" = "#CCCCCC",
          "Up-regulated" = "#FF6B6B",
          "Down-regulated" = "#4ECDC4"
        )
      } else if (color_by == "experiment") {
        if ("experiment" %in% names(plot_data)) {
          plot_data$color_group <- plot_data$experiment
          color_scale <- NULL # Use default plotly colors
        } else {
          plot_data$color_group <- ifelse(plot_data$significant, "Significant", "Not Significant")
          color_scale <- c("Significant" = "#FF6B6B", "Not Significant" = "#CCCCCC")
        }
      } else {  # color by gene/mutation
        if ("gene" %in% names(plot_data)) {
          plot_data$color_group <- plot_data$gene
          color_scale <- NULL # Use default plotly colors
        } else {
          plot_data$color_group <- ifelse(plot_data$significant, "Significant", "Not Significant")
          color_scale <- c("Significant" = "#FF6B6B", "Not Significant" = "#CCCCCC")
        }
      }
      
      # Create volcano plot
      p <- plot_ly(
        data = plot_data,
        x = ~log2FC,
        y = ~negLog10p,
        color = ~color_group,
        colors = color_scale,
        type = 'scatter',
        mode = 'markers',
        marker = list(size = 5, opacity = 0.7),
        text = ~paste("Gene:", gene_name,
                     "<br>Log2FC:", round(log2FC, 3),
                     "<br>P-value:", format(pvalue, digits = 3),
                     "<br>Experiment:", experiment,
                     "<br>Mutation/Perturbation:", gene),
        hoverinfo = "text"
      ) %>%
        layout(
          title = list(text = plot_title, 
                      font = list(size = 14)),
          xaxis = list(title = "Log2 Fold Change", zeroline = TRUE),
          yaxis = list(title = "-Log10 P-value", zeroline = FALSE),
          showlegend = TRUE,
          legend = list(
            orientation = "h", 
            x = 0.5, 
            y = -0.2,  # Moved slightly lower to ensure consistent spacing
            xanchor = "center",
            bgcolor = "rgba(255, 255, 255, 0.9)",
            bordercolor = "rgba(0, 0, 0, 0.1)",
            borderwidth = 1,
            font = list(size = 9),  # Slightly smaller font for compactness
            itemsizing = "constant"  # Consistent legend item sizes
          ),
          margin = list(b = 100, l = 50, r = 50, t = 60),  # Consistent margins for all plots
          autosize = TRUE,  # Enable responsive sizing
          height = 350  # Fixed height to match container
        )
      
      # Add threshold lines using layout shapes instead of add_trace to avoid vector length issues
      p <- p %>%
        layout(
          shapes = list(
            # Vertical line at log2FC = -1
            list(
              type = "line",
              x0 = -1, x1 = -1,
              y0 = 0, y1 = max(plot_data$negLog10p, na.rm = TRUE),
              line = list(color = "gray", dash = "dash", width = 1)
            ),
            # Vertical line at log2FC = 1
            list(
              type = "line", 
              x0 = 1, x1 = 1,
              y0 = 0, y1 = max(plot_data$negLog10p, na.rm = TRUE),
              line = list(color = "gray", dash = "dash", width = 1)
            ),
            # Horizontal line at p = 0.05
            list(
              type = "line",
              x0 = min(plot_data$log2FC, na.rm = TRUE), 
              x1 = max(plot_data$log2FC, na.rm = TRUE),
              y0 = -log10(0.05), y1 = -log10(0.05),
              line = list(color = "gray", dash = "dash", width = 1)
            )
          )
        )
      
      return(p)
      
      }, error = function(e) {
        # BULLETPROOF: If any error occurs, return a simple error plot
        plot_ly() %>%
          layout(
            title = paste(analysis_type, "- Error occurred"),
            xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
            yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
            annotations = list(
              x = 0,
              y = 5,
              text = paste("Error in plot generation:", e$message),
              showarrow = FALSE,
              font = list(size = 14, color = "red")
            )
          )
      })
    }
    
    # Generate static publication-quality volcano plot using EnhancedVolcano
    generate_static_volcano_plot <- function(de_data, analysis_type, selected_cluster, current_gene = NULL, experiment_info = NULL) {
      # Check if EnhancedVolcano is available
      if (!enhancedvolcano_available) {
        # Return installation message plot
        return(
          ggplot() +
            annotate("text", x = 0, y = 0, 
                    label = "EnhancedVolcano package required for static plots\nPlease install with: BiocManager::install('EnhancedVolcano')", 
                    size = 5, color = "red", hjust = 0.5, vjust = 0.5) +
            theme_void() +
            labs(title = paste(analysis_type, "Volcano Plot - Package Missing"))
        )
      }
      
      tryCatch({
        if (is.null(de_data) || nrow(de_data) == 0) {
          return(
            ggplot() +
              annotate("text", x = 0, y = 0, 
                      label = "No differential expression data available", 
                      size = 5, color = "gray50", hjust = 0.5, vjust = 0.5) +
              theme_void() +
              labs(title = paste(analysis_type, "- No Data Available"))
          )
        }
        
        # Filter by cluster if selected
        if (!is.null(selected_cluster) && selected_cluster != "All") {
          plot_data <- de_data[de_data$cluster == selected_cluster, ]
          cluster_label <- gsub("cluster_", "", selected_cluster)
        } else {
          plot_data <- de_data
          cluster_label <- "All Clusters"
        }
        
        if (nrow(plot_data) == 0) {
          return(
            ggplot() +
              annotate("text", x = 0, y = 0, 
                      label = paste("No data for", selected_cluster), 
                      size = 5, color = "gray50", hjust = 0.5, vjust = 0.5) +
              theme_void() +
              labs(title = paste(analysis_type, "- No Data for Selected Cluster"))
          )
        }
        
        # Prepare data for EnhancedVolcano
        # CRITICAL FIX: Ensure gene names are proper, not row numbers
        if ("gene_name" %in% colnames(plot_data)) {
          # Check if gene_name contains actual names or just numbers
          if (all(grepl("^[0-9]+$", plot_data$gene_name[1:min(10, nrow(plot_data))]))) {
            cat("[STATIC VOLCANO DEBUG] gene_name appears to be row numbers, attempting to fix\n")
            # If gene_name contains only numbers, create meaningful names
            plot_data$gene_name <- paste0("Gene_", seq_len(nrow(plot_data)))
          }
          rownames(plot_data) <- make.unique(plot_data$gene_name)  # Ensure unique rownames
        } else {
          plot_data$gene_name <- paste0("Gene_", seq_len(nrow(plot_data)))
          rownames(plot_data) <- plot_data$gene_name
        }
        
        # Generate highly descriptive title with multi-line format
        if (analysis_type == "MAST") {
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            plot_title <- paste0(current_gene, " mutation vs isogenic eWT controls\n", 
                                "MAST Analysis\n",
                                "Cluster ", cluster_label)
          } else {
            plot_title <- paste0("MAST mutation analysis vs isogenic eWT controls\n",
                                "All Genes\n", 
                                "Cluster ", cluster_label)
          }
        } else if (analysis_type == "MixScale") {
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            if (!is.null(experiment_info) && experiment_info != "" && experiment_info != "default") {
              plot_title <- paste0(current_gene, " CRISPRi knockdown vs Non-Targeting\n",
                                  "MixScale Analysis (", experiment_info, ")\n",
                                  "Cluster ", cluster_label)
            } else {
              plot_title <- paste0(current_gene, " CRISPRi knockdown vs Non-Targeting controls\n",
                                  "MixScale Analysis\n",
                                  "Cluster ", cluster_label)
            }
          } else {
            if (!is.null(experiment_info) && experiment_info != "" && experiment_info != "default") {
              plot_title <- paste0("CRISPRi knockdown vs Non-Targeting\n",
                                  "MixScale Analysis (", experiment_info, ")\n", 
                                  "Cluster ", cluster_label)
            } else {
              plot_title <- paste0("CRISPRi knockdown vs Non-Targeting controls\n",
                                  "MixScale Analysis\n",
                                  "Cluster ", cluster_label)
            }
          }
        } else {
          plot_title <- paste(analysis_type, "\nCluster", cluster_label)
        }
        
        # Select top significant genes for labeling (up to 20)
        plot_data$significant <- plot_data$pvalue < 0.05 & abs(plot_data$log2FC) > 0.25
        significant_genes <- plot_data[plot_data$significant, ]
        
        # Sort by combined significance and fold change for best labels
        if (nrow(significant_genes) > 0) {
          significant_genes$rank_score <- abs(significant_genes$log2FC) * (-log10(significant_genes$pvalue + 1e-300))
          significant_genes <- significant_genes[order(significant_genes$rank_score, decreasing = TRUE), ]
          
          # Select top 15-20 genes for labeling
          n_labels <- min(20, nrow(significant_genes))
          genes_to_label <- rownames(significant_genes)[1:n_labels]
        } else {
          # If no significant genes, label top genes by absolute fold change
          top_genes <- plot_data[order(abs(plot_data$log2FC), decreasing = TRUE), ]
          n_labels <- min(15, nrow(top_genes))
          genes_to_label <- rownames(top_genes)[1:n_labels]
        }
        
        # Create the enhanced volcano plot
        p <- EnhancedVolcano(
          plot_data,
          lab = rownames(plot_data),
          x = 'log2FC',
          y = 'pvalue',
          title = plot_title,
          subtitle = paste("Total genes:", nrow(plot_data), "| Significant:", sum(plot_data$significant)),
          caption = paste("Significance: p < 0.05 & |log2FC| > 0.25"),
          pCutoff = 0.05,
          FCcutoff = 0.25,
          selectLab = genes_to_label,
          xlab = bquote(~Log[2]~ 'fold change'),
          ylab = bquote(~-Log[10]~'P-value'),
          axisLabSize = 12,
          titleLabSize = 14,
          subtitleLabSize = 11,
          captionLabSize = 10,
          labSize = 3.5,
          labCol = 'black',
          labFace = 'bold',
          boxedLabels = TRUE,
          colAlpha = 0.7,
          legendPosition = 'bottom',
          legendLabSize = 8,
          legendIconSize = 3.0,
          drawConnectors = TRUE,
          widthConnectors = 0.3,
          colConnectors = 'black',
          max.overlaps = Inf,
          # Color scheme
          col = c('grey30', 'forestgreen', 'royalblue', 'red2'),
          colCustom = NULL,
          # Grid and background
          gridlines.major = TRUE,
          gridlines.minor = FALSE,
          border = 'partial',
          borderWidth = 0.8,
          borderColour = 'black'
        )
        
        # Additional styling for publication quality
        p <- p + 
          theme(
            plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 11, hjust = 0.5),
            plot.caption = element_text(size = 10, hjust = 0.5),
            panel.background = element_rect(fill = "white", color = "black"),
            plot.background = element_rect(fill = "white", color = NA),
            legend.background = element_rect(fill = "white", color = "black"),
            legend.key = element_rect(fill = "white"),
            text = element_text(family = "Arial", color = "black")
          )
        
        return(p)
        
      }, error = function(e) {
        # Return error plot
        ggplot() +
          annotate("text", x = 0, y = 0, 
                  label = paste("Error generating plot:", e$message), 
                  size = 4, color = "red", hjust = 0.5, vjust = 0.5) +
          theme_void() +
          labs(title = paste(analysis_type, "- Plot Generation Error"))
      })
    }
    
    # REMOVED DUPLICATE CONTAINER DEFINITIONS - Using the first implementations with proper width fixes
    
    # Render MAST static volcano plot
    output$mast_volcano_static <- renderPlot({
      tryCatch({
        cat("[STATIC PLOT] Starting MAST volcano plot generation\n")
        
        # Load DE data on demand when volcano plot is first rendered
        load_de_data()
        
        # Check if cluster is selected
        if (is.null(values$selected_cluster) || values$selected_cluster == "") {
          cat("[STATIC PLOT] No cluster selected\n")
          return(
            ggplot() +
              annotate("text", x = 0, y = 0, 
                      label = "Click on a cluster in the UMAP\nto view its differential expression results", 
                      size = 5, color = "#3c8dbc", hjust = 0.5, vjust = 0.5) +
              theme_void() +
              labs(title = "MAST Volcano Plot") +
              theme(plot.background = element_rect(fill = "white", color = NA))
          )
        } else if (is.null(values$de_data_mast)) {
          cat("[STATIC PLOT] No MAST data available\n")
          return(
            ggplot() +
              annotate("text", x = 0, y = 0, 
                      label = "MAST DE results not available.\nPlease ensure full_DE_results.rds is present.", 
                      size = 5, color = "gray50", hjust = 0.5, vjust = 0.5) +
              theme_void() +
              labs(title = "MAST Volcano Plot - No Data") +
              theme(plot.background = element_rect(fill = "white", color = NA))
          )
        } else {
          # Filter by global gene selection
          current_gene <- global_selection()$gene
          current_experiment <- global_selection()$experiment
          
          cat("[STATIC PLOT] Current gene:", current_gene, "| Cluster:", values$selected_cluster, "\n")
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            filtered_data <- values$de_data_mast[values$de_data_mast$gene == current_gene, ]
          } else {
            filtered_data <- values$de_data_mast
          }
          
          cat("[STATIC PLOT] Filtered data rows:", nrow(filtered_data), "\n")
          
          # Test with simple plot first if EnhancedVolcano fails
          if (!enhancedvolcano_available) {
            cat("[STATIC PLOT] EnhancedVolcano not available, creating simple plot\n")
            return(
              ggplot() +
                annotate("text", x = 0, y = 0, 
                        label = "EnhancedVolcano package required for static plots\nPlease install with: BiocManager::install('EnhancedVolcano')", 
                        size = 4, color = "red", hjust = 0.5, vjust = 0.5) +
                theme_void() +
                labs(title = "MAST Volcano Plot - Package Missing") +
                theme(plot.background = element_rect(fill = "white", color = NA))
            )
          }
          
          result_plot <- generate_static_volcano_plot(filtered_data, "MAST", values$selected_cluster, current_gene, current_experiment)
          cat("[STATIC PLOT] Generated static plot, class:", class(result_plot), "\n")
          return(result_plot)
        }
      }, error = function(e) {
        cat("[STATIC PLOT] Error:", e$message, "\n")
        ggplot() +
          annotate("text", x = 0, y = 0, 
                  label = paste("Error generating plot:", e$message), 
                  size = 4, color = "red", hjust = 0.5, vjust = 0.5) +
          theme_void() +
          labs(title = "MAST Volcano Plot - Error") +
          theme(plot.background = element_rect(fill = "white", color = NA))
      })
    }, width = 600, height = 400)
    
    # Render MixScale static volcano plot
    output$mixscale_volcano_static <- renderPlot({
      tryCatch({
        # Load DE data on demand when volcano plot is first rendered
        load_de_data()
        
        # Check if cluster is selected
        if (is.null(values$selected_cluster) || values$selected_cluster == "") {
          ggplot() +
            annotate("text", x = 0, y = 0, 
                    label = "Click on a cluster in the UMAP\nto view its differential expression results", 
                    size = 5, color = "#3c8dbc", hjust = 0.5, vjust = 0.5) +
            theme_void() +
            labs(title = "MixScale Volcano Plot")
        } else if (is.null(values$de_data_mixscale)) {
          ggplot() +
            annotate("text", x = 0, y = 0, 
                    label = "MixScale DE results not available.\nPlease ensure full_DE_results.rds is present.", 
                    size = 5, color = "gray50", hjust = 0.5, vjust = 0.5) +
            theme_void() +
            labs(title = "MixScale Volcano Plot - No Data")
        } else {
          # Filter by global gene selection
          current_gene <- global_selection()$gene
          current_experiment <- global_selection()$experiment
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            # CRITICAL FIX: Apply gene harmonization for MixScale volcano plot
            # Map MAST gene names to corresponding MixScale names
            mixscale_gene <- current_gene
            if (current_gene == "PRKN") {
              mixscale_gene <- "PARK2"
            } else if (current_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
              mixscale_gene <- "SNCA" 
            } else if (current_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
              mixscale_gene <- "VPS13C"
            }
            
            filtered_data <- values$de_data_mixscale[values$de_data_mixscale$gene == mixscale_gene, ]
            cat("[MixScale Volcano] Gene harmonization applied - MAST:", current_gene, "→ MixScale:", mixscale_gene, "\n")
          } else {
            filtered_data <- values$de_data_mixscale
          }
          
          # Get experiment info for title
          experiment_info <- if (!is.null(current_experiment) && current_experiment != "default" && current_experiment != "") {
            current_experiment
          } else if ("experiment" %in% colnames(filtered_data) && length(unique(filtered_data$experiment)) == 1) {
            unique(filtered_data$experiment)[1]
          } else {
            NULL
          }
          
          generate_static_volcano_plot(filtered_data, "MixScale", values$selected_cluster, current_gene, experiment_info)
        }
      }, error = function(e) {
        ggplot() +
          annotate("text", x = 0, y = 0, 
                  label = paste("Error generating plot:", e$message), 
                  size = 4, color = "red", hjust = 0.5, vjust = 0.5) +
          theme_void() +
          labs(title = "MixScale Volcano Plot - Error")
      })
    })
    
    # Render MAST interactive volcano plot (renamed from original)
    output$mast_volcano_interactive <- renderPlotly({
      tryCatch({
        
        # Load DE data on demand when volcano plot is first rendered
        load_de_data()
        
        # Check if cluster is selected
        if (is.null(values$selected_cluster) || values$selected_cluster == "") {
          plot_ly() %>%
            layout(
              title = "MAST Volcano Plot",
              xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
              yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
              annotations = list(
                x = 0,
                y = 5,
                text = "Click on a cluster in the UMAP to view its differential expression results",
                showarrow = FALSE,
                font = list(size = 16, color = "#3c8dbc")
              )
            )
        } else if (is.null(values$de_data_mast)) {
          # No data available - show empty plot with message
          plot_ly() %>%
            layout(
              title = "MAST Volcano Plot - No Data",
              xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
              yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
              annotations = list(
                x = 0,
                y = 5,
                text = "MAST DE results not available.\nPlease ensure full_DE_results.rds is present.",
                showarrow = FALSE,
                font = list(size = 16, color = "gray")
              )
            )
        } else {
          # Filter by global gene selection - CLEANED: Removed cat() statements
          current_gene <- global_selection()$gene
          current_experiment <- global_selection()$experiment
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            filtered_data <- values$de_data_mast[values$de_data_mast$gene == current_gene, ]
          } else {
            filtered_data <- values$de_data_mast
          }
          
          generate_volcano_plot(filtered_data, "MAST", values$selected_cluster, input$color_by, current_gene, current_experiment)
        }
      }, error = function(e) {
        # BULLETPROOF: Return a simple error plot instead of plotly_empty()
        showNotification("Error rendering MAST volcano plot", type = "error")
        plot_ly() %>%
          layout(
            title = "MAST Volcano Plot - Error",
            xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
            yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
            annotations = list(
              x = 0,
              y = 5,
              text = paste("Error:", e$message),
              showarrow = FALSE,
              font = list(size = 14, color = "red")
            )
          )
      })
    })
    
    # Render MixScale interactive volcano plot (renamed from original)
    output$mixscale_volcano_interactive <- renderPlotly({
      tryCatch({
        
        # Load DE data on demand when volcano plot is first rendered
        load_de_data()
        
        # Check if cluster is selected
        if (is.null(values$selected_cluster) || values$selected_cluster == "") {
          plot_ly() %>%
            layout(
              title = "MixScale Volcano Plot",
              xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
              yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
              annotations = list(
                x = 0,
                y = 5,
                text = "Click on a cluster in the UMAP to view its differential expression results",
                showarrow = FALSE,
                font = list(size = 16, color = "#3c8dbc")
              )
            )
        } else if (is.null(values$de_data_mixscale)) {
          # No data available - show empty plot with message
          plot_ly() %>%
            layout(
              title = "MixScale Volcano Plot - No Data",
              xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
              yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
              annotations = list(
                x = 0,
                y = 5,
                text = "MixScale DE results not available.\nPlease ensure full_DE_results.rds is present.",
                showarrow = FALSE,
                font = list(size = 16, color = "gray")
              )
            )
        } else {
          # Filter by global gene selection - CLEANED: Removed cat() statements
          current_gene <- global_selection()$gene
          current_experiment <- global_selection()$experiment
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            # CRITICAL FIX: Apply gene harmonization for MixScale interactive volcano plot
            # Map MAST gene names to corresponding MixScale names
            mixscale_gene <- current_gene
            if (current_gene == "PRKN") {
              mixscale_gene <- "PARK2"
            } else if (current_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
              mixscale_gene <- "SNCA" 
            } else if (current_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
              mixscale_gene <- "VPS13C"
            }
            
            filtered_data <- values$de_data_mixscale[values$de_data_mixscale$gene == mixscale_gene, ]
          } else {
            filtered_data <- values$de_data_mixscale
          }
          
          # Get experiment info for title
          experiment_info <- if (!is.null(current_experiment) && current_experiment != "default" && current_experiment != "") {
            current_experiment
          } else if ("experiment" %in% colnames(filtered_data) && length(unique(filtered_data$experiment)) == 1) {
            unique(filtered_data$experiment)[1]
          } else {
            NULL
          }
          
          generate_volcano_plot(filtered_data, "MixScale", values$selected_cluster, input$color_by, current_gene, experiment_info)
        }
      }, error = function(e) {
        # BULLETPROOF: Return a simple error plot instead of plotly_empty()
        showNotification("Error rendering MixScale volcano plot", type = "error")
        plot_ly() %>%
          layout(
            title = "MixScale Volcano Plot - Error",
            xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
            yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
            annotations = list(
              x = 0,
              y = 5,
              text = paste("Error:", e$message),
              showarrow = FALSE,
              font = list(size = 14, color = "red")
            )
          )
      })
    })
    
    # REACTIVE LOOP FIX: Move side effects to separate observer
    observeEvent(global_selection(), {
      # Clear overlap data values when global selection changes to prevent stale data
      # This ensures when switching to genes with no MixScale data (like GBA), 
      # the overlap section properly shows "no data" instead of stale results
      values$mast_sig_data <- NULL
      values$mixscale_sig_data <- NULL
      values$overlap_data <- NULL
    }, priority = 100)  # High priority to clear data before other reactions
    
    # REACTIVE LOOP FIX: Pure reactive for statistics computation (no side effects)
    stats_data <- reactive({
      # Ensure proper reactivity by accessing all reactive dependencies
      req(global_selection())
      current_selection <- global_selection()
      
      # Also ensure reactivity to local cluster selection changes
      # Note: values$selected_cluster can be NULL initially, so don't require it
      
      # Get display text for current settings
      cluster_text <- if (is.null(values$selected_cluster) || values$selected_cluster == "All") {
        "all clusters"
      } else {
        values$selected_cluster
      }
      
      gene_text <- if (is.null(current_selection$gene) || current_selection$gene == "All" || current_selection$gene == "") {
        "All genes"
      } else {
        current_selection$gene
      }
      
      # Calculate stats from actual data
      mast_sig <- 0
      mixscale_sig <- 0
      overlap <- 0
      
      # Get current global selection for filtering
      current_gene <- current_selection$gene
      
      # Calculate MAST significant genes
      if (!is.null(values$de_data_mast) && nrow(values$de_data_mast) > 0) {
        mast_filtered <- values$de_data_mast
        
        # Filter by selected cluster
        if (!is.null(values$selected_cluster) && values$selected_cluster != "All" && values$selected_cluster != "") {
          mast_filtered <- mast_filtered[mast_filtered$cluster == values$selected_cluster, ]
        }
        
        # Filter by global gene selection (MAST uses original gene names)
        if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
          mast_filtered <- mast_filtered[mast_filtered$gene == current_gene, ]
        }
        
        # Apply significance criteria
        if (nrow(mast_filtered) > 0) {
          mast_sig <- sum(!is.na(mast_filtered$pvalue) & !is.na(mast_filtered$log2FC) &
                         mast_filtered$pvalue < 0.05 & abs(mast_filtered$log2FC) > 0.25)
        }
      }
      
      # Calculate MixScale significant genes
      if (!is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0) {
        mixscale_filtered <- values$de_data_mixscale
        
        # Filter by selected cluster
        if (!is.null(values$selected_cluster) && values$selected_cluster != "All" && values$selected_cluster != "") {
          mixscale_filtered <- mixscale_filtered[mixscale_filtered$cluster == values$selected_cluster, ]
        }
        
        # Filter by global gene selection with gene harmonization
        if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
          # Apply gene harmonization for MixScale data
          mixscale_gene <- current_gene
          if (current_gene == "PRKN") {
            mixscale_gene <- "PARK2"
          } else if (current_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
            mixscale_gene <- "SNCA" 
          } else if (current_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
            mixscale_gene <- "VPS13C"
          }
          
          mixscale_filtered <- mixscale_filtered[mixscale_filtered$gene == mixscale_gene, ]
        }
        
        # Apply significance criteria
        if (nrow(mixscale_filtered) > 0) {
          mixscale_sig <- sum(!is.na(mixscale_filtered$pvalue) & !is.na(mixscale_filtered$log2FC) &
                             mixscale_filtered$pvalue < 0.05 & abs(mixscale_filtered$log2FC) > 0.25)
        }
      }
      
      # Calculate overlap between MAST and MixScale significant genes WITH DIRECTION INFO
      mast_sig_genes <- character(0)
      mixscale_sig_genes <- character(0)
      mast_sig_data <- NULL
      mixscale_sig_data <- NULL
      
      if (mast_sig > 0 && mixscale_sig > 0 && !is.null(values$de_data_mast) && !is.null(values$de_data_mixscale)) {
        # Get gene names that are significant in both
        mast_filtered <- values$de_data_mast
        mixscale_filtered <- values$de_data_mixscale
        
        # Apply same filtering as above
        if (!is.null(values$selected_cluster) && values$selected_cluster != "All" && values$selected_cluster != "") {
          mast_filtered <- mast_filtered[mast_filtered$cluster == values$selected_cluster, ]
          mixscale_filtered <- mixscale_filtered[mixscale_filtered$cluster == values$selected_cluster, ]
        }
        
        if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
          # Apply gene harmonization for cross-method compatibility
          mast_gene <- current_gene
          mixscale_gene <- current_gene
          
          # Map MAST gene names to corresponding MixScale names  
          if (current_gene == "PRKN") {
            mixscale_gene <- "PARK2"
          } else if (current_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
            mixscale_gene <- "SNCA" 
          } else if (current_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
            mixscale_gene <- "VPS13C"
          }
          
          mast_filtered <- mast_filtered[mast_filtered$gene == mast_gene, ]
          mixscale_filtered <- mixscale_filtered[mixscale_filtered$gene == mixscale_gene, ]
        }
        
        # Get significant gene data with direction info
        if (nrow(mast_filtered) > 0) {
          mast_sig_idx <- !is.na(mast_filtered$pvalue) & !is.na(mast_filtered$log2FC) &
                         mast_filtered$pvalue < 0.05 & abs(mast_filtered$log2FC) > 0.25
          mast_sig_data <- mast_filtered[mast_sig_idx, c("gene_name", "log2FC", "pvalue")]
          mast_sig_genes <- mast_sig_data$gene_name
        }
        
        if (nrow(mixscale_filtered) > 0) {
          mixscale_sig_idx <- !is.na(mixscale_filtered$pvalue) & !is.na(mixscale_filtered$log2FC) &
                             mixscale_filtered$pvalue < 0.05 & abs(mixscale_filtered$log2FC) > 0.25
          mixscale_sig_data <- mixscale_filtered[mixscale_sig_idx, c("gene_name", "log2FC", "pvalue")]
          mixscale_sig_genes <- mixscale_sig_data$gene_name
        }
        
        # Calculate overlap
        overlap <- length(intersect(mast_sig_genes, mixscale_sig_genes))
        
        # Return significant gene data (no side effects)
        # Side effects will be handled by separate observer
      }
      
      # UPDATED: Statistical significance with proper background gene handling (intersection/union approaches)
      fisher_results <- list(
        intersection_approach = list(fisher_p = NA, fisher_or = NA, background_size = 0),
        union_approach = list(fisher_p = NA, fisher_or = NA, background_size = 0)
      )
      
      # DEBUG: Log DE data availability for user troubleshooting
      cat("[DE Results] DEBUG - DE data status:\n")
      cat("  MAST DE data loaded:", !is.null(values$de_data_mast), "\n")
      if (!is.null(values$de_data_mast)) {
        cat("  MAST DE data rows:", nrow(values$de_data_mast), "\n")
      }
      cat("  MixScale DE data loaded:", !is.null(values$de_data_mixscale), "\n") 
      if (!is.null(values$de_data_mixscale)) {
        cat("  MixScale DE data rows:", nrow(values$de_data_mixscale), "\n")
      }
      cat("  Current overlap conditions: MAST sig >0:", mast_sig > 0, ", MixScale sig >0:", mixscale_sig > 0, ", overlap >=0:", overlap >= 0, "\n")
      
      if (mast_sig > 0 && mixscale_sig > 0 && overlap >= 0) {
        # Get proper background genes from all tested genes (not just significant ones)
        if (!is.null(values$de_data_mast) && !is.null(values$de_data_mixscale)) {
          # CRITICAL FIX: Apply same filtering logic as significance calculation 
          # to get genes tested in CURRENT selection only (not all comparisons)
          mast_background_data <- values$de_data_mast
          mixscale_background_data <- values$de_data_mixscale
          
          # Filter by selected cluster (same logic as significance calculation)
          if (!is.null(values$selected_cluster) && values$selected_cluster != "All" && values$selected_cluster != "") {
            mast_background_data <- mast_background_data[mast_background_data$cluster == values$selected_cluster, ]
            mixscale_background_data <- mixscale_background_data[mixscale_background_data$cluster == values$selected_cluster, ]
          }
          
          # Filter by global gene selection with gene harmonization
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            # Apply gene harmonization for cross-method compatibility
            mast_gene <- current_gene
            mixscale_gene <- current_gene
            
            # Map MAST gene names to corresponding MixScale names  
            if (current_gene == "PRKN") {
              mixscale_gene <- "PARK2"
            } else if (current_gene %in% c("SNCA_A30P", "SNCA_A53T")) {
              mixscale_gene <- "SNCA" 
            } else if (current_gene %in% c("VPS13C_A444P", "VPS13C_W395C")) {
              mixscale_gene <- "VPS13C"
            }
            
            mast_background_data <- mast_background_data[mast_background_data$gene == mast_gene, ]
            mixscale_background_data <- mixscale_background_data[mixscale_background_data$gene == mixscale_gene, ]
            
            cat("[DE Results] Gene harmonization applied - MAST:", mast_gene, ", MixScale:", mixscale_gene, "\n")
          }
          
          # NOW extract genes tested in the current specific comparison
          mast_all_genes <- unique(mast_background_data$gene_name)
          mixscale_all_genes <- unique(mixscale_background_data$gene_name)
          
          cat("[DE Results] Filtered background data - MAST rows:", nrow(mast_background_data), ", MixScale rows:", nrow(mixscale_background_data), "\n")
          
          # INTERSECTION APPROACH (Conservative): Genes tested in BOTH methods
          intersection_background <- intersect(mast_all_genes, mixscale_all_genes)
          intersection_size <- length(intersection_background)
          
          # UNION APPROACH (Liberal): Genes tested in EITHER method  
          union_background <- unique(c(mast_all_genes, mixscale_all_genes))
          union_size <- length(union_background)
          
          cat("[DE Results] Background genes - Intersection:", intersection_size, ", Union:", union_size, "\n")
          cat("[DE Results] Attempting Fisher's exact test calculations...\n")
          
          # Calculate Fisher's exact test for both approaches
          if (intersection_size > max(mast_sig, mixscale_sig) && intersection_size > 0) {
            # Filter significant genes to intersection background
            mast_sig_filtered <- intersect(mast_sig_genes, intersection_background)
            mixscale_sig_filtered <- intersect(mixscale_sig_genes, intersection_background)
            overlap_filtered <- length(intersect(mast_sig_filtered, mixscale_sig_filtered))
            
            # Contingency table for intersection approach
            mast_only <- length(mast_sig_filtered) - overlap_filtered
            mixscale_only <- length(mixscale_sig_filtered) - overlap_filtered  
            neither <- intersection_size - length(mast_sig_filtered) - mixscale_only
            
            if (mast_only >= 0 && mixscale_only >= 0 && neither >= 0) {
              contingency_matrix <- matrix(c(overlap_filtered, mast_only, mixscale_only, neither), nrow = 2, byrow = TRUE)
              tryCatch({
                fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
                fisher_results$intersection_approach$fisher_p <- fisher_result$p.value
                fisher_results$intersection_approach$fisher_or <- fisher_result$estimate
                fisher_results$intersection_approach$background_size <- intersection_size
              }, error = function(e) {
                cat("[DE Results] Intersection Fisher's test failed:", e$message, "\n")
              })
            }
          }
          
          # Calculate Fisher's exact test for union approach
          if (union_size > max(mast_sig, mixscale_sig) && union_size > 0) {
            # Use original overlap (genes can be in either background)
            mast_only <- mast_sig - overlap
            mixscale_only <- mixscale_sig - overlap
            neither <- union_size - mast_sig - mixscale_only
            
            if (mast_only >= 0 && mixscale_only >= 0 && neither >= 0) {
              contingency_matrix <- matrix(c(overlap, mast_only, mixscale_only, neither), nrow = 2, byrow = TRUE)
              tryCatch({
                fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
                fisher_results$union_approach$fisher_p <- fisher_result$p.value
                fisher_results$union_approach$fisher_or <- fisher_result$estimate  
                fisher_results$union_approach$background_size <- union_size
              }, error = function(e) {
                cat("[DE Results] Union Fisher's test failed:", e$message, "\n")
              })
            }
          }
        } else {
          cat("[DE Results] No DE data available for proper background calculation\n")
        }
      }
      
      # DEBUG: Log final Fisher's test results
      cat("[DE Results] Final Fisher's test results:\n")
      cat("  Intersection p-value:", fisher_results$intersection_approach$fisher_p, "\n")
      cat("  Union p-value:", fisher_results$union_approach$fisher_p, "\n")
      cat("  Will show dual Fisher's test?:", (!is.na(fisher_results$intersection_approach$fisher_p) || !is.na(fisher_results$union_approach$fisher_p)), "\n")
      
      # Return computed data (no UI rendering in reactive)
      return(list(
        gene_text = gene_text,
        cluster_text = cluster_text,
        mast_sig = mast_sig,
        mixscale_sig = mixscale_sig,
        overlap = overlap,
        fisher_results = fisher_results,
        mast_sig_data = mast_sig_data,
        mixscale_sig_data = mixscale_sig_data
      ))
    })
    
    # Observer to handle side effects (updating values$ based on computed stats)
    observeEvent(stats_data(), {
      stats <- stats_data()
      if (!is.null(stats)) {
        values$mast_sig_data <- stats$mast_sig_data
        values$mixscale_sig_data <- stats$mixscale_sig_data
      }
    })
    
    # Pure UI renderer (no side effects)
    output$stats_content <- renderUI({
      stats <- stats_data()
      req(stats)
      
      tagList(
        # Dynamic title showing current settings
        div(class = "text-center", style = "margin-bottom: 20px;",
          h4(paste("Summary Statistics:", stats$gene_text, "-", stats$cluster_text), 
             style = "color: #2c3e50; margin-bottom: 5px;")
        ),
        
        # Explanatory note about analysis scope
        div(class = "alert alert-info", style = "margin-bottom: 15px; font-size: 0.9em;",
          icon("info-circle"),
          strong(" Analysis Scope: "), 
          "This shows ", tags$strong("direct overlap of differentially expressed genes"), " meeting significance criteria (p < 0.05, |log2FC| > 0.25). ",
          "This differs from the pathway-based gene overlaps shown on the Signature Nomination page, which analyze genes contributing to enriched pathways. ",
          "For genes appearing in enriched pathways only (smaller subset), see Signature Heatmap in Signature Nomination."
        ),
        
        fluidRow(
          column(4,
            div(class = "text-center",
              h5("MAST Significant"),
              h6("(All DE Genes)", style = "color: #7f8c8d; margin-top: -5px;"),
              h3(stats$mast_sig, style = "color: #3c8dbc;"),
              p("(p < 0.05, |log2FC| > 0.25)")
            )
          ),
          column(4,
            div(class = "text-center",
              h5("MixScale Significant"),
              h6("(All DE Genes)", style = "color: #7f8c8d; margin-top: -5px;"),
              h3(stats$mixscale_sig, style = "color: #5cb85c;"),
              p("(p < 0.05, |log2FC| > 0.25)")
            )
          ),
          column(4,
            div(class = "text-center",
              h5("Overlapping Genes"),
              h6("(All DE Genes)", style = "color: #7f8c8d; margin-top: -5px;"),
              h3(stats$overlap, style = "color: #f0ad4e;"),
              p("(significant in both)")
            )
          )
        ),
        
        # UPDATED: Statistical significance test results with dual approach
        if (!is.na(stats$fisher_results$intersection_approach$fisher_p) || !is.na(stats$fisher_results$union_approach$fisher_p)) {
          tagList(
            hr(),
            div(class = "text-center",
              h5("Overlap Significance Test", style = "margin-bottom: 15px;"),
              div(class = "alert alert-info", style = "margin-bottom: 15px; font-size: 0.9em;",
                icon("info-circle"),
                strong(" Statistical Approaches: "), 
                "Intersection (conservative) uses genes tested in BOTH methods. ",
                "Union (liberal) uses genes tested in EITHER method. ",
                "Both approaches use proper background genes from DE analysis."
              ),
              fluidRow(
                # Intersection approach column
                if (!is.na(stats$fisher_results$intersection_approach$fisher_p)) {
                  column(6,
                    div(
                      strong("Fisher's Test (Intersection)"),
                      br(),
                      span("p-value: ", style = "color: #666;"),
                      span(format_pvalue(stats$fisher_results$intersection_approach$fisher_p), 
                           style = paste0("color: ", if (stats$fisher_results$intersection_approach$fisher_p < 0.05) "#d9534f" else "#333", ";")),
                      br(),
                      span("Background: ", style = "color: #666;"),
                      span(format(stats$fisher_results$intersection_approach$background_size, big.mark = ",")),
                      br(),
                      span("Odds Ratio: ", style = "color: #666;"),
                      span(format(stats$fisher_results$intersection_approach$fisher_or, digits = 2)),
                      br(),
                      span("(Conservative)", style = "color: #5cb85c; font-size: 0.9em;")
                    )
                  )
                } else {
                  column(6,
                    div(
                      strong("Fisher's Test (Intersection)"),
                      br(),
                      span("Not available", style = "color: #999;"),
                      br(),
                      span("(Insufficient data)", style = "color: #666; font-size: 0.9em;")
                    )
                  )
                },
                
                # Union approach column
                if (!is.na(stats$fisher_results$union_approach$fisher_p)) {
                  column(6,
                    div(
                      strong("Fisher's Test (Union)"),
                      br(),
                      span("p-value: ", style = "color: #666;"),
                      span(format_pvalue(stats$fisher_results$union_approach$fisher_p), 
                           style = paste0("color: ", if (stats$fisher_results$union_approach$fisher_p < 0.05) "#d9534f" else "#333", ";")),
                      br(),
                      span("Background: ", style = "color: #666;"),
                      span(format(stats$fisher_results$union_approach$background_size, big.mark = ",")),
                      br(),
                      span("Odds Ratio: ", style = "color: #666;"),
                      span(format(stats$fisher_results$union_approach$fisher_or, digits = 2)),
                      br(),
                      span("(Liberal)", style = "color: #f0ad4e; font-size: 0.9em;")
                    )
                  )
                } else {
                  column(6,
                    div(
                      strong("Fisher's Test (Union)"),
                      br(),
                      span("Not available", style = "color: #999;"),
                      br(),
                      span("(Insufficient data)", style = "color: #666; font-size: 0.9em;")
                    )
                  )
                }
              )
            )
          )
        } else {
          tagList(
            hr(),
            div(class = "text-center text-muted",
              p("Statistical test requires both MAST and MixScale data")
            )
          )
        },
        
        hr(),
        p(paste("Statistics for", stats$cluster_text), class = "text-muted text-center")
      )
    })
    
    # Render overlapping genes content
    output$overlap_content <- renderUI({
      if (is.null(values$mast_sig_data) || is.null(values$mixscale_sig_data)) {
        return(
          div(class = "text-center text-muted",
            p("No overlapping genes found"),
            tags$small("Requires significant genes in both MAST and MixScale results")
          )
        )
      }
      
      # Get overlapping genes with direction information
      mast_genes <- values$mast_sig_data
      mixscale_genes <- values$mixscale_sig_data
      
      # Find overlapping genes
      overlap_genes <- intersect(mast_genes$gene_name, mixscale_genes$gene_name)
      
      if (length(overlap_genes) == 0) {
        return(
          div(class = "text-center text-muted",
            p("No overlapping genes found"),
            tags$small("Different significant genes between MAST and MixScale")
          )
        )
      }
      
      # Categorize overlapping genes by direction
      same_direction_up <- character(0)
      same_direction_down <- character(0)
      opposite_direction <- character(0)
      
      for (gene in overlap_genes) {
        mast_fc <- mast_genes$log2FC[mast_genes$gene_name == gene][1]
        mixscale_fc <- mixscale_genes$log2FC[mixscale_genes$gene_name == gene][1]
        
        if (!is.na(mast_fc) && !is.na(mixscale_fc)) {
          if (mast_fc > 0 && mixscale_fc > 0) {
            same_direction_up <- c(same_direction_up, gene)
          } else if (mast_fc < 0 && mixscale_fc < 0) {
            same_direction_down <- c(same_direction_down, gene)
          } else {
            opposite_direction <- c(opposite_direction, gene)
          }
        }
      }
      
      # Create downloadable gene lists
      all_overlap_data <- data.frame(
        Gene = overlap_genes,
        MAST_log2FC = sapply(overlap_genes, function(g) {
          idx <- which(mast_genes$gene_name == g)[1]
          if (length(idx) > 0 && !is.na(idx)) mast_genes$log2FC[idx] else NA
        }),
        MAST_pvalue = sapply(overlap_genes, function(g) {
          idx <- which(mast_genes$gene_name == g)[1]
          if (length(idx) > 0 && !is.na(idx)) mast_genes$pvalue[idx] else NA
        }),
        MixScale_log2FC = sapply(overlap_genes, function(g) {
          idx <- which(mixscale_genes$gene_name == g)[1]
          if (length(idx) > 0 && !is.na(idx)) mixscale_genes$log2FC[idx] else NA
        }),
        MixScale_pvalue = sapply(overlap_genes, function(g) {
          idx <- which(mixscale_genes$gene_name == g)[1]
          if (length(idx) > 0 && !is.na(idx)) mixscale_genes$pvalue[idx] else NA
        }),
        stringsAsFactors = FALSE
      )
      
      # Store for download
      values$overlap_data <- all_overlap_data
      
      # Create UI
      tagList(
        # Summary counts
        div(style = "margin-bottom: 15px;",
          fluidRow(
            column(4,
              div(class = "text-center",
                h6("Same ↑", style = "color: #d9534f; margin-bottom: 5px;"),
                h4(length(same_direction_up), style = "color: #d9534f; margin: 0;"),
                tags$small("Both up-regulated")
              )
            ),
            column(4,
              div(class = "text-center",
                h6("Same ↓", style = "color: #5bc0de; margin-bottom: 5px;"),
                h4(length(same_direction_down), style = "color: #5bc0de; margin: 0;"),
                tags$small("Both down-regulated")
              )
            ),
            column(4,
              div(class = "text-center",
                h6("Opposite", style = "color: #f0ad4e; margin-bottom: 5px;"),
                h4(length(opposite_direction), style = "color: #f0ad4e; margin: 0;"),
                tags$small("Different directions")
              )
            )
          )
        ),
        
        hr(style = "margin: 10px 0;"),
        
        # Gene lists with expand/collapse
        if (length(same_direction_up) > 0) {
          tagList(
            h6("Same Direction Up-regulated:", style = "color: #d9534f;"),
            div(style = "max-height: 80px; overflow-y: auto; background-color: #f9f9f9; padding: 8px; margin-bottom: 10px; border-radius: 3px;",
              paste(same_direction_up, collapse = ", ")
            )
          )
        },
        
        if (length(same_direction_down) > 0) {
          tagList(
            h6("Same Direction Down-regulated:", style = "color: #5bc0de;"),
            div(style = "max-height: 80px; overflow-y: auto; background-color: #f9f9f9; padding: 8px; margin-bottom: 10px; border-radius: 3px;",
              paste(same_direction_down, collapse = ", ")
            )
          )
        },
        
        if (length(opposite_direction) > 0) {
          tagList(
            h6("Opposite Directions:", style = "color: #f0ad4e;"),
            div(style = "max-height: 80px; overflow-y: auto; background-color: #f9f9f9; padding: 8px; margin-bottom: 10px; border-radius: 3px;",
              paste(opposite_direction, collapse = ", ")
            )
          )
        },
        
        # Download button
        div(class = "text-center", style = "margin-top: 15px;",
          downloadButton(
            session$ns("download_overlap"),
            "Download Gene List",
            class = "btn-sm btn-primary",
            icon = icon("download")
          )
        )
      )
    })
    
    # Download handler for overlapping genes
    output$download_overlap <- downloadHandler(
      filename = function() {
        cluster_suffix <- if (is.null(values$selected_cluster) || values$selected_cluster == "All") {
          "all_clusters"
        } else {
          values$selected_cluster
        }
        
        gene_suffix <- if (is.null(global_selection()$gene) || global_selection()$gene == "All") {
          "all_genes"
        } else {
          global_selection()$gene
        }
        
        paste0("overlapping_DE_genes_", gene_suffix, "_", cluster_suffix, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        if (!is.null(values$overlap_data)) {
          write.csv(values$overlap_data, file, row.names = FALSE)
        } else {
          # Create empty file with headers
          empty_df <- data.frame(
            Gene = character(0),
            MAST_log2FC = numeric(0),
            MAST_pvalue = numeric(0),
            MixScale_log2FC = numeric(0),
            MixScale_pvalue = numeric(0)
          )
          write.csv(empty_df, file, row.names = FALSE)
        }
      }
    )
    
    # === CROSS-CLUSTER ANALYSIS LOGIC ===
    
    # Populate gene choices for cross-cluster analysis
    observe({
      if (!is.null(values$de_data_mast) || !is.null(values$de_data_mixscale)) {
        gene_choices <- c("Choose a gene..." = "")
        
        # Get available genes from MAST data
        if (!is.null(values$de_data_mast)) {
          mast_genes <- unique(values$de_data_mast$gene)
          gene_choices <- c(gene_choices, setNames(mast_genes, paste0(mast_genes, " (MAST)")))
        }
        
        # Get available genes from MixScale data
        if (!is.null(values$de_data_mixscale)) {
          mixscale_genes <- unique(values$de_data_mixscale$gene)
          gene_choices <- c(gene_choices, setNames(mixscale_genes, paste0(mixscale_genes, " (MixScale)")))
        }
        
        updateSelectInput(session, "cross_cluster_gene", choices = gene_choices)
      }
    })
    
    # Cross-cluster analysis table
    output$cross_cluster_table <- DT::renderDataTable({
      req(input$cross_cluster_gene)
      req(input$cross_cluster_gene != "")
      
      # Determine if this is MAST or MixScale gene
      selected_gene <- input$cross_cluster_gene
      
      cross_cluster_data <- NULL
      
      # Process MAST data
      if (!is.null(values$de_data_mast)) {
        mast_subset <- values$de_data_mast[values$de_data_mast$gene == selected_gene, ]
        
        if (nrow(mast_subset) > 0) {
          cross_cluster_data <- mast_subset %>%
            filter(p_val_adj < input$pvalue_threshold) %>%
            group_by(gene_name) %>%
            summarise(
              clusters_affected = n(),
              clusters_list = paste(unique(cluster), collapse = ", "),
              mean_log2FC = round(mean(log2FC, na.rm = TRUE), 3),
              min_pvalue = min(p_val_adj, na.rm = TRUE),
              method = "MAST",
              .groups = "drop"
            ) %>%
            filter(clusters_affected >= input$min_clusters) %>%
            arrange(desc(clusters_affected), min_pvalue)
        }
      }
      
      # Process MixScale data
      if (!is.null(values$de_data_mixscale) && is.null(cross_cluster_data)) {
        mixscale_subset <- values$de_data_mixscale[values$de_data_mixscale$gene == selected_gene, ]
        
        if (nrow(mixscale_subset) > 0) {
          cross_cluster_data <- mixscale_subset %>%
            filter(p_val_adj < input$pvalue_threshold) %>%
            group_by(gene_name) %>%
            summarise(
              clusters_affected = n(),
              clusters_list = paste(unique(cluster), collapse = ", "),
              mean_log2FC = round(mean(log2FC, na.rm = TRUE), 3),
              min_pvalue = min(p_val_adj, na.rm = TRUE),
              method = "MixScale",
              .groups = "drop"
            ) %>%
            filter(clusters_affected >= input$min_clusters) %>%
            arrange(desc(clusters_affected), min_pvalue)
        }
      }
      
      if (is.null(cross_cluster_data) || nrow(cross_cluster_data) == 0) {
        return(DT::datatable(
          data.frame(Message = "No genes found meeting the criteria"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      # Format for display (with safe column selection)
      available_cols <- c("gene_name", "clusters_affected", "mean_log2FC", "min_pvalue", "clusters_list")
      if ("method" %in% colnames(cross_cluster_data)) {
        available_cols <- c(available_cols, "method")
      }
      
      display_data <- cross_cluster_data %>%
        select(all_of(available_cols))
      
      # Set appropriate column names
      if ("method" %in% colnames(display_data)) {
        colnames(display_data) <- c("Gene", "# Clusters", "Mean log2FC", "Min p-value", "Clusters", "Method")
      } else {
        colnames(display_data) <- c("Gene", "# Clusters", "Mean log2FC", "Min p-value", "Clusters")
      }
      
      DT::datatable(display_data,
                   options = list(
                     pageLength = 15,
                     scrollX = TRUE,
                     order = list(list(1, 'desc'))  # Sort by # Clusters descending
                   ),
                   rownames = FALSE) %>%
        DT::formatRound(c("Mean log2FC"), digits = 3) %>%
        DT::formatSignif(c("Min p-value"), digits = 3)
    })
    
    # Cross-cluster heatmap
    output$cross_cluster_heatmap <- renderPlot({
      req(input$cross_cluster_gene)
      req(input$cross_cluster_gene != "")
      
      selected_gene <- input$cross_cluster_gene
      
      # Get data for selected gene
      plot_data <- NULL
      
      # Process MAST data
      if (!is.null(values$de_data_mast)) {
        mast_subset <- values$de_data_mast[values$de_data_mast$gene == selected_gene, ]
        
        if (nrow(mast_subset) > 0) {
          plot_data <- mast_subset %>%
            filter(p_val_adj < input$pvalue_threshold) %>%
            select(gene_name, cluster, log2FC, p_val_adj) %>%
            mutate(
              significant = p_val_adj < input$pvalue_threshold,
              log_p = -log10(p_val_adj)
            )
        }
      }
      
      # Process MixScale data
      if (!is.null(values$de_data_mixscale) && is.null(plot_data)) {
        mixscale_subset <- values$de_data_mixscale[values$de_data_mixscale$gene == selected_gene, ]
        
        if (nrow(mixscale_subset) > 0) {
          plot_data <- mixscale_subset %>%
            filter(p_val_adj < input$pvalue_threshold) %>%
            select(gene_name, cluster, log2FC, p_val_adj) %>%
            mutate(
              significant = p_val_adj < input$pvalue_threshold,
              log_p = -log10(p_val_adj)
            )
        }
      }
      
      if (is.null(plot_data) || nrow(plot_data) == 0) {
        # Create empty plot with message
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, label = "No significant genes found", size = 6) +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
      } else {
        # Filter to genes that appear in multiple clusters
        genes_multi_cluster <- plot_data %>%
          group_by(gene_name) %>%
          summarise(n_clusters = n_distinct(cluster), .groups = "drop") %>%
          filter(n_clusters >= input$min_clusters) %>%
          pull(gene_name)
        
        plot_data_filtered <- plot_data %>%
          filter(gene_name %in% genes_multi_cluster)
        
        if (nrow(plot_data_filtered) == 0) {
          ggplot() +
            annotate("text", x = 0.5, y = 0.5, label = "No genes meet minimum cluster criteria", size = 6) +
            theme_void() +
            xlim(0, 1) + ylim(0, 1)
        } else {
          # Create heatmap
          ggplot(plot_data_filtered, aes(x = cluster, y = gene_name, fill = log2FC)) +
            geom_tile(color = "white", size = 0.5) +
            scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                               name = "log2FC") +
            theme_minimal() +
            theme(
              axis.text.x = element_text(angle = 45, hjust = 1),
              axis.text.y = element_text(size = 8),
              panel.grid = element_blank()
            ) +
            labs(
              title = paste("Cross-Cluster DE Status:", selected_gene),
              x = "Cluster",
              y = "Gene",
              subtitle = paste("Genes significant in ≥", input$min_clusters, "clusters")
            )
        }
      }
    })
    
    # === DE GENE OVERLAP HEATMAP LOGIC ===
    
    # Populate cluster choices for heatmap
    observe({
      if (!is.null(values$de_data_mast) || !is.null(values$de_data_mixscale)) {
        cluster_choices <- c("Choose a cluster..." = "")
        
        # Get available clusters
        if (!is.null(values$de_data_mast)) {
          mast_clusters <- unique(values$de_data_mast$cluster)
          cluster_choices <- c(cluster_choices, mast_clusters)
        }
        
        if (!is.null(values$de_data_mixscale)) {
          mixscale_clusters <- unique(values$de_data_mixscale$cluster)
          cluster_choices <- c(cluster_choices, unique(c(cluster_choices, mixscale_clusters)))
        }
        
        updateSelectInput(session, "heatmap_cluster", choices = cluster_choices)
      }
    })
    
    # DE gene overlap heatmap
    output$de_overlap_heatmap <- renderPlot({
      req(input$heatmap_cluster)
      req(input$heatmap_cluster != "")
      
      selected_cluster <- input$heatmap_cluster
      selected_method <- input$heatmap_method
      selected_direction <- input$heatmap_direction
      min_genes <- input$min_genes_overlap
      
      # Get data for selected cluster and method with direction filtering
      gene_lists <- list()
      
      if (selected_method %in% c("both", "mast") && !is.null(values$de_data_mast)) {
        mast_cluster_data <- values$de_data_mast[values$de_data_mast$cluster == selected_cluster, ]
        
        if (nrow(mast_cluster_data) > 0) {
          for (gene in unique(mast_cluster_data$gene)) {
            gene_subset <- mast_cluster_data[mast_cluster_data$gene == gene, ]
            
            # Apply direction filtering to prevent inflation
            if (selected_direction == "UP") {
              significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC > 0, ]$gene_name
            } else if (selected_direction == "DOWN") {
              significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC < 0, ]$gene_name
            } else {  # ALL - deduplicated
              significant_genes <- unique(gene_subset[gene_subset$p_val_adj < 0.05, ]$gene_name)
            }
            
            if (length(significant_genes) >= min_genes) {
              gene_lists[[paste0(gene, " (MAST)")]] <- significant_genes
            }
          }
        }
      }
      
      if (selected_method %in% c("both", "mixscale") && !is.null(values$de_data_mixscale)) {
        mixscale_cluster_data <- values$de_data_mixscale[values$de_data_mixscale$cluster == selected_cluster, ]
        
        if (nrow(mixscale_cluster_data) > 0) {
          for (gene in unique(mixscale_cluster_data$gene)) {
            gene_subset <- mixscale_cluster_data[mixscale_cluster_data$gene == gene, ]
            
            # Apply direction filtering to prevent inflation
            if (selected_direction == "UP") {
              significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC > 0, ]$gene_name
            } else if (selected_direction == "DOWN") {
              significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC < 0, ]$gene_name
            } else {  # ALL - deduplicated
              significant_genes <- unique(gene_subset[gene_subset$p_val_adj < 0.05, ]$gene_name)
            }
            
            if (length(significant_genes) >= min_genes) {
              gene_lists[[paste0(gene, " (MixScale)")]] <- significant_genes
            }
          }
        }
      }
      
      if (length(gene_lists) < 2) {
        # Provide better user feedback about available data
        feedback_text <- paste0(
          "Need at least 2 gene sets with ≥", min_genes, " DE genes\n",
          "Found: ", length(gene_lists), " gene set(s)\n\n",
          "Try:\n• Lower minimum genes threshold\n",
          "• Change direction filter\n",
          "• Select 'both' methods\n",
          "• Choose a different cluster"
        )
        
        ggplot() +
          annotate("text", x = 0.5, y = 0.5, 
                  label = feedback_text, size = 4, hjust = 0.5, vjust = 0.5) +
          theme_void() +
          xlim(0, 1) + ylim(0, 1)
      } else {
        # Calculate Fisher's exact test matrix
        gene_names <- names(gene_lists)
        n_genes <- length(gene_names)
        
        # Get all unique genes for background
        all_genes <- unique(unlist(gene_lists))
        background_size <- length(all_genes)
        
        # Initialize matrices
        p_matrix <- matrix(1, nrow = n_genes, ncol = n_genes)
        overlap_matrix <- matrix(0, nrow = n_genes, ncol = n_genes)
        
        rownames(p_matrix) <- colnames(p_matrix) <- gene_names
        rownames(overlap_matrix) <- colnames(overlap_matrix) <- gene_names
        
        # Calculate pairwise Fisher's tests
        for (i in 1:n_genes) {
          for (j in 1:n_genes) {
            if (i != j) {
              genes1 <- gene_lists[[i]]
              genes2 <- gene_lists[[j]]
              
              overlap <- intersect(genes1, genes2)
              overlap_count <- length(overlap)
              
              overlap_matrix[i, j] <- overlap_count
              
              if (overlap_count > 0) {
                # Fisher's exact test
                genes1_only <- length(genes1) - overlap_count
                genes2_only <- length(genes2) - overlap_count
                neither <- background_size - length(genes1) - length(genes2) + overlap_count
                
                if (genes1_only >= 0 && genes2_only >= 0 && neither >= 0) {
                  contingency_matrix <- matrix(c(overlap_count, genes1_only, genes2_only, neither), 
                                             nrow = 2, byrow = TRUE)
                  
                  tryCatch({
                    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
                    p_matrix[i, j] <- fisher_result$p.value
                  }, error = function(e) {
                    p_matrix[i, j] <- 1
                  })
                }
              }
            } else {
              overlap_matrix[i, j] <- length(gene_lists[[i]])
            }
          }
        }
        
        # Convert to data frame for ggplot
        p_matrix_long <- expand.grid(Gene1 = gene_names, Gene2 = gene_names)
        p_matrix_long$p_value <- as.vector(p_matrix)
        p_matrix_long$overlap <- as.vector(overlap_matrix)
        p_matrix_long$neg_log_p <- -log10(p_matrix_long$p_value)
        
        # Cap extremely high values for visualization
        p_matrix_long$neg_log_p[p_matrix_long$neg_log_p > 10] <- 10
        
        # Create enhanced heatmap with significance threshold visualization
        significance_threshold <- -log10(0.05)  # 1.301
        
        ggplot(p_matrix_long, aes(x = Gene1, y = Gene2, fill = neg_log_p)) +
          geom_tile(color = "white", size = 0.5) +
          geom_text(aes(label = overlap), color = "black", size = 3) +
          # Enhanced color scale with discrete bins for significance
          scale_fill_gradientn(
            colors = c("#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
            values = c(0, 0.13, 0.3, 0.5, 0.7, 1),  # Custom breaks emphasizing significance threshold
            name = "-log10(p)", 
            limits = c(0, 10),
            breaks = c(0, round(significance_threshold, 1), 2, 5, 10),
            labels = c("0 (p=1)", paste0(round(significance_threshold, 1), " (p=0.05)"), "2", "5", "10+")
          ) +
          theme_minimal() +
          theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_blank(),
            panel.grid = element_blank(),
            legend.position = "right"
          ) +
          labs(
            title = paste("DE Gene Overlap Significance -", selected_cluster, "-", selected_direction),
            subtitle = paste0("Fisher's exact test p-values (numbers show overlap counts) - Direction-filtered to prevent inflation\n",
                             "Significance threshold: p=0.05 (-log10 = ", round(significance_threshold, 2), ")")
          ) +
          coord_equal()
      }
    })
    
    # Download handler for overlap matrix
    output$download_overlap_matrix <- downloadHandler(
      filename = function() {
        paste0("de_overlap_matrix_", input$heatmap_cluster, "_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(input$heatmap_cluster)
        req(input$heatmap_cluster != "")
        
        # Generate the same matrix data as in the plot
        selected_cluster <- input$heatmap_cluster
        selected_method <- input$heatmap_method
        selected_direction <- input$heatmap_direction
        min_genes <- input$min_genes_overlap
        
        # Recreate gene lists (same logic as plot with direction filtering)
        gene_lists <- list()
        
        if (selected_method %in% c("both", "mast") && !is.null(values$de_data_mast)) {
          mast_cluster_data <- values$de_data_mast[values$de_data_mast$cluster == selected_cluster, ]
          
          if (nrow(mast_cluster_data) > 0) {
            for (gene in unique(mast_cluster_data$gene)) {
              gene_subset <- mast_cluster_data[mast_cluster_data$gene == gene, ]
              
              # Apply direction filtering to prevent inflation
              if (selected_direction == "UP") {
                significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC > 0, ]$gene_name
              } else if (selected_direction == "DOWN") {
                significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC < 0, ]$gene_name
              } else {  # ALL - deduplicated
                significant_genes <- unique(gene_subset[gene_subset$p_val_adj < 0.05, ]$gene_name)
              }
              
              if (length(significant_genes) >= min_genes) {
                gene_lists[[paste0(gene, " (MAST)")]] <- significant_genes
              }
            }
          }
        }
        
        if (selected_method %in% c("both", "mixscale") && !is.null(values$de_data_mixscale)) {
          mixscale_cluster_data <- values$de_data_mixscale[values$de_data_mixscale$cluster == selected_cluster, ]
          
          if (nrow(mixscale_cluster_data) > 0) {
            for (gene in unique(mixscale_cluster_data$gene)) {
              gene_subset <- mixscale_cluster_data[mixscale_cluster_data$gene == gene, ]
              
              # Apply direction filtering to prevent inflation
              if (selected_direction == "UP") {
                significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC > 0, ]$gene_name
              } else if (selected_direction == "DOWN") {
                significant_genes <- gene_subset[gene_subset$p_val_adj < 0.05 & gene_subset$log2FC < 0, ]$gene_name
              } else {  # ALL - deduplicated
                significant_genes <- unique(gene_subset[gene_subset$p_val_adj < 0.05, ]$gene_name)
              }
              
              if (length(significant_genes) >= min_genes) {
                gene_lists[[paste0(gene, " (MixScale)")]] <- significant_genes
              }
            }
          }
        }
        
        if (length(gene_lists) >= 2) {
          # Create the overlap matrix
          gene_names <- names(gene_lists)
          n_genes <- length(gene_names)
          all_genes <- unique(unlist(gene_lists))
          background_size <- length(all_genes)
          
          overlap_data <- data.frame()
          
          for (i in 1:n_genes) {
            for (j in 1:n_genes) {
              genes1 <- gene_lists[[i]]
              genes2 <- gene_lists[[j]]
              overlap <- intersect(genes1, genes2)
              overlap_count <- length(overlap)
              
              # Calculate Fisher's test if different genes
              p_value <- 1
              if (i != j && overlap_count > 0) {
                genes1_only <- length(genes1) - overlap_count
                genes2_only <- length(genes2) - overlap_count
                neither <- background_size - length(genes1) - length(genes2) + overlap_count
                
                if (genes1_only >= 0 && genes2_only >= 0 && neither >= 0) {
                  contingency_matrix <- matrix(c(overlap_count, genes1_only, genes2_only, neither), 
                                             nrow = 2, byrow = TRUE)
                  
                  tryCatch({
                    fisher_result <- fisher.test(contingency_matrix, alternative = "greater")
                    p_value <- fisher_result$p.value
                  }, error = function(e) {
                    p_value <- 1
                  })
                }
              } else if (i == j) {
                p_value <- NA  # Self-comparison
              }
              
              overlap_data <- rbind(overlap_data, data.frame(
                Gene1 = gene_names[i],
                Gene2 = gene_names[j],
                Overlap_Count = overlap_count,
                P_Value = p_value,
                Gene1_Total = length(genes1),
                Gene2_Total = length(genes2)
              ))
            }
          }
          
          write.csv(overlap_data, file, row.names = FALSE)
        } else {
          # Write empty file with message
          write.csv(data.frame(Message = "Insufficient data for analysis"), file, row.names = FALSE)
        }
      }
    )
    
    # Return values for potential use by other modules
    return(list(
      selected_cluster = reactive({ values$selected_cluster }),
      de_data_mast = reactive({ values$de_data_mast }),
      de_data_mixscale = reactive({ values$de_data_mixscale })
    ))
    
  })
}