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

process_mixscale_for_volcano <- function(mixscale_data) {
  # Convert MixScale data structure to volcano plot format
  volcano_data <- data.frame()
  
  for (gene in names(mixscale_data)) {
    for (cluster in names(mixscale_data[[gene]])) {
      if (!is.null(mixscale_data[[gene]][[cluster]]$results)) {
        de_results <- mixscale_data[[gene]][[cluster]]$results
        
        # Find log2FC and p-value columns
        log2fc_cols <- grep("^log2FC_", names(de_results), value = TRUE)
        
        if (length(log2fc_cols) > 0) {
          # Use the first log2FC column and corresponding p-value
          log2fc_col <- log2fc_cols[1]
          # Extract experiment name from column
          exp <- gsub("^log2FC_", "", log2fc_col)
          pval_col <- paste0("p_cell_type", exp, ":weight")
          
          if (pval_col %in% colnames(de_results)) {
            cluster_data <- data.frame(
              gene = gene,
              cluster = cluster,
              gene_name = rownames(de_results),
              log2FC = de_results[[log2fc_col]],
              pvalue = de_results[[pval_col]],
              experiment = exp,
              stringsAsFactors = FALSE
            )
            volcano_data <- rbind(volcano_data, cluster_data)
          }
        }
      }
    }
  }
  
  return(volcano_data)
}

# UI function
mod_de_results_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      # Left panel: UMAP with cluster selection
      column(6,
        wellPanel(
          h3("Cell Cluster Analysis", icon("object-group")),
          
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
          
          # Color options
          div(style = "margin-bottom: 15px;",
            radioButtons(ns("color_by"),
                        "Color Points By:",
                        choices = c("Significance" = "significance",
                                  "Experiment" = "experiment",
                                  "Gene/Mutation" = "gene"),
                        selected = "significance",
                        inline = TRUE)
          ),
          
          # MAST volcano plot
          div(style = "margin-bottom: 20px;",
            h4("MAST Results"),
            shinycssloaders::withSpinner(
              plotlyOutput(ns("mast_volcano"), height = "350px"),
              type = 6,
              color = "#3c8dbc"
            )
          ),
          
          # MixScale volcano plot
          div(
            h4("MixScale Results"),
            shinycssloaders::withSpinner(
              plotlyOutput(ns("mixscale_volcano"), height = "350px"),
              type = 6,
              color = "#3c8dbc"
            )
          )
        )
      )
    ),
    
    # Bottom panel: Summary statistics
    fluidRow(
      column(12,
        wellPanel(
          h4("Summary Statistics"),
          div(id = ns("summary_stats"),
            uiOutput(ns("stats_content"))
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
    
    # Load UMAP data using the same approach as the landing page
    observe({
      cat("[DE Results] UMAP data observe block triggered\n")
      cat("[DE Results] app_data$data_loaded =", app_data$data_loaded, "\n")
      
      req(app_data$data_loaded)
      cat("[DE Results] app_data is loaded, attempting to load UMAP data...\n")
      
      # Determine which dataset to load based on app data
      has_crispri <- any(grepl("MixScale", app_data$consolidated_data$method))
      has_mutations <- any(grepl("MAST", app_data$consolidated_data$method))
      
      if (has_crispri && has_mutations) {
        dataset_to_load <- "Full_Dataset"
      } else if (has_crispri) {
        dataset_to_load <- "iSCORE_PD_CRISPRi"
      } else {
        dataset_to_load <- "iSCORE_PD"
      }
      
      cat("[DE Results] Determined dataset to load:", dataset_to_load, "\n")
      
      # Try to load the appropriate UMAP data
      possible_paths <- c(
        system.file("extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds"), 
                    package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds")),
        paste0("../../inst/extdata/umap_data/", dataset_to_load, "_umap_data.rds")
      )
      
      data_loaded <- FALSE
      for (path in possible_paths) {
        cat("[DE Results] Checking path:", path, "\n")
        if (file.exists(path)) {
          tryCatch({
            sce <- readRDS(path)
            cat("[DE Results] Successfully loaded UMAP SCE data from:", path, "\n")
            
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
              data_loaded <- TRUE
            }
            break
          }, error = function(e) {
            cat("[DE Results] Failed to load from", path, ":", e$message, "\n")
          })
        }
      }
      
      if (!data_loaded) {
        showNotification("UMAP data not found. Please run extract_umap_data.R first.", 
                       type = "error",
                       duration = NULL)
      }
      
      cat("[DE Results] UMAP data population complete. values$umap_data is", 
          ifelse(is.null(values$umap_data), "NULL", "populated"), "\n")
    })
    
    # Load DE results data from full_DE_results.rds
    observe({
      req(app_data$data_loaded)
      cat("[DE Results] Loading DE results data...\n")
      
      # Get the dataset directory from environment
      data_dir <- dirname(Sys.getenv("ISCORE_DATA_FILE", ""))
      
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
            }
            
            if ("CRISPRi_Mixscale" %in% names(de_results)) {
              # Convert MixScale data to volcano plot format  
              mixscale_data <- de_results$CRISPRi_Mixscale
              values$de_data_mixscale <- process_mixscale_for_volcano(mixscale_data)
              cat("[DE Results] Processed MixScale data:", nrow(values$de_data_mixscale), "rows\n")
            }
            
            de_loaded <- TRUE
            break
          }, error = function(e) {
            cat("[DE Results] Failed to load DE results from", path, ":", e$message, "\n")
          })
        }
      }
      
      if (!de_loaded) {
        showNotification("DE results file not found. Volcano plots will not be available.", 
                       type = "warning",
                       duration = 5)
      }
    })
    
    # Initialize with global cluster or "All"
    observe({
      if (is.null(values$selected_cluster)) {
        values$selected_cluster <- isolate(global_selection()$cluster)
        if (is.null(values$selected_cluster) || values$selected_cluster == "") {
          values$selected_cluster <- "All"
        }
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
      
      clusters <- sort(unique(values$umap_data$cluster))
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
        values$selected_cluster <- input$cluster_selector
      } else {
        values$selected_cluster <- NULL
      }
    })
    
    # Render UMAP plot with ggplot2
    output$umap_plot <- renderPlot({
      req(values$umap_data)
      
      # Get cluster colors
      clusters <- sort(unique(values$umap_data$cluster))
      n_clusters <- length(clusters)
      ditto_colors <- get_ditto_colors(n_clusters)
      names(ditto_colors) <- clusters
      
      # Create display data with highlighting
      plot_data <- values$umap_data
      
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
          legend.position = "none",
          panel.grid = element_blank(),
          axis.text = element_blank(),
          axis.title = element_text(size = 12),
          plot.background = element_rect(fill = "white", color = NA)
        ) +
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
        # Show all labels if not too many clusters
        p <- p + 
          geom_text(
            data = cluster_centers,
            aes(x = x, y = y, label = gsub("cluster_", "", label)),
            size = 4,
            fontface = "bold"
          )
      }
      
      p
    }, height = 600, width = 600)
    
    # Render cluster information
    output$cluster_info <- renderUI({
      req(input$cluster_selector)
      req(values$umap_data)
      
      cluster_cells <- values$umap_data %>%
        filter(cluster == input$cluster_selector)
      
      n_cells <- nrow(cluster_cells)
      pct_cells <- round(100 * n_cells / nrow(values$umap_data), 1)
      
      # Get DE summary if available
      de_summary <- "Calculating..."
      if (!is.null(values$de_data_mast) || !is.null(values$de_data_mixscale)) {
        mast_de <- if (!is.null(values$de_data_mast)) {
          sum(values$de_data_mast$cluster == input$cluster_selector & 
              values$de_data_mast$pvalue < 0.05)
        } else 0
        
        mixscale_de <- if (!is.null(values$de_data_mixscale)) {
          sum(values$de_data_mixscale$cluster == input$cluster_selector & 
              values$de_data_mixscale$pvalue < 0.05)
        } else 0
        
        de_summary <- paste(mast_de + mixscale_de, "DE genes")
      }
      
      tagList(
        tags$div(style = "display: flex; justify-content: space-between;",
          tags$div(
            tags$strong("Cells: "),
            tags$span(format(n_cells, big.mark = ","), 
                     paste0(" (", pct_cells, "%)")
            )
          ),
          tags$div(
            tags$strong("DE genes: "),
            tags$span(de_summary)
          )
        )
      )
    })
    
    # Generate volcano plot function
    generate_volcano_plot <- function(de_data, analysis_type, selected_cluster, color_by) {
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
      
      # Filter by cluster if selected
      if (!is.null(selected_cluster) && selected_cluster != "All") {
        plot_data <- de_data[de_data$cluster == selected_cluster, ]
        title_suffix <- paste("- Cluster", gsub("cluster_", "", selected_cluster))
      } else {
        plot_data <- de_data
        title_suffix <- "- All Clusters"
      }
      
      if (nrow(plot_data) == 0) {
        # No data for this cluster
        p <- plot_ly() %>%
          layout(
            title = paste(analysis_type, title_suffix, "- No data"),
            xaxis = list(title = "Log2 Fold Change", range = c(-5, 5)),
            yaxis = list(title = "-Log10 P-value", range = c(0, 10)),
            annotations = list(
              x = 0,
              y = 5,
              text = "No data for selected cluster",
              showarrow = FALSE,
              font = list(size = 16, color = "gray")
            )
          )
        return(p)
      }
      
      # Calculate -log10 p-value
      plot_data$negLog10p <- -log10(plot_data$pvalue + 1e-300)  # Add small value to avoid log(0)
      
      # Determine significance
      plot_data$significant <- plot_data$pvalue < 0.05 & abs(plot_data$log2FC) > 1
      
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
          title = list(text = paste(analysis_type, "Volcano Plot", title_suffix), 
                      font = list(size = 14)),
          xaxis = list(title = "Log2 Fold Change", zeroline = TRUE),
          yaxis = list(title = "-Log10 P-value", zeroline = FALSE),
          showlegend = TRUE,
          legend = list(orientation = "v", x = 1.02, y = 0.5)
        ) %>%
        # Add threshold lines
        add_trace(
          x = c(-1, -1),
          y = c(0, max(plot_data$negLog10p, na.rm = TRUE)),
          type = "scatter",
          mode = "lines",
          line = list(color = "gray", dash = "dash", width = 1),
          showlegend = FALSE,
          hoverinfo = "skip"
        ) %>%
        add_trace(
          x = c(1, 1),
          y = c(0, max(plot_data$negLog10p, na.rm = TRUE)),
          type = "scatter",
          mode = "lines",
          line = list(color = "gray", dash = "dash", width = 1),
          showlegend = FALSE,
          hoverinfo = "skip"
        ) %>%
        add_trace(
          x = c(min(plot_data$log2FC, na.rm = TRUE), max(plot_data$log2FC, na.rm = TRUE)),
          y = c(-log10(0.05), -log10(0.05)),
          type = "scatter",
          mode = "lines",
          line = list(color = "gray", dash = "dash", width = 1),
          showlegend = FALSE,
          hoverinfo = "skip"
        )
      
      p
    }
    
    # Render MAST volcano plot
    output$mast_volcano <- renderPlotly({
      tryCatch({
        cat("[DE Results] Attempting to render MAST volcano plot...\n")
        cat("[DE Results] values$selected_cluster =", values$selected_cluster, "\n")
        cat("[DE Results] values$de_data_mast is", ifelse(is.null(values$de_data_mast), "NULL", "populated"), "\n")
        
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
                text = "Please select a cluster to view differential expression results",
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
          generate_volcano_plot(values$de_data_mast, "MAST", values$selected_cluster, input$color_by)
        }
      }, error = function(e) {
        cat("[DE Results] Error rendering MAST volcano plot:", e$message, "\n")
        showNotification("Error rendering MAST volcano plot", type = "error")
        plotly::plotly_empty()
      })
    })
    
    # Render MixScale volcano plot
    output$mixscale_volcano <- renderPlotly({
      tryCatch({
        cat("[DE Results] Attempting to render MixScale volcano plot...\n")
        cat("[DE Results] values$selected_cluster =", values$selected_cluster, "\n")
        cat("[DE Results] values$de_data_mixscale is", ifelse(is.null(values$de_data_mixscale), "NULL", "populated"), "\n")
        
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
                text = "Please select a cluster to view differential expression results",
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
          generate_volcano_plot(values$de_data_mixscale, "MixScale", values$selected_cluster, input$color_by)
        }
      }, error = function(e) {
        cat("[DE Results] Error rendering MixScale volcano plot:", e$message, "\n")
        showNotification("Error rendering MixScale volcano plot", type = "error")
        plotly::plotly_empty()
      })
    })
    
    # Render summary statistics
    output$stats_content <- renderUI({
      cluster_text <- if (is.null(values$selected_cluster) || values$selected_cluster == "All") {
        "all clusters"
      } else {
        values$selected_cluster
      }
      
      # Calculate stats (mock for now)
      mast_sig <- 156  # Would calculate from actual data
      mixscale_sig <- 203
      overlap <- 47
      
      tagList(
        fluidRow(
          column(4,
            div(class = "text-center",
              h5("MAST Significant"),
              h3(mast_sig, style = "color: #3c8dbc;"),
              p("(p < 0.05, |log2FC| > 1)")
            )
          ),
          column(4,
            div(class = "text-center",
              h5("MixScale Significant"),
              h3(mixscale_sig, style = "color: #5cb85c;"),
              p("(p < 0.05, |log2FC| > 1)")
            )
          ),
          column(4,
            div(class = "text-center",
              h5("Overlapping Genes"),
              h3(overlap, style = "color: #f0ad4e;"),
              p("(significant in both)")
            )
          )
        ),
        hr(),
        p(paste("Statistics for", cluster_text), class = "text-muted text-center")
      )
    })
    
    # Return values for potential use by other modules
    return(list(
      selected_cluster = reactive({ values$selected_cluster }),
      de_data_mast = reactive({ values$de_data_mast }),
      de_data_mixscale = reactive({ values$de_data_mixscale })
    ))
    
  })
}