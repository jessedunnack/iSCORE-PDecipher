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
      
      if (has_crispri && has_mutations) {
        dataset_to_load <- "Full_Dataset"
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
          clusters <- sort(unique(values$umap_data$cluster))
          cluster_choices <- setNames(clusters, paste("Cluster", gsub("cluster_", "", clusters)))
          
          updateSelectInput(session, "cluster_selector",
                           choices = c("Choose a cluster to analyze..." = "", cluster_choices))
        }
      }
    })
    
    # DISABLED: Automatic DE results loading to prevent interference with other modules
    # Users must manually load volcano plot data when on the DE Results tab
    
    # # Look for full_DE_results.rds in the dataset directory
    # possible_de_paths <- c(
    #   file.path(data_dir, "full_DE_results.rds"),
    #   Sys.getenv("ISCORE_DE_FILE", ""),
    #   file.path(dirname(Sys.getenv("ISCORE_ENRICHMENT_DIR", "")), "full_DE_results.rds")
    # )
    # 
    # # Remove empty paths
    # possible_de_paths <- possible_de_paths[possible_de_paths != ""]
    # 
    # de_loaded <- FALSE
    # for (path in possible_de_paths) {
    #   cat("[DE Results] Checking DE path:", path, "\n")
    #   if (file.exists(path)) {
    #     tryCatch({
    #       de_results <- readRDS(path)
    #       cat("[DE Results] Successfully loaded DE results from:", path, "\n")
    #       cat("[DE Results] DE results structure:", paste(names(de_results), collapse=", "), "\n")
    #       
    #       # Extract MAST and MixScale data
    #       if ("iSCORE_PD_MAST" %in% names(de_results)) {
    #         # Convert MAST data to volcano plot format
    #         mast_data <- de_results$iSCORE_PD_MAST
    #         values$de_data_mast <- process_mast_for_volcano(mast_data)
    #         cat("[DE Results] Processed MAST data:", nrow(values$de_data_mast), "rows\n")
    #       cat("[DE Results] Available MAST genes:", paste(unique(values$de_data_mast$gene), collapse=", "), "\n")
    #       }
    #       
    #       if ("CRISPRi_Mixscale" %in% names(de_results)) {
    #         # Convert MixScale data to volcano plot format  
    #         mixscale_data <- de_results$CRISPRi_Mixscale
    #         values$de_data_mixscale <- process_mixscale_for_volcano(mixscale_data)
    #         cat("[DE Results] Processed MixScale data:", nrow(values$de_data_mixscale), "rows\n")
    #       cat("[DE Results] Available MixScale genes:", paste(unique(values$de_data_mixscale$gene), collapse=", "), "\n")
    #       }
    #       
    #       de_loaded <- TRUE
    #       break
    #     }, error = function(e) {
    #       cat("[DE Results] Failed to load DE results from", path, ":", e$message, "\n")
    #     })
    #   }
    # }
    # 
    # if (!de_loaded) {
    #   showNotification("DE results file not found. Volcano plots will not be available.", 
    #                  type = "warning",
    #                  duration = 5)
    # }
    # })
    
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
    
    # Dynamic MAST volcano plot container with adaptive height
    output$mast_volcano_container <- renderUI({
      ns <- session$ns
      has_mixscale <- !is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0
      
      # If no MixScale data, expand MAST plot to full height (700px)
      # If MixScale data available, use standard height (350px)
      plot_height <- if (has_mixscale) "350px" else "700px"
      
      plotlyOutput(ns("mast_volcano"), height = plot_height)
    })
    
    # Dynamic MixScale volcano plot container - ONLY renders if data available
    output$mixscale_volcano_container <- renderUI({
      ns <- session$ns
      has_mixscale <- !is.null(values$de_data_mixscale) && nrow(values$de_data_mixscale) > 0
      
      # Only render MixScale plot if data is available
      if (has_mixscale) {
        div(
          h4("MixScale Results"),
          shinycssloaders::withSpinner(
            plotlyOutput(ns("mixscale_volcano"), height = "350px"),
            type = 6,
            color = "#3c8dbc"
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
    
    # Generate volcano plot function - BULLETPROOF: Always returns valid plotly object
    generate_volcano_plot <- function(de_data, analysis_type, selected_cluster, color_by) {
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
        title_suffix <- paste("- Cluster", gsub("cluster_", "", selected_cluster))
      } else {
        plot_data <- de_data
        title_suffix <- "- All Clusters"
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
    
    # Render MAST volcano plot - CLEANED: Removed interfering cat() statements
    output$mast_volcano <- renderPlotly({
      tryCatch({
        
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
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            filtered_data <- values$de_data_mast[values$de_data_mast$gene == current_gene, ]
          } else {
            filtered_data <- values$de_data_mast
          }
          
          generate_volcano_plot(filtered_data, "MAST", values$selected_cluster, input$color_by)
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
    
    # Render MixScale volcano plot - CLEANED: Removed interfering cat() statements
    output$mixscale_volcano <- renderPlotly({
      tryCatch({
        
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
          
          if (!is.null(current_gene) && current_gene != "" && current_gene != "All") {
            filtered_data <- values$de_data_mixscale[values$de_data_mixscale$gene == current_gene, ]
          } else {
            filtered_data <- values$de_data_mixscale
          }
          
          generate_volcano_plot(filtered_data, "MixScale", values$selected_cluster, input$color_by)
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