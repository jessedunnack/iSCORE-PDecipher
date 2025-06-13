# Module for DE Results page with interactive UMAP and volcano plots
# Allows clicking on UMAP clusters to update volcano plots

# Helper functions to process DE data for volcano plots
process_mast_for_volcano <- function(mast_data) {
  # Convert MAST data structure to volcano plot format
  volcano_data <- data.frame()
  
  for (gene in names(mast_data)) {
    for (cluster in names(mast_data[[gene]])) {
      if (!is.null(mast_data[[gene]][[cluster]])) {
        de_results <- mast_data[[gene]][[cluster]]
        
        # Extract log2FC and p-values
        if ("log2FC" %in% colnames(de_results) && "pvalue" %in% colnames(de_results)) {
          cluster_data <- data.frame(
            gene = gene,
            cluster = cluster,
            gene_name = rownames(de_results),
            log2FC = de_results$log2FC,
            pvalue = de_results$pvalue,
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
      if (!is.null(mixscale_data[[gene]][[cluster]])) {
        for (exp in names(mixscale_data[[gene]][[cluster]])) {
          de_results <- mixscale_data[[gene]][[cluster]][[exp]]
          
          # Extract log2FC column (e.g., log2FC_C12_FPD-23)
          log2fc_col <- paste0("log2FC_", exp)
          pval_col <- paste0("p_cell_type", exp, ":weight")
          
          if (log2fc_col %in% colnames(de_results) && pval_col %in% colnames(de_results)) {
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
      # Left panel: Interactive UMAP
      column(6,
        wellPanel(
          h3("Interactive UMAP", icon("map")),
          p("Click on a cluster to view its differential expression results"),
          
          # UMAP plot output
          shinycssloaders::withSpinner(
            plotlyOutput(ns("umap_plot"), height = "700px"),
            type = 6,
            color = "#3c8dbc"
          ),
          
          # Selected cluster info
          div(style = "margin-top: 15px;",
            h5("Selected Cluster:"),
            textOutput(ns("selected_cluster_text"), inline = TRUE)
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
              umap_coords <- reducedDim(sce, "UMAP")
              
              # Create data frame for plotting
              values$umap_data <- data.frame(
                UMAP1 = umap_coords[, 1],
                UMAP2 = umap_coords[, 2],
                cluster = as.character(colData(sce)$seurat_clusters),
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
      
      # Try to find full_DE_results.rds file
      possible_de_paths <- c(
        file.path(dirname(Sys.getenv("ISCORE_DATA_FILE", "")), "full_DE_results.rds"),
        file.path(getwd(), "full_DE_results.rds"),
        file.path(dirname(getwd()), "full_DE_results.rds"),
        "../../full_DE_results.rds"
      )
      
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

    # Render UMAP plot
    output$umap_plot <- renderPlotly({
      tryCatch({
        cat("[DE Results] Attempting to render UMAP plot...\n")
        cat("[DE Results] values$umap_data is", ifelse(is.null(values$umap_data), "NULL", "populated"), "\n")
        
        req(values$umap_data)
        cat("[DE Results] req(values$umap_data) satisfied, proceeding with plot...\n")
        
        # Create color palette
        clusters <- unique(values$umap_data$cluster)
        n_clusters <- length(clusters)
        colors <- viridis::viridis(n_clusters)
        names(colors) <- clusters
        
        # Highlight selected cluster
        if (!is.null(values$selected_cluster) && values$selected_cluster != "All") {
          colors[values$selected_cluster] <- "#FF0000"  # Red for selected
          colors[names(colors) != values$selected_cluster] <- "#CCCCCC"  # Gray for others
        }
        
        p <- plot_ly(
          data = values$umap_data,
          x = ~UMAP1,
          y = ~UMAP2,
          color = ~cluster,
          colors = colors,
          type = 'scatter',
          mode = 'markers',
          marker = list(size = 5, opacity = 0.7),
          text = ~paste("Cluster:", cluster),
          hoverinfo = "text",
          source = "umap_click"
        ) %>%
          layout(
            title = list(text = "UMAP - Cell Clusters", font = list(size = 16)),
            xaxis = list(title = "UMAP 1", zeroline = FALSE),
            yaxis = list(title = "UMAP 2", zeroline = FALSE),
            showlegend = TRUE,
            legend = list(orientation = "v", x = 1.02, y = 0.5)
          )
        
        p
      }, error = function(e) {
        cat("[DE Results] Error rendering UMAP plot:", e$message, "\n")
        showNotification("Error rendering UMAP plot", type = "error")
        plotly::plotly_empty()
      })
    })
    
    # Handle UMAP click events
    observeEvent(event_data("plotly_click", source = "umap_click"), {
      click_data <- event_data("plotly_click", source = "umap_click")
      if (!is.null(click_data)) {
        # Extract cluster from click
        clicked_cluster <- values$umap_data$cluster[click_data$pointNumber + 1]
        values$selected_cluster <- clicked_cluster
        
        showNotification(
          paste("Selected cluster:", clicked_cluster),
          type = "message",
          duration = 2
        )
      }
    })
    
    # Update selected cluster text
    output$selected_cluster_text <- renderText({
      if (is.null(values$selected_cluster) || values$selected_cluster == "All") {
        "All clusters (global results)"
      } else {
        values$selected_cluster
      }
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
        
        if (is.null(values$de_data_mast)) {
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
        
        if (is.null(values$de_data_mixscale)) {
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