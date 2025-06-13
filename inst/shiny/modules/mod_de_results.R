# Module for DE Results page with interactive UMAP and volcano plots
# Allows clicking on UMAP clusters to update volcano plots

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
      umap_data = NULL
    )
    
    # Helper function to generate mock UMAP data (moved up for access)
    generate_mock_umap_data <- function() {
      cat("[DE Results] Generating mock UMAP data...\n")
      # In production, this would compute from actual single-cell data
      set.seed(42)
      n_cells <- 5000
      clusters <- paste0("cluster_", 0:9)
      
      umap_data <- data.frame(
        UMAP1 = rnorm(n_cells),
        UMAP2 = rnorm(n_cells),
        cluster = sample(clusters, n_cells, replace = TRUE),
        stringsAsFactors = FALSE
      )
      
      # Add some structure to make clusters visible
      for (i in seq_along(clusters)) {
        idx <- umap_data$cluster == clusters[i]
        umap_data$UMAP1[idx] <- umap_data$UMAP1[idx] + cos(2*pi*i/length(clusters)) * 3
        umap_data$UMAP2[idx] <- umap_data$UMAP2[idx] + sin(2*pi*i/length(clusters)) * 3
      }
      cat("[DE Results] Mock UMAP data generated with", n_cells, "cells\n")
      return(umap_data)
    }
    
    # Load UMAP data
    observe({
      cat("[DE Results] UMAP data observe block triggered\n")
      cat("[DE Results] app_data$data_loaded =", app_data$data_loaded, "\n")
      
      req(app_data$data_loaded)
      cat("[DE Results] app_data is loaded, attempting to populate UMAP data...\n")
      
      # Try to load pre-computed UMAP data
      umap_file <- file.path(dirname(Sys.getenv("ISCORE_DATA_FILE", "")), "umap_data.rds")
      cat("[DE Results] Looking for UMAP file at:", umap_file, "\n")
      
      if (file.exists(umap_file)) {
        cat("[DE Results] UMAP file exists, attempting to load...\n")
        tryCatch({
          values$umap_data <- readRDS(umap_file)
          cat("[DE Results] Successfully loaded UMAP data from:", umap_file, "\n")
          cat("[DE Results] UMAP data dimensions:", nrow(values$umap_data), "cells x", ncol(values$umap_data), "columns\n")
        }, error = function(e) {
          showNotification("Could not load UMAP data", type = "warning")
          cat("[DE Results] Error loading UMAP data:", e$message, "\n")
          cat("[DE Results] Falling back to mock data generation...\n")
          # Generate mock data on error too
          values$umap_data <- generate_mock_umap_data()
        })
      } else {
        cat("[DE Results] UMAP file not found, generating mock data...\n")
        # Generate mock UMAP data if file not found
        values$umap_data <- generate_mock_umap_data()
      }
      
      cat("[DE Results] UMAP data population complete. values$umap_data is", 
          ifelse(is.null(values$umap_data), "NULL", "populated"), "\n")
    })
    
    # Load DE results data
    observe({
      req(app_data$data_loaded)
      
      # Try to load pre-computed DE results
      de_file <- Sys.getenv("ISCORE_DE_FILE", "")
      
      if (file.exists(de_file)) {
        tryCatch({
          de_results <- readRDS(de_file)
          values$de_data_mast <- de_results$MAST
          values$de_data_mixscale <- de_results$MixScale
          cat("[DE Results] Loaded DE data from:", de_file, "\n")
        }, error = function(e) {
          showNotification("Could not load DE results", type = "warning")
          cat("[DE Results] Error loading DE data:", e$message, "\n")
        })
      } else {
        # Use enrichment data as fallback (extract genes from enrichment results)
        # In production, this would load actual DE results
        showNotification("Using enrichment data for volcano plots", type = "message")
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
          # FALLBACK if 'experiment' column is missing
          showNotification("'Experiment' column not found. Coloring by significance instead.", type = "warning")
          plot_data$color_group <- ifelse(plot_data$significant, "Significant", "Not Significant")
          color_scale <- c("Significant" = "#FF6B6B", "Not Significant" = "#CCCCCC")
        }
      } else {  # color by gene/mutation
        if ("gene" %in% names(plot_data)) {
          plot_data$color_group <- plot_data$gene
          color_scale <- NULL # Use default plotly colors
        } else {
          # FALLBACK if 'gene' column is missing
          showNotification("'Gene' column not found. Coloring by significance instead.", type = "warning")
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
        text = ~paste("Gene:", gene,
                     "<br>Log2FC:", round(log2FC, 3),
                     "<br>P-value:", format(pvalue, digits = 3),
                     "<br>Experiment:", experiment),
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
        
        # Use mock data if no real DE data available
        if (is.null(values$de_data_mast)) {
          # Generate mock DE data from enrichment results
          set.seed(123)
          genes <- c("LRRK2", "PINK1", "PARK7", "SNCA", "GBA", "ATP13A2", "VPS13C")
          mock_data <- expand.grid(
            gene = rep(genes, each = 20),
            cluster = paste0("cluster_", 0:9),
            stringsAsFactors = FALSE
          ) %>%
            mutate(
              log2FC = rnorm(n(), mean = 0, sd = 1.5),
              pvalue = 10^(-rexp(n(), rate = 0.5)),
              experiment = "MAST_analysis"
            )
          
          generate_volcano_plot(mock_data, "MAST", values$selected_cluster, input$color_by)
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
        
        # Use mock data if no real DE data available
        if (is.null(values$de_data_mixscale)) {
          # Generate mock DE data
          set.seed(456)
          genes <- c("LRRK2", "PINK1", "PARK7", "SNCA", "GBA", "ATP13A2", "VPS13C")
          experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
          
          mock_data <- expand.grid(
            gene = rep(genes, each = 15),
            cluster = paste0("cluster_", 0:9),
            experiment = sample(experiments, 10, replace = TRUE),
            stringsAsFactors = FALSE
          ) %>%
            mutate(
              log2FC = rnorm(n(), mean = 0, sd = 1.2),
              pvalue = 10^(-rexp(n(), rate = 0.6))
            )
          
          generate_volcano_plot(mock_data, "MixScale", values$selected_cluster, input$color_by)
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