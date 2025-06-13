# Landing Page Module with UMAP Visualization (Version 2)
# Compact layout with UMAP on left, summary boxes on right
# UMAP dataset determined by loaded data, not user selection

# Source the UMAP viewer module
source("modules/mod_umap_viewer.R")

#' Landing Page UI with UMAP (Version 2)
#' 
#' @param id Module namespace
landingPageWithUmapUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Add compact styling for the markers section
    tags$style(HTML(paste0("
      #", ns(""), " .form-group {
        margin-bottom: 8px !important;
      }
      #", ns(""), " .control-label {
        margin-bottom: 3px !important;
        font-size: 12px !important;
      }
      #", ns("markers_table"), " .dataTables_wrapper {
        padding: 0px !important;
      }
      #", ns("markers_table"), " table {
        font-size: 11px !important;
      }
    "))),
    # Main content area with two columns - MAXIMIZED FOR SCREEN SPACE
    fluidRow(class = "landing-row-maximized", style = "min-height: 80vh;", # Use 80% of viewport height
      # Left column - UMAP visualization (expanded and taller)
      column(8, class = "umap-column-maximized",  # Increased from 7 to 8 for more width
        div(class = "box box-primary", style = "margin-top: 0; height: calc(80vh - 40px);", # Full height minus margins
          div(class = "box-header with-border",
            h3(class = "box-title", 
               icon("chart-line"),
               "Dataset UMAP Visualization")
          ),
          div(class = "box-body umap-container", style = "padding: 15px; height: calc(100% - 60px);", # Simplified styling
            div(class = "umap-maximized",
              withSpinner(
                plotOutput(ns("umap_plot"), 
                          height = "100%",  # Fill available height
                          width = "100%"),  # Fill available width
                type = 4,
                color = "#3c8dbc"
              )
            )
          )
        )
      ),
      
      # Right column - Summary statistics + Cluster Markers (reduced width but same height)
      column(4, class = "markers-column-optimized",  # Reduced from 5 to 4 to give UMAP more space
        # Summary statistics cards in a grid (top half) - Compact layout
        div(class = "value-box-compact",
          fluidRow(
            column(6,
              uiOutput(ns("total_cells_box"))
            ),
            column(6,
              uiOutput(ns("total_clusters_box"))
            )
          ),
          fluidRow(
            column(6,
              uiOutput(ns("total_results_box"))
            ),
            column(6,
              uiOutput(ns("total_genes_box"))
            )
          ),
          fluidRow(
            column(6,
              uiOutput(ns("total_experiments_box"))
            ),
            column(6,
              uiOutput(ns("enrichment_types_box"))
            )
          )
        ),
        
        # Cluster Markers section (fill remaining height)
        div(class = "box box-info", style = "margin-top: 20px; height: calc(80vh - 260px);", # Fill remaining space after metrics boxes
          div(class = "box-header with-border",
            h3(class = "box-title", 
               icon("dna"),
               "Cluster Marker Genes")
          ),
          div(class = "box-body", style = "padding: 10px; height: calc(100% - 60px); display: flex; flex-direction: column;",
            # Compact controls
            fluidRow(
              column(6,
                selectInput(ns("selected_cluster"),
                           "Select Cluster:",
                           choices = NULL,
                           width = "100%")
              ),
              column(6,
                numericInput(ns("max_markers"),
                            "Max Markers:",
                            value = 15,
                            min = 5,
                            max = 50,
                            step = 5,
                            width = "100%")
              )
            ),
            # Markers table - use flex-grow to fill remaining space
            div(class = "markers-table-optimized", style = "flex-grow: 1; overflow-y: auto; margin-top: 10px;",
              withSpinner(
                DT::dataTableOutput(ns("markers_table"), height = "100%"),
                type = 1,
                color = "#3c8dbc"
              )
            )
          )
        )
      )
    ),
    
    # Detailed breakdown - full width below
    fluidRow(
      # By Analysis Type
      column(4,
        div(class = "box box-primary",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Analysis Type")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("analysis_type_plot"), height = "350px"))
          )
        )
      ),
      
      # By Enrichment Type
      column(4,
        div(class = "box box-success",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Enrichment Database")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("enrichment_type_plot"), height = "350px"))
          )
        )
      ),
      
      # By Direction
      column(4,
        div(class = "box box-info",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Direction")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("direction_plot"), height = "350px"))
          )
        )
      )
    ),
    
    # Detailed tables
    fluidRow(
      column(12,
        div(class = "nav-tabs-custom",
          h3("Detailed Result Counts", style = "margin-top: 20px; margin-bottom: 20px;"),
          tabsetPanel(
            id = ns("detail_tabs"),
            
            # By Gene/Mutation
            tabPanel(
              "By Gene/Mutation",
              div(style = "margin-top: 15px;",
                DT::dataTableOutput(ns("gene_table"))
              )
            ),
            
            # By Cluster
            tabPanel(
              "By Cluster",
              div(style = "margin-top: 15px;",
                DT::dataTableOutput(ns("cluster_table"))
              )
            ),
            
            # Complete Matrix
            tabPanel(
              "Complete Matrix",
              div(style = "margin-top: 15px;",
                p("This table shows the number of significant terms for each combination of parameters."),
                p("Use the filters to explore specific combinations."),
                DT::dataTableOutput(ns("matrix_table"))
              )
            ),
            
            # Top Terms
            tabPanel(
              "Top Enriched Terms",
              div(style = "margin-top: 15px;",
                fluidRow(
                  column(4,
                    selectInput(ns("top_terms_method"),
                               "Filter by Method:",
                               choices = c("All Methods" = "all",
                                         "MAST" = "MAST",
                                         "MixScale" = "MixScale"),
                               selected = "all")
                  ),
                  column(4,
                    selectInput(ns("top_terms_gene"),
                               "Filter by Gene:",
                               choices = c("All Genes" = "all"),
                               selected = "all")
                  ),
                  column(4,
                    numericInput(ns("top_terms_n"),
                                "Number of top terms:",
                                value = 20,
                                min = 10,
                                max = 100,
                                step = 10)
                  )
                ),
                DT::dataTableOutput(ns("top_terms_table"))
              )
            )
          )
        )
      )
    )
  )
}

#' Landing Page Server with UMAP (Version 2)
#' 
#' @param id Module namespace
#' @param data Reactive data object from app
landingPageWithUmapServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for UMAP data and markers
    umap_data <- reactiveValues(
      sce = NULL,
      dataset_name = NULL,
      loaded = FALSE,
      markers = NULL
    )
    
    # Determine which UMAP dataset to load based on app data
    observe({
      req(data$data_loaded)
      
      # Determine dataset based on loaded data
      # Check for specific markers in the data
      has_crispri <- any(grepl("MixScale", data$consolidated_data$method))
      has_mutations <- any(grepl("MAST", data$consolidated_data$method))
      
      # Determine which dataset to load
      if (has_crispri && has_mutations) {
        dataset_to_load <- "Full_Dataset"
      } else if (has_crispri) {
        dataset_to_load <- "iSCORE_PD_CRISPRi"
      } else {
        dataset_to_load <- "iSCORE_PD"
      }
      
      # Try to load the appropriate UMAP data
      possible_paths <- c(
        system.file("extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds"), 
                    package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", paste0(dataset_to_load, "_umap_data.rds")),
        paste0("../../inst/extdata/umap_data/", dataset_to_load, "_umap_data.rds")
      )
      
      for (path in possible_paths) {
        if (file.exists(path)) {
          tryCatch({
            umap_data$sce <- readRDS(path)
            umap_data$dataset_name <- dataset_to_load
            umap_data$loaded <- TRUE
            
            # Try to load corresponding markers
            markers_path <- file.path(dirname(path), paste0(dataset_to_load, "_cluster_markers.rds"))
            if (file.exists(markers_path)) {
              umap_data$markers <- readRDS(markers_path)
              message("Loaded markers for ", dataset_to_load)
            } else {
              message("No markers found for ", dataset_to_load, " at ", markers_path)
            }
            
            break
          }, error = function(e) {
            message("Failed to load UMAP from ", path, ": ", e$message)
          })
        }
      }
    })
    
    # Render UMAP plot
    output$umap_plot <- renderPlot({
      if (!umap_data$loaded || is.null(umap_data$sce)) {
        # Placeholder when no data
        plot.new()
        text(0.5, 0.5, "UMAP data not available\nPlease run extract_umap_data.R", 
             cex = 1.2, col = "gray60")
        return()
      }
      
      # Check if dittoSeq is available
      if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        plot.new()
        text(0.5, 0.5, "dittoSeq package required\nInstall with: BiocManager::install('dittoSeq')", 
             cex = 1.2, col = "red")
        return()
      }
      
      library(dittoSeq)
      
      # Create UMAP plot colored by clusters - OPTIMIZED FOR LARGER DISPLAY
      tryCatch({
        p <- dittoDimPlot(
          umap_data$sce,
          var = "seurat_clusters",
          reduction.use = "UMAP",
          size = 0.7,  # Increased point size for better visibility
          do.label = TRUE,
          labels.size = 6,  # DOUBLED label size as requested
          legend.show = TRUE,
          main = ""  # No title to save space
        )
        
        # Additional customization to maximize plot area and improve readability
        p <- p + 
          theme(
            # Maximize plot area by minimizing margins
            plot.margin = margin(2, 2, 2, 2, "pt"),  # Minimal margins
            # Larger legend text for better readability
            legend.text = element_text(size = 14),
            legend.title = element_text(size = 16),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 10, "pt"),
            # Ensure axis text is readable but compact
            axis.text = element_text(size = 12),
            axis.title = element_text(size = 14),
            # Remove unnecessary background elements
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            # Optimize panel grid
            panel.grid.major = element_line(color = "grey95", size = 0.2),
            panel.grid.minor = element_blank()
          ) +
          # Use coord_fixed but allow some flexibility for the container
          coord_fixed(ratio = 1, expand = FALSE)
        
        return(p)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating UMAP:\n", e$message), 
             cex = 1.2, col = "red")
      })
    })
    
    # Update cluster choices when UMAP data is loaded
    observe({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        clusters <- sort(unique(SummarizedExperiment::colData(umap_data$sce)$seurat_clusters))
        cluster_choices <- setNames(as.character(clusters), paste("Cluster", clusters))
        
        updateSelectInput(session, "selected_cluster",
                         choices = cluster_choices,
                         selected = cluster_choices[1])
      }
    })
    
    # Render markers table
    output$markers_table <- DT::renderDataTable({
      req(input$selected_cluster)
      req(umap_data$markers)
      
      # Filter markers for selected cluster
      cluster_markers <- umap_data$markers %>%
        filter(cluster == input$selected_cluster) %>%
        arrange(desc(avg_log2FC)) %>%
        head(input$max_markers) %>%
        select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
        mutate(
          avg_log2FC = round(avg_log2FC, 3),
          p_val_adj = formatC(p_val_adj, format = "e", digits = 2),
          pct.1 = round(pct.1, 3),
          pct.2 = round(pct.2, 3)
        )
      
      # Create DataTable with settings optimized for compact right column
      DT::datatable(
        cluster_markers,
        options = list(
          pageLength = 12,  # Show fewer rows to fit in compact space
          scrollY = "240px",  # Fit in smaller container height
          scrollCollapse = TRUE,
          dom = 't',  # Only show table (no search/pagination)
          autoWidth = FALSE,  # Control column widths
          columnDefs = list(
            list(width = '60px', targets = 0),  # Gene column
            list(width = '50px', targets = 1),  # Log2FC
            list(width = '60px', targets = 2),  # P-val
            list(width = '40px', targets = 3),  # % in cluster
            list(width = '40px', targets = 4),  # % in other
            list(className = 'dt-center', targets = 1:4)  # Center align numeric columns
          )
        ),
        rownames = FALSE,
        colnames = c('Gene', 'Log2FC', 'P-val', '% in', '% out')  # Short names for compact space
      ) %>%
        DT::formatStyle(
          'avg_log2FC',
          background = DT::styleColorBar(cluster_markers$avg_log2FC, 'lightblue'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    })
    
    # Cell count box
    output$total_cells_box <- renderUI({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        n_cells <- ncol(umap_data$sce)
        valueBox(
          value = format(n_cells, big.mark = ","),
          subtitle = "Total Cells",
          icon = icon("microscope"),
          color = "aqua"
        )
      } else {
        valueBox(
          value = "N/A",
          subtitle = "Total Cells",
          icon = icon("microscope"),
          color = "gray"
        )
      }
    })
    
    # Cluster count box
    output$total_clusters_box <- renderUI({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        n_clusters <- length(unique(SummarizedExperiment::colData(umap_data$sce)$seurat_clusters))
        valueBox(
          value = n_clusters,
          subtitle = "Cell Clusters",
          icon = icon("object-group"),
          color = "yellow"
        )
      } else {
        n_clusters <- length(unique(data$consolidated_data$cluster))
        valueBox(
          value = n_clusters,
          subtitle = "Cell Clusters",
          icon = icon("object-group"),
          color = "yellow"
        )
      }
    })
    
    # Results count box
    output$total_results_box <- renderUI({
      valueBox(
        value = format(nrow(data$consolidated_data), big.mark = ","),
        subtitle = "Enrichment Results",
        icon = icon("chart-bar"),
        color = "blue"
      )
    })
    
    # Genes count box
    output$total_genes_box <- renderUI({
      n_genes <- length(unique(data$consolidated_data$gene))
      valueBox(
        value = n_genes,
        subtitle = "Genes Analyzed",
        icon = icon("dna"),
        color = "green"
      )
    })
    
    # Experiments count box
    output$total_experiments_box <- renderUI({
      n_exp <- length(unique(data$consolidated_data$method))
      valueBox(
        value = n_exp,
        subtitle = "Analysis Methods",
        icon = icon("flask"),
        color = "purple"
      )
    })
    
    # Enrichment types box
    output$enrichment_types_box <- renderUI({
      n_types <- length(unique(data$consolidated_data$enrichment_type))
      valueBox(
        value = n_types,
        subtitle = "Enrichment Types",
        icon = icon("database"),
        color = "orange"
      )
    })
    
    # Dataset info panel
    output$dataset_info <- renderUI({
      if (umap_data$loaded) {
        dataset_label <- switch(umap_data$dataset_name,
          "iSCORE_PD" = "iSCORE-PD (Mutations Only)",
          "iSCORE_PD_CRISPRi" = "iSCORE-PD + CRISPRi",
          "Full_Dataset" = "Full Dataset (CRISPRi + CRISPRa)",
          umap_data$dataset_name
        )
        
        tagList(
          p(strong("Active Dataset:"), dataset_label),
          p(strong("Analysis Types:"), paste(unique(data$consolidated_data$method), collapse = ", ")),
          p(strong("Data Loaded:"), format(Sys.time(), "%Y-%m-%d %H:%M"))
        )
      } else {
        p("Loading dataset information...", style = "color: gray;")
      }
    })
    
    # Analysis type plot
    output$analysis_type_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        group_by(method) %>%
        summarise(count = n(), .groups = 'drop')
      
      plot_ly(summary_data, 
              x = ~method, 
              y = ~count, 
              type = 'bar',
              width = 0.8,
              marker = list(color = c('#374E55', '#DF8F44'))) %>%
        layout(title = NULL,
               xaxis = list(title = "", automargin = TRUE),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE,
               autosize = TRUE,
               margin = list(l = 50, r = 20, t = 20, b = 40),
               bargap = 0.3)
    })
    
    # Enrichment type plot
    output$enrichment_type_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        group_by(enrichment_type) %>%
        summarise(count = n(), .groups = 'drop') %>%
        arrange(desc(count))
      
      # Define colors for each enrichment type
      colors <- c(
        "GO_BP" = "#8dd3c7",
        "GO_CC" = "#ffffb3", 
        "GO_MF" = "#bebada",
        "KEGG" = "#fb8072",
        "Reactome" = "#80b1d3",
        "WikiPathways" = "#fdb462",
        "STRING" = "#b3de69",
        "GSEA" = "#fccde5"
      )
      
      plot_ly(summary_data,
              x = ~enrichment_type,
              y = ~count,
              type = 'bar',
              width = 0.8,
              marker = list(color = unname(colors[summary_data$enrichment_type]))) %>%
        layout(title = NULL,
               xaxis = list(title = "", tickangle = -45, automargin = TRUE),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE,
               autosize = TRUE,
               margin = list(l = 50, r = 20, t = 20, b = 80),
               bargap = 0.2)
    })
    
    # Direction plot
    output$direction_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        filter(direction != "RANKED") %>%
        group_by(direction) %>%
        summarise(count = n(), .groups = 'drop')
      
      colors <- c("UP" = "#DC3220", "DOWN" = "#005AB5", "ALL" = "#7C4DFF")
      
      plot_ly(summary_data,
              labels = ~direction,
              values = ~count,
              type = 'pie',
              hole = 0.3,
              marker = list(colors = colors[summary_data$direction])) %>%
        layout(title = NULL,
               showlegend = TRUE,
               autosize = TRUE,
               margin = list(l = 20, r = 20, t = 20, b = 20))
    })
    
    # Gene table
    output$gene_table <- DT::renderDataTable({
      gene_summary <- data$consolidated_data %>%
        group_by(gene, method) %>%
        summarise(
          total_terms = n(),
          enrichment_types = n_distinct(enrichment_type),
          clusters = n_distinct(cluster),
          .groups = 'drop'
        ) %>%
        arrange(desc(total_terms))
      
      DT::datatable(gene_summary,
                    options = list(pageLength = 15),
                    rownames = FALSE) %>%
        formatStyle('total_terms',
                   background = styleColorBar(gene_summary$total_terms, 'lightblue'),
                   backgroundSize = '100% 90%',
                   backgroundRepeat = 'no-repeat',
                   backgroundPosition = 'center')
    })
    
    # Cluster table
    output$cluster_table <- DT::renderDataTable({
      cluster_summary <- data$consolidated_data %>%
        group_by(cluster) %>%
        summarise(
          total_terms = n(),
          genes = n_distinct(gene),
          methods = paste(unique(method), collapse = ", "),
          .groups = 'drop'
        ) %>%
        arrange(cluster)
      
      DT::datatable(cluster_summary,
                    options = list(pageLength = 15),
                    rownames = FALSE)
    })
    
    # Complete matrix table
    output$matrix_table <- DT::renderDataTable({
      matrix_data <- data$consolidated_data %>%
        group_by(gene, cluster, method, enrichment_type, direction) %>%
        summarise(n_terms = n(), .groups = 'drop') %>%
        arrange(gene, cluster, method)
      
      DT::datatable(matrix_data,
                    filter = 'top',
                    options = list(
                      pageLength = 20,
                      scrollX = TRUE
                    ),
                    rownames = FALSE)
    })
    
    # Update gene choices for top terms
    observe({
      gene_choices <- c("All Genes" = "all")
      unique_genes <- sort(unique(data$consolidated_data$gene))
      gene_choices <- c(gene_choices, setNames(unique_genes, unique_genes))
      
      updateSelectInput(session, "top_terms_gene",
                       choices = gene_choices)
    })
    
    # Top terms table
    output$top_terms_table <- DT::renderDataTable({
      top_data <- data$consolidated_data
      
      # Apply filters
      if (input$top_terms_method != "all") {
        top_data <- top_data %>%
          filter(method == input$top_terms_method)
      }
      
      if (!is.null(input$top_terms_gene) && input$top_terms_gene != "all") {
        top_data <- top_data %>%
          filter(gene == input$top_terms_gene)
      }
      
      # Get top terms by p.adjust
      top_terms <- top_data %>%
        arrange(p.adjust) %>%
        head(input$top_terms_n) %>%
        select(gene, cluster, method, enrichment_type, direction, 
               Description, p.adjust, Count) %>%
        mutate(p.adjust = format(p.adjust, scientific = TRUE, digits = 3))
      
      DT::datatable(top_terms,
                    options = list(
                      pageLength = 20,
                      scrollX = TRUE
                    ),
                    rownames = FALSE)
    })
    
  })
}

# Helper function for value boxes (if not already defined)
valueBox <- function(value, subtitle, icon = NULL, color = "blue", width = 12) {
  div(class = paste0("small-box bg-", color), style = "margin-bottom: 15px;",
    div(class = "inner",
      h3(value, style = "margin-bottom: 5px;"),
      p(subtitle, style = "margin: 0;")
    ),
    if (!is.null(icon)) {
      div(class = "icon",
        icon
      )
    }
  )
}