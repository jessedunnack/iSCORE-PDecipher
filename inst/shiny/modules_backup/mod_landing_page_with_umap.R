# Landing Page Module with UMAP Visualization
# Comprehensive overview of all available enrichment results with dataset UMAP plots

# Source the UMAP viewer module
source("modules/mod_umap_viewer.R")

#' Landing Page UI with UMAP
#' 
#' @param id Module namespace
landingPageWithUmapUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Title and UMAP section
    fluidRow(
      column(12,
        h2("PD Enrichment Explorer - Interactive Data Overview", 
           class = "text-center", 
           style = "margin-bottom: 30px;")
      )
    ),
    
    # UMAP Visualization Section
    fluidRow(
      column(12,
        div(class = "box box-primary",
          div(class = "box-header with-border",
            h3(class = "box-title", 
               icon("chart-scatter", lib = "font-awesome"),
               "Dataset UMAP Visualization"),
            div(class = "box-tools pull-right",
              actionButton(ns("toggle_umap"), 
                          "Toggle UMAP", 
                          class = "btn-sm btn-default",
                          icon = icon("chevron-up"))
            )
          ),
          div(class = "box-body", 
              id = ns("umap_container"),
              mod_umap_viewer_ui(ns("umap_viewer"))
          )
        )
      )
    ),
    
    # Summary statistics cards
    fluidRow(
      column(12,
        h3("Enrichment Analysis Summary", 
           style = "margin-top: 30px; margin-bottom: 20px;"),
        
        fluidRow(
          uiOutput(ns("total_results_box")),
          uiOutput(ns("total_genes_box")),
          uiOutput(ns("total_clusters_box")),
          uiOutput(ns("total_experiments_box"))
        )
      )
    ),
    
    # Detailed breakdown
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
          h3("Detailed Result Counts", style = "margin-top: 0; margin-bottom: 20px;"),
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
    ),
    
    # Instructions
    fluidRow(
      column(12,
        div(class = "box box-default collapsed-box",
          div(class = "box-header with-border",
            h3(class = "box-title", "Quick Start Guide"),
            div(class = "box-tools pull-right",
              tags$button(
                class = "btn btn-box-tool",
                `data-widget` = "collapse",
                tags$i(class = "fa fa-plus")
              )
            )
          ),
          div(class = "box-body", style = "display: none;",
            p("Welcome to the PD Enrichment Explorer! This application provides interactive visualization and analysis of functional enrichment results from Parkinson's Disease studies."),
            h4("Getting Started:"),
            tags$ol(
              tags$li("Explore the UMAP visualization above to see the cellular landscape of our datasets"),
              tags$li("Review the summary statistics to understand the scope of available data"),
              tags$li("Use the sidebar on the left to navigate to specific analyses"),
              tags$li("Click on any chart or table for interactive exploration")
            ),
            h4("Available Analyses:"),
            tags$ul(
              tags$li(strong("Precomputed Results:"), " Browse and visualize enrichment results by gene, cluster, and analysis type"),
              tags$li(strong("Comparison Analysis:"), " Compare enrichment patterns across different conditions"),
              tags$li(strong("Heatmap Explorer:"), " Create customizable heatmaps of enrichment results"),
              tags$li(strong("Pathway Viewer:"), " Visualize enriched pathways in detail")
            )
          )
        )
      )
    ),
    
    # Add JavaScript for toggle functionality
    tags$script(HTML(sprintf("
      $(document).ready(function() {
        $('#%s').click(function() {
          $('#%s').slideToggle();
          var icon = $(this).find('i');
          if (icon.hasClass('fa-chevron-up')) {
            icon.removeClass('fa-chevron-up').addClass('fa-chevron-down');
          } else {
            icon.removeClass('fa-chevron-down').addClass('fa-chevron-up');
          }
        });
      });
    ", ns("toggle_umap"), ns("umap_container"))))
  )
}

#' Landing Page Server with UMAP
#' 
#' @param id Module namespace
#' @param data Reactive data object from app
landingPageWithUmapServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    
    # Call the UMAP viewer module
    umap_data <- mod_umap_viewer_server("umap_viewer", data)
    
    # Summary statistics
    output$total_results_box <- renderUI({
      valueBox(
        value = format(nrow(data$consolidated_data), big.mark = ","),
        subtitle = "Total Enrichment Results",
        icon = icon("chart-bar"),
        color = "blue"
      )
    })
    
    output$total_genes_box <- renderUI({
      n_genes <- length(unique(data$consolidated_data$gene))
      valueBox(
        value = n_genes,
        subtitle = "Genes/Mutations Analyzed",
        icon = icon("dna"),
        color = "green"
      )
    })
    
    output$total_clusters_box <- renderUI({
      n_clusters <- length(unique(data$consolidated_data$cluster))
      valueBox(
        value = n_clusters,
        subtitle = "Cell Clusters",
        icon = icon("object-group"),
        color = "yellow"
      )
    })
    
    output$total_experiments_box <- renderUI({
      n_exp <- length(unique(data$consolidated_data$method))
      valueBox(
        value = n_exp,
        subtitle = "Analysis Methods",
        icon = icon("flask"),
        color = "purple"
      )
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
              marker = list(color = c('#374E55', '#DF8F44'))) %>%
        layout(title = NULL,
               xaxis = list(title = ""),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE)
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
        "STRING" = "#b3de69"
      )
      
      plot_ly(summary_data,
              x = ~enrichment_type,
              y = ~count,
              type = 'bar',
              marker = list(color = colors[summary_data$enrichment_type])) %>%
        layout(title = NULL,
               xaxis = list(title = "", tickangle = -45),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE)
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
              marker = list(colors = colors[summary_data$direction])) %>%
        layout(title = NULL,
               showlegend = TRUE)
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
valueBox <- function(value, subtitle, icon = NULL, color = "blue", width = 3) {
  column(width,
    div(class = paste0("small-box bg-", color),
      div(class = "inner",
        h3(value),
        p(subtitle)
      ),
      if (!is.null(icon)) {
        div(class = "icon",
          icon
        )
      }
    )
  )
}