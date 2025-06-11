# Landing Page Module
# Comprehensive overview of all available enrichment results

#' Landing Page UI
#' 
#' @param id Module namespace
landingPageUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Overview cards
    fluidRow(
      column(12,
        h2("PD Enrichment Explorer - Data Overview", 
           class = "text-center", 
           style = "margin-bottom: 30px;"),
        
        # Summary statistics cards
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
            
            # Data Validation
            tabPanel(
              "Data Validation",
              div(style = "margin-top: 15px;",
                h4("Data Quality Checks"),
                uiOutput(ns("validation_results")),
                br(),
                h4("Missing Combinations"),
                p("These parameter combinations have no significant results:"),
                DT::dataTableOutput(ns("missing_table"))
              )
            )
          )
        )
      )
    ),
    
    # Navigation helper
    fluidRow(
      column(12,
        div(class = "box box-warning collapsed-box",
          div(class = "box-header with-border",
            h3(class = "box-title", "Quick Navigation"),
            div(class = "box-tools pull-right",
              tags$button(
                class = "btn btn-box-tool",
                `data-widget` = "collapse",
                icon("plus")
              )
            )
          ),
          div(class = "box-body", style = "display: none;",
            p("Click on any row in the tables above to navigate directly to those results in the analysis tabs."),
            
            fluidRow(
              column(4,
                h5("Top Genes by Result Count:"),
                uiOutput(ns("top_genes"))
              ),
              column(4,
                h5("Most Active Clusters:"),
                uiOutput(ns("top_clusters"))
              ),
              column(4,
                h5("Enrichment Databases:"),
                uiOutput(ns("enrichment_summary"))
              )
            )
          )
        )
      )
    )
  )
}

#' Landing Page Server
#' 
#' @param id Module namespace
#' @param app_data Reactive values containing consolidated data
landingPageServer <- function(id, app_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive: Calculate comprehensive statistics
    result_stats <- reactive({
      req(app_data$consolidated_data)
      data <- app_data$consolidated_data
      
      # Overall statistics
      list(
        total_results = nrow(data),
        unique_genes = length(unique(data$gene)),
        unique_clusters = length(unique(data$cluster)),
        unique_experiments = length(unique(data$experiment)),
        unique_methods = length(unique(data$method)),
        unique_enrichments = length(unique(data$enrichment_type)),
        unique_directions = length(unique(data$direction))
      )
    })
    
    # Custom value boxes (avoiding shinydashboard dependency)
    output$total_results_box <- renderUI({
      div(class = "col-sm-3",
        div(class = "small-box bg-blue",
          div(class = "inner",
            h3(format(result_stats()$total_results, big.mark = ",")),
            p("Total Significant Terms")
          ),
          div(class = "icon",
            icon("database")
          )
        )
      )
    })
    
    output$total_genes_box <- renderUI({
      div(class = "col-sm-3",
        div(class = "small-box bg-green",
          div(class = "inner",
            h3(result_stats()$unique_genes),
            p("Genes/Mutations")
          ),
          div(class = "icon",
            icon("dna")
          )
        )
      )
    })
    
    output$total_clusters_box <- renderUI({
      div(class = "col-sm-3",
        div(class = "small-box bg-yellow",
          div(class = "inner",
            h3(result_stats()$unique_clusters),
            p("Cell Clusters")
          ),
          div(class = "icon",
            icon("circle-nodes")
          )
        )
      )
    })
    
    output$total_experiments_box <- renderUI({
      div(class = "col-sm-3",
        div(class = "small-box bg-red",
          div(class = "inner",
            h3(result_stats()$unique_experiments),
            p("Experiments")
          ),
          div(class = "icon",
            icon("flask")
          )
        )
      )
    })
    
    # Analysis type breakdown
    output$analysis_type_plot <- renderPlotly({
      data <- app_data$consolidated_data
      
      summary_data <- data %>%
        group_by(method) %>%
        summarise(
          count = n(),
          genes = n_distinct(gene),
          clusters = n_distinct(cluster)
        ) %>%
        arrange(desc(count))
      
      p <- plot_ly(summary_data, 
                   x = ~method, 
                   y = ~count,
                   type = 'bar',
                   name = 'Total Terms',
                   text = ~paste("Terms:", format(count, big.mark = ","),
                               "<br>Genes:", genes,
                               "<br>Clusters:", clusters),
                   hovertemplate = "%{text}<extra></extra>",
                   marker = list(color = '#3c8dbc'))
      
      p %>% layout(
        xaxis = list(title = "Analysis Method"),
        yaxis = list(title = "Number of Significant Terms"),
        showlegend = FALSE
      )
    })
    
    # Enrichment type breakdown
    output$enrichment_type_plot <- renderPlotly({
      data <- app_data$consolidated_data
      
      # Simplify enrichment types
      summary_data <- data %>%
        mutate(enrichment_simple = gsub("_ALL|_BP|_CC|_MF", "", enrichment_type)) %>%
        group_by(enrichment_simple) %>%
        summarise(
          count = n(),
          subtypes = n_distinct(enrichment_type)
        ) %>%
        arrange(desc(count))
      
      p <- plot_ly(summary_data, 
                   x = ~reorder(enrichment_simple, count), 
                   y = ~count,
                   type = 'bar',
                   orientation = 'v',
                   text = ~paste("Terms:", format(count, big.mark = ","),
                               "<br>Subtypes:", subtypes),
                   hovertemplate = "%{text}<extra></extra>",
                   marker = list(color = '#00a65a'))
      
      p %>% layout(
        xaxis = list(title = "Enrichment Database"),
        yaxis = list(title = "Number of Significant Terms"),
        showlegend = FALSE
      )
    })
    
    # Direction breakdown
    output$direction_plot <- renderPlotly({
      data <- app_data$consolidated_data
      
      summary_data <- data %>%
        group_by(direction) %>%
        summarise(count = n()) %>%
        mutate(
          percentage = round(100 * count / sum(count), 1),
          label = paste0(direction, "\n", format(count, big.mark = ","), " terms\n", percentage, "%")
        )
      
      colors <- c("UP" = "#e74c3c", "DOWN" = "#3498db", "ALL" = "#2ecc71", "RANKED" = "#9b59b6")
      
      p <- plot_ly(summary_data, 
                   labels = ~direction, 
                   values = ~count,
                   type = 'pie',
                   text = ~label,
                   textposition = 'inside',
                   textinfo = 'text',
                   hovertemplate = "%{label}<extra></extra>",
                   marker = list(colors = colors[summary_data$direction]))
      
      p %>% layout(showlegend = TRUE)
    })
    
    # Detailed gene table
    output$gene_table <- DT::renderDataTable({
      data <- app_data$consolidated_data
      
      gene_summary <- data %>%
        group_by(gene, method) %>%
        summarise(
          total_terms = n(),
          clusters = n_distinct(cluster),
          experiments = n_distinct(experiment),
          enrichment_types = n_distinct(enrichment_type),
          up_terms = sum(direction == "UP"),
          down_terms = sum(direction == "DOWN"),
          all_terms = sum(direction == "ALL"),
          .groups = "drop"
        ) %>%
        arrange(desc(total_terms))
      
      DT::datatable(
        gene_summary,
        options = list(
          pageLength = 15,
          scrollX = TRUE,
          order = list(list(2, 'desc')),
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        rownames = FALSE,
        selection = 'single',
        class = 'display compact'
      ) %>%
        formatStyle(
          'total_terms',
          background = styleColorBar(gene_summary$total_terms, 'lightblue'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
    }, server = FALSE)
    
    # Detailed cluster table
    output$cluster_table <- DT::renderDataTable({
      data <- app_data$consolidated_data
      
      cluster_summary <- data %>%
        group_by(cluster) %>%
        summarise(
          total_terms = n(),
          genes = n_distinct(gene),
          methods = paste(unique(method), collapse = ", "),
          enrichment_types = n_distinct(enrichment_type),
          up_terms = sum(direction == "UP"),
          down_terms = sum(direction == "DOWN"),
          all_terms = sum(direction == "ALL"),
          .groups = "drop"
        ) %>%
        arrange(desc(total_terms))
      
      DT::datatable(
        cluster_summary,
        options = list(
          pageLength = 12,
          scrollX = TRUE,
          order = list(list(1, 'desc')),
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel')
        ),
        rownames = FALSE,
        selection = 'single',
        class = 'display compact'
      )
    }, server = FALSE)
    
    # Complete matrix table
    output$matrix_table <- DT::renderDataTable({
      data <- app_data$consolidated_data
      
      matrix_summary <- data %>%
        group_by(method, gene, cluster, experiment, enrichment_type, direction) %>%
        summarise(
          term_count = n(),
          .groups = "drop"
        ) %>%
        arrange(method, gene, cluster, experiment, enrichment_type, direction)
      
      DT::datatable(
        matrix_summary,
        options = list(
          pageLength = 25,
          scrollX = TRUE,
          scrollY = "400px",
          dom = 'Bfrtip',
          buttons = c('copy', 'csv', 'excel'),
          columnDefs = list(
            list(width = '80px', targets = c(0, 5, 6))
          )
        ),
        rownames = FALSE,
        filter = 'top',
        class = 'display compact'
      ) %>%
        formatStyle(
          'term_count',
          backgroundColor = styleInterval(
            c(10, 50, 100, 500),
            c('white', '#e8f4f8', '#b8e0ea', '#7fc9d9', '#3c8dbc')
          )
        )
    }, server = FALSE)
    
    # Data validation checks
    output$validation_results <- renderUI({
      data <- app_data$consolidated_data
      
      checks <- list()
      
      # Check for p-value range
      if ("p.adjust" %in% names(data)) {
        max_pval <- max(data$p.adjust, na.rm = TRUE)
        if (max_pval > 0.05) {
          checks[[length(checks) + 1]] <- tags$div(
            class = "alert alert-warning",
            icon("exclamation-triangle"),
            sprintf("Warning: Maximum p.adjust value is %.4f (expected ≤ 0.05)", max_pval)
          )
        } else {
          checks[[length(checks) + 1]] <- tags$div(
            class = "alert alert-success",
            icon("check-circle"),
            "✓ All p.adjust values ≤ 0.05"
          )
        }
      }
      
      # Check for missing values
      missing_cols <- c("gene", "cluster", "experiment", "enrichment_type", "direction")
      for (col in missing_cols) {
        if (any(is.na(data[[col]]))) {
          n_missing <- sum(is.na(data[[col]]))
          checks[[length(checks) + 1]] <- tags$div(
            class = "alert alert-warning",
            icon("exclamation-triangle"),
            sprintf("Warning: %d missing values in %s column", n_missing, col)
          )
        }
      }
      
      # Check for expected enrichment types
      expected_enrichments <- c("GO_BP", "GO_CC", "GO_MF", "GO_ALL", "KEGG", "Reactome", 
                               "WikiPathways", "STRING", "GSEA")
      actual_enrichments <- unique(data$enrichment_type)
      missing_enrichments <- setdiff(expected_enrichments, actual_enrichments)
      
      if (length(missing_enrichments) > 0) {
        checks[[length(checks) + 1]] <- tags$div(
          class = "alert alert-info",
          icon("info-circle"),
          paste("Missing enrichment types:", paste(missing_enrichments, collapse = ", "))
        )
      }
      
      # Check STRING dominance
      string_count <- sum(data$enrichment_type == "STRING")
      string_percentage <- round(100 * string_count / nrow(data), 1)
      if (string_percentage > 80) {
        checks[[length(checks) + 1]] <- tags$div(
          class = "alert alert-info",
          icon("info-circle"),
          sprintf("STRING results comprise %.1f%% of all data (%s terms). This is expected due to comprehensive protein interaction analysis.",
                  string_percentage, format(string_count, big.mark = ","))
        )
      }
      
      # Check for balanced directions
      direction_counts <- table(data$direction)
      if (length(direction_counts) < 3) {
        checks[[length(checks) + 1]] <- tags$div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          paste("Missing directions:", paste(setdiff(c("UP", "DOWN", "ALL"), names(direction_counts)), collapse = ", "))
        )
      }
      
      if (length(checks) == 0) {
        checks[[1]] <- tags$div(
          class = "alert alert-success",
          icon("check-circle"),
          "✓ All data validation checks passed!"
        )
      }
      
      do.call(tagList, checks)
    })
    
    # Missing combinations table
    output$missing_table <- DT::renderDataTable({
      data <- app_data$consolidated_data
      
      # Get all possible combinations
      all_combinations <- expand.grid(
        method = unique(data$method),
        gene = unique(data$gene),
        cluster = unique(data$cluster),
        enrichment_type = unique(data$enrichment_type),
        direction = unique(data$direction),
        stringsAsFactors = FALSE
      )
      
      # Get existing combinations
      existing <- data %>%
        group_by(method, gene, cluster, enrichment_type, direction) %>%
        summarise(exists = TRUE, .groups = "drop")
      
      # Find missing combinations
      missing <- all_combinations %>%
        left_join(existing, by = c("method", "gene", "cluster", "enrichment_type", "direction")) %>%
        filter(is.na(exists)) %>%
        select(-exists) %>%
        arrange(method, gene, cluster, enrichment_type, direction)
      
      if (nrow(missing) > 1000) {
        missing <- missing[1:1000, ]
        showNotification(
          "Showing first 1000 missing combinations only",
          type = "info",
          duration = 5
        )
      }
      
      DT::datatable(
        missing,
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          dom = 'frtip'
        ),
        rownames = FALSE,
        class = 'display compact'
      )
    }, server = FALSE)
    
    # Top genes list
    output$top_genes <- renderUI({
      data <- app_data$consolidated_data
      
      top_genes <- data %>%
        group_by(gene) %>%
        summarise(count = n()) %>%
        arrange(desc(count)) %>%
        head(5)
      
      gene_links <- lapply(1:nrow(top_genes), function(i) {
        tags$div(
          actionLink(
            session$ns(paste0("gene_link_", i)),
            paste0(top_genes$gene[i], " (", format(top_genes$count[i], big.mark = ","), " terms)"),
            style = "cursor: pointer;"
          )
        )
      })
      
      do.call(tagList, gene_links)
    })
    
    # Top clusters list
    output$top_clusters <- renderUI({
      data <- app_data$consolidated_data
      
      top_clusters <- data %>%
        group_by(cluster) %>%
        summarise(
          count = n(),
          genes = n_distinct(gene)
        ) %>%
        arrange(desc(count)) %>%
        head(5)
      
      cluster_text <- lapply(1:nrow(top_clusters), function(i) {
        tags$div(
          paste0(top_clusters$cluster[i], 
                 " (", format(top_clusters$count[i], big.mark = ","), " terms, ",
                 top_clusters$genes[i], " genes)")
        )
      })
      
      do.call(tagList, cluster_text)
    })
    
    # Enrichment summary
    output$enrichment_summary <- renderUI({
      data <- app_data$consolidated_data
      
      enrichment_counts <- data %>%
        mutate(enrichment_base = gsub("_ALL|_BP|_CC|_MF", "", enrichment_type)) %>%
        group_by(enrichment_base) %>%
        summarise(count = n()) %>%
        arrange(desc(count))
      
      enrichment_text <- lapply(1:nrow(enrichment_counts), function(i) {
        tags$div(
          paste0(enrichment_counts$enrichment_base[i], 
                 " (", format(enrichment_counts$count[i], big.mark = ","), " terms)")
        )
      })
      
      do.call(tagList, enrichment_text)
    })
    
    # Handle table row clicks for navigation
    observeEvent(input$gene_table_rows_selected, {
      selected_row <- input$gene_table_rows_selected
      if (!is.null(selected_row)) {
        data <- app_data$consolidated_data
        gene_summary <- data %>%
          group_by(gene, method) %>%
          summarise(total_terms = n(), .groups = "drop") %>%
          arrange(desc(total_terms))
        
        selected_gene <- gene_summary$gene[selected_row]
        selected_method <- gene_summary$method[selected_row]
        
        showNotification(
          paste("Navigate to:", selected_gene, "-", selected_method),
          type = "info",
          duration = 3
        )
        
        # Update global selections to navigate
        updateSelectInput(session$parent, "global_analysis", selected = selected_method)
        updateSelectInput(session$parent, "global_gene", selected = selected_gene)
        
        # Switch to pre-computed results tab
        updateNavbarPage(session$parent, "navbar", selected = "precomputed")
      }
    })
    
  })
}