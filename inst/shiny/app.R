# Full-Featured PerturbSeq Analysis App with Heatmaps
# Includes all original visualization functionality

library(shiny)
library(shinyWidgets)
library(shinycssloaders)
library(shinyjs)
library(DT)
library(ggplot2)
library(plotly)
library(dplyr)
library(tidyr)

# Source heatmap functions
source("R/heatmap_functions.R")

# Load the complete dataset
cat("Loading complete enrichment dataset...\n")

# Check for environment variable first (set by launch_iscore_app)
data_file <- Sys.getenv("ISCORE_DATA_FILE", unset = "")

# If no environment variable, try default location
if (data_file == "" || !file.exists(data_file)) {
  data_file <- "../all_enrichment_padj005_complete_with_direction.rds"
}

# If still not found, try the data_files directory
if (!file.exists(data_file)) {
  data_file <- "../data_files/all_enrichment_padj005_complete_with_direction.rds"
}

if (!file.exists(data_file)) {
  stop("Dataset file not found. Tried:\n", 
       "  - Environment variable ISCORE_DATA_FILE\n",
       "  - ../all_enrichment_padj005_complete_with_direction.rds\n",
       "  - ../data_files/all_enrichment_padj005_complete_with_direction.rds\n",
       "Please provide the correct path to launch_iscore_app()")
}

# Load data
full_data <- readRDS(data_file)
cat("Loaded", nrow(full_data), "enrichment terms\n")

# Ensure required columns exist
if (!"gene" %in% names(full_data) && "mutation_perturbation" %in% names(full_data)) {
  names(full_data)[names(full_data) == "mutation_perturbation"] <- "gene"
}

# Add modality based on method column
if (!"modality" %in% names(full_data)) {
  full_data$modality <- ifelse(full_data$method == "MAST", "MAST", 
                              ifelse(full_data$method == "MixScale_CRISPRa", "CRISPRa", "CRISPRi"))
}

# Filter to significant results
full_data <- full_data[!is.na(full_data$p.adjust) & full_data$p.adjust <= 0.05, ]
cat("Filtered to", nrow(full_data), "significant results\n")

# Get unique values for UI
all_genes <- sort(unique(full_data$gene))
all_clusters <- sort(unique(full_data$cluster))
all_enrichment_types <- sort(unique(full_data$enrichment_type))

# UI
ui <- fluidPage(
  title = "PerturbSeq Analysis - Full Featured",
  useShinyjs(),
  
  tags$head(
    tags$style(HTML("
      .content-header { 
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 20px;
        margin-bottom: 20px;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }
      .filter-section { 
        background: #ffffff;
        padding: 15px;
        border: 1px solid #dee2e6;
        border-radius: 5px;
        margin-bottom: 15px;
        box-shadow: 0 1px 3px rgba(0,0,0,0.05);
      }
      .stats-box {
        background: #f8f9fa;
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 10px;
        border-left: 4px solid #667eea;
      }
      .tab-content { padding-top: 20px; }
      .heatmap-controls {
        background: #e9ecef;
        padding: 15px;
        border-radius: 5px;
        margin-bottom: 20px;
      }
    "))
  ),
  
  div(class = "content-header",
    h1("🧬 PerturbSeq Enrichment Analysis", style = "margin: 0;"),
    h4("Full-Featured Local Version with Heatmap Visualizations", style = "margin: 5px 0; opacity: 0.9;"),
    p(format(nrow(full_data), big.mark = ","), "enrichment terms • ",
      length(unique(full_data$gene)), "genes • ",
      length(unique(full_data$cluster)), "clusters",
      style = "margin: 0; opacity: 0.8;")
  ),
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      
      div(class = "filter-section",
        h4("📊 Dataset Overview", style = "margin-top: 0;"),
        div(class = "stats-box",
          h5("Modality Distribution", style = "margin-top: 0;"),
          tableOutput("modality_stats")
        )
      ),
      
      div(class = "filter-section",
        h4("🔍 Global Filters", style = "margin-top: 0;"),
        
        pickerInput("selected_modalities", "Modalities:",
                   choices = c("MAST", "CRISPRi", "CRISPRa"),
                   selected = c("MAST", "CRISPRi", "CRISPRa"),
                   multiple = TRUE,
                   options = list(`actions-box` = TRUE)),
        
        pickerInput("selected_genes", "Genes:",
                   choices = all_genes,
                   selected = NULL,
                   multiple = TRUE,
                   options = list(
                     `actions-box` = TRUE,
                     `live-search` = TRUE,
                     `selected-text-format` = "count > 3"
                   )),
        
        selectInput("enrichment_type_filter", "Enrichment Type:",
                   choices = c("All" = "ALL", all_enrichment_types),
                   selected = "GO_BP"),
        
        radioButtons("direction_filter", "Direction:",
                    choices = c("All" = "ALL", "Up" = "UP", "Down" = "DOWN"),
                    selected = "ALL"),
        
        sliderInput("p_cutoff", "P-value cutoff:",
                   min = 0, max = 0.05, value = 0.01, step = 0.001),
        
        actionButton("update_filters", "Update Visualizations",
                    class = "btn-primary btn-block")
      )
    ),
    
    mainPanel(
      width = 9,
      
      tabsetPanel(
        id = "main_tabs",
        
        # Heatmap Tab
        tabPanel("🔥 Enrichment Heatmaps",
          value = "heatmaps",
          div(class = "tab-content",
            
            div(class = "heatmap-controls",
              fluidRow(
                column(4,
                  radioButtons("heatmap_type", "Heatmap Type:",
                             choices = c(
                               "Cross-condition comparison" = "cross_condition",
                               "Modality comparison" = "modality",
                               "Gene-centric view" = "gene_centric"
                             ),
                             selected = "cross_condition")
                ),
                column(4,
                  numericInput("max_terms", "Max terms to show:",
                             value = 30, min = 10, max = 100, step = 5)
                ),
                column(4,
                  numericInput("min_frequency", "Min frequency (%):",
                             value = 20, min = 0, max = 100, step = 5)
                )
              )
            ),
            
            conditionalPanel(
              condition = "input.heatmap_type == 'gene_centric'",
              fluidRow(
                column(6,
                  selectInput("heatmap_cluster", "Focus on cluster:",
                            choices = c("All clusters" = "ALL", all_clusters))
                ),
                column(6,
                  p("This view shows enrichment patterns across genes for terms that appear in multiple genes",
                    style = "color: #6c757d; margin-top: 8px;")
                )
              )
            ),
            
            withSpinner(
              plotlyOutput("interactive_heatmap", height = "auto"),
              type = 6,
              color = "#667eea"
            ),
            
            br(),
            
            fluidRow(
              column(6,
                downloadButton("download_heatmap_data", "📊 Download Heatmap Data", 
                             class = "btn-secondary")
              ),
              column(6, align = "right",
                downloadButton("download_heatmap_plot", "🖼️ Download Plot (PNG)", 
                             class = "btn-secondary")
              )
            )
          )
        ),
        
        # Filtered Results Tab
        tabPanel("📋 Filtered Results",
          value = "results",
          div(class = "tab-content",
            div(style = "margin-bottom: 15px;",
              downloadButton("download_results", "📥 Download Filtered Results", 
                           class = "btn-primary"),
              span(style = "margin-left: 20px;",
                   textOutput("result_summary", inline = TRUE))
            ),
            withSpinner(
              DT::dataTableOutput("results_table"),
              type = 6
            )
          )
        ),
        
        # Term Frequency Analysis Tab
        tabPanel("📊 Term Frequency Analysis",
          value = "frequency",
          div(class = "tab-content",
            h4("Most Frequently Enriched Terms Across Conditions"),
            p("Terms that appear significantly enriched in multiple genes/clusters"),
            
            fluidRow(
              column(6,
                withSpinner(
                  plotlyOutput("frequency_plot", height = "500px"),
                  type = 6
                )
              ),
              column(6,
                h5("Top Frequent Terms"),
                withSpinner(
                  DT::dataTableOutput("frequency_table"),
                  type = 6
                )
              )
            )
          )
        ),
        
        # Comparison Tab
        tabPanel("🔄 Modality Comparison",
          value = "comparison",
          div(class = "tab-content",
            h4("Compare Enrichment Between Modalities"),
            p("Identify terms that are consistently enriched across MAST, CRISPRi, and CRISPRa experiments"),
            
            fluidRow(
              column(12,
                div(class = "heatmap-controls",
                  selectInput("comparison_genes", "Select genes to compare:",
                            choices = all_genes,
                            selected = all_genes[1:min(5, length(all_genes))],
                            multiple = TRUE)
                )
              )
            ),
            
            withSpinner(
              plotlyOutput("modality_comparison_plot", height = "600px"),
              type = 6
            ),
            
            br(),
            
            h5("Shared Terms Across Modalities"),
            withSpinner(
              DT::dataTableOutput("shared_terms_table"),
              type = 6
            )
          )
        ),
        
        # Summary Statistics Tab
        tabPanel("📈 Summary Statistics",
          value = "summary",
          div(class = "tab-content",
            h4("Dataset Summary Statistics"),
            
            fluidRow(
              column(6,
                h5("Enrichment by Gene"),
                withSpinner(
                  plotlyOutput("gene_summary_plot", height = "400px"),
                  type = 6
                )
              ),
              column(6,
                h5("Enrichment by Type"),
                withSpinner(
                  plotlyOutput("enrichment_type_plot", height = "400px"),
                  type = 6
                )
              )
            ),
            
            br(),
            
            h5("Detailed Statistics"),
            verbatimTextOutput("detailed_stats")
          )
        )
      )
    )
  )
)

# Server
server <- function(input, output, session) {
  
  # Reactive values
  filtered_data <- reactiveVal(NULL)
  heatmap_data <- reactiveVal(NULL)
  
  # Initial data filtering
  observe({
    data <- full_data
    
    # Filter by modality
    if (!is.null(input$selected_modalities) && length(input$selected_modalities) > 0) {
      data <- data[data$modality %in% input$selected_modalities, ]
    }
    
    # Filter by genes
    if (!is.null(input$selected_genes) && length(input$selected_genes) > 0) {
      data <- data[data$gene %in% input$selected_genes, ]
    }
    
    # Filter by enrichment type
    if (!is.null(input$enrichment_type_filter) && input$enrichment_type_filter != "ALL") {
      data <- data[data$enrichment_type == input$enrichment_type_filter, ]
    }
    
    # Filter by direction
    if (input$direction_filter != "ALL") {
      data <- data[data$direction == input$direction_filter, ]
    }
    
    # Filter by p-value
    data <- data[data$p.adjust <= input$p_cutoff, ]
    
    filtered_data(data)
  })
  
  # Update filters button
  observeEvent(input$update_filters, {
    updateTabsetPanel(session, "main_tabs", selected = "heatmaps")
  })
  
  # Modality statistics
  output$modality_stats <- renderTable({
    stats <- table(full_data$modality)
    data.frame(
      Modality = names(stats),
      Terms = as.numeric(stats),
      Percent = paste0(round(100 * as.numeric(stats) / sum(stats), 1), "%")
    )
  }, striped = TRUE, hover = TRUE, width = "100%")
  
  # Interactive heatmap
  output$interactive_heatmap <- renderPlotly({
    req(filtered_data())
    
    if (input$heatmap_type == "cross_condition") {
      # Cross-condition heatmap
      hmap_data <- prepare_enrichment_heatmap(
        filtered_data(),
        genes = input$selected_genes,
        enrichment_types = if(input$enrichment_type_filter == "ALL") {
          c("GO_BP", "KEGG", "Reactome")
        } else {
          input$enrichment_type_filter
        },
        direction = input$direction_filter,
        max_terms = input$max_terms,
        min_frequency = input$min_frequency / 100,
        p_cutoff = input$p_cutoff
      )
      
      heatmap_data(hmap_data)
      
      if (nrow(hmap_data$matrix) > 0) {
        create_interactive_heatmap(hmap_data, 
                                 title = "Enrichment Terms Across Conditions")
      } else {
        plotly_empty() %>%
          layout(title = "No enriched terms found with current filters")
      }
      
    } else if (input$heatmap_type == "modality") {
      # Modality comparison heatmap
      hmap_data <- create_modality_comparison_heatmap(
        filtered_data(),
        genes = input$selected_genes,
        enrichment_type = if(input$enrichment_type_filter == "ALL") "GO_BP" else input$enrichment_type_filter,
        max_terms = input$max_terms,
        p_cutoff = input$p_cutoff
      )
      
      heatmap_data(hmap_data)
      
      if (nrow(hmap_data$matrix) > 0) {
        plot_ly(
          z = hmap_data$matrix,
          x = colnames(hmap_data$matrix),
          y = rownames(hmap_data$matrix),
          type = "heatmap",
          colorscale = "RdBu",
          reversescale = TRUE
        ) %>%
          layout(
            title = "Terms Enriched Across Modalities",
            xaxis = list(title = "Modality"),
            yaxis = list(title = "Enrichment Term"),
            height = 600 + nrow(hmap_data$matrix) * 15,
            margin = list(l = 300)
          )
      } else {
        plotly_empty() %>%
          layout(title = "No shared terms found across modalities")
      }
      
    } else {
      # Gene-centric heatmap
      hmap_data <- create_gene_enrichment_heatmap(
        filtered_data(),
        enrichment_type = if(input$enrichment_type_filter == "ALL") "GO_BP" else input$enrichment_type_filter,
        cluster = if(input$heatmap_cluster == "ALL") NULL else input$heatmap_cluster,
        direction = input$direction_filter,
        max_terms = input$max_terms,
        p_cutoff = input$p_cutoff
      )
      
      heatmap_data(hmap_data)
      
      if (nrow(hmap_data$matrix) > 0) {
        plot_ly(
          z = hmap_data$matrix,
          x = colnames(hmap_data$matrix),
          y = rownames(hmap_data$matrix),
          type = "heatmap",
          colorscale = "Viridis"
        ) %>%
          layout(
            title = "Terms Enriched Across Genes",
            xaxis = list(title = "Gene", tickangle = -45),
            yaxis = list(title = "Enrichment Term"),
            height = 600 + nrow(hmap_data$matrix) * 15,
            margin = list(l = 300, b = 100)
          )
      } else {
        plotly_empty() %>%
          layout(title = "No terms found enriched in multiple genes")
      }
    }
  })
  
  # Results table
  output$results_table <- DT::renderDataTable({
    req(filtered_data())
    
    display_cols <- c("gene", "cluster", "modality", "enrichment_type", 
                     "Description", "p.adjust", "direction")
    optional_cols <- c("count", "fold_enrichment", "Count")
    
    for (col in optional_cols) {
      if (col %in% names(filtered_data())) {
        display_cols <- c(display_cols, col)
      }
    }
    
    display_data <- filtered_data()[, intersect(display_cols, names(filtered_data())), drop = FALSE]
    
    DT::datatable(
      display_data,
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        dom = 'Bfrtip'
      ),
      filter = 'top',
      rownames = FALSE
    ) %>%
    DT::formatRound(columns = intersect(c("p.adjust", "fold_enrichment"), names(display_data)), digits = 4)
  })
  
  # Result summary
  output$result_summary <- renderText({
    req(filtered_data())
    paste(format(nrow(filtered_data()), big.mark = ","), "enrichment terms")
  })
  
  # Term frequency plot
  output$frequency_plot <- renderPlotly({
    req(filtered_data())
    
    # Calculate term frequency
    term_freq <- filtered_data() %>%
      group_by(Description) %>%
      summarise(
        n_conditions = n_distinct(paste(gene, cluster, direction)),
        mean_p = mean(p.adjust)
      ) %>%
      filter(n_conditions >= 3) %>%
      arrange(desc(n_conditions)) %>%
      head(20)
    
    if (nrow(term_freq) > 0) {
      plot_ly(
        data = term_freq,
        x = ~n_conditions,
        y = ~reorder(Description, n_conditions),
        type = 'bar',
        orientation = 'h',
        marker = list(color = ~-log10(mean_p), colorscale = 'Viridis'),
        text = ~paste("Count:", n_conditions, "<br>Mean p:", round(mean_p, 4)),
        hoverinfo = 'text'
      ) %>%
        layout(
          title = "Most Frequently Enriched Terms",
          xaxis = list(title = "Number of Conditions"),
          yaxis = list(title = ""),
          margin = list(l = 300)
        )
    } else {
      plotly_empty() %>%
        layout(title = "No frequently enriched terms found")
    }
  })
  
  # Frequency table
  output$frequency_table <- DT::renderDataTable({
    req(filtered_data())
    
    term_freq <- filtered_data() %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        Conditions = n_distinct(paste(gene, cluster, direction)),
        Genes = n_distinct(gene),
        Mean_P = mean(p.adjust),
        .groups = 'drop'
      ) %>%
      filter(Conditions >= 2) %>%
      arrange(desc(Conditions), Mean_P)
    
    DT::datatable(
      term_freq,
      options = list(pageLength = 15, dom = 'ft'),
      rownames = FALSE
    ) %>%
    DT::formatRound(columns = "Mean_P", digits = 4)
  })
  
  # Modality comparison plot
  output$modality_comparison_plot <- renderPlotly({
    req(filtered_data(), input$comparison_genes)
    
    comp_data <- filtered_data() %>%
      filter(gene %in% input$comparison_genes)
    
    if (nrow(comp_data) > 0) {
      hmap_data <- create_modality_comparison_heatmap(
        comp_data,
        genes = input$comparison_genes,
        enrichment_type = if(input$enrichment_type_filter == "ALL") "GO_BP" else input$enrichment_type_filter,
        max_terms = 30,
        p_cutoff = input$p_cutoff
      )
      
      if (nrow(hmap_data$matrix) > 0) {
        plot_ly(
          z = hmap_data$matrix,
          x = colnames(hmap_data$matrix),
          y = rownames(hmap_data$matrix),
          type = "heatmap",
          colorscale = "RdBu",
          reversescale = TRUE
        ) %>%
          layout(
            title = paste("Modality Comparison:", paste(input$comparison_genes, collapse = ", ")),
            xaxis = list(title = "Modality"),
            yaxis = list(title = "Enrichment Term"),
            height = max(400, 200 + nrow(hmap_data$matrix) * 20),
            margin = list(l = 300)
          )
      } else {
        plotly_empty() %>%
          layout(title = "No shared terms found")
      }
    }
  })
  
  # Shared terms table
  output$shared_terms_table <- renderDataTable({
    req(filtered_data())
    
    shared_terms <- filtered_data() %>%
      group_by(Description, enrichment_type) %>%
      summarise(
        MAST = sum(modality == "MAST"),
        CRISPRi = sum(modality == "CRISPRi"),
        CRISPRa = sum(modality == "CRISPRa"),
        Total_Modalities = n_distinct(modality),
        .groups = 'drop'
      ) %>%
      filter(Total_Modalities >= 2) %>%
      arrange(desc(Total_Modalities), desc(MAST + CRISPRi + CRISPRa))
    
    DT::datatable(
      shared_terms,
      options = list(pageLength = 15),
      rownames = FALSE
    )
  })
  
  # Gene summary plot
  output$gene_summary_plot <- renderPlotly({
    req(filtered_data())
    
    gene_summary <- filtered_data() %>%
      group_by(gene, modality) %>%
      summarise(n_terms = n(), .groups = 'drop')
    
    plot_ly(
      data = gene_summary,
      x = ~gene,
      y = ~n_terms,
      color = ~modality,
      type = 'bar'
    ) %>%
      layout(
        title = "Enrichment Terms by Gene",
        xaxis = list(title = "Gene", tickangle = -45),
        yaxis = list(title = "Number of Terms"),
        barmode = 'stack',
        margin = list(b = 100)
      )
  })
  
  # Enrichment type plot
  output$enrichment_type_plot <- renderPlotly({
    req(filtered_data())
    
    type_summary <- filtered_data() %>%
      group_by(enrichment_type) %>%
      summarise(n_terms = n(), .groups = 'drop') %>%
      arrange(desc(n_terms))
    
    plot_ly(
      data = type_summary,
      labels = ~enrichment_type,
      values = ~n_terms,
      type = 'pie',
      textposition = 'inside',
      textinfo = 'label+percent'
    ) %>%
      layout(title = "Distribution of Enrichment Types")
  })
  
  # Detailed statistics
  output$detailed_stats <- renderText({
    req(filtered_data())
    
    data <- filtered_data()
    
    paste(
      "=== Current Filter Statistics ===",
      paste("Total enrichment terms:", format(nrow(data), big.mark = ",")),
      paste("Unique genes:", length(unique(data$gene))),
      paste("Unique clusters:", length(unique(data$cluster))),
      paste("Enrichment types:", paste(unique(data$enrichment_type), collapse = ", ")),
      "",
      "=== Modality Breakdown ===",
      paste(capture.output(print(table(data$modality))), collapse = "\n"),
      "",
      "=== Direction Breakdown ===",
      paste(capture.output(print(table(data$direction))), collapse = "\n"),
      "",
      "=== P-value Statistics ===",
      paste("Min p.adjust:", format(min(data$p.adjust), scientific = TRUE)),
      paste("Median p.adjust:", format(median(data$p.adjust), scientific = TRUE)),
      paste("Max p.adjust:", format(max(data$p.adjust), scientific = TRUE)),
      sep = "\n"
    )
  })
  
  # Download handlers
  output$download_results <- downloadHandler(
    filename = function() {
      paste0("enrichment_results_filtered_", Sys.Date(), ".csv")
    },
    content = function(file) {
      write.csv(filtered_data(), file, row.names = FALSE)
    }
  )
  
  output$download_heatmap_data <- downloadHandler(
    filename = function() {
      paste0("heatmap_data_", input$heatmap_type, "_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(heatmap_data())
      write.csv(heatmap_data()$matrix, file)
    }
  )
  
  output$download_heatmap_plot <- downloadHandler(
    filename = function() {
      paste0("heatmap_", input$heatmap_type, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      req(heatmap_data())
      
      # Create static heatmap
      p <- create_static_heatmap(
        heatmap_data(),
        title = paste("Enrichment Heatmap -", input$heatmap_type)
      )
      
      # Save to file
      png(file, width = 1200, height = 800 + nrow(heatmap_data()$matrix) * 15)
      print(p)
      dev.off()
    }
  )
}

# Run the app
cat("🚀 Starting full-featured PerturbSeq app with heatmap visualizations...\n")
shinyApp(ui = ui, server = server)