# Signature Trends Analysis Module
# Data-driven discovery of signature patterns with interactive explanations

#' Signature Trends Analysis UI
#'
#' @param id Module ID
#' @return Shiny UI
mod_signature_trends_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Include CSS for tooltips
    tags$head(
      tags$style(HTML("
        .help-icon {
          color: #337ab7;
          margin-left: 5px;
          cursor: help;
        }
        .help-icon:hover {
          color: #23527c;
        }
        .metric-explanation {
          background-color: #f8f9fa;
          border: 1px solid #dee2e6;
          border-radius: 4px;
          padding: 10px;
          margin: 10px 0;
        }
        .trend-summary-box {
          background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
          color: white;
          padding: 20px;
          border-radius: 8px;
          margin-bottom: 20px;
        }
      "))
    ),
    
    fluidRow(
      # Control Panel
      column(4,
        wellPanel(
          h4("Signature Trends Analysis", icon("chart-line")),
          p("Data-driven discovery of signature patterns without manual curation."),
          
          # Analysis Controls
          div(style = "margin-bottom: 15px;",
            h5("Analysis Parameters"),
            
            numericInput(ns("min_frequency"), 
                        label = tagList("Minimum Frequency", 
                                      actionLink(ns("help_frequency"), "", 
                                               icon = icon("info-circle"), 
                                               class = "help-icon")),
                        value = 2, min = 1, max = 10, step = 1),
            
            numericInput(ns("top_n_results"),
                        label = tagList("Top Results to Display",
                                      actionLink(ns("help_top_n"), "", 
                                               icon = icon("info-circle"), 
                                               class = "help-icon")),
                        value = 50, min = 10, max = 200, step = 10)
          ),
          
          # Run Analysis Button
          div(class = "text-center", style = "margin: 20px 0;",
            actionButton(ns("run_trends_analysis"), 
                        "Run Trends Analysis", 
                        class = "btn-primary btn-lg",
                        icon = icon("chart-line"))
          ),
          
          # Progress indicator
          conditionalPanel(
            condition = "input.run_trends_analysis > 0",
            ns = ns,
            div(
              h5("Analysis Progress"),
              shinycssloaders::withSpinner(
                textOutput(ns("analysis_progress")),
                type = 6, color = "#337ab7"
              )
            )
          )
        ),
        
        # Metric Explanations Panel
        wellPanel(
          h4("Metric Definitions", icon("question-circle")),
          
          div(class = "metric-explanation",
            h5("Signature Strength", style = "color: #2c3e50;"),
            p("Composite score combining gene overlap count, pathway overlap count, and statistical significance (Fisher's p-value). Higher values indicate stronger cross-method agreement."),
            tags$small("Range: 0-10+ | Interpretation: >2.0 = Strong, >1.0 = Moderate, <1.0 = Weak")
          ),
          
          div(class = "metric-explanation", 
            h5("Frequency Score", style = "color: #2c3e50;"),
            p("How often a signature pattern appears across different cluster/condition combinations. Normalized to 0-1 scale."),
            tags$small("Range: 0-1 | Interpretation: >0.8 = Very Common, >0.5 = Common, <0.5 = Rare")
          ),
          
          div(class = "metric-explanation",
            h5("Impact Score", style = "color: #2c3e50;"),
            p("Statistical strength normalized to the strongest signature in the dataset. Based on significance and effect size."),
            tags$small("Range: 0-1 | Interpretation: >0.9 = Highest Impact, >0.7 = High Impact, <0.5 = Moderate Impact")
          ),
          
          div(class = "metric-explanation",
            h5("Pattern Categories", style = "color: #2c3e50;"),
            tags$ul(
              tags$li(tags$strong("Mitochondrial:"), " Terms related to mitochondrial function, respiratory chain, oxidative phosphorylation"),
              tags$li(tags$strong("Autophagy:"), " Terms related to autophagy, lysosomal function, cellular cleanup mechanisms"),
              tags$li(tags$strong("Protein Quality:"), " Terms related to protein folding, proteasome, chaperone functions"),
              tags$li(tags$strong("Neurotransmission:"), " Terms related to dopamine, neurotransmitter metabolism"),
              tags$li(tags$strong("Synaptic:"), " Terms related to synaptic function, vesicle transport")
            )
          )
        )
      ),
      
      # Results Panel
      column(8,
        # Summary Statistics
        conditionalPanel(
          condition = "output.trends_available == true",
          ns = ns,
          div(class = "trend-summary-box",
            h4("Trends Analysis Summary", style = "margin-top: 0;"),
            uiOutput(ns("trends_summary_content"))
          )
        ),
        
        # Results Tabs
        conditionalPanel(
          condition = "output.trends_available == true",
          ns = ns,
          tabsetPanel(id = ns("trends_tabs"),
            
            # Frequency Analysis Tab
            tabPanel("Signature Frequency", value = "frequency",
              wellPanel(
                h4("Most Frequently Occurring Signatures"),
                p("Signatures that appear consistently across multiple conditions, indicating robust cross-method patterns."),
                
                DT::dataTableOutput(ns("frequency_table")),
                
                div(style = "margin-top: 20px;",
                  h5("Frequency Distribution"),
                  plotlyOutput(ns("frequency_distribution_plot"), height = "300px")
                )
              )
            ),
            
            # Impact Analysis Tab
            tabPanel("Impact Rankings", value = "impact",
              wellPanel(
                h4("Highest Impact Signatures"),
                p("Signatures ranked by statistical strength and effect size, showing the most significant cross-method agreements."),
                
                DT::dataTableOutput(ns("impact_table")),
                
                div(style = "margin-top: 20px;",
                  h5("Impact Distribution"),
                  plotlyOutput(ns("impact_distribution_plot"), height = "300px")
                )
              )
            ),
            
            # Term Patterns Tab
            tabPanel("Term Patterns", value = "terms",
              wellPanel(
                h4("Most Frequent Enrichment Terms"),
                p("Enrichment terms that appear most commonly across signatures, revealing dominant biological themes."),
                
                DT::dataTableOutput(ns("terms_table")),
                
                div(style = "margin-top: 20px;",
                  h5("Pattern Categories"),
                  plotlyOutput(ns("pattern_categories_plot"), height = "300px")
                )
              )
            ),
            
            # Trend Visualizations Tab
            tabPanel("Trend Visualizations", value = "visualizations",
              wellPanel(
                h4("Interactive Trend Visualizations"),
                
                fluidRow(
                  column(6,
                    h5("Frequency vs Impact Scatter"),
                    plotlyOutput(ns("frequency_vs_impact_plot"), height = "400px")
                  ),
                  column(6,
                    h5("Signature Strength Heatmap"),
                    plotlyOutput(ns("signature_heatmap"), height = "400px")
                  )
                )
              )
            )
          )
        ),
        
        # No Data Message
        conditionalPanel(
          condition = "output.trends_available == false",
          ns = ns,
          div(class = "alert alert-info",
            h4("No Trends Analysis Available"),
            p("Run the trends analysis to discover data-driven signature patterns. This analysis will:"),
            tags$ul(
              tags$li("Identify the most frequently occurring signatures across conditions"),
              tags$li("Rank signatures by statistical impact and significance"),
              tags$li("Discover the most common enrichment terms without manual curation"),
              tags$li("Provide objective, data-driven insights for manuscript development")
            )
          )
        )
      )
    )
  )
}

#' Signature Trends Analysis Server
#'
#' @param id Module ID
#' @param analysis_results Reactive containing signature analysis results
#' @param enrichment_data Reactive containing enrichment data
#' @return Module server function
mod_signature_trends_server <- function(id, analysis_results, enrichment_data) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for trends analysis
    values <- reactiveValues(
      trends_results = NULL,
      analysis_running = FALSE
    )
    
    # Check if trends analysis is available
    output$trends_available <- reactive({
      !is.null(values$trends_results) && !is.null(values$trends_results$trends_summary)
    })
    outputOptions(output, "trends_available", suspendWhenHidden = FALSE)
    
    # Help tooltips using shinyBS if available, otherwise simple alerts
    observeEvent(input$help_frequency, {
      showModal(modalDialog(
        title = "Minimum Frequency Explanation",
        p("Sets the minimum number of times a signature pattern must appear across different conditions to be included in the frequency analysis."),
        tags$ul(
          tags$li(tags$strong("Low values (1-2):"), " Include rare signatures, comprehensive analysis"),
          tags$li(tags$strong("Medium values (3-5):"), " Focus on moderately common patterns"),
          tags$li(tags$strong("High values (6+):"), " Only very frequent, robust patterns")
        ),
        p(tags$strong("Recommendation:"), " Start with 2 for comprehensive discovery, increase to focus on robust patterns."),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    observeEvent(input$help_top_n, {
      showModal(modalDialog(
        title = "Top Results Explanation", 
        p("Number of top-ranked results to display in each analysis category (frequency, impact, terms)."),
        tags$ul(
          tags$li(tags$strong("Small values (10-25):"), " Focus on most significant findings"),
          tags$li(tags$strong("Medium values (26-75):"), " Balanced overview of patterns"),
          tags$li(tags$strong("Large values (76+):"), " Comprehensive analysis for exploration")
        ),
        p(tags$strong("Recommendation:"), " Use 50 for manuscript development, 100+ for comprehensive exploration."),
        easyClose = TRUE,
        footer = modalButton("Close")
      ))
    })
    
    # Run trends analysis
    observeEvent(input$run_trends_analysis, {
      req(analysis_results(), enrichment_data())
      
      values$analysis_running <- TRUE
      
      # Progress updates
      progress_text <- reactive({
        if (values$analysis_running) {
          "Running comprehensive trends analysis..."
        } else {
          "Analysis complete!"
        }
      })
      
      output$analysis_progress <- renderText({ progress_text() })
      
      # Run analysis with error handling
      tryCatch({
        # Convert signature results to expected format
        # analysis_results() returns list(analysis_results = reactive, analysis_running = reactive)
        # We need to get the actual analysis results from the reactive
        signature_data <- analysis_results()$analysis_results()
        
        # Debug: Check what we're actually getting
        cat("[TRENDS DEBUG] Signature data structure:", class(signature_data), "\n")
        if (is.list(signature_data)) {
          cat("[TRENDS DEBUG] List elements:", names(signature_data), "\n")
        }
        
        converted_results <- convert_signature_results_for_trends(signature_data)
        
        trends_result <- analyze_signature_trends(
          signature_results = converted_results,
          enrichment_data = enrichment_data(),
          min_frequency = input$min_frequency,
          top_n = input$top_n_results
        )
        
        values$trends_results <- trends_result
        values$analysis_running <- FALSE
        
        # Show success notification
        showNotification(
          "Trends analysis completed successfully!",
          type = "message",
          duration = 3
        )
        
      }, error = function(e) {
        values$analysis_running <- FALSE
        showNotification(
          paste("Trends analysis failed:", e$message),
          type = "error",
          duration = 5
        )
      })
    })
    
    # Render trends summary
    output$trends_summary_content <- renderUI({
      req(values$trends_results)
      
      summary_stats <- values$trends_results$trends_summary
      metadata <- values$trends_results$analysis_metadata
      
      tagList(
        fluidRow(
          column(4,
            div(style = "text-align: center;",
              h5("Most Frequent", style = "margin-bottom: 5px;"),
              p(summary_stats$most_frequent_signature, style = "font-size: 14px; font-weight: bold;")
            )
          ),
          column(4,
            div(style = "text-align: center;", 
              h5("Highest Impact", style = "margin-bottom: 5px;"),
              p(summary_stats$highest_impact_signature, style = "font-size: 14px; font-weight: bold;")
            )
          ),
          column(4,
            div(style = "text-align: center;",
              h5("Patterns Analyzed", style = "margin-bottom: 5px;"),
              p(summary_stats$total_patterns_analyzed, style = "font-size: 14px; font-weight: bold;")
            )
          )
        ),
        hr(style = "border-color: rgba(255,255,255,0.3);"),
        p(paste("Analysis completed:", format(metadata$analysis_timestamp, "%Y-%m-%d %H:%M")),
          style = "font-size: 12px; text-align: center; margin-bottom: 0;")
      )
    })
    
    # Render frequency table
    output$frequency_table <- DT::renderDataTable({
      req(values$trends_results)
      
      freq_data <- values$trends_results$frequency_analysis$top_frequent_signatures
      
      if (nrow(freq_data) == 0) {
        return(DT::datatable(
          data.frame(Message = "No frequent signatures found with current parameters"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      # Format data for display
      display_data <- freq_data
      colnames(display_data) <- c("Gene Pair", "Frequency Count", "Mean Strength", 
                                 "Max Strength", "Frequency Score", "Cluster Breadth")
      
      DT::datatable(display_data,
                   options = list(pageLength = 15, scrollX = TRUE,
                                 columnDefs = list(list(className = 'dt-center', targets = 1:5))),
                   rownames = FALSE) %>%
        DT::formatRound(c("Mean Strength", "Max Strength", "Frequency Score"), digits = 3)
    })
    
    # Render impact table  
    output$impact_table <- DT::renderDataTable({
      req(values$trends_results)
      
      impact_data <- values$trends_results$impact_analysis$top_impact_signatures
      
      if (nrow(impact_data) == 0) {
        return(DT::datatable(
          data.frame(Message = "No impact signatures available"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      # Select relevant columns for display
      display_cols <- intersect(c("gene_pair", "signature_strength", "impact_score", "cluster", 
                                 "gene_overlap_count", "pathway_overlap_count"), 
                               colnames(impact_data))
      display_data <- impact_data[, display_cols, drop = FALSE]
      
      # Clean column names
      colnames(display_data) <- gsub("_", " ", stringr::str_to_title(colnames(display_data)))
      
      DT::datatable(display_data,
                   options = list(pageLength = 15, scrollX = TRUE,
                                 columnDefs = list(list(className = 'dt-center', targets = 1:ncol(display_data)))),
                   rownames = FALSE) %>%
        DT::formatRound(which(sapply(display_data, is.numeric)), digits = 3)
    })
    
    # Render terms table
    output$terms_table <- DT::renderDataTable({
      req(values$trends_results)
      
      terms_data <- values$trends_results$term_patterns$top_frequent_terms
      
      if (nrow(terms_data) == 0) {
        return(DT::datatable(
          data.frame(Message = "No term patterns available"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      # Format for display
      display_data <- terms_data
      colnames(display_data) <- c("Enrichment Term", "Frequency", "Pattern Category")
      
      DT::datatable(display_data,
                   options = list(pageLength = 15, scrollX = TRUE),
                   rownames = FALSE) %>%
        DT::formatStyle("Pattern Category",
                       backgroundColor = DT::styleEqual(
                         c("Mitochondrial", "Autophagy", "Protein Quality", "Neurotransmission", "Synaptic"),
                         c("#ff9999", "#99ccff", "#99ff99", "#ffcc99", "#cc99ff")
                       ))
    })
    
    # Create visualizations (simplified for now)
    output$frequency_distribution_plot <- renderPlotly({
      req(values$trends_results)
      
      freq_dist <- values$trends_results$frequency_analysis$frequency_distribution
      
      if (nrow(freq_dist) == 0) {
        return(plotly_empty("No frequency data available"))
      }
      
      p <- ggplot2::ggplot(freq_dist, ggplot2::aes(x = frequency_bin, y = signature_count)) +
        ggplot2::geom_col(fill = "#337ab7", alpha = 0.8) +
        ggplot2::labs(title = "Distribution of Signature Frequencies",
                     x = "Frequency (Number of Occurrences)", 
                     y = "Number of Signatures") +
        ggplot2::theme_minimal()
      
      plotly::ggplotly(p, tooltip = c("x", "y"))
    })
    
    output$impact_distribution_plot <- renderPlotly({
      req(values$trends_results)
      
      impact_dist <- values$trends_results$impact_analysis$impact_distribution
      
      if (nrow(impact_dist) == 0) {
        return(plotly_empty("No impact data available"))
      }
      
      p <- ggplot2::ggplot(impact_dist, ggplot2::aes(x = strength_bin, y = signature_count)) +
        ggplot2::geom_col(fill = "#d73027", alpha = 0.8) +
        ggplot2::labs(title = "Distribution of Signature Strengths",
                     x = "Signature Strength", 
                     y = "Number of Signatures") +
        ggplot2::theme_minimal()
      
      plotly::ggplotly(p, tooltip = c("x", "y"))
    })
    
    output$pattern_categories_plot <- renderPlotly({
      req(values$trends_results)
      
      pattern_cats <- values$trends_results$term_patterns$pattern_categories
      
      if (length(pattern_cats) == 0) {
        return(plotly_empty("No pattern categories available"))
      }
      
      cat_data <- data.frame(
        Category = names(pattern_cats),
        Count = as.numeric(pattern_cats),
        stringsAsFactors = FALSE
      )
      
      p <- ggplot2::ggplot(cat_data, ggplot2::aes(x = reorder(Category, Count), y = Count)) +
        ggplot2::geom_col(fill = "#2166ac", alpha = 0.8) +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Enrichment Term Pattern Categories",
                     x = "Category", 
                     y = "Term Count") +
        ggplot2::theme_minimal()
      
      plotly::ggplotly(p, tooltip = c("x", "y"))
    })
    
    # MISSING VISUALIZATIONS FOR TREND VISUALIZATIONS TAB
    
    # Frequency vs Impact scatter plot
    output$frequency_vs_impact_plot <- renderPlotly({
      req(values$trends_results)
      
      freq_data <- values$trends_results$frequency_analysis$top_frequent_signatures
      impact_data <- values$trends_results$impact_analysis$top_impact_signatures
      
      if (is.null(freq_data) || nrow(freq_data) == 0 || 
          is.null(impact_data) || nrow(impact_data) == 0) {
        return(plotly_empty("No frequency vs impact data available"))
      }
      
      # Merge frequency and impact data
      tryCatch({
        # Assume both have gene_pair column for merging
        if ("gene_pair" %in% colnames(freq_data) && "gene_pair" %in% colnames(impact_data)) {
          merged_data <- merge(freq_data, impact_data, by = "gene_pair", all = FALSE)
          
          if (nrow(merged_data) == 0) {
            return(plotly_empty("No overlapping signatures between frequency and impact analysis"))
          }
          
          # Debug merged data structure
          cat("[TRENDS DEBUG] Merged data columns:", paste(names(merged_data), collapse = ", "), "\n")
          cat("[TRENDS DEBUG] Merged data rows:", nrow(merged_data), "\n")
          
          # Find appropriate columns for plotting
          freq_col <- if ("frequency_score" %in% names(merged_data)) "frequency_score" else if ("frequency_count" %in% names(merged_data)) "frequency_count" else names(merged_data)[2]
          impact_col <- if ("impact_score" %in% names(merged_data)) "impact_score" else if ("signature_strength" %in% names(merged_data)) "signature_strength" else names(merged_data)[3]
          size_col <- if ("mean_strength" %in% names(merged_data)) "mean_strength" else if ("signature_strength" %in% names(merged_data)) "signature_strength" else NULL
          
          # Create scatter plot with safe column references
          if (!is.null(size_col)) {
            p <- ggplot2::ggplot(merged_data, ggplot2::aes_string(x = freq_col, y = impact_col)) +
              ggplot2::geom_point(ggplot2::aes_string(size = size_col), alpha = 0.7, color = "#2166ac") +
              ggplot2::labs(title = "Signature Frequency vs Impact",
                           x = "Frequency Score", 
                           y = "Impact Score",
                           size = "Signature Strength") +
              ggplot2::theme_minimal()
          } else {
            p <- ggplot2::ggplot(merged_data, ggplot2::aes_string(x = freq_col, y = impact_col)) +
              ggplot2::geom_point(alpha = 0.7, color = "#2166ac") +
              ggplot2::labs(title = "Signature Frequency vs Impact",
                           x = "Frequency Score", 
                           y = "Impact Score") +
              ggplot2::theme_minimal()
          }
          
          plotly::ggplotly(p)
        } else {
          return(plotly_empty("Cannot merge frequency and impact data - missing gene_pair column"))
        }
      }, error = function(e) {
        return(plotly_empty(paste("Error creating frequency vs impact plot:", e$message)))
      })
    })
    
    # Signature heatmap for trends
    output$signature_heatmap <- renderPlotly({
      req(values$trends_results)
      
      freq_data <- values$trends_results$frequency_analysis$top_frequent_signatures
      
      if (is.null(freq_data) || nrow(freq_data) == 0) {
        return(plotly_empty("No signature data available for heatmap"))
      }
      
      tryCatch({
        # Debug available columns
        cat("[TRENDS DEBUG] freq_data columns for heatmap:", paste(names(freq_data), collapse = ", "), "\n")
        cat("[TRENDS DEBUG] freq_data rows:", nrow(freq_data), "\n")
        
        # Find available numeric columns for heatmap
        numeric_cols <- sapply(freq_data, is.numeric)
        available_metrics <- names(freq_data)[numeric_cols]
        available_metrics <- available_metrics[available_metrics != "gene_pair"]  # Exclude gene_pair
        
        cat("[TRENDS DEBUG] Available numeric columns:", paste(available_metrics, collapse = ", "), "\n")
        
        if (length(available_metrics) >= 2 && "gene_pair" %in% names(freq_data)) {
          
          # Take top 20 signatures for readability
          plot_data <- head(freq_data, 20)
          
          # Use first two available numeric columns
          metric_cols <- head(available_metrics, 2)
          cat("[TRENDS DEBUG] Using metrics for heatmap:", paste(metric_cols, collapse = ", "), "\n")
          
          # Reshape data for heatmap
          library(tidyr)
          heatmap_data <- plot_data %>%
            dplyr::select(gene_pair, all_of(metric_cols)) %>%
            tidyr::pivot_longer(cols = all_of(metric_cols), 
                               names_to = "metric", values_to = "value")
          
          # Create heatmap
          p <- ggplot2::ggplot(heatmap_data, ggplot2::aes(x = metric, y = gene_pair, fill = value)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_viridis_c() +
            ggplot2::labs(title = "Top Signatures Heatmap",
                         x = "Metric", 
                         y = "Gene Pair",
                         fill = "Value") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8))
          
          plotly::ggplotly(p, tooltip = c("x", "y", "fill"))
        } else {
          return(plotly_empty(paste("Insufficient data for heatmap. Available columns:", paste(names(freq_data), collapse = ", "), ". Need gene_pair + 2 numeric columns")))
        }
      }, error = function(e) {
        return(plotly_empty(paste("Error creating signature heatmap:", e$message)))
      })
    })
    
    # Helper function for empty plots
    plotly_empty <- function(message) {
      plotly::plot_ly() %>%
        plotly::add_text(x = 0.5, y = 0.5, text = message, 
                        textfont = list(size = 16, color = "gray"),
                        showlegend = FALSE) %>%
        plotly::layout(
          xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
          yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
          plot_bgcolor = "rgba(0,0,0,0)",
          paper_bgcolor = "rgba(0,0,0,0)"
        )
    }
    
    # Return reactive values for potential use by other modules
    return(list(
      trends_results = reactive({ values$trends_results }),
      analysis_running = reactive({ values$analysis_running })
    ))
  })
}