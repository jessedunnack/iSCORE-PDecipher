# Signature Nomination Module for Cross-Method Comparison
# Interactive interface for discovering shared signatures between MAST and CRISPRi

# Source required R functions for signature analysis
# Find package root directory
pkg_root <- system.file(package = "iSCORE.PDecipher")
if (pkg_root == "") {
  # Development mode - try to find R directory
  if (file.exists("R/gene_harmonization.R")) {
    source("R/gene_harmonization.R")
    source("R/signature_analysis.R")
    source("R/manuscript_signature_discovery.R")
  } else if (file.exists("../../R/gene_harmonization.R")) {
    source("../../R/gene_harmonization.R")
    source("../../R/signature_analysis.R")
    source("../../R/manuscript_signature_discovery.R")
  } else {
    # Try from working directory parent
    base_dir <- getwd()
    while (!file.exists(file.path(base_dir, "R", "gene_harmonization.R")) && base_dir != dirname(base_dir)) {
      base_dir <- dirname(base_dir)
    }
    if (file.exists(file.path(base_dir, "R", "gene_harmonization.R"))) {
      source(file.path(base_dir, "R", "gene_harmonization.R"))
      source(file.path(base_dir, "R", "signature_analysis.R"))
      source(file.path(base_dir, "R", "manuscript_signature_discovery.R"))
    }
  }
} else {
  # Installed package mode - functions should be available via namespace
  # If not, this indicates a package installation issue
}

# UI function
mod_signature_nomination_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      # Left panel: Analysis constraints and controls
      column(4,
        wellPanel(
          h4("Signature Nomination Analysis", icon("search")),
          
          # Simple analysis panel (always visible)
          div(id = ns("simple_panel"),
            h5("Quick Analysis", style = "color: #3c8dbc;"),
            
            # Gene selection with priority indicators
            selectInput(ns("gene_selection"),
                       "Gene Pairs to Analyze:",
                       choices = c("All available pairs" = "all"),
                       selected = "all",
                       multiple = TRUE),
            
            # Cluster selection with dopaminergic priority
            checkboxGroupInput(ns("cluster_selection"),
                              "Clusters to Include:",
                              choices = c(),  # Will be populated by server
                              selected = NULL),
            
            # Variant handling
            div(style = "margin-bottom: 15px;",
              checkboxInput(ns("combine_snca"), "Combine SNCA variants (A30P + A53T)", value = TRUE),
              checkboxInput(ns("combine_vps13c"), "Combine VPS13C variants (A444P + W395C)", value = TRUE)
            ),
            
            # Run analysis button
            div(class = "text-center", style = "margin: 20px 0;",
              actionButton(ns("run_analysis"), 
                          "Discover Signatures", 
                          class = "btn-primary btn-lg",
                          icon = icon("rocket")),
              br(), br(),
              actionButton(ns("run_quick_test"), 
                          "Quick Test (2 clusters)", 
                          class = "btn-warning btn-sm",
                          icon = icon("zap"))
            )
          ),
          
          # Advanced settings (collapsible)
          div(style = "margin-top: 15px;",
            actionButton(ns("toggle_advanced"), "Advanced Settings", 
                        class = "btn-link btn-sm",
                        icon = icon("cog"))
          ),
          
          conditionalPanel(
            condition = "input.show_advanced == true",
            ns = ns,
            div(id = ns("advanced_panel"), style = "margin-top: 15px;",
              h5("Advanced Parameters", style = "color: #f39c12;"),
              
              # Statistical thresholds
              numericInput(ns("min_overlap"), "Minimum Gene Overlap:", 
                          value = 2, min = 1, max = 50),
              numericInput(ns("fisher_threshold"), "Fisher's p-value threshold:", 
                          value = 0.05, min = 0.001, max = 0.1, step = 0.001),
              
              # Clustering parameters
              selectInput(ns("clustering_method"), "Clustering Method:",
                         choices = c("Ward.D2" = "ward.D2", "Complete" = "complete", 
                                   "Average" = "average", "Single" = "single"),
                         selected = "ward.D2"),
              
              # Analysis scope
              checkboxGroupInput(ns("analysis_scope"),
                               "Analysis Components:",
                               choices = c("Gene overlaps" = "genes", 
                                         "Pathway overlaps" = "pathways",
                                         "PD-relevant focus" = "pd_focus"),
                               selected = c("genes", "pathways", "pd_focus")),
              
              # Performance settings
              numericInput(ns("top_signatures"), "Top signatures to display:",
                          value = 20, min = 5, max = 100),
              numericInput(ns("min_cluster_breadth"), "Min clusters for pan-cluster signatures:",
                          value = 8, min = 3, max = 14)
            )
          ),
          
          # Progress and status
          div(id = ns("progress_section"), style = "margin-top: 20px;",
            conditionalPanel(
              condition = "input.run_analysis > 0",
              ns = ns,
              div(
                h5("Analysis Progress"),
                shinycssloaders::withSpinner(
                  textOutput(ns("analysis_progress")),
                  type = 6, color = "#3c8dbc"
                )
              )
            )
          )
        )
      ),
      
      # Right panel: Results and visualizations
      column(8,
        # Results tabs
        tabsetPanel(id = ns("results_tabs"),
          
          # Overview tab
          tabPanel("Signature Overview", value = "overview",
            wellPanel(
              h4("Signature Discovery Results"),
              uiOutput(ns("analysis_summary")),
              
              # Quick action buttons
              div(class = "text-center", style = "margin: 15px 0;",
                downloadButton(ns("download_summary"), "Download Summary (CSV)", 
                              class = "btn-success btn-sm"),
                downloadButton(ns("download_full"), "Download Full Results (Excel)", 
                              class = "btn-info btn-sm", style = "margin-left: 10px;")
              )
            )
          ),
          
          # Pan-cluster signatures tab
          tabPanel("Pan-Cluster Signatures", value = "pan_cluster",
            wellPanel(
              h4("Signatures Across Multiple Clusters"),
              p("These signatures appear consistently across many cell clusters, suggesting core biological effects."),
              
              div(style = "margin-bottom: 15px;",
                DT::dataTableOutput(ns("pan_cluster_table"))
              ),
              
              # Heatmap for pan-cluster signatures
              div(style = "margin-top: 20px;",
                h5("Pan-Cluster Signature Heatmap"),
                plotlyOutput(ns("pan_cluster_heatmap"), height = "500px")
              )
            )
          ),
          
          # Cluster-specific signatures tab  
          tabPanel("Cluster-Specific Signatures", value = "cluster_specific",
            wellPanel(
              h4("Cell Type-Specific Signatures"),
              p("These signatures are unique to specific cell clusters, revealing cell type-specific effects."),
              
              # Cluster selector for detailed view
              div(style = "margin-bottom: 15px;",
                selectInput(ns("selected_cluster_detail"), "View Cluster:",
                           choices = c(), selected = NULL)
              ),
              
              div(style = "margin-bottom: 15px;",
                DT::dataTableOutput(ns("cluster_specific_table"))
              ),
              
              # Visualization for selected cluster
              div(style = "margin-top: 20px;",
                h5("Cluster-Specific Signature Visualization"),
                plotlyOutput(ns("cluster_specific_plot"), height = "400px")
              )
            )
          ),
          
          # Gene pair details tab
          tabPanel("Gene Pair Analysis", value = "gene_pairs", 
            wellPanel(
              h4("Detailed Gene Pair Comparisons"),
              
              # Gene pair selector
              div(style = "margin-bottom: 15px;",
                selectInput(ns("selected_gene_pair"), "Select Gene Pair:",
                           choices = c(), selected = NULL)
              ),
              
              # Detailed results for selected pair
              div(
                uiOutput(ns("gene_pair_details"))
              ),
              
              # Correlation plot for selected pair
              div(style = "margin-top: 20px;",
                h5("Effect Size Correlation"),
                plotlyOutput(ns("gene_pair_correlation"), height = "400px")
              )
            )
          ),
          
          # Interactive heatmap tab
          tabPanel("Signature Heatmap", value = "heatmap",
            wellPanel(
              h4("Interactive Signature Strength Heatmap"),
              
              # Heatmap controls
              fluidRow(
                column(3,
                  selectInput(ns("heatmap_metric"), "Display Metric:",
                             choices = c("Signature Strength" = "signature_strength",
                                       "Gene Overlap Count" = "gene_overlap_count", 
                                       "Fisher p-value" = "gene_fisher_p",
                                       "Jaccard Index" = "gene_jaccard"),
                             selected = "signature_strength")
                ),
                column(3,
                  selectInput(ns("heatmap_clustering"), "Clustering:",
                             choices = c("Both" = "both", "Rows" = "row", 
                                       "Columns" = "column", "None" = "none"),
                             selected = "both")
                ),
                column(3,
                  selectInput(ns("color_scale"), "Color Scale:",
                             choices = c("Viridis" = "viridis", "RdBu" = "RdBu", 
                                       "Red" = "Reds", "Blue" = "Blues"),
                             selected = "viridis")
                ),
                column(3,
                  downloadButton(ns("download_heatmap_html"), "Download Interactive",
                                class = "btn-primary btn-sm")
                )
              ),
              
              # Main heatmap
              div(style = "margin-top: 15px;",
                shinycssloaders::withSpinner(
                  plotlyOutput(ns("signature_heatmap"), height = "600px"),
                  type = 6, color = "#3c8dbc"
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
mod_signature_nomination_server <- function(id, global_selection, app_data) {
  moduleServer(id, function(input, output, session) {
    
    cat("[Signature Nomination] Module server starting...\n")
    
    # Reactive values for storing analysis results
    values <- reactiveValues(
      analysis_results = NULL,
      gene_pairs = NULL,
      available_clusters = NULL,
      analysis_running = FALSE
    )
    
    # Initialize available options when data is loaded
    observe({
      req(app_data$data_loaded)
      cat("[Signature Nomination] Initializing with loaded data...\n")
      
      # Get available clusters from the data
      if (!is.null(app_data$consolidated_data)) {
        unique_clusters <- unique(app_data$consolidated_data$cluster)
        unique_clusters <- sort(unique_clusters[!is.na(unique_clusters)])
        values$available_clusters <- unique_clusters
        
        # Update cluster selection UI
        cluster_choices <- setNames(unique_clusters, paste("Cluster", gsub("cluster_", "", unique_clusters)))
        
        # Mark dopaminergic clusters (typically lower numbered clusters in PD data)
        # You may need to adjust this based on your specific cluster annotations
        dopaminergic_clusters <- unique_clusters[1:min(4, length(unique_clusters))]  # First 4 clusters often dopaminergic
        
        updateCheckboxGroupInput(session, "cluster_selection",
                                choices = cluster_choices,
                                selected = dopaminergic_clusters)
        
        # Initialize gene pair options
        gene_pairs <- get_comparable_gene_pairs(combine_snca_variants = TRUE,
                                               combine_vps13c_variants = TRUE,
                                               include_mast_only = FALSE)
        values$gene_pairs <- gene_pairs
        
        # Update gene selection UI
        gene_choices <- c("All available pairs" = "all")
        for (i in seq_len(nrow(gene_pairs))) {
          pair_name <- paste(gene_pairs$mast_gene[i], "vs", gene_pairs$crispri_gene[i])
          gene_choices[[pair_name]] <- paste0(gene_pairs$mast_gene[i], "_vs_", gene_pairs$crispri_gene[i])
        }
        
        updateSelectInput(session, "gene_selection",
                         choices = gene_choices,
                         selected = "all")
      }
    })
    
    # Toggle advanced settings
    observeEvent(input$toggle_advanced, {
      shinyjs::toggle("advanced_panel")
      current_state <- if(is.null(input$show_advanced)) FALSE else input$show_advanced
      updateCheckboxInput(session, "show_advanced", value = !current_state)
    })
    
    # Main analysis function
    perform_signature_analysis <- function() {
      
      # Validate inputs
      req(app_data$consolidated_data)
      req(input$cluster_selection)
      req(values$gene_pairs)
      
      # Set up progress tracking
      progress <- shiny::Progress$new(session)
      progress$set(message = "Initializing signature analysis...", value = 0)
      on.exit(progress$close())
      
      # Initialize timing
      start_time <- Sys.time()
      last_update <- start_time
      
      # Filter data based on user selections
      selected_clusters <- input$cluster_selection
      
      cat("[SIGNATURE] Starting analysis at", format(start_time), "\n")
      cat("[SIGNATURE] Selected clusters:", paste(selected_clusters, collapse = ", "), "\n")
      
      # Filter enrichment data
      filtered_data <- app_data$consolidated_data
      
      # Filter by clusters
      if (length(selected_clusters) > 0) {
        filtered_data <- filtered_data[filtered_data$cluster %in% selected_clusters, ]
      }
      
      # Filter by methods (MAST vs CRISPRi only - no CRISPRa)
      filtered_data <- filtered_data[filtered_data$method %in% c("MAST", "MixScale"), ]
      
      cat("[SIGNATURE] Filtered data:", nrow(filtered_data), "enrichment terms\n")
      
      progress$set(message = paste("Processing", nrow(filtered_data), "enrichment terms..."), value = 0.1)
      
      # Calculate total work for accurate progress
      total_gene_pairs <- nrow(values$gene_pairs)
      total_clusters <- length(selected_clusters)
      total_combinations <- total_gene_pairs * total_clusters
      
      cat("[SIGNATURE] Total work:", total_gene_pairs, "gene pairs ×", total_clusters, "clusters =", total_combinations, "combinations\n")
      
      progress$set(message = paste("Analyzing", total_gene_pairs, "gene pairs across", total_clusters, "clusters..."), value = 0.2)
      
      # Run signature discovery
      tryCatch({
        signature_results <- discover_top_signatures(
          enrichment_data = filtered_data,
          top_n = input$top_signatures %||% 20,
          min_cluster_breadth = input$min_cluster_breadth %||% 8,
          combine_variants = input$combine_snca && input$combine_vps13c,
          progress_callback = function(msg, value = NULL, detail = NULL) {
            current_time <- Sys.time()
            elapsed <- as.numeric(difftime(current_time, start_time, units = "mins"))
            
            # Enhanced progress message with timing
            if (!is.null(value)) {
              remaining_work <- (1 - value) / max(value, 0.01) * elapsed
              full_msg <- paste0(msg, 
                                sprintf(" (%.1f min elapsed", elapsed),
                                if(!is.null(detail)) paste0(", ", detail) else "",
                                if(value > 0.1) sprintf(", ~%.1f min remaining", remaining_work) else "",
                                ")")
              progress$set(message = full_msg, value = 0.2 + (value * 0.7))
              
              cat("[SIGNATURE]", format(current_time, "%H:%M:%S"), "-", msg, 
                  sprintf("(%.1f%% complete)\n", value * 100))
            } else {
              full_msg <- paste0(msg, sprintf(" (%.1f min elapsed)", elapsed))
              progress$set(message = full_msg)
              cat("[SIGNATURE]", format(current_time, "%H:%M:%S"), "-", msg, "\n")
            }
            
            # Force UI update every 30 seconds
            if (difftime(current_time, last_update, units = "secs") >= 30) {
              last_update <<- current_time
              Sys.sleep(0.1)  # Brief pause to allow UI update
            }
          }
        )
        
        progress$set(message = "Analysis complete!", value = 1.0)
        
        return(signature_results)
        
      }, error = function(e) {
        progress$set(message = paste("Error:", e$message), value = 1.0)
        showNotification(paste("Analysis failed:", e$message), type = "error")
        return(NULL)
      })
    }
    
    # Run analysis when button is clicked
    observeEvent(input$run_analysis, {
      values$analysis_running <- TRUE
      
      # Show progress
      output$analysis_progress <- renderText({
        "Analysis in progress..."
      })
      
      # Perform analysis
      results <- perform_signature_analysis()
      
      if (!is.null(results)) {
        values$analysis_results <- results
        values$analysis_running <- FALSE
        
        # Update UI elements with results
        update_results_ui()
        
        showNotification("Signature analysis completed successfully!", type = "message")
      } else {
        values$analysis_running <- FALSE
        showNotification("Analysis failed. Please check your settings and try again.", type = "error")
      }
    })
    
    # Quick test with limited clusters
    observeEvent(input$run_quick_test, {
      # Temporarily override cluster selection for quick test
      original_clusters <- input$cluster_selection
      quick_clusters <- head(values$available_clusters, 2)  # Just first 2 clusters
      
      values$analysis_running <- TRUE
      
      # Show progress
      output$analysis_progress <- renderText({
        "Running quick test with 2 clusters..."
      })
      
      # Temporarily update cluster selection for this analysis
      updateCheckboxGroupInput(session, "cluster_selection", selected = quick_clusters)
      
      # Perform analysis
      results <- perform_signature_analysis()
      
      # Restore original cluster selection
      updateCheckboxGroupInput(session, "cluster_selection", selected = original_clusters)
      
      if (!is.null(results)) {
        values$analysis_results <- results
        values$analysis_running <- FALSE
        
        # Update UI elements with results
        update_results_ui()
        
        showNotification(paste("Quick test completed! Used clusters:", paste(quick_clusters, collapse = ", ")), type = "message")
      } else {
        values$analysis_running <- FALSE
        showNotification("Quick test failed. Please check your settings and try again.", type = "error")
      }
    })
    
    # Function to update results UI elements
    update_results_ui <- function() {
      req(values$analysis_results)
      
      results <- values$analysis_results
      
      # Update cluster-specific selector
      if (length(results$cluster_specific_signatures) > 0) {
        cluster_choices <- names(results$cluster_specific_signatures)
        updateSelectInput(session, "selected_cluster_detail",
                         choices = setNames(cluster_choices, paste("Cluster", gsub("cluster_", "", cluster_choices))),
                         selected = cluster_choices[1])
      }
      
      # Update gene pair selector
      if (nrow(results$top_signatures) > 0) {
        gene_pair_choices <- unique(results$top_signatures$gene_pair)
        updateSelectInput(session, "selected_gene_pair",
                         choices = gene_pair_choices,
                         selected = gene_pair_choices[1])
      }
    }
    
    # Render analysis summary
    output$analysis_summary <- renderUI({
      req(values$analysis_results)
      
      results <- values$analysis_results
      summary_stats <- results$analysis_summary
      
      # Defensive programming - check if summary_stats exists and has required fields
      if (is.null(summary_stats)) {
        return(div(
          h4("Analysis Error"),
          p("Summary statistics not available. Please try running the analysis again.")
        ))
      }
      
      # Provide default values for missing fields
      total_signatures <- summary_stats$total_signatures %||% 0
      pan_cluster_count <- summary_stats$pan_cluster_count %||% 0
      total_gene_pairs <- summary_stats$total_gene_pairs %||% 0
      cluster_specific_count <- summary_stats$cluster_specific_count %||% 0
      strongest_gene_pair <- summary_stats$strongest_gene_pair %||% "None"
      
      tagList(
        div(class = "row",
          div(class = "col-md-3",
            div(class = "info-box bg-blue",
              div(class = "info-box-content",
                span(class = "info-box-text", "Total Signatures"),
                span(class = "info-box-number", total_signatures)
              )
            )
          ),
          div(class = "col-md-3", 
            div(class = "info-box bg-green",
              div(class = "info-box-content",
                span(class = "info-box-text", "Pan-Cluster"),
                span(class = "info-box-number", pan_cluster_count)
              )
            )
          ),
          div(class = "col-md-3",
            div(class = "info-box bg-yellow",
              div(class = "info-box-content",
                span(class = "info-box-text", "Gene Pairs"),
                span(class = "info-box-number", total_gene_pairs)
              )
            )
          ),
          div(class = "col-md-3",
            div(class = "info-box bg-red",
              div(class = "info-box-content",
                span(class = "info-box-text", "Clusters"),
                span(class = "info-box-number", cluster_specific_count)
              )
            )
          )
        ),
        
        div(style = "margin-top: 20px;",
          h5("Top Finding:"),
          p(if(!is.null(strongest_gene_pair) && 
               length(strongest_gene_pair) > 0 && 
               strongest_gene_pair != "None") {
            paste("Strongest signature:", strongest_gene_pair)
          } else {
            "No significant signatures found with current parameters."
          })
        )
      )
    })
    
    # Render pan-cluster table
    output$pan_cluster_table <- DT::renderDataTable({
      req(values$analysis_results)
      req(values$analysis_results$pan_cluster_signatures)
      
      pan_cluster_data <- values$analysis_results$pan_cluster_signatures
      
      if (nrow(pan_cluster_data) > 0) {
        display_data <- pan_cluster_data[, c("gene_pair", "cluster_count", "mean_signature_strength", 
                                           "total_gene_overlaps", "total_pathway_overlaps")]
        colnames(display_data) <- c("Gene Pair", "Clusters", "Mean Strength", "Gene Overlaps", "Pathway Overlaps")
        
        DT::datatable(display_data, 
                     options = list(pageLength = 10, scrollX = TRUE),
                     rownames = FALSE) %>%
          DT::formatRound(c("Mean Strength"), digits = 2)
      } else {
        DT::datatable(data.frame(Message = "No pan-cluster signatures found"),
                     options = list(dom = 't'), rownames = FALSE)
      }
    })
    
    # Helper function for null coalescing
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    # Return reactive values for potential use by other modules
    return(list(
      analysis_results = reactive({ values$analysis_results }),
      analysis_running = reactive({ values$analysis_running })
    ))
  })
}