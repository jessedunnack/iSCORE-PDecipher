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
    source("R/pd_signature_interpretation.R")
  } else if (file.exists("../../R/gene_harmonization.R")) {
    source("../../R/gene_harmonization.R")
    source("../../R/signature_analysis.R")
    source("../../R/manuscript_signature_discovery.R")
    source("../../R/pd_signature_interpretation.R")
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
      source(file.path(base_dir, "R", "pd_signature_interpretation.R"))
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
          h4("Cross-Method Signature Discovery", icon("search")),
          p("Find shared effects between genetic mutations (MAST) and gene knockdowns (CRISPRi)", 
            style = "font-size: 12px; color: #666;"),
          
          # Statistical explanation (collapsible)
          actionLink(ns("show_stats_info"), "Statistical Methods Explained", icon = icon("info-circle")),
          conditionalPanel(
            condition = "input.show_stats_info % 2 == 1",
            ns = ns,
            div(class = "alert alert-info", style = "margin-top: 10px; font-size: 11px;",
              h5("Cross-Method Comparison Statistics", style = "margin-top: 0; color: #2c3e50;"),
              
              div(style = "margin-bottom: 10px;",
                h6(tags$b("Gene Overlap vs Pathway Overlap"), style = "color: #8b4513;"),
                tags$ul(
                  tags$li(tags$b("Gene Overlap:"), " Often non-significant - mutations vs knockdowns affect different specific genes"),
                  tags$li(tags$b("Pathway Overlap:"), " More meaningful - different genes can converge on same biological functions"),
                  tags$li(tags$b("Biological Insight:"), " Functional convergence matters more than exact gene matches")
                )
              ),
              
              div(style = "margin-bottom: 10px;",
                h6(tags$b("Statistical Tests"), style = "color: #1f77b4;"),
                tags$ul(
                  tags$li(tags$b("Fisher's Exact Test:"), " Tests if overlap exceeds chance (2x2 contingency table)"),
                  tags$li(tags$b("Jaccard Index:"), " Overlap/(Union) - size-independent similarity measure"),
                  tags$li(tags$b("Signature Score:"), " Weighted combination of gene + pathway + statistical evidence")
                )
              ),
              
              div(
                h6(tags$b("Interpretation Guide"), style = "color: #d73027;"),
                tags$ul(
                  tags$li(tags$b("Non-significant gene overlap:"), " Expected and normal for cross-method comparisons"),
                  tags$li(tags$b("Significant pathway overlap:"), " Strong evidence of functional convergence"),
                  tags$li(tags$b("High Jaccard index:"), " Methods targeting similar gene sets (>0.1 is good)")
                )
              )
            )
          ),
          
          # Simple analysis panel (always visible)
          div(id = ns("simple_panel"), style = "margin-top: 15px;",
            h5("Analysis Settings", style = "color: #3c8dbc;"),
            
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
              checkboxInput(ns("combine_snca"), "Combine SNCA variants (A30P + A53T)", value = FALSE),
              checkboxInput(ns("combine_vps13c"), "Combine VPS13C variants (A444P + W395C)", value = FALSE)
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
          
          # Enhanced Analysis Transparency (v0.2.6)
          div(style = "margin-top: 15px;",
            actionLink(ns("show_enhanced_info"), "Enhanced Analysis Information", 
                      icon = icon("microscope"),
                      style = "font-size: 12px; color: #17a2b8;")
          ),
          conditionalPanel(
            condition = "input.show_enhanced_info % 2 == 1",
            ns = ns,
            div(class = "alert alert-info", style = "margin-top: 10px; font-size: 11px;",
              h5("Enhanced Direction-Aware Analysis (v0.2.6)", style = "margin-top: 0; color: #2c3e50;"),
              
              div(style = "margin-bottom: 10px;",
                h6(tags$b("Direction Analysis"), style = "color: #8b4513;"),
                tags$ul(
                  tags$li(tags$b("LRRK2:"), " Opposing effects (gain-of-function mutation vs loss-of-function CRISPRi)"),
                  tags$li(tags$b("SNCA variants:"), " Same-direction effects (aggregation-related pathways)"),
                  tags$li(tags$b("Other mutations:"), " Same-direction effects (loss-of-function patterns)")
                )
              ),
              
              div(style = "margin-bottom: 10px;",
                h6(tags$b("Experiment Weighting"), style = "color: #1f77b4;"),
                tags$ul(
                  tags$li(tags$b("Primary:"), " C12_FPD-24 (highest cell count)"),
                  tags$li(tags$b("Secondary:"), " C12_FPD-23, C18_FPD-23"),
                  tags$li(tags$b("Meta-analysis:"), " Weighted by perturbed cell counts per cluster")
                )
              ),
              
              div(
                h6(tags$b("Statistical Enhancements"), style = "color: #d73027;"),
                tags$ul(
                  tags$li(tags$b("Hierarchical FDR:"), " Benjamini-Yekutieli correction for dependent tests"),
                  tags$li(tags$b("Fisher's Method:"), " Weighted p-value combination across experiments"),
                  tags$li(tags$b("Background Genes:"), " Intersection approach (conservative) & Union approach (liberal)")
                )
              )
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
              
              # Enhanced Analysis Control (v0.2.6)
              div(style = "margin-bottom: 15px;",
                checkboxInput(ns("use_enhanced_analysis"), 
                             "Enhanced Direction-Aware Analysis (v0.2.6)", 
                             value = TRUE),
                div(style = "font-size: 11px; color: #666; margin-top: 5px;",
                  "Includes experiment weighting, direction analysis, and hierarchical FDR correction"
                )
              ),
              
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
        # Results tabs - PD Biology Focus set as default tab
        tabsetPanel(id = ns("results_tabs"), selected = "pd_biology",
          
          # PD Biological Interpretation tab (NOW DEFAULT)
          tabPanel("PD Biology Focus", value = "pd_biology",
            wellPanel(
              h4("Parkinson's Disease Biological Interpretation", icon("brain")),
              p("Focused analysis of signatures for Parkinson's disease-relevant biological processes."),
              
              # PD analysis summary
              div(style = "margin-bottom: 20px;",
                h5("Key Biological Findings", style = "color: #d73027;"),
                uiOutput(ns("pd_biological_summary"))
              ),
              
              # Top PD-relevant signatures table
              div(style = "margin-bottom: 20px;",
                h5("Top PD-Relevant Signatures"),
                DT::dataTableOutput(ns("pd_signatures_table"))
              ),
              
              # PD pathway frequency analysis
              div(style = "margin-bottom: 20px;",
                h5("Most Frequently Disrupted PD Pathways"),
                plotlyOutput(ns("pd_pathway_frequency"), height = "400px")
              ),
              
              # Biological category breakdown
              div(style = "margin-bottom: 20px;",
                h5("Biological Process Categories"),
                plotlyOutput(ns("pd_categories_plot"), height = "350px")
              ),
              
              # Download PD analysis results
              div(class = "text-center", style = "margin: 15px 0;",
                downloadButton(ns("download_pd_analysis"), "Download PD Analysis (CSV)", 
                              class = "btn-warning btn-sm")
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
              
              # Contextual overview explaining the analysis (collapsible)
              div(style = "margin-bottom: 15px;",
                actionButton(ns("toggle_help"), "Understanding Cross-Method Analysis", 
                           icon = icon("info-circle"), class = "btn-info btn-sm",
                           style = "margin-bottom: 10px;")
              ),
              
              conditionalPanel(
                condition = "input.toggle_help % 2 == 1",
                ns = ns,
                div(class = "alert alert-info", style = "margin-bottom: 15px; font-size: 0.9em;",
                p(
                  strong("What this analysis shows:"), " For the selected gene pair (e.g., LRRK2 mutation vs LRRK2 knockdown), ",
                  "we compare differentially expressed (DE) genes identified in MAST (genetic mutations) with those from CRISPRi (gene knockdown) ",
                  "within each cell cluster."
                ),
                p(
                  strong("Statistical Testing:"), " We use Fisher's exact test to determine if the overlap in DE genes between methods ",
                  "is greater than expected by chance. Two approaches are shown:"
                ),
                tags$ul(
                  tags$li(tags$strong("Intersection Approach (Conservative):"), " Tests overlap among genes that both methods could detect (more stringent)"),
                  tags$li(tags$strong("Union Approach (Liberal):"), " Tests overlap among all genes either method could detect (more inclusive)")
                ),
                p(
                  strong("FDR Correction (NEW):"), " P-values are corrected for multiple testing using ", 
                  tags$strong("hierarchical Benjamini-Hochberg FDR correction"), ": ",
                  "first within each gene pair (across clusters and test types), then across all gene pairs. ",
                  "When available, FDR-corrected p-values are displayed by default with clear '(FDR)' labeling."
                ),
                p(
                  strong("Significance Levels:"), " ",
                  tags$span("*", style = "color: #28a745; font-weight: bold;"), " = FDR < 0.05 (suggestive), ",
                  tags$span("**", style = "color: #17a2b8; font-weight: bold;"), " = FDR < 0.01 (significant), ",
                  tags$span("***", style = "color: #dc3545; font-weight: bold;"), " = FDR < 0.001 (highly significant), ",
                  tags$span("ns", style = "color: #6c757d;"), " = not significant"
                ),
                p(
                  tags$em("Note: Gene overlap may be non-significant even when pathway overlap is significant, ",
                  "as different genes can converge on the same biological functions.")
                ),
                div(class = "alert alert-warning", style = "margin-top: 15px; font-size: 0.85em;",
                  icon("exclamation-triangle"),
                  h6("Gene Overlap Definition", style = "margin-top: 0; color: #856404;"),
                  p(
                    tags$strong("Important distinction:"), " The gene overlaps shown here are ", 
                    tags$strong("genes contributing to significantly enriched pathways"), 
                    ", not direct differential expression overlaps. This differs from the 'Shared DE Genes' ",
                    "count on the DE Genes page, which shows direct overlap of genes meeting ",
                    "DE criteria (p < 0.05, |log2FC| > 0.25).", style = "margin-bottom: 8px;"
                  ),
                  p(
                    tags$strong("Why they differ:"), " A gene can be differentially expressed but not contribute ",
                    "to enriched pathways, or contribute to pathways without being the most significantly DE. ",
                    "This pathway-based overlap analysis reveals functional convergence between methods.",
                    style = "margin-bottom: 0;"
                  )
                )
                )
              ),
              
              # Analysis approach toggle (intersection vs union)
              div(style = "margin-bottom: 15px;",
                h5("Analysis Approach", style = "color: #3c8dbc; margin-bottom: 10px;"),
                radioButtons(ns("analysis_approach"), 
                           label = "Statistical Approach:",
                           choices = list(
                             "Intersection (Conservative)" = "intersection",
                             "Union (Liberal)" = "union"
                           ),
                           selected = "intersection",
                           inline = TRUE),
                div(class = "help-text", style = "font-size: 0.85em; color: #666; margin-top: 5px;",
                  HTML("<strong>Intersection:</strong> Tests overlap among genes both methods can detect (more stringent)<br>",
                       "<strong>Union:</strong> Tests overlap among all genes either method can detect (more inclusive)")
                )
              ),
              
              # Dynamic summary of significant clusters
              div(style = "margin-bottom: 15px;",
                uiOutput(ns("gene_pair_summary"))
              ),
              
              # Detailed results for selected pair
              div(
                DT::dataTableOutput(ns("gene_pair_table"))
              ),
              
              # Explanatory footnote for table
              div(class = "alert alert-info", style = "margin-top: 10px; font-size: 0.85em; padding: 8px;",
                icon("info-circle", style = "margin-right: 5px;"),
                tags$strong("Table Explanation:"), tags$br(),
                tags$span(style = "margin-left: 15px;",
                  "• ", tags$strong("Shared DE Genes*:"), " Differentially expressed genes with statistically significant overlap between MAST (genetic mutation) and CRISPRi (gene knockdown), determined by Fisher's exact test.", tags$br(),
                  "• ", tags$strong("p-values:"), " Test if gene overlap exceeds chance expectation. FDR correction controls for multiple testing.", tags$br(),
                  "• ", tags$strong("Background Genes:"), " Total number of genes tested in the selected approach (intersection=conservative, union=liberal)."
                )
              ),
              
              # Correlation plot for selected pair
              div(style = "margin-top: 20px;",
                h5("Effect Size Correlation"),
                div(class = "alert alert-info", style = "margin-bottom: 15px; font-size: 0.85em;",
                  icon("info-circle"),
                  tags$strong(" Purpose: "), 
                  "Compare log2 fold changes between MAST (mutations) and CRISPRi (knockdowns) for overlapping genes. ",
                  "High correlation indicates similar effect sizes between methods."
                ),
                
                # Correlation plot controls
                fluidRow(
                  column(4,
                    selectInput(ns("crispri_experiment"), 
                               "CRISPRi Experiment:",
                               choices = c("C12_FPD-24" = "C12_FPD-24",
                                         "C12_FPD-23" = "C12_FPD-23", 
                                         "C18_FPD-23" = "C18_FPD-23"),
                               selected = "C12_FPD-24")
                  ),
                  column(4,
                    selectInput(ns("gene_filter_approach"), 
                               "Gene Selection:",
                               choices = c(
                                 "Top 200 genes (recommended)" = "top_200",
                                 "Top 100 genes (strongest)" = "top_100", 
                                 "Top 500 genes (broader)" = "top_500",
                                 "All genes (current)" = "all_genes"
                               ),
                               selected = "top_200")
                  ),
                  column(4,
                    div(style = "margin-top: 25px;",
                      checkboxInput(ns("show_all_experiments"), 
                                   "Show all experiments", 
                                   value = TRUE)
                    )
                  ),
                  column(4,
                    div(style = "margin-top: 25px;",
                      checkboxInput(ns("cluster_grid_view"), 
                                   "Cluster-specific view (vertical)", 
                                   value = TRUE)
                    )
                  )
                ),
                
                # Gene filtering explanation
                div(class = "alert alert-info", style = "margin-top: 10px; margin-bottom: 15px; font-size: 0.9em;",
                  icon("info-circle"),
                  strong(" Gene Selection Impact: "), 
                  "Filtering to top N most changed genes dramatically improves correlations ",
                  "(6-11x stronger) by focusing on genes most affected by each method. ",
                  "'All genes' includes thousands of unchanged genes that dilute the signal."
                ),
                
                # Correlation plot explanation
                div(class = "alert alert-success", style = "margin-bottom: 15px; font-size: 0.9em;",
                  icon("chart-line"),
                  strong(" What you're seeing: "), 
                  "Each point = one gene. X-axis = MAST log2FC (mutation effect). ",
                  "Y-axis = CRISPRi log2FC (knockdown effect). ",
                  "Hover over points to see gene names. Strong correlations show similar effects ",
                  "between genetic mutation and CRISPR knockdown."
                ),
                plotlyOutput(ns("gene_pair_correlation"), height = "700px")
              )
            )
          ),
          
          # Interactive heatmap tab
          tabPanel("Signature Heatmap", value = "heatmap",
            wellPanel(
              h4("Interactive Signature Strength Heatmap"),
              
              # Explanatory note about pathway-focused analysis
              div(class = "alert alert-info", style = "margin-bottom: 15px; font-size: 0.9em;",
                icon("info-circle"),
                strong(" Analysis Focus: "), 
                "This shows genes that appear in ANY enriched pathway for MAST ", tags$em("AND"), " any enriched pathway for CRISPRi (same cluster). ",
                "Includes all enrichment databases (GO, KEGG, Reactome, etc.). ",
                "For ALL DE gene overlaps (larger numbers), see Summary Statistics in DE Genes page."
              ),
              
              # Heatmap controls
              fluidRow(
                column(3,
                  selectInput(ns("heatmap_metric"), "Display Metric:",
                             choices = c("Signature Strength" = "signature_strength",
                                       "Gene Overlap Count" = "gene_overlap_count", 
                                       "Fisher p-value (FDR)" = "intersection_fisher_p_fdr_enhanced_hierarchical",
                                       "Fisher p-value (raw)" = "gene_fisher_p",
                                       "Jaccard Index" = "gene_jaccard"),
                             selected = "intersection_fisher_p_fdr_enhanced_hierarchical")
                ),
                column(3,
                  selectInput(ns("heatmap_clustering"), "Clustering:",
                             choices = c("Both" = "both", "Rows" = "row", 
                                       "Columns" = "column", "None" = "none"),
                             selected = "both")
                ),
                column(3,
                  selectInput(ns("color_scale"), "Color Scale:",
                             choices = c("Red" = "Reds", "Viridis" = "viridis", "RdBu" = "RdBu", 
                                       "Blue" = "Blues"),
                             selected = "Reds")
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
    
    # Capture namespace function for use in reactive contexts
    ns <- session$ns
    
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
        unique_clusters <- natural_sort_clusters(unique_clusters[!is.na(unique_clusters)])
        values$available_clusters <- unique_clusters
        
        # Update cluster selection UI
        cluster_choices <- setNames(unique_clusters, paste("Cluster", gsub("cluster_", "", unique_clusters)))
        
        # Select ALL clusters by default for comprehensive analysis
        all_clusters <- unique_clusters  # All available clusters selected by default
        
        updateCheckboxGroupInput(session, "cluster_selection",
                                choices = cluster_choices,
                                selected = all_clusters)
        
        # Initialize gene pair options - analyze variants independently by default
        gene_pairs <- get_comparable_gene_pairs(combine_snca_variants = FALSE,
                                               combine_vps13c_variants = FALSE,
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
    
    # Update gene pairs when variant combination checkboxes change
    observeEvent(c(input$combine_snca, input$combine_vps13c), {
      req(app_data$data_loaded)
      
      # Regenerate gene pairs based on current checkbox states
      gene_pairs <- get_comparable_gene_pairs(
        combine_snca_variants = input$combine_snca %||% FALSE,
        combine_vps13c_variants = input$combine_vps13c %||% FALSE,
        include_mast_only = FALSE
      )
      values$gene_pairs <- gene_pairs
      
      # Update gene selection UI with new pairs
      gene_choices <- c("All available pairs" = "all")
      for (i in seq_len(nrow(gene_pairs))) {
        pair_name <- paste(gene_pairs$mast_gene[i], "vs", gene_pairs$crispri_gene[i])
        gene_choices[[pair_name]] <- paste0(gene_pairs$mast_gene[i], "_vs_", gene_pairs$crispri_gene[i])
      }
      
      updateSelectInput(session, "gene_selection",
                       choices = gene_choices,
                       selected = "all")
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
      
      # Load DE data for proper Fisher's exact test background genes
      de_data <- NULL
      data_dir <- Sys.getenv("ISCORE_DATA_DIR", "")
      de_file_path <- file.path(data_dir, "full_DE_results.rds")
      
      if (file.exists(de_file_path)) {
        tryCatch({
          cat("[SIGNATURE ANALYSIS] Loading DE data for proper background genes...\n")
          de_data <- readRDS(de_file_path)
          cat("[SIGNATURE ANALYSIS] ✓ Loaded DE data - proper Fisher's exact tests enabled\n")
        }, error = function(e) {
          cat("[SIGNATURE ANALYSIS] ⚠ Warning: Could not load DE data:", e$message, "\n")
          de_data <- NULL
        })
      } else {
        cat("[SIGNATURE ANALYSIS] ⚠ Warning: DE data file not found, using legacy background calculation\n")
      }
      
      # Run signature discovery
      tryCatch({
        signature_results <- discover_top_signatures(
          enrichment_data = filtered_data,
          de_data = de_data,
          top_n = input$top_signatures %||% 20,
          min_cluster_breadth = input$min_cluster_breadth %||% 8,
          combine_variants = input$combine_snca && input$combine_vps13c,
          use_enhanced_analysis = input$use_enhanced_analysis %||% TRUE,
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
        
        progress$set(message = "Running PD biological analysis...", value = 0.9)
        
        # Run PD-focused biological interpretation
        pd_analysis <- NULL
        cat("[PD ANALYSIS] Checking analysis scope...\n")
        cat("[PD ANALYSIS] Analysis scope:", paste(input$analysis_scope, collapse = ", "), "\n")
        
        if ("pd_focus" %in% input$analysis_scope) {
          cat("[PD ANALYSIS] PD focus enabled, starting analysis...\n")
          tryCatch({
            # Check if function exists
            if (!exists("analyze_pd_signatures", mode = "function")) {
              cat("[PD ANALYSIS] ERROR: analyze_pd_signatures function not found!\n")
              pd_analysis <- NULL
            } else {
              cat("[PD ANALYSIS] Function found, running analysis...\n")
              pd_analysis <- analyze_pd_signatures(
                signature_results = signature_results,
                enrichment_data = filtered_data,
                focus_on_pan_cluster = TRUE
              )
              cat("[PD ANALYSIS] Analysis completed, checking results...\n")
              if (!is.null(pd_analysis)) {
                cat("[PD ANALYSIS] Success! Enhanced signatures:", length(pd_analysis$enhanced_signatures), "\n")
              } else {
                cat("[PD ANALYSIS] Warning: Analysis returned NULL\n")
              }
            }
          }, error = function(e) {
            cat("[PD ANALYSIS] Error in PD analysis:", e$message, "\n")
            pd_analysis <<- NULL
          })
        } else {
          cat("[PD ANALYSIS] PD focus not in analysis scope, skipping\n")
        }
        
        # Combine results
        final_results <- signature_results
        final_results$pd_analysis <- pd_analysis
        final_results$de_data <- de_data  # Add DE data for correlation analysis
        
        progress$set(message = "Analysis complete!", value = 1.0)
        
        return(final_results)
        
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
        
        # Debug: Log what results we got
        cat("[RESULTS DEBUG] Analysis completed successfully\n")
        cat("[RESULTS DEBUG] Results structure:", paste(names(results), collapse = ", "), "\n")
        
        if ("analysis_summary" %in% names(results)) {
          summary_stats <- results$analysis_summary
          cat("[RESULTS DEBUG] Summary stats:", paste(names(summary_stats), collapse = ", "), "\n")
          cat("[RESULTS DEBUG] Total signatures:", summary_stats$total_signatures %||% "NULL", "\n")
          cat("[RESULTS DEBUG] Pan-cluster count:", summary_stats$pan_cluster_count %||% "NULL", "\n")
          cat("[RESULTS DEBUG] Strongest gene pair:", summary_stats$strongest_gene_pair %||% "NULL", "\n")
        }
        
        if ("top_signatures" %in% names(results)) {
          cat("[RESULTS DEBUG] Top signatures data frame rows:", nrow(results$top_signatures), "\n")
          if (nrow(results$top_signatures) > 0) {
            cat("[RESULTS DEBUG] Top signature columns:", paste(names(results$top_signatures), collapse = ", "), "\n")
            cat("[RESULTS DEBUG] Top signature strength:", safe_max_signature_strength(results$top_signatures), "\n")
          }
        }
        
        # Update UI elements with results
        update_results_ui()
        
        showNotification("Signature analysis completed successfully!", type = "message")
      } else {
        values$analysis_running <- FALSE
        cat("[RESULTS DEBUG] Analysis returned NULL results\n")
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
        # Set PRKN_vs_PARK2 as default if available (best statistical correlation)
        default_pair <- if ("PRKN_vs_PARK2" %in% gene_pair_choices) {
          "PRKN_vs_PARK2"
        } else {
          gene_pair_choices[1]
        }
        updateSelectInput(session, "selected_gene_pair",
                         choices = gene_pair_choices,
                         selected = default_pair)
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
    
    # Render pan-cluster table with improved column names and tooltips
    output$pan_cluster_table <- DT::renderDataTable({
      req(values$analysis_results)
      req(values$analysis_results$pan_cluster_signatures)
      
      pan_cluster_data <- values$analysis_results$pan_cluster_signatures
      
      if (nrow(pan_cluster_data) > 0) {
        # Build display columns dynamically based on available data
        display_cols <- c("gene_pair", "cluster_count", "mean_signature_strength", 
                         "total_gene_overlaps", "total_pathway_overlaps")
        new_names <- c("Gene Pair (MAST vs CRISPRi)", 
                      "Shared Clusters", 
                      "Avg Signature Score", 
                      "Overlapping DE Genes",
                      "Shared Pathways")
        
        # Add gene and pathway Fisher's p-values if available (from all_signatures data)
        all_sigs <- values$analysis_results$all_signatures
        if (!is.null(all_sigs) && nrow(all_sigs) > 0) {
          # Calculate average Fisher's p-values for each gene pair
          pan_cluster_data$avg_gene_fisher_p <- sapply(pan_cluster_data$gene_pair, function(pair) {
            pair_sigs <- all_sigs[all_sigs$gene_pair == pair, ]
            if (nrow(pair_sigs) > 0 && "gene_fisher_p" %in% colnames(pair_sigs)) {
              mean(pair_sigs$gene_fisher_p, na.rm = TRUE)
            } else {
              NA
            }
          })
          
          pan_cluster_data$avg_pathway_fisher_p <- sapply(pan_cluster_data$gene_pair, function(pair) {
            pair_sigs <- all_sigs[all_sigs$gene_pair == pair, ]
            if (nrow(pair_sigs) > 0 && "pathway_fisher_p" %in% colnames(pair_sigs)) {
              mean(pair_sigs$pathway_fisher_p, na.rm = TRUE)
            } else {
              NA
            }
          })
          
          pan_cluster_data$avg_jaccard <- sapply(pan_cluster_data$gene_pair, function(pair) {
            pair_sigs <- all_sigs[all_sigs$gene_pair == pair, ]
            if (nrow(pair_sigs) > 0 && "gene_jaccard" %in% colnames(pair_sigs)) {
              mean(pair_sigs$gene_jaccard, na.rm = TRUE)
            } else {
              NA
            }
          })
          
          display_cols <- c(display_cols, "avg_gene_fisher_p", "avg_pathway_fisher_p", "avg_jaccard")
          new_names <- c(new_names, "Gene p-value", "Pathway p-value", "Avg Jaccard")
        }
        
        display_data <- pan_cluster_data[, display_cols, drop = FALSE]
        colnames(display_data) <- new_names
        
        # Add interpretation columns for both gene and pathway significance
        if ("Gene p-value" %in% colnames(display_data)) {
          display_data$`Gene Sig` <- ifelse(
            is.na(display_data$`Gene p-value`), "n/a",
            ifelse(display_data$`Gene p-value` < 0.001, "***",
                   ifelse(display_data$`Gene p-value` < 0.01, "**",
                          ifelse(display_data$`Gene p-value` < 0.05, "*", "ns")))
          )
        }
        
        if ("Pathway p-value" %in% colnames(display_data)) {
          display_data$`Pathway Sig` <- ifelse(
            is.na(display_data$`Pathway p-value`), "n/a",
            ifelse(display_data$`Pathway p-value` < 0.001, "***",
                   ifelse(display_data$`Pathway p-value` < 0.01, "**",
                          ifelse(display_data$`Pathway p-value` < 0.05, "*", "ns")))
          )
        }
        
        # Create the base datatable
        dt_table <- DT::datatable(display_data, 
                     options = list(
                       pageLength = 10, 
                       scrollX = TRUE,
                       columnDefs = list(
                         list(targets = "_all", className = 'dt-center'),
                         list(targets = c(0), className = 'dt-left')  # Gene pair column left-aligned
                       )
                     ),
                     rownames = FALSE,
                     caption = "Pan-cluster signatures: Focus on 'Pathway Sig' column - pathway overlap is more biologically meaningful than exact gene overlap for cross-method comparisons. *** p<0.001, ** p<0.01, * p<0.05, ns = not significant")
        
        # Apply formatting only for columns that exist
        round_cols <- intersect(c("Avg Signature Score", "Avg Jaccard"), colnames(display_data))
        if (length(round_cols) > 0) {
          dt_table <- dt_table %>% DT::formatRound(round_cols, digits = 2)
        }
        
        signif_cols <- intersect(c("Gene p-value", "Pathway p-value"), colnames(display_data))
        if (length(signif_cols) > 0) {
          dt_table <- dt_table %>% DT::formatSignif(signif_cols, digits = 3)
        }
        
        dt_table %>%
          DT::formatStyle("Pathway Sig", 
                         backgroundColor = DT::styleEqual(
                           c("***", "**", "*", "ns", "n/a"),
                           c("#d4edda", "#d1ecf1", "#fff3cd", "#f8d7da", "#e2e3e5")
                         )) %>%
          DT::formatStyle("Gene Sig",
                         backgroundColor = DT::styleEqual(
                           c("***", "**", "*", "ns", "n/a"),
                           c("#e8f5e8", "#e6f3ff", "#fffacd", "#ffe6e6", "#f5f5f5")
                         ))
      } else {
        DT::datatable(data.frame(Message = "No pan-cluster signatures found"),
                     options = list(dom = 't'), rownames = FALSE)
      }
    })
    
    # === PD BIOLOGICAL ANALYSIS OUTPUTS (NEW) ===
    
    # PD biological summary
    output$pd_biological_summary <- renderUI({
      req(values$analysis_results)
      req(values$analysis_results$pd_analysis)
      
      pd_analysis <- values$analysis_results$pd_analysis
      
      if (is.null(pd_analysis) || is.null(pd_analysis$pd_summary)) {
        return(div(
          class = "alert alert-info",
          h5("PD Analysis Not Available"),
          p("PD-focused analysis was not performed or failed. Ensure 'PD-relevant focus' is selected in analysis scope.")
        ))
      }
      
      summary_stats <- pd_analysis$pd_summary$summary_stats
      biological_categories <- pd_analysis$pd_summary$biological_categories
      
      # Get top 3 biological categories
      category_ranking <- biological_categories[order(unlist(biological_categories), decreasing = TRUE)]
      top_categories <- head(names(category_ranking)[unlist(category_ranking) > 0], 3)
      
      tagList(
        div(class = "row",
          div(class = "col-md-6",
            h6("Analysis Overview", style = "color: #2c3e50; font-weight: bold;"),
            tags$ul(
              tags$li(paste("Analysis Type:", str_to_title(summary_stats$analysis_type))),
              tags$li(paste("Signatures Analyzed:", summary_stats$total_signatures)),
              tags$li(paste("Mean PD Relevance:", round(summary_stats$mean_pd_relevance, 2)))
            )
          ),
          div(class = "col-md-6",
            h6("Top Biological Categories", style = "color: #d73027; font-weight: bold;"),
            if (length(top_categories) > 0) {
              tags$ol(
                lapply(seq_along(top_categories), function(i) {
                  cat_name <- top_categories[i]
                  count <- biological_categories[[cat_name]]
                  display_name <- str_to_title(gsub("_", " ", cat_name))
                  tags$li(paste0(display_name, " (", count, " occurrences)"))
                })
              )
            } else {
              p("No major biological categories identified", style = "color: #666;")
            }
          )
        ),
        div(style = "margin-top: 15px;",
          h6("Key Insights", style = "color: #27ae60; font-weight: bold;"),
          div(class = "well well-sm", style = "background-color: #f1f8e9;",
            if (!is.null(summary_stats$most_relevant_signature) && summary_stats$most_relevant_signature != "None") {
              p(paste("🔬 Most PD-relevant signature:", summary_stats$most_relevant_signature))
            },
            if (length(top_categories) > 0) {
              p(paste("🧠 Dominant biological theme:", str_to_title(gsub("_", " ", summary_stats$top_biological_category))))
            }
          )
        )
      )
    })
    
    # PD signatures table
    output$pd_signatures_table <- DT::renderDataTable({
      req(values$analysis_results)
      req(values$analysis_results$pd_analysis)
      
      pd_analysis <- values$analysis_results$pd_analysis
      
      if (is.null(pd_analysis$enhanced_signatures) || length(pd_analysis$enhanced_signatures) == 0) {
        return(DT::datatable(
          data.frame(Message = "No PD-relevant signatures found"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      # Create summary table from enhanced signatures
      table_data <- data.frame()
      for (i in seq_along(pd_analysis$enhanced_signatures)) {
        sig <- pd_analysis$enhanced_signatures[[i]]
        
        # Handle different column names for signature strength
        sig_strength <- if (!is.null(sig$signature$signature_strength)) {
          sig$signature$signature_strength
        } else if (!is.null(sig$signature$mean_signature_strength)) {
          sig$signature$mean_signature_strength
        } else if (!is.null(sig$signature$max_signature_strength)) {
          sig$signature$max_signature_strength
        } else {
          NA
        }
        
        # Handle missing cluster info for pan-cluster signatures
        cluster_info <- if (!is.null(sig$signature$cluster)) {
          sig$signature$cluster
        } else if (!is.null(sig$signature$cluster_count)) {
          paste0("Pan-cluster (", sig$signature$cluster_count, " clusters)")
        } else {
          "Pan-cluster"
        }
        
        row_data <- data.frame(
          Rank = i,
          Gene_Pair = sig$signature$gene_pair,
          Cluster = cluster_info,
          Signature_Strength = if(!is.na(sig_strength)) round(sig_strength, 2) else NA,
          PD_Relevance = round(sig$pd_relevance_score, 2),
          Shared_PD_Pathways = nrow(sig$shared_pd_pathways),
          Mitochondrial = sig$biological_categories$mitochondrial,
          Protein_Quality = sig$biological_categories$protein_quality,
          Autophagy = sig$biological_categories$autophagy,
          Dopamine = sig$biological_categories$dopamine,
          stringsAsFactors = FALSE
        )
        table_data <- rbind(table_data, row_data)
      }
      
      # Sort by PD relevance
      table_data <- table_data[order(-table_data$PD_Relevance, -table_data$Signature_Strength), ]
      
      DT::datatable(table_data,
                   options = list(pageLength = 10, scrollX = TRUE, autoWidth = TRUE),
                   rownames = FALSE) %>%
        DT::formatRound(c("Signature_Strength", "PD_Relevance"), digits = 2) %>%
        DT::formatStyle("PD_Relevance",
                       background = DT::styleColorBar(range(table_data$PD_Relevance), "#fee0d2"),
                       backgroundSize = "100% 90%",
                       backgroundRepeat = "no-repeat",
                       backgroundPosition = "center")
    })
    
    # PD pathway frequency plot
    output$pd_pathway_frequency <- renderPlotly({
      req(values$analysis_results)
      req(values$analysis_results$pd_analysis)
      
      pd_analysis <- values$analysis_results$pd_analysis
      
      if (is.null(pd_analysis$enhanced_signatures) || length(pd_analysis$enhanced_signatures) == 0) {
        return(plotly_empty("No PD pathway data available"))
      }
      
      # Aggregate pathway information across all signatures
      all_pathways <- data.frame()
      for (sig in pd_analysis$enhanced_signatures) {
        if (nrow(sig$shared_pd_pathways) > 0) {
          pathway_data <- sig$shared_pd_pathways
          pathway_data$gene_pair <- sig$signature$gene_pair
          all_pathways <- rbind(all_pathways, pathway_data)
        }
      }
      
      if (nrow(all_pathways) == 0) {
        return(plotly_empty("No shared PD pathways found"))
      }
      
      # Create frequency table
      pathway_frequency <- as.data.frame(table(all_pathways$pathway))
      names(pathway_frequency) <- c("Pathway", "Frequency")
      pathway_frequency <- pathway_frequency[order(-pathway_frequency$Frequency), ]
      pathway_frequency <- head(pathway_frequency, 15)  # Top 15 pathways
      
      # Create plot
      p <- ggplot(pathway_frequency, aes(x = reorder(Pathway, Frequency), y = Frequency)) +
        geom_col(fill = "#d73027", alpha = 0.8) +
        coord_flip() +
        labs(title = "Most Frequently Disrupted PD-Relevant Pathways",
             x = "Pathway", y = "Frequency Across Signatures") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 10))
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    # PD categories plot
    output$pd_categories_plot <- renderPlotly({
      req(values$analysis_results)
      req(values$analysis_results$pd_analysis)
      
      pd_analysis <- values$analysis_results$pd_analysis
      
      if (is.null(pd_analysis$pd_summary$biological_categories)) {
        return(plotly_empty("No biological category data available"))
      }
      
      categories <- pd_analysis$pd_summary$biological_categories
      
      # Convert to data frame
      cat_data <- data.frame(
        Category = names(categories),
        Count = unlist(categories),
        stringsAsFactors = FALSE
      )
      
      # Remove categories with 0 count
      cat_data <- cat_data[cat_data$Count > 0, ]
      
      if (nrow(cat_data) == 0) {
        return(plotly_empty("No biological categories with counts > 0"))
      }
      
      # Clean up category names for display
      cat_data$Category_Display <- str_to_title(gsub("_", " ", cat_data$Category))
      cat_data <- cat_data[order(-cat_data$Count), ]
      
      # Create plot
      p <- ggplot(cat_data, aes(x = reorder(Category_Display, Count), y = Count)) +
        geom_col(fill = "#2166ac", alpha = 0.8) +
        coord_flip() +
        labs(title = "Biological Process Categories in PD Signatures",
             x = "Biological Process", y = "Pathway Count") +
        theme_minimal() +
        theme(axis.text.y = element_text(size = 11))
      
      ggplotly(p, tooltip = c("x", "y"))
    })
    
    
    # === CLUSTER-SPECIFIC SIGNATURES OUTPUTS ===
    
    # Cluster-specific signatures table
    output$cluster_specific_table <- DT::renderDataTable({
      req(values$analysis_results)
      req(input$selected_cluster_detail)
      
      cluster_sigs <- values$analysis_results$cluster_specific_signatures
      selected_cluster <- input$selected_cluster_detail
      
      if (is.null(cluster_sigs) || length(cluster_sigs) == 0 || 
          !selected_cluster %in% names(cluster_sigs)) {
        return(DT::datatable(
          data.frame(Message = "No cluster-specific signatures found for this cluster"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      cluster_data <- cluster_sigs[[selected_cluster]]
      
      if (nrow(cluster_data) > 0) {
        # Build columns dynamically based on available data
        display_cols <- c("gene_pair")
        col_names <- c("Gene Pair (MAST vs CRISPRi)")
        
        # Add signature strength if available
        strength_col <- if ("signature_strength" %in% colnames(cluster_data)) {
          "signature_strength"
        } else if ("mean_signature_strength" %in% colnames(cluster_data)) {
          "mean_signature_strength"
        } else {
          NULL
        }
        
        if (!is.null(strength_col)) {
          display_cols <- c(display_cols, strength_col)
          col_names <- c(col_names, "Signature Score")
        }
        
        # Add overlap counts and statistics
        display_cols <- c(display_cols, "gene_overlap_count", "pathway_overlap_count")
        col_names <- c(col_names, "DE Genes Shared", "Pathways Shared")
        
        # Add Fisher's tests for both genes and pathways
        if ("gene_fisher_p" %in% colnames(cluster_data)) {
          display_cols <- c(display_cols, "gene_fisher_p")
          col_names <- c(col_names, "Gene p-value")
        }
        
        if ("pathway_fisher_p" %in% colnames(cluster_data)) {
          display_cols <- c(display_cols, "pathway_fisher_p")
          col_names <- c(col_names, "Pathway p-value")
        }
        
        # Add Jaccard index if available
        if ("gene_jaccard" %in% colnames(cluster_data)) {
          display_cols <- c(display_cols, "gene_jaccard")
          col_names <- c(col_names, "Jaccard Index")
        }
        
        display_data <- cluster_data[, display_cols, drop = FALSE]
        colnames(display_data) <- col_names
        
        # Add significance interpretations
        if ("Gene p-value" %in% colnames(display_data)) {
          display_data$`Gene Sig` <- ifelse(
            is.na(display_data$`Gene p-value`), "n/a",
            ifelse(display_data$`Gene p-value` < 0.001, "***",
                   ifelse(display_data$`Gene p-value` < 0.01, "**",
                          ifelse(display_data$`Gene p-value` < 0.05, "*", "ns")))
          )
        }
        
        if ("Pathway p-value" %in% colnames(display_data)) {
          display_data$`Pathway Sig` <- ifelse(
            is.na(display_data$`Pathway p-value`), "n/a",
            ifelse(display_data$`Pathway p-value` < 0.001, "***",
                   ifelse(display_data$`Pathway p-value` < 0.01, "**",
                          ifelse(display_data$`Pathway p-value` < 0.05, "*", "ns")))
          )
        }
        
        DT::datatable(display_data,
                     options = list(
                       pageLength = 10, 
                       scrollX = TRUE,
                       order = list(list(1, 'desc')),  # Sort by signature score if available
                       columnDefs = list(
                         list(targets = "_all", className = 'dt-center'),
                         list(targets = c(0), className = 'dt-left')
                       )
                     ),
                     rownames = FALSE,
                     caption = paste("Cluster", gsub("cluster_", "", selected_cluster), 
                                   "signatures. Gene overlap often non-significant (expected). Focus on pathway overlap for biological convergence.")) %>%
          DT::formatRound(c("Signature Score", "Jaccard Index"), digits = 2) %>%
          DT::formatSignif(c("Gene p-value", "Pathway p-value"), digits = 3) %>%
          DT::formatStyle("Pathway Sig", 
                         backgroundColor = DT::styleEqual(
                           c("***", "**", "*", "ns", "n/a"),
                           c("#d4edda", "#d1ecf1", "#fff3cd", "#f8d7da", "#e2e3e5")
                         )) %>%
          DT::formatStyle("Gene Sig",
                         backgroundColor = DT::styleEqual(
                           c("***", "**", "*", "ns", "n/a"),
                           c("#e8f5e8", "#e6f3ff", "#fffacd", "#ffe6e6", "#f5f5f5")
                         ))
      } else {
        DT::datatable(data.frame(Message = "No signatures found for this cluster"),
                     options = list(dom = 't'), rownames = FALSE)
      }
    })
    
    # === GENE PAIR ANALYSIS OUTPUTS ===
    
    # Gene pair analysis table  
    output$gene_pair_table <- DT::renderDataTable({
      req(values$analysis_results)
      req(input$selected_gene_pair)
      
      # Get analysis results with better validation
      analysis_results <- values$analysis_results
      
      # Debug what we actually have
      cat("[GENE PAIR DEBUG] Analysis results class:", class(analysis_results), "\n")
      if (is.list(analysis_results)) {
        cat("[GENE PAIR DEBUG] Analysis results elements:", names(analysis_results), "\n")
      }
      
      # Get all_signatures with robust error checking
      all_sigs <- NULL
      if (is.list(analysis_results) && "all_signatures" %in% names(analysis_results)) {
        all_sigs <- analysis_results$all_signatures
      }
      
      selected_pair <- input$selected_gene_pair
      
      # Better error handling for data structure
      if (is.null(all_sigs)) {
        cat("[GENE PAIR DEBUG] all_signatures is NULL\n")
        return(DT::datatable(
          data.frame(Message = "No signature data available - all_signatures is NULL"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      if (!is.data.frame(all_sigs)) {
        cat("[GENE PAIR DEBUG] all_signatures is not a data frame, class:", class(all_sigs), "\n")
        return(DT::datatable(
          data.frame(Message = paste("Signature data is not a data frame, got:", class(all_sigs))),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      if (nrow(all_sigs) == 0) {
        cat("[GENE PAIR DEBUG] all_signatures data frame is empty\n")
        return(DT::datatable(
          data.frame(Message = "No signature data available - empty data frame"),
          options = list(dom = 't'), rownames = FALSE
        ))
      }
      
      cat("[GENE PAIR DEBUG] all_signatures is valid data frame with", nrow(all_sigs), "rows and columns:", paste(names(all_sigs), collapse = ", "), "\n")
      cat("[GENE PAIR DEBUG] Selected gene pair:", selected_pair, "\n")
      cat("[GENE PAIR DEBUG] Available gene pairs:", paste(unique(all_sigs$gene_pair), collapse = ", "), "\n")
      
      pair_data <- all_sigs[all_sigs$gene_pair == selected_pair, ]
      cat("[GENE PAIR DEBUG] Filtered pair_data rows:", nrow(pair_data), "\n")
      
      if (nrow(pair_data) > 0) {
        # Handle different column names for signature strength
        strength_col <- if ("signature_strength" %in% colnames(pair_data)) {
          "signature_strength"
        } else if ("mean_signature_strength" %in% colnames(pair_data)) {
          "mean_signature_strength"
        } else if ("max_signature_strength" %in% colnames(pair_data)) {
          "max_signature_strength"
        } else {
          NULL
        }
        
        # Build display columns dynamically
        display_cols <- c("cluster")
        col_names <- c("Cluster")
        
        if (!is.null(strength_col)) {
          display_cols <- c(display_cols, strength_col)
          col_names <- c(col_names, "Signature Score")
        }
        
        # Check if we have new intersection/union Fisher's test results
        has_intersection_union <- all(c("intersection_fisher_p", "union_fisher_p", 
                                       "intersection_background_size", "union_background_size") %in% colnames(pair_data))
        
        # Determine whether to show FDR-corrected or raw p-values (check for enhanced hierarchical columns)
        has_fdr_correction <- any(c("intersection_fisher_p_fdr_hierarchical", "intersection_fisher_p_fdr_enhanced_hierarchical") %in% names(pair_data))
        
        # Check if we have enhanced direction analysis results (v0.2.6)
        has_enhanced_analysis <- any(c("biological_expectation", "primary_direction_pattern", 
                                     "same_direction_count", "opposite_direction_count") %in% colnames(pair_data))
        
        # Get user's selected approach (intersection vs union)
        selected_approach <- input$analysis_approach %||% "intersection"
        
        if (has_intersection_union) {
          # PRIORITIZE SIGNIFICANCE COLUMNS FIRST (after cluster)
          
          if (selected_approach == "intersection") {
            # Show intersection approach - prioritize FDR when available
            if ("intersection_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
              display_cols <- c(display_cols, "intersection_fisher_p_fdr_enhanced_hierarchical", "intersection_background_size")
              col_names <- c(col_names, "p-value (FDR)", "Background Genes")
            } else if ("intersection_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
              display_cols <- c(display_cols, "intersection_fisher_p_fdr_hierarchical", "intersection_background_size")
              col_names <- c(col_names, "p-value (FDR)", "Background Genes")
            } else {
              display_cols <- c(display_cols, "intersection_fisher_p", "intersection_background_size")
              col_names <- c(col_names, "p-value", "Background Genes")
            }
          } else {
            # Show union approach - prioritize FDR when available  
            if ("union_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
              display_cols <- c(display_cols, "union_fisher_p_fdr_enhanced_hierarchical", "union_background_size")  
              col_names <- c(col_names, "p-value (FDR)", "Background Genes")
            } else if ("union_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
              display_cols <- c(display_cols, "union_fisher_p_fdr_hierarchical", "union_background_size")  
              col_names <- c(col_names, "p-value (FDR)", "Background Genes")
            } else {
              display_cols <- c(display_cols, "union_fisher_p", "union_background_size")  
              col_names <- c(col_names, "p-value", "Background Genes")
            }
          }
          
          # Add supporting data columns after significance
          display_cols <- c(display_cols, "gene_overlap_count", "pathway_overlap_count", "gene_jaccard")
          col_names <- c(col_names, "Shared DE Genes*", "Enriched Pathways", "Jaccard Index")
          
          # Add enhanced direction analysis indicators (v0.2.6)
          if (has_enhanced_analysis) {
            if ("biological_expectation" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "biological_expectation")
              col_names <- c(col_names, "Direction Expectation")
            }
            if ("primary_direction_pattern" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "primary_direction_pattern") 
              col_names <- c(col_names, "Detected Pattern")
            }
            if ("same_direction_count" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "same_direction_count")
              col_names <- c(col_names, "Same Direction")
            }
            if ("opposite_direction_count" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "opposite_direction_count")
              col_names <- c(col_names, "Opposite Direction")
            }
          }
          
        } else {
          # LEGACY: Display old single approach for backwards compatibility
          
          # Prioritize significance column first
          if (has_fdr_correction && "gene_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
            display_cols <- c(display_cols, "gene_fisher_p_fdr_hierarchical")
            col_names <- c(col_names, "Significance (FDR)")
          } else {
            display_cols <- c(display_cols, "gene_fisher_p")
            col_names <- c(col_names, "Significance (raw)")
          }
          
          # Add supporting data columns
          display_cols <- c(display_cols, "gene_overlap_count", "pathway_overlap_count", "gene_jaccard")
          col_names <- c(col_names, "Pathway Genes", "Shared Pathways", "Jaccard Index")
          
          # Add enhanced direction analysis indicators (v0.2.6) - Legacy section
          if (has_enhanced_analysis) {
            if ("biological_expectation" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "biological_expectation")
              col_names <- c(col_names, "Direction Expectation")
            }
            if ("primary_direction_pattern" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "primary_direction_pattern") 
              col_names <- c(col_names, "Detected Pattern")
            }
            if ("same_direction_count" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "same_direction_count")
              col_names <- c(col_names, "Same Direction")
            }
            if ("opposite_direction_count" %in% colnames(pair_data)) {
              display_cols <- c(display_cols, "opposite_direction_count")
              col_names <- c(col_names, "Opposite Direction")
            }
          }
          
          # Add background information if available (legacy)
          if ("background_size" %in% colnames(pair_data) && "background_type" %in% colnames(pair_data)) {
            display_cols <- c(display_cols, "background_size", "background_type")
            col_names <- c(col_names, "Background Genes", "Test Type")
          }
        }
        
        # Debug before subsetting
        cat("[GENE PAIR DEBUG] About to subset pair_data with display_cols:", paste(display_cols, collapse = ", "), "\n")
        cat("[GENE PAIR DEBUG] pair_data class before subsetting:", class(pair_data), "\n")
        cat("[GENE PAIR DEBUG] pair_data columns:", paste(names(pair_data), collapse = ", "), "\n")
        
        # Ensure we have a proper data frame before subsetting
        if (!is.data.frame(pair_data)) {
          cat("[GENE PAIR DEBUG] pair_data is not a data frame, converting from:", class(pair_data), "\n")
          pair_data <- as.data.frame(pair_data)
        }
        
        # Check if display_cols exist in pair_data
        missing_cols <- setdiff(display_cols, names(pair_data))
        if (length(missing_cols) > 0) {
          cat("[GENE PAIR DEBUG] Missing columns:", paste(missing_cols, collapse = ", "), "\n")
          return(DT::datatable(
            data.frame(Message = paste("Missing required columns:", paste(missing_cols, collapse = ", "))),
            options = list(dom = 't'), rownames = FALSE
          ))
        }
        
        display_data <- pair_data[, display_cols, drop = FALSE]
        cat("[GENE PAIR DEBUG] display_data created successfully, class:", class(display_data), "\n")
        
        colnames(display_data) <- col_names
        
        # CRITICAL: Validate display_data right before DT operations
        cat("[GENE PAIR DEBUG] display_data class after colnames:", class(display_data), "\n")
        cat("[GENE PAIR DEBUG] display_data dimensions:", dim(display_data), "\n")
        cat("[GENE PAIR DEBUG] display_data is.data.frame:", is.data.frame(display_data), "\n")
        
        # Add significance interpretation based on selected approach and FDR availability
        if (has_intersection_union) {
          # Create significance column based on selected approach
          cat("[FDR DEBUG] Selected approach:", selected_approach, "\n")
          cat("[FDR DEBUG] FDR correction available:", has_fdr_correction, "\n")
          
          if (selected_approach == "intersection") {
            if (has_fdr_correction && "intersection_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
              selected_p_values <- pair_data$intersection_fisher_p_fdr_enhanced_hierarchical[match(display_data$Cluster, pair_data$cluster)]
              cat("[FDR DEBUG] Using intersection FDR-corrected p-values\n")
            } else {
              selected_p_values <- pair_data$intersection_fisher_p[match(display_data$Cluster, pair_data$cluster)]
              cat("[FDR DEBUG] Using intersection raw p-values\n")
            }
          } else {
            if (has_fdr_correction && "union_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
              selected_p_values <- pair_data$union_fisher_p_fdr_enhanced_hierarchical[match(display_data$Cluster, pair_data$cluster)]
              cat("[FDR DEBUG] Using union FDR-corrected p-values\n")
            } else {
              selected_p_values <- pair_data$union_fisher_p[match(display_data$Cluster, pair_data$cluster)]
              cat("[FDR DEBUG] Using union raw p-values\n")
            }
          }
          
          # Create single significance column for selected approach
          display_data$`Sig Level` <- ifelse(
            is.na(selected_p_values), "n/a",
            ifelse(selected_p_values < 0.001, "***",
                   ifelse(selected_p_values < 0.01, "**",
                          ifelse(selected_p_values < 0.05, "*", "ns")))
          )
          
        } else {
          # Legacy significance interpretation
          if (has_fdr_correction && "gene_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
            cat("[FDR DEBUG] Using hierarchical FDR-corrected p-values for legacy significance\n")
            fdr_gene_p <- pair_data$gene_fisher_p_fdr_hierarchical[match(display_data$Cluster, pair_data$cluster)]
            display_data$`Sig Level` <- ifelse(
              is.na(fdr_gene_p), "n/a",
              ifelse(fdr_gene_p < 0.001, "***",
                     ifelse(fdr_gene_p < 0.01, "**",
                            ifelse(fdr_gene_p < 0.05, "*", "ns")))
            )
          } else {
            # Use original data column for raw p-values
            raw_gene_p <- pair_data$gene_fisher_p[match(display_data$Cluster, pair_data$cluster)]
            display_data$`Sig Level` <- ifelse(
              is.na(raw_gene_p), "n/a",
              ifelse(raw_gene_p < 0.001, "***",
                     ifelse(raw_gene_p < 0.01, "**",
                            ifelse(raw_gene_p < 0.05, "*", "ns")))
            )
          }
        }
        
        # FINAL VALIDATION: Right before DT::datatable call
        cat("[GENE PAIR DEBUG] FINAL CHECK - display_data class:", class(display_data), "\n")
        cat("[GENE PAIR DEBUG] FINAL CHECK - display_data is.data.frame:", is.data.frame(display_data), "\n")
        cat("[GENE PAIR DEBUG] FINAL CHECK - display_data nrow:", nrow(display_data), "\n")
        cat("[GENE PAIR DEBUG] FINAL CHECK - display_data ncol:", ncol(display_data), "\n")
        
        # Force display_data to be a proper data.frame
        display_data <- as.data.frame(display_data, stringsAsFactors = FALSE)
        
        # Add a Details button column
        display_data$Details <- paste0(
          '<button class="btn btn-primary btn-sm" onclick="Shiny.setInputValue(\'', 
          ns('show_details'), '\', {cluster: ', 1:nrow(display_data), 
          ', random: Math.random()}); event.stopPropagation();">',
          '<i class="fa fa-info-circle"></i> Details</button>'
        )
        
        cat("[GENE PAIR DEBUG] AFTER as.data.frame - class:", class(display_data), "\n")
        
        # SIMPLIFIED DT call to isolate the issue
        tryCatch({
          # Simple caption to avoid DT::datatable issues
          simple_caption <- paste("Cross-Method Analysis:", selected_pair)
          
          dt_result <- DT::datatable(display_data,
                           options = list(
                             pageLength = 15, 
                             scrollX = TRUE,
                             columnDefs = list(
                               list(targets = which(names(display_data) == "Details") - 1,
                                    orderable = FALSE)
                             )
                           ),
                           rownames = FALSE,
                           caption = simple_caption,
                           escape = FALSE)
          
          cat("[GENE PAIR DEBUG] DT::datatable created successfully\n")
          
          # Apply formatting only if base datatable succeeds
          final_dt <- dt_result %>%
            DT::formatRound(intersect(c("Signature Score", "Jaccard Index"), names(display_data)), digits = 2)
          
          # Add color coding for FDR-corrected significance columns (PRIORITIZE FDR)
          if ("Intersection Sig (FDR)" %in% names(display_data)) {
            final_dt <- final_dt %>%
              DT::formatStyle("Intersection Sig (FDR)",
                            backgroundColor = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("#dc3545", "#17a2b8", "#28a745", "#e9ecef", "#f8f9fa")  # Red, Blue, Green, Light gray, Very light gray
                            ),
                            color = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("white", "white", "white", "#6c757d", "#6c757d")  # White text for colored cells, gray for others
                            ),
                            fontWeight = styleEqual(
                              c("***", "**", "*"),
                              c("bold", "bold", "bold")
                            ))
          } else if ("Intersection Sig" %in% names(display_data)) {
            final_dt <- final_dt %>%
              DT::formatStyle("Intersection Sig",
                            backgroundColor = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("#dc3545", "#17a2b8", "#28a745", "#e9ecef", "#f8f9fa")  # Red, Blue, Green, Light gray, Very light gray
                            ),
                            color = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("white", "white", "white", "#6c757d", "#6c757d")  # White text for colored cells, gray for others
                            ),
                            fontWeight = styleEqual(
                              c("***", "**", "*"),
                              c("bold", "bold", "bold")
                            ))
          }
          
          # Also add color coding for Union Sig (FDR) if available
          if ("Union Sig (FDR)" %in% names(display_data)) {
            final_dt <- final_dt %>%
              DT::formatStyle("Union Sig (FDR)",
                            backgroundColor = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("#dc3545", "#17a2b8", "#28a745", "#e9ecef", "#f8f9fa")  
                            ),
                            color = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("white", "white", "white", "#6c757d", "#6c757d")  
                            ),
                            fontWeight = styleEqual(
                              c("***", "**", "*"),
                              c("bold", "bold", "bold")
                            ))
          }
          
          # Color coding for legacy FDR significance
          if ("Significance (FDR)" %in% names(display_data)) {
            final_dt <- final_dt %>%
              DT::formatStyle("Significance (FDR)",
                            backgroundColor = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("#dc3545", "#17a2b8", "#28a745", "#e9ecef", "#f8f9fa")  
                            ),
                            color = styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("white", "white", "white", "#6c757d", "#6c757d")  
                            ),
                            fontWeight = styleEqual(
                              c("***", "**", "*"),
                              c("bold", "bold", "bold")
                            ))
          }
          
          # Color coding for unified Sig Level column (CRITICAL FIX - this was missing!)
          if ("Sig Level" %in% names(display_data)) {
            final_dt <- final_dt %>%
              DT::formatStyle("Sig Level",
                            backgroundColor = DT::styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("#dc3545", "#17a2b8", "#28a745", "#e9ecef", "#f8f9fa")  
                            ),
                            color = DT::styleEqual(
                              c("***", "**", "*", "ns", "n/a"),
                              c("white", "white", "white", "#6c757d", "#6c757d")  
                            ),
                            fontWeight = DT::styleEqual(
                              c("***", "**", "*"),
                              c("bold", "bold", "bold")
                            ))
          }
          
          return(final_dt)
          
        }, error = function(e) {
          cat("[GENE PAIR DEBUG] ERROR in DT::datatable:", e$message, "\n")
          # Return a simple fallback datatable
          return(DT::datatable(
            data.frame(Error = paste("Failed to create datatable:", e$message), 
                      DataClass = paste(class(display_data), collapse = ", "),
                      Dimensions = paste(dim(display_data), collapse = " x ")),
            options = list(dom = 't'), rownames = FALSE
          ))
        })
      } else {
        DT::datatable(data.frame(Message = "No data found for this gene pair"),
                     options = list(dom = 't'), rownames = FALSE)
      }
    })
    
    # Gene pair summary output - dynamic summary of significant clusters
    output$gene_pair_summary <- renderUI({
      req(values$analysis_results)
      req(input$selected_gene_pair)
      
      all_sigs <- values$analysis_results$all_signatures
      selected_pair <- input$selected_gene_pair
      
      if (is.null(all_sigs) || !is.data.frame(all_sigs) || nrow(all_sigs) == 0) {
        return(NULL)
      }
      
      pair_data <- all_sigs[all_sigs$gene_pair == selected_pair, ]
      
      if (nrow(pair_data) == 0) {
        return(NULL)
      }
      
      # Check if we have intersection/union data
      has_intersection_union <- all(c("intersection_fisher_p", "union_fisher_p") %in% colnames(pair_data))
      
      if (has_intersection_union) {
        # Determine which approach to use based on user selection
        approach <- input$analysis_approach %||% "intersection"
        
        # Select appropriate p-value column based on approach
        if (approach == "intersection") {
          if ("intersection_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
            p_col <- "intersection_fisher_p_fdr_enhanced_hierarchical"
            p_col_raw <- "intersection_fisher_p"
            cat("[FDR SUMMARY] Using enhanced FDR-corrected intersection p-values for summary\n")
          } else if ("intersection_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
            p_col <- "intersection_fisher_p_fdr_hierarchical"
            p_col_raw <- "intersection_fisher_p"
            cat("[FDR SUMMARY] Using FDR-corrected intersection p-values for summary\n")
          } else {
            p_col <- "intersection_fisher_p"
            p_col_raw <- "intersection_fisher_p"
            cat("[FDR SUMMARY] FDR-corrected p-values not available, using raw intersection p-values\n")
          }
        } else {  # union approach
          if ("union_fisher_p_fdr_enhanced_hierarchical" %in% names(pair_data)) {
            p_col <- "union_fisher_p_fdr_enhanced_hierarchical"
            p_col_raw <- "union_fisher_p"
            cat("[FDR SUMMARY] Using enhanced FDR-corrected union p-values for summary\n")
          } else if ("union_fisher_p_fdr_hierarchical" %in% names(pair_data)) {
            p_col <- "union_fisher_p_fdr_hierarchical"
            p_col_raw <- "union_fisher_p"
            cat("[FDR SUMMARY] Using FDR-corrected union p-values for summary\n")
          } else {
            p_col <- "union_fisher_p"
            p_col_raw <- "union_fisher_p"
            cat("[FDR SUMMARY] FDR-corrected p-values not available, using raw union p-values\n")
          }
        }
        
        # Filter for significant clusters using selected approach
        sig_clusters_strict <- pair_data[!is.na(pair_data[[p_col]]) & 
                                        pair_data[[p_col]] < 0.05, ]
        highly_sig_clusters <- pair_data[!is.na(pair_data[[p_col]]) & 
                                        pair_data[[p_col]] < 0.001, ]
        
        # Parse gene names from the pair
        genes <- strsplit(selected_pair, "_vs_")[[1]]
        mast_gene <- genes[1]
        crispri_gene <- genes[2]
        
        # Create summary message
        if (nrow(sig_clusters_strict) == 0) {
          summary_html <- div(
            class = "alert alert-warning",
            h5(icon("chart-line"), " Overlap Analysis Summary", style = "margin-top: 0;"),
            p(strong("No significant gene overlap detected"), " between ", 
              tags$span(mast_gene, style = "color: #d73027; font-weight: bold;"), " (MAST mutation) and ",
              tags$span(crispri_gene, style = "color: #4575b4; font-weight: bold;"), " (CRISPRi knockdown) ",
              "in any cluster using the ", ifelse(approach == "intersection", "conservative intersection", "liberal union"), " approach."),
            p(tags$em("This is common and expected - check the pathway overlap results which may still show biological convergence."))
          )
        } else {
          # Build cluster summary
          cluster_summary <- lapply(1:nrow(sig_clusters_strict), function(i) {
            row <- sig_clusters_strict[i, ]
            cluster_name <- gsub("cluster_", "", row$cluster)
            p_val <- row[[p_col_raw]]  # Use the raw p-value for display (not FDR-corrected)
            overlap_count <- row$gene_overlap_count
            
            sig_level <- if (p_val < 0.001) "***" else if (p_val < 0.01) "**" else "*"
            sig_color <- if (p_val < 0.001) "#dc3545" else if (p_val < 0.01) "#17a2b8" else "#28a745"
            
            tags$li(
              tags$strong(paste0("Cluster ", cluster_name, ":")),
              " ", overlap_count, " overlapping DE genes ",
              tags$span(paste0("(p = ", formatC(p_val, format = "e", digits = 2), " ", sig_level, ")"),
                       style = paste0("color: ", sig_color, "; font-weight: bold;"))
            )
          })
          
          summary_html <- div(
            class = "alert alert-success",
            h5(icon("check-circle"), " Significant Gene Overlap Detected!", style = "margin-top: 0;"),
            p(strong(nrow(sig_clusters_strict), " cluster(s)"), " show significant overlap between ",
              tags$span(mast_gene, style = "color: #d73027; font-weight: bold;"), " (MAST mutation) and ",
              tags$span(crispri_gene, style = "color: #4575b4; font-weight: bold;"), " (CRISPRi knockdown):"),
            tags$ul(cluster_summary),
            if (nrow(highly_sig_clusters) > 0) {
              p(tags$strong(nrow(highly_sig_clusters), " cluster(s)"), 
                " show ", tags$span("highly significant (***)", style = "color: #dc3545; font-weight: bold;"), 
                " overlap, suggesting robust cross-method agreement.")
            } else NULL
          )
        }
      } else {
        # Legacy data format
        summary_html <- div(
          class = "alert alert-info",
          p("Detailed overlap statistics not available for this analysis. Re-run the analysis for enhanced statistics.")
        )
      }
      
      return(summary_html)
    })
    
    # === HEATMAP OUTPUTS ===
    
    # Interactive signature heatmap (MISSING IMPLEMENTATION - This was causing the loading wheel!)
    output$signature_heatmap <- renderPlotly({
      req(values$analysis_results)
      req(values$analysis_results$all_signatures)
      
      cat("[SIGNATURE HEATMAP] Generating heatmap with metric:", input$heatmap_metric, "\n")
      
      # Get signature data
      signature_data <- values$analysis_results$all_signatures
      
      if (nrow(signature_data) == 0) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, text = "No signature data available for heatmap",
                               textfont = list(size = 16, color = "gray")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      tryCatch({
        # Create heatmap using the visualization function with all UI parameters
        heatmap_plot <- create_interactive_signature_heatmap_enhanced(
          signature_data = signature_data,
          metric = input$heatmap_metric %||% "signature_strength",
          cluster_filter = NULL,  # Include all clusters
          clustering = input$heatmap_clustering %||% "both",
          color_scale = input$color_scale %||% "Reds"
        )
        
        cat("[SIGNATURE HEATMAP] Heatmap generated successfully\n")
        return(heatmap_plot)
        
      }, error = function(e) {
        cat("[SIGNATURE HEATMAP] Error generating heatmap:", e$message, "\n")
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = paste("Error generating heatmap:", e$message),
                               textfont = list(size = 14, color = "red")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      })
    })
    
    # Download handler for interactive heatmap
    output$download_heatmap_html <- downloadHandler(
      filename = function() {
        paste0("signature_heatmap_", Sys.Date(), ".html")
      },
      content = function(file) {
        req(values$analysis_results$all_signatures)
        
        # Create heatmap with all UI parameters
        heatmap_plot <- create_interactive_signature_heatmap_enhanced(
          signature_data = values$analysis_results$all_signatures,
          metric = input$heatmap_metric %||% "signature_strength",
          cluster_filter = NULL,
          clustering = input$heatmap_clustering %||% "both",
          color_scale = input$color_scale %||% "Reds"
        )
        
        # Save as HTML
        htmlwidgets::saveWidget(heatmap_plot, file, selfcontained = TRUE)
      }
    )
    
    # Pan-cluster visualization - replace confusing heatmap with informative bar chart
    output$pan_cluster_heatmap <- renderPlotly({
      req(values$analysis_results)
      
      pan_cluster_data <- values$analysis_results$pan_cluster_signatures
      
      if (is.null(pan_cluster_data) || nrow(pan_cluster_data) == 0) {
        return(plotly_empty("No pan-cluster signatures available for visualization"))
      }
      
      # Create stacked bar chart showing cluster distribution
      # First, get detailed cluster information if available
      all_sigs <- values$analysis_results$all_signatures
      
      if (!is.null(all_sigs) && nrow(all_sigs) > 0) {
        # Count signatures per cluster for each gene pair
        filtered_sigs <- all_sigs %>%
          dplyr::filter(gene_pair %in% pan_cluster_data$gene_pair)
        
        # Add signature strength safely
        filtered_sigs$strength_safe <- get_signature_strength(filtered_sigs)
        
        cluster_counts <- filtered_sigs %>%
          dplyr::group_by(gene_pair, cluster) %>%
          dplyr::summarise(sig_strength = mean(strength_safe, na.rm = TRUE), .groups = "drop") %>%
          dplyr::arrange(gene_pair, cluster)
        
        # Create stacked bar chart
        plot_ly(cluster_counts, 
                x = ~gene_pair, 
                y = ~sig_strength,
                color = ~cluster,
                type = 'bar',
                text = ~paste("Cluster:", cluster, "<br>Avg Strength:", round(sig_strength, 2)),
                hovertemplate = "%{text}<extra></extra>") %>%
          layout(
            title = "Signature Strength Distribution Across Clusters",
            xaxis = list(title = "Gene Pair (MAST vs CRISPRi)", tickangle = -45),
            yaxis = list(title = "Average Signature Strength"),
            barmode = 'stack',
            height = 500,
            margin = list(b = 100)
          )
      } else {
        # Fallback: Simple bar chart of average strengths
        pan_cluster_data <- pan_cluster_data[order(pan_cluster_data$mean_signature_strength, decreasing = TRUE), ]
        
        plot_ly(pan_cluster_data,
                x = ~gene_pair,
                y = ~mean_signature_strength,
                type = 'bar',
                marker = list(color = ~cluster_count, 
                            colorscale = 'Blues',
                            showscale = TRUE,
                            colorbar = list(title = "Shared<br>Clusters")),
                text = ~paste("Shared Clusters:", cluster_count,
                            "<br>Avg Strength:", round(mean_signature_strength, 2),
                            "<br>DE Genes:", total_gene_overlaps),
                hovertemplate = "%{text}<extra></extra>") %>%
          layout(
            title = "Pan-Cluster Signatures Ranked by Average Strength",
            xaxis = list(title = "Gene Pair (MAST vs CRISPRi)", tickangle = -45),
            yaxis = list(title = "Average Signature Strength"),
            height = 500,
            margin = list(b = 100)
          )
      }
    })
    
    # Download handlers for PD analysis
    output$download_pd_analysis <- downloadHandler(
      filename = function() {
        paste0("pd_signature_analysis_", Sys.Date(), ".csv")
      },
      content = function(file) {
        req(values$analysis_results$pd_analysis)
        
        # Create detailed table from enhanced signatures
        table_data <- data.frame()
        for (i in seq_along(values$analysis_results$pd_analysis$enhanced_signatures)) {
          sig <- values$analysis_results$pd_analysis$enhanced_signatures[[i]]
          
          row_data <- data.frame(
            rank = i,
            gene_pair = sig$signature$gene_pair,
            mast_gene = sig$signature$mast_gene,
            crispri_gene = sig$signature$crispri_gene,
            cluster = safe_signature_access(sig$signature, "cluster"),
            signature_strength = safe_signature_access(sig$signature, "signature_strength"),
            pd_relevance_score = sig$pd_relevance_score,
            shared_pd_pathways = nrow(sig$shared_pd_pathways),
            mitochondrial_pathways = sig$biological_categories$mitochondrial,
            protein_quality_pathways = sig$biological_categories$protein_quality,
            autophagy_pathways = sig$biological_categories$autophagy,
            dopamine_pathways = sig$biological_categories$dopamine,
            synaptic_pathways = sig$biological_categories$synaptic,
            oxidative_stress_pathways = sig$biological_categories$oxidative_stress,
            neuronal_pathways = sig$biological_categories$neuronal,
            stringsAsFactors = FALSE
          )
          table_data <- rbind(table_data, row_data)
        }
        
        write.csv(table_data, file, row.names = FALSE)
      }
    )
    
    
    # Helper function for empty plotly plots
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
    
    # === HELPER FUNCTIONS ===
    
    # Helper function to filter genes based on selected approach
    filter_genes_for_correlation <- function(gene_data, approach) {
      if (approach == "all_genes" || nrow(gene_data) == 0) {
        return(gene_data)
      }
      
      # Filter to top N genes by absolute log2FC
      top_n <- switch(approach,
                     "top_100" = 100,
                     "top_200" = 200, 
                     "top_500" = 500,
                     "top_1000" = 1000,
                     200)  # fallback
      
      # Order by absolute log2FC and take top N
      gene_data_ordered <- gene_data[order(abs(gene_data$log2FC), decreasing = TRUE), ]
      top_genes <- head(gene_data_ordered, min(top_n, nrow(gene_data_ordered)))
      
      return(top_genes)
    }
    
    # === EFFECT SIZE CORRELATION PLOT ===
    
    # Effect size correlation plot with interactive filtering
    output$gene_pair_correlation <- renderPlotly({
      req(values$analysis_results)
      req(input$selected_gene_pair)
      
      # Check if DE data is available in analysis results
      analysis_results <- values$analysis_results
      if (!"de_data" %in% names(analysis_results) || is.null(analysis_results$de_data)) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = "DE data not available for correlation analysis",
                               textfont = list(size = 16, color = "gray")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      de_data <- analysis_results$de_data
      selected_pair <- input$selected_gene_pair
      
      # Parse gene pair (format: "MAST_GENE_vs_CRISPRI_GENE")
      genes <- strsplit(selected_pair, "_vs_")[[1]]
      if (length(genes) != 2) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = "Invalid gene pair format",
                               textfont = list(size = 16, color = "red")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      mast_gene <- genes[1]
      crispri_gene <- genes[2]
      
      cat("[CORRELATION] Analyzing", mast_gene, "vs", crispri_gene, "\n")
      
      # Get available clusters from current analysis
      all_sigs <- values$analysis_results$all_signatures
      if (is.null(all_sigs) || nrow(all_sigs) == 0) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = "No signature analysis results available",
                               textfont = list(size = 16, color = "gray")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      # Get clusters for this gene pair
      pair_clusters <- unique(all_sigs[all_sigs$gene_pair == selected_pair, "cluster"])
      if (length(pair_clusters) == 0) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = "No clusters found for this gene pair",
                               textfont = list(size = 16, color = "gray")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      cat("[CORRELATION] Found", length(pair_clusters), "clusters:", paste(pair_clusters, collapse = ", "), "\n")
      
      # Extract log2FC data across all clusters
      all_correlation_data <- list()
      
      for (cluster in pair_clusters) {
        cat("[CORRELATION] Processing cluster:", cluster, "\n")
        
        # Extract MAST data
        mast_data <- NULL
        if (mast_gene %in% names(de_data$iSCORE_PD_MAST) && 
            cluster %in% names(de_data$iSCORE_PD_MAST[[mast_gene]])) {
          mast_results <- de_data$iSCORE_PD_MAST[[mast_gene]][[cluster]]$results
          if (!is.null(mast_results) && "avg_log2FC" %in% colnames(mast_results)) {
            mast_data <- data.frame(
              gene_name = rownames(mast_results),
              log2FC = mast_results$avg_log2FC,
              stringsAsFactors = FALSE
            )
            mast_data <- mast_data[!is.na(mast_data$log2FC), ]
            
            # Apply gene filtering based on user selection
            gene_filter <- input$gene_filter_approach %||% "top_200"
            mast_data <- filter_genes_for_correlation(mast_data, gene_filter)
            
            cat("[CORRELATION] MAST data after filtering (", gene_filter, "): ", nrow(mast_data), "genes\n")
          }
        }
        
        # Extract CRISPRi data
        if (!is.null(mast_data) && nrow(mast_data) > 0 && 
            crispri_gene %in% names(de_data$CRISPRi_Mixscale) && 
            cluster %in% names(de_data$CRISPRi_Mixscale[[crispri_gene]])) {
          
          crispri_results <- de_data$CRISPRi_Mixscale[[crispri_gene]][[cluster]]$results
          if (!is.null(crispri_results)) {
            
            # Get available experiments
            available_experiments <- c("C12_FPD-23", "C12_FPD-24", "C18_FPD-23")
            
            if (input$show_all_experiments) {
              # Multi-experiment mode: include all experiments
              for (exp in available_experiments) {
                log2fc_col <- paste0("log2FC_", exp)
                if (log2fc_col %in% colnames(crispri_results)) {
                  crispri_data <- data.frame(
                    gene_name = crispri_results$gene_ID,
                    log2FC = crispri_results[[log2fc_col]],
                    experiment = exp,
                    stringsAsFactors = FALSE
                  )
                  crispri_data <- crispri_data[!is.na(crispri_data$log2FC), ]
                  
                  # Apply same gene filtering as MAST data
                  crispri_data <- filter_genes_for_correlation(crispri_data, gene_filter)
                  
                  if (nrow(crispri_data) > 0) {
                    # Calculate correlation for this cluster and experiment
                    cor_result <- calculate_effect_size_correlation(mast_data, crispri_data[, c("gene_name", "log2FC")])
                    if (!is.null(cor_result) && !is.na(cor_result$correlation)) {
                      all_correlation_data[[paste(cluster, exp, sep = "_")]] <- list(
                        cluster = cluster,
                        experiment = exp,
                        correlation = cor_result,
                        mast_data = mast_data,
                        crispri_data = crispri_data
                      )
                    }
                  }
                }
              }
            } else {
              # Single experiment mode
              selected_exp <- input$crispri_experiment %||% "C12_FPD-23"
              log2fc_col <- paste0("log2FC_", selected_exp)
              
              if (log2fc_col %in% colnames(crispri_results)) {
                crispri_data <- data.frame(
                  gene_name = crispri_results$gene_ID,
                  log2FC = crispri_results[[log2fc_col]],
                  stringsAsFactors = FALSE
                )
                crispri_data <- crispri_data[!is.na(crispri_data$log2FC), ]
                
                # Apply same gene filtering as MAST data
                crispri_data <- filter_genes_for_correlation(crispri_data, gene_filter)
                
                if (nrow(crispri_data) > 0) {
                  # Calculate correlation for this cluster
                  cor_result <- calculate_effect_size_correlation(mast_data, crispri_data)
                  if (!is.null(cor_result) && !is.na(cor_result$correlation)) {
                    all_correlation_data[[cluster]] <- list(
                      cluster = cluster,
                      experiment = selected_exp,
                      correlation = cor_result,
                      mast_data = mast_data,
                      crispri_data = crispri_data
                    )
                  }
                }
              }
            }
          }
        }
      }
      
      if (length(all_correlation_data) == 0) {
        return(plotly::plot_ly() %>% 
               plotly::add_text(x = 0.5, y = 0.5, 
                               text = "No overlapping genes found between MAST and CRISPRi",
                               textfont = list(size = 16, color = "orange")) %>%
               plotly::layout(
                 xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                 yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
               ))
      }
      
      cat("[CORRELATION] Generated", length(all_correlation_data), "correlation datasets\n")
      
      # Helper function to identify gene of interest for highlighting
      identify_gene_of_interest <- function(gene_pair) {
        # Extract the actual gene names for highlighting
        if (grepl("PRKN_vs_PARK2", gene_pair)) {
          return(c("PRKN", "PARK2"))  # Both genes are relevant
        } else if (grepl("VPS13C.*_vs_VPS13C", gene_pair)) {
          return("VPS13C")
        } else if (grepl("SNCA.*_vs_SNCA", gene_pair)) {
          return("SNCA")
        } else {
          # For cases like PARK7_vs_PARK7, DNAJC6_vs_DNAJC6, etc.
          gene_name <- strsplit(gene_pair, "_vs_")[[1]][1]
          return(gene_name)
        }
      }
      
      # Get gene of interest for current pair
      selected_pair <- input$selected_gene_pair
      mast_gene <- strsplit(selected_pair, "_vs_")[[1]][1]
      crispri_gene <- strsplit(selected_pair, "_vs_")[[1]][2]
      genes_of_interest <- identify_gene_of_interest(selected_pair)
      
      # Create plotly visualization
      if (input$cluster_grid_view %||% TRUE) {
        # CLUSTER VERTICAL VIEW: Create vertically stacked plots for each cluster
        unique_clusters <- unique(sapply(all_correlation_data, function(x) x$cluster))
        unique_clusters <- sort(unique_clusters)
        
        # Vertical stacking - no grid calculations needed
        n_clusters <- length(unique_clusters)
        
        # Create individual plots for each cluster
        cluster_plots <- list()
        experiment_colors <- c("C12_FPD-23" = "#1f77b4", "C12_FPD-24" = "#ff7f0e", "C18_FPD-23" = "#2ca02c")
        
        for (i in seq_along(unique_clusters)) {
          cluster <- unique_clusters[i]
          
          # Get data for this cluster
          cluster_data <- data.frame()
          for (data_key in names(all_correlation_data)) {
            corr_data <- all_correlation_data[[data_key]]
            if (corr_data$cluster == cluster) {
              merged_data <- corr_data$correlation$merged_data
              merged_data$cluster <- corr_data$cluster
              merged_data$experiment <- corr_data$experiment
              cluster_data <- rbind(cluster_data, merged_data)
            }
          }
          
          if (nrow(cluster_data) > 0) {
            # Clean data
            cluster_data_clean <- cluster_data[complete.cases(cluster_data[c("log2FC_mast", "log2FC_crispri")]), ]
            
            # Calculate overall correlation for this cluster (across all experiments)
            cluster_cor_result <- NULL
            if (nrow(cluster_data_clean) >= 3) {
              tryCatch({
                cor_test <- cor.test(cluster_data_clean$log2FC_mast, cluster_data_clean$log2FC_crispri, method = "pearson")
                cluster_cor_result <- list(
                  correlation = cor_test$estimate,
                  p_value = cor_test$p.value,
                  n_genes = nrow(cluster_data_clean)
                )
              }, error = function(e) {
                cluster_cor_result <- NULL
              })
            }
            
            # Create plot title with correlation statistics
            cluster_title <- paste("Cluster", cluster)
            if (!is.null(cluster_cor_result)) {
              r_val <- round(cluster_cor_result$correlation, 3)
              p_val <- cluster_cor_result$p_value
              n_genes <- cluster_cor_result$n_genes
              
              # Format p-value
              if (p_val < 0.001) {
                p_text <- "p < 0.001"
              } else {
                p_text <- paste("p =", round(p_val, 3))
              }
              
              cluster_title <- paste0(cluster_title, "\n", "r = ", r_val, ", ", p_text, " (n = ", n_genes, ")")
            }
            
            # Add gene highlighting
            cluster_data_clean$is_gene_of_interest <- cluster_data_clean$gene_name %in% genes_of_interest
            cluster_data_clean$point_color <- ifelse(cluster_data_clean$is_gene_of_interest, "red", 
                                                     experiment_colors[cluster_data_clean$experiment])
            cluster_data_clean$point_symbol <- ifelse(cluster_data_clean$is_gene_of_interest, "star", "circle")
            
            # Create hover text
            cluster_data_clean$hover_text <- paste("Gene:", cluster_data_clean$gene_name, 
                                                  "<br>Cluster:", cluster_data_clean$cluster,
                                                  "<br>Experiment:", cluster_data_clean$experiment,
                                                  "<br>MAST log2FC:", round(cluster_data_clean$log2FC_mast, 3),
                                                  "<br>CRISPRi log2FC:", round(cluster_data_clean$log2FC_crispri, 3),
                                                  ifelse(cluster_data_clean$is_gene_of_interest, "<br><b>Gene of Interest</b>", ""))
            
            # Create ggplot for this cluster
            p_cluster <- ggplot2::ggplot(cluster_data_clean, ggplot2::aes(x = log2FC_mast, y = log2FC_crispri)) +
              # Add bold reference lines at x=0 and y=0
              ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
              ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
              ggplot2::geom_point(ggplot2::aes(text = hover_text, color = experiment, 
                                              shape = ifelse(is_gene_of_interest, "Gene of Interest", "Other")), 
                                 size = ifelse(cluster_data_clean$is_gene_of_interest, 3, 2), 
                                 alpha = ifelse(cluster_data_clean$is_gene_of_interest, 1, 0.7)) +
              ggplot2::scale_color_manual(values = experiment_colors) +
              ggplot2::scale_shape_manual(values = c("Gene of Interest" = 17, "Other" = 16)) +
              ggplot2::labs(title = cluster_title,
                           x = paste("MAST log2FC (", mast_gene, ")"),
                           y = paste("CRISPRi log2FC (", crispri_gene, ")")) +
              ggplot2::theme_minimal() +
              ggplot2::theme(legend.position = "none")  # Remove individual legends
            
            # Add trend line if enough points
            if (nrow(cluster_data_clean) >= 3) {
              p_cluster <- p_cluster + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.5)
            }
            
            cluster_plots[[i]] <- p_cluster
          } else {
            # Create empty plot for clusters with no data
            cluster_plots[[i]] <- ggplot2::ggplot() + 
              ggplot2::annotate("text", x = 0, y = 0, label = "No data", size = 4, color = "gray") +
              ggplot2::labs(title = paste("Cluster", cluster)) +
              ggplot2::theme_void()
          }
        }
        
        # Stack plots vertically for better readability
        if (length(cluster_plots) > 1) {
          # Convert plots to plotly, handling empty plots carefully
          plotly_plots <- list()
          for (i in seq_along(cluster_plots)) {
            tryCatch({
              plotly_plots[[i]] <- plotly::ggplotly(cluster_plots[[i]], tooltip = "text")
            }, error = function(e) {
              # Create a simple plotly object for empty clusters
              plotly_plots[[i]] <- plotly::plot_ly() %>%
                plotly::add_annotations(
                  text = paste("Cluster", unique_clusters[i], "- No data"),
                  x = 0.5, y = 0.5,
                  xref = "paper", yref = "paper",
                  showarrow = FALSE,
                  font = list(size = 14, color = "gray")
                ) %>%
                plotly::layout(
                  xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
                  yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
                )
            })
          }
          
          # Remove any NULL plots
          plotly_plots <- plotly_plots[!sapply(plotly_plots, is.null)]
          
          if (length(plotly_plots) > 0) {
            p <- plotly::subplot(
              plotly_plots,
              nrows = length(plotly_plots), 
              ncols = 1,
              shareX = FALSE,  # Don't share X axis so each plot can have its own scale
              shareY = FALSE,  # Don't share Y axis so each plot can have its own scale
              titleX = TRUE,   # Show X axis titles for each plot
              titleY = TRUE,   # Show Y axis titles for each plot
              margin = 0.08    # Add some margin between plots
            ) %>%
            plotly::layout(
              title = list(text = paste0("Cluster-Specific Correlations: ", mast_gene, " vs ", crispri_gene, 
                                        " (", ifelse(input$show_all_experiments, "All Experiments", "Single Experiment"), 
                                        ") - Scroll to view all clusters"),
                          font = list(size = 16)),
              showlegend = TRUE,
              height = 400 * length(plotly_plots)  # Dynamic height based on number of clusters
            )
          } else {
            p <- plotly::plot_ly() %>%
              plotly::add_annotations(
                text = "No correlation data available",
                x = 0.5, y = 0.5,
                xref = "paper", yref = "paper",
                showarrow = FALSE
              )
          }
        } else if (length(cluster_plots) == 1) {
          tryCatch({
            p <- plotly::ggplotly(cluster_plots[[1]], tooltip = "text")
          }, error = function(e) {
            p <- plotly::plot_ly() %>%
              plotly::add_annotations(
                text = "No correlation data available",
                x = 0.5, y = 0.5,
                xref = "paper", yref = "paper",
                showarrow = FALSE
              )
          })
        } else {
          p <- plotly::plot_ly() %>%
            plotly::add_annotations(
              text = "No correlation data available",
              x = 0.5, y = 0.5,
              xref = "paper", yref = "paper",
              showarrow = FALSE
            )
        }
        
        return(p)
        
      } else if (input$show_all_experiments) {
        # PAN-CLUSTER VIEW: Multi-experiment plot with different colors
        plot_data <- data.frame()
        experiment_colors <- c("C12_FPD-23" = "#1f77b4", "C12_FPD-24" = "#ff7f0e", "C18_FPD-23" = "#2ca02c")
        
        for (data_key in names(all_correlation_data)) {
          corr_data <- all_correlation_data[[data_key]]
          merged_data <- corr_data$correlation$merged_data
          cat("[CORRELATION DEBUG] Adding data for", data_key, "- rows:", nrow(merged_data), "\n")
          merged_data$cluster <- corr_data$cluster
          merged_data$experiment <- corr_data$experiment
          plot_data <- rbind(plot_data, merged_data)
        }
        
        cat("[CORRELATION DEBUG] Total plot_data before cleaning - rows:", nrow(plot_data), "\n")
        
        if (nrow(plot_data) > 0) {
          # Clean data for plotly to prevent tibble column size mismatches
          plot_data_clean <- plot_data[complete.cases(plot_data[c("log2FC_mast", "log2FC_crispri")]), ]
          cat("[CORRELATION DEBUG] Total plot_data_clean after filtering - rows:", nrow(plot_data_clean), "\n")
          
          # Calculate overall correlation for pan-cluster view
          pan_cluster_cor_result <- NULL
          if (nrow(plot_data_clean) >= 3) {
            tryCatch({
              cor_test <- cor.test(plot_data_clean$log2FC_mast, plot_data_clean$log2FC_crispri, method = "pearson")
              pan_cluster_cor_result <- list(
                correlation = cor_test$estimate,
                p_value = cor_test$p.value,
                n_genes = nrow(plot_data_clean)
              )
            }, error = function(e) {
              pan_cluster_cor_result <- NULL
            })
          }
          
          # Add gene highlighting for pan-cluster view
          plot_data_clean$is_gene_of_interest <- plot_data_clean$gene_name %in% genes_of_interest
          
          # Create hover text AFTER filtering to ensure size matches
          plot_data_clean$hover_text <- paste("Gene:", plot_data_clean$gene_name, 
                                             "<br>Cluster:", plot_data_clean$cluster,
                                             "<br>Experiment:", plot_data_clean$experiment,
                                             "<br>MAST log2FC:", round(plot_data_clean$log2FC_mast, 3),
                                             "<br>CRISPRi log2FC:", round(plot_data_clean$log2FC_crispri, 3),
                                             ifelse(plot_data_clean$is_gene_of_interest, "<br><b>Gene of Interest</b>", ""))
          
          cat("[MULTI EXP DEBUG] Final plot_data_clean before plotly - rows:", nrow(plot_data_clean), 
              "cols:", ncol(plot_data_clean), "\n")
          cat("[MULTI EXP DEBUG] hover_text length:", length(plot_data_clean$hover_text), "\n")
          cat("[MULTI EXP DEBUG] x data length:", length(plot_data_clean$log2FC_mast), "\n")
          cat("[MULTI EXP DEBUG] y data length:", length(plot_data_clean$log2FC_crispri), "\n")
          
          # Sample data if too large to prevent plotly internal filtering issues
          if (nrow(plot_data_clean) > 5000) {
            cat("[MULTI EXP DEBUG] Sampling large dataset from", nrow(plot_data_clean), "to 5000 points\n")
            sample_indices <- sample(nrow(plot_data_clean), 5000)
            plot_data_sampled <- plot_data_clean[sample_indices, ]
          } else {
            plot_data_sampled <- plot_data_clean
          }
          
          cat("[MULTI EXP DEBUG] Final sampled data - rows:", nrow(plot_data_sampled), "\n")
          
          # Create interactive visualization using ggplot2 + ggplotly for gene highlighting
          gene_filter_text <- switch(input$gene_filter_approach %||% "top_200",
                                   "top_100" = "top 100",
                                   "top_200" = "top 200", 
                                   "top_500" = "top 500",
                                   "top_1000" = "top 1000",
                                   "all_genes" = "all",
                                   "top 200")
          
          # Create plot title with correlation statistics
          plot_title <- paste0("Pan-Cluster Correlation: ", mast_gene, " vs ", crispri_gene, " (All Experiments)")
          if (!is.null(pan_cluster_cor_result)) {
            r_val <- round(pan_cluster_cor_result$correlation, 3)
            p_val <- pan_cluster_cor_result$p_value
            n_genes <- pan_cluster_cor_result$n_genes
            
            # Format p-value
            if (p_val < 0.001) {
              p_text <- "p < 0.001"
            } else {
              p_text <- paste("p =", round(p_val, 3))
            }
            
            plot_title <- paste0(plot_title, "\nr = ", r_val, ", ", p_text, " (n = ", n_genes, ") | ",
                               "Based on ", gene_filter_text, " most changed genes per method")
          } else {
            plot_title <- paste0(plot_title, "\nBased on ", gene_filter_text, " most changed genes per method")
          }
          
          # Create ggplot2 scatter plot with gene highlighting
          p <- ggplot2::ggplot(plot_data_sampled, ggplot2::aes(x = log2FC_mast, y = log2FC_crispri)) +
            # Add bold reference lines at x=0 and y=0
            ggplot2::geom_hline(yintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
            ggplot2::geom_vline(xintercept = 0, color = "black", linewidth = 0.8, alpha = 0.7) +
            ggplot2::geom_point(ggplot2::aes(text = hover_text, color = experiment, 
                                            shape = ifelse(is_gene_of_interest, "Gene of Interest", "Other")), 
                               size = ifelse(plot_data_sampled$is_gene_of_interest, 3, 2), 
                               alpha = ifelse(plot_data_sampled$is_gene_of_interest, 1, 0.7)) +
            ggplot2::scale_color_manual(values = c("C12_FPD-23" = "#1f77b4", "C12_FPD-24" = "#ff7f0e", "C18_FPD-23" = "#2ca02c")) +
            ggplot2::scale_shape_manual(values = c("Gene of Interest" = 17, "Other" = 16)) +
            ggplot2::labs(
              title = plot_title,
              x = paste("MAST log2FC (", mast_gene, "mutation)"),
              y = paste("CRISPRi log2FC (", crispri_gene, "knockdown)"),
              color = "Experiment",
              shape = "Gene Type"
            ) +
            ggplot2::theme_minimal()
          
          # Add trend line if enough points
          if (nrow(plot_data_sampled) >= 3) {
            p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "gray", linetype = "dashed", linewidth = 0.8)
          }
          
          # Convert to interactive plotly
          p <- plotly::ggplotly(p, tooltip = "text")
          
          return(p)
        }
      } else {
        # Single experiment plot
        # Combine data from all clusters for overall correlation
        combined_plot_data <- data.frame()
        
        for (data_key in names(all_correlation_data)) {
          corr_data <- all_correlation_data[[data_key]]
          merged_data <- corr_data$correlation$merged_data
          cat("[SINGLE EXP DEBUG] Adding data for", data_key, "- rows:", nrow(merged_data), "\n")
          merged_data$cluster <- corr_data$cluster
          combined_plot_data <- rbind(combined_plot_data, merged_data)
        }
        
        cat("[SINGLE EXP DEBUG] Total combined_plot_data before cleaning - rows:", nrow(combined_plot_data), "\n")
        
        if (nrow(combined_plot_data) > 0) {
          # Calculate overall correlation with error handling
          overall_cor <- tryCatch({
            calculate_effect_size_correlation(
              data.frame(gene_name = combined_plot_data$gene_name, log2FC = combined_plot_data$log2FC_mast),
              data.frame(gene_name = combined_plot_data$gene_name, log2FC = combined_plot_data$log2FC_crispri)
            )
          }, error = function(e) {
            list(correlation = NA, p_value = NA, n_genes = 0, error = e$message)
          })
          
          selected_exp <- input$crispri_experiment %||% "C12_FPD-23"
          
          # Safe correlation display - extract only scalar values to avoid serialization issues
          cor_text <- if (!is.null(overall_cor) && !is.na(overall_cor$correlation)) {
            # Extract scalar values to prevent "[object Object]" errors in plotly
            cor_val <- as.numeric(overall_cor$correlation)
            p_val <- as.numeric(overall_cor$p_value) 
            n_genes <- as.numeric(overall_cor$n_genes)
            
            paste0("r = ", round(cor_val, 3), 
                   ", p = ", format.pval(p_val, digits = 3),
                   ", n = ", n_genes, " genes")
          } else {
            "Correlation calculation failed"
          }
          
          # Clean data for plotly to prevent tibble column size mismatches
          plot_data_clean <- combined_plot_data[, c("gene_name", "log2FC_mast", "log2FC_crispri", "cluster")]
          plot_data_clean <- plot_data_clean[complete.cases(plot_data_clean), ]
          
          # Convert to data frame to ensure proper structure
          plot_data_clean <- as.data.frame(plot_data_clean)
          
          # Create hover text AFTER filtering to ensure size matches
          plot_data_clean$hover_text <- paste("Gene:", plot_data_clean$gene_name, 
                                             "<br>Cluster:", plot_data_clean$cluster,
                                             "<br>MAST log2FC:", round(plot_data_clean$log2FC_mast, 3),
                                             "<br>CRISPRi log2FC:", round(plot_data_clean$log2FC_crispri, 3))
          
          cat("[SINGLE EXP DEBUG] Final plot_data_clean before plotly - rows:", nrow(plot_data_clean), 
              "cols:", ncol(plot_data_clean), "\n")
          cat("[SINGLE EXP DEBUG] Column names:", paste(colnames(plot_data_clean), collapse = ", "), "\n")
          cat("[SINGLE EXP DEBUG] hover_text length:", length(plot_data_clean$hover_text), "\n")
          cat("[SINGLE EXP DEBUG] x data length:", length(plot_data_clean$log2FC_mast), "\n")
          cat("[SINGLE EXP DEBUG] y data length:", length(plot_data_clean$log2FC_crispri), "\n")
          
          # Sample data if too large to prevent plotly internal filtering issues
          if (nrow(plot_data_clean) > 5000) {
            cat("[SINGLE EXP DEBUG] Sampling large dataset from", nrow(plot_data_clean), "to 5000 points\n")
            sample_indices <- sample(nrow(plot_data_clean), 5000)
            plot_data_sampled <- plot_data_clean[sample_indices, ]
          } else {
            plot_data_sampled <- plot_data_clean
          }
          
          # Ensure sampled data is also a proper data frame
          plot_data_sampled <- as.data.frame(plot_data_sampled)
          
          cat("[SINGLE EXP DEBUG] Final sampled data - rows:", nrow(plot_data_sampled), "\n")
          cat("[SINGLE EXP DEBUG] Sampled data columns:", paste(colnames(plot_data_sampled), collapse = ", "), "\n")
          cat("[SINGLE EXP DEBUG] Sampled hover_text exists:", "hover_text" %in% colnames(plot_data_sampled), "\n")
          
          # Create gene filter text for subtitle
          gene_filter_text <- switch(input$gene_filter_approach %||% "top_200",
                                   "top_100" = "top 100",
                                   "top_200" = "top 200", 
                                   "top_500" = "top 500",
                                   "top_1000" = "top 1000",
                                   "all_genes" = "all",
                                   "top 200")
          
          # Create ggplot2 scatter plot (more stable than pure plotly)
          p <- ggplot2::ggplot(plot_data_sampled, ggplot2::aes(x = log2FC_mast, y = log2FC_crispri)) +
            ggplot2::geom_point(ggplot2::aes(text = hover_text), size = 2, alpha = 0.7, color = "#1f77b4") +
            ggplot2::labs(
              title = paste0("Effect Size Correlation: ", mast_gene, " vs ", crispri_gene, " (", selected_exp, ")\n",
                           "Based on ", gene_filter_text, " most changed genes per method | ", cor_text),
              x = paste("MAST log2FC (", mast_gene, "mutation)"),
              y = paste("CRISPRi log2FC (", crispri_gene, "knockdown)")
            ) +
            ggplot2::theme_minimal()
          
          # Add trend line if we have enough data points
          if (nrow(plot_data_sampled) >= 3) {
            p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.8)
          }
          
          # Convert to interactive plotly
          p <- plotly::ggplotly(p, tooltip = "text")
          
          # STATIC PLOT FALLBACK (uncomment if ggplotly still fails):
          # p <- ggplot2::ggplot(plot_data_sampled, ggplot2::aes(x = log2FC_mast, y = log2FC_crispri)) +
          #   ggplot2::geom_point(size = 2, alpha = 0.7, color = "#1f77b4") +
          #   ggplot2::labs(
          #     title = paste0("Effect Size Correlation: ", mast_gene, " vs ", crispri_gene, " (", selected_exp, ")\n",
          #                  "Based on ", gene_filter_text, " most changed genes per method | ", cor_text),
          #     x = paste("MAST log2FC (", mast_gene, "mutation)"),
          #     y = paste("CRISPRi log2FC (", crispri_gene, "knockdown)")
          #   ) +
          #   ggplot2::theme_minimal()
          # if (nrow(plot_data_sampled) >= 3) {
          #   p <- p + ggplot2::geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed", linewidth = 0.8)
          # }
          
          return(p)
        }
      }
      
      # Fallback empty plot
      return(plotly::plot_ly() %>% 
             plotly::add_text(x = 0.5, y = 0.5, 
                             text = "Unable to generate correlation plot",
                             textfont = list(size = 16, color = "red")) %>%
             plotly::layout(
               xaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE),
               yaxis = list(showgrid = FALSE, showticklabels = FALSE, zeroline = FALSE)
             ))
    })
    
    # === MODAL HANDLER FOR GENE DETAILS ===
    observeEvent(input$show_details, {
      req(input$show_details)
      req(values$analysis_results)
      req(input$selected_gene_pair)
      
      # Extract cluster index from the button click
      cluster_idx <- input$show_details$cluster
      
      # Get the cluster name from the current table data
      all_sigs <- values$analysis_results$all_signatures
      selected_pair <- input$selected_gene_pair
      pair_data <- all_sigs[all_sigs$gene_pair == selected_pair, ]
      
      if (nrow(pair_data) >= cluster_idx) {
        selected_cluster <- pair_data$cluster[cluster_idx]
        
        # Extract shared genes and pathways
        shared_genes_data <- extract_shared_genes(values$analysis_results, selected_pair, selected_cluster)
        shared_pathways_data <- extract_shared_pathways(values$analysis_results, selected_pair, selected_cluster)
        
        # Show modal with details
        showModal(modalDialog(
          title = paste("Details for", selected_pair, "- Cluster", selected_cluster),
          size = "l",
          easyClose = TRUE,
          footer = modalButton("Close"),
          
          tabsetPanel(
            tabPanel("Shared DE Genes",
              br(),
              if (nrow(shared_genes_data) > 0) {
                DT::dataTableOutput(ns("modal_genes_table"))
              } else {
                p("No shared DE genes found for this cluster.")
              }
            ),
            tabPanel("Enriched Pathways",
              br(),
              if (nrow(shared_pathways_data) > 0) {
                DT::dataTableOutput(ns("modal_pathways_table"))
              } else {
                p("No enriched pathways found for this cluster.")
              }
            )
          )
        ))
      }
    })
    
    # Modal output for genes table
    output$modal_genes_table <- DT::renderDataTable({
      req(input$show_details)
      cluster_idx <- input$show_details$cluster
      all_sigs <- values$analysis_results$all_signatures
      selected_pair <- input$selected_gene_pair
      pair_data <- all_sigs[all_sigs$gene_pair == selected_pair, ]
      
      if (nrow(pair_data) >= cluster_idx) {
        selected_cluster <- pair_data$cluster[cluster_idx]
        shared_genes_data <- extract_shared_genes(values$analysis_results, selected_pair, selected_cluster)
        
        DT::datatable(shared_genes_data,
                     options = list(pageLength = 10, scrollX = TRUE),
                     rownames = FALSE)
      }
    })
    
    # Modal output for pathways table
    output$modal_pathways_table <- DT::renderDataTable({
      req(input$show_details)
      cluster_idx <- input$show_details$cluster
      all_sigs <- values$analysis_results$all_signatures
      selected_pair <- input$selected_gene_pair
      pair_data <- all_sigs[all_sigs$gene_pair == selected_pair, ]
      
      if (nrow(pair_data) >= cluster_idx) {
        selected_cluster <- pair_data$cluster[cluster_idx]
        shared_pathways_data <- extract_shared_pathways(values$analysis_results, selected_pair, selected_cluster)
        
        DT::datatable(shared_pathways_data,
                     options = list(pageLength = 10, scrollX = TRUE),
                     rownames = FALSE)
      }
    })
    
    # === HELPER FUNCTIONS FOR MODAL DATA EXTRACTION ===
    
    # Extract shared genes from analysis results
    extract_shared_genes <- function(analysis_results, selected_pair, selected_cluster) {
      tryCatch({
        cat("[DETAILS DEBUG] === EXTRACTING SHARED GENES FROM ENRICHMENT TERMS ===\n")
        cat("[DETAILS DEBUG] Selected pair:", selected_pair, "\n")
        cat("[DETAILS DEBUG] Selected cluster:", selected_cluster, "\n")
        
        # Parse gene names from the pair
        genes <- strsplit(selected_pair, "_vs_")[[1]]
        mast_gene <- genes[1]
        crispri_gene <- genes[2]
        cat("[DETAILS DEBUG] Parsed genes - MAST:", mast_gene, "CRISPRi:", crispri_gene, "\n")
        
        # Gene name harmonization (PRKN = PARK2)
        # Apply the same harmonization used in the Fisher's exact test
        if (mast_gene == "PRKN") mast_gene <- "PARK2"
        if (crispri_gene == "PRKN") crispri_gene <- "PARK2"
        cat("[DETAILS DEBUG] After harmonization - MAST:", mast_gene, "CRISPRi:", crispri_gene, "\n")
        
        # Get the enrichment data using the data manager
        enrichment_data <- NULL
        
        # First try the app_data if available
        if (!is.null(app_data) && !is.null(app_data$consolidated_data)) {
          enrichment_data <- app_data$consolidated_data
          cat("[DETAILS DEBUG] Got enrichment data from app_data$consolidated_data\n")
        } else {
          # Use the data manager directly
          cat("[DETAILS DEBUG] app_data not available, using data manager\n")
          if (exists("get_enrichment_data", mode = "function")) {
            enrichment_data <- get_enrichment_data()
            cat("[DETAILS DEBUG] Got enrichment data from data manager\n")
          } else {
            cat("[DETAILS DEBUG] Data manager function not available\n")
            return(data.frame(Message = "Enrichment data not loaded. Please ensure the enrichment data file is loaded."))
          }
        }
        
        if (is.null(enrichment_data)) {
          cat("[DETAILS DEBUG] Enrichment data is NULL\n")
          return(data.frame(Message = "Enrichment data not available"))
        }
        
        cat("[DETAILS DEBUG] Enrichment data available with", nrow(enrichment_data), "rows\n")
        if (nrow(enrichment_data) > 0) {
          cat("[DETAILS DEBUG] First few mutation values:", paste(head(enrichment_data$mutation_perturbation, 10), collapse=", "), "\n")
          cat("[DETAILS DEBUG] First few cluster values:", paste(head(enrichment_data$cluster, 10), collapse=", "), "\n")
          cat("[DETAILS DEBUG] Looking for MAST:", mast_gene, "cluster:", selected_cluster, "\n")
          cat("[DETAILS DEBUG] Looking for CRISPRi:", crispri_gene, "cluster:", selected_cluster, "\n")
          
          # Check what mutations are available for PRKN/PARK2
          prkn_variants <- unique(enrichment_data$mutation_perturbation[grep("PRKN|PARK", enrichment_data$mutation_perturbation)])
          cat("[DETAILS DEBUG] Available PRKN/PARK variants in data:", paste(prkn_variants, collapse=", "), "\n")
          
          # Check what's available for this cluster
          cluster_mutations <- unique(enrichment_data$mutation_perturbation[enrichment_data$cluster == selected_cluster])
          cat("[DETAILS DEBUG] Mutations available in", selected_cluster, ":", paste(head(cluster_mutations, 10), collapse=", "), "\n")
          
          # Check methods available
          cat("[DETAILS DEBUG] Available methods:", paste(unique(enrichment_data$method), collapse=", "), "\n")
        }
        
        # Get enrichment terms for MAST
        mast_terms <- enrichment_data[
          enrichment_data$mutation_perturbation == mast_gene & 
          enrichment_data$cluster == selected_cluster &
          enrichment_data$method == "MAST", 
        ]
        
        # Get enrichment terms for CRISPRi
        crispri_terms <- enrichment_data[
          enrichment_data$mutation_perturbation == crispri_gene & 
          enrichment_data$cluster == selected_cluster &
          enrichment_data$method == "MixScale", 
        ]
        
        cat("[DETAILS DEBUG] Found", nrow(mast_terms), "MAST enrichment terms\n")
        cat("[DETAILS DEBUG] Found", nrow(crispri_terms), "CRISPRi enrichment terms\n")
        
        if (nrow(mast_terms) == 0 || nrow(crispri_terms) == 0) {
          return(data.frame(Message = "No enrichment data found for one or both methods"))
        }
        
        # Find shared enrichment terms (same Description)
        shared_terms <- intersect(mast_terms$Description, crispri_terms$Description)
        cat("[DETAILS DEBUG] Found", length(shared_terms), "shared enrichment terms\n")
        
        if (length(shared_terms) == 0) {
          return(data.frame(Message = "No shared enrichment terms found between methods"))
        }
        
        # Get genes associated with these shared terms
        shared_mast_terms <- mast_terms[mast_terms$Description %in% shared_terms, ]
        shared_crispri_terms <- crispri_terms[crispri_terms$Description %in% shared_terms, ]
        
        # Extract genes from geneID column for shared terms
        mast_genes <- unique(unlist(strsplit(as.character(shared_mast_terms$geneID), "/")))
        mast_genes <- mast_genes[mast_genes != "" & !is.na(mast_genes)]
        
        crispri_genes <- unique(unlist(strsplit(as.character(shared_crispri_terms$geneID), "/")))
        crispri_genes <- crispri_genes[crispri_genes != "" & !is.na(crispri_genes)]
        
        # Get genes that appear in shared terms from BOTH methods
        shared_genes <- intersect(mast_genes, crispri_genes)
        
        cat("[DETAILS DEBUG] Genes from MAST shared terms:", length(mast_genes), "\n")
        cat("[DETAILS DEBUG] Genes from CRISPRi shared terms:", length(crispri_genes), "\n")
        cat("[DETAILS DEBUG] Overlapping genes:", length(shared_genes), "\n")
        
        if (length(shared_genes) == 0) {
          return(data.frame(
            Message = paste0("Found ", length(shared_terms), " shared enrichment terms, ",
                           "but no overlapping genes between the methods")
          ))
        }
        
        # Return the shared genes
        return(data.frame(
          Gene = sort(shared_genes),
          stringsAsFactors = FALSE
        ))
      }, error = function(e) {
        cat("[DETAILS DEBUG] ERROR in extract_shared_genes:", e$message, "\n")
        return(data.frame(Error = paste("Failed to extract genes:", e$message)))
      })
    }
    
    # Extract shared pathways from analysis results
    extract_shared_pathways <- function(analysis_results, selected_pair, selected_cluster) {
      tryCatch({
        cat("[PATHWAYS DEBUG] === EXTRACTING SHARED PATHWAYS ===\n")
        cat("[PATHWAYS DEBUG] Selected pair:", selected_pair, "\n")
        cat("[PATHWAYS DEBUG] Selected cluster:", selected_cluster, "\n")
        
        # Parse gene names from the pair
        genes <- strsplit(selected_pair, "_vs_")[[1]]
        mast_gene <- genes[1]
        crispri_gene <- genes[2]
        
        # Get the enrichment data using the data manager
        enrichment_data <- NULL
        
        # First try the app_data if available
        if (!is.null(app_data) && !is.null(app_data$consolidated_data)) {
          enrichment_data <- app_data$consolidated_data
          cat("[PATHWAYS DEBUG] Got enrichment data from app_data$consolidated_data\n")
        } else {
          # Use the data manager directly
          cat("[PATHWAYS DEBUG] app_data not available, using data manager\n")
          if (exists("get_enrichment_data", mode = "function")) {
            enrichment_data <- get_enrichment_data()
            cat("[PATHWAYS DEBUG] Got enrichment data from data manager\n")
          } else {
            cat("[PATHWAYS DEBUG] Data manager function not available\n")
            return(data.frame(Message = "Enrichment data not loaded."))
          }
        }
        
        if (is.null(enrichment_data)) {
          cat("[PATHWAYS DEBUG] Enrichment data is NULL\n")
          return(data.frame(Message = "Enrichment data not available"))
        }
        
        # Get enrichment terms for MAST
        mast_terms <- enrichment_data[
          enrichment_data$mutation_perturbation == mast_gene & 
          enrichment_data$cluster == selected_cluster &
          enrichment_data$method == "MAST", 
        ]
        
        # Get enrichment terms for CRISPRi
        crispri_terms <- enrichment_data[
          enrichment_data$mutation_perturbation == crispri_gene & 
          enrichment_data$cluster == selected_cluster &
          enrichment_data$method == "MixScale", 
        ]
        
        cat("[PATHWAYS DEBUG] Found", nrow(mast_terms), "MAST terms,", nrow(crispri_terms), "CRISPRi terms\n")
        
        if (nrow(mast_terms) == 0 || nrow(crispri_terms) == 0) {
          return(data.frame(Message = "No enrichment data found for one or both methods"))
        }
        
        # Find shared enrichment terms (same Description)
        shared_descriptions <- intersect(mast_terms$Description, crispri_terms$Description)
        cat("[PATHWAYS DEBUG] Found", length(shared_descriptions), "shared pathway descriptions\n")
        
        if (length(shared_descriptions) == 0) {
          return(data.frame(Message = "No shared pathways found between methods"))
        }
        
        # Get the shared terms with their details
        shared_pathways <- enrichment_data[
          enrichment_data$Description %in% shared_descriptions &
          enrichment_data$mutation_perturbation %in% c(mast_gene, crispri_gene) &
          enrichment_data$cluster == selected_cluster,
        ]
        
        # Create display dataframe with unique pathways
        unique_pathways <- unique(shared_pathways[, c("ID", "Description", "enrichment_type")])
        
        if (nrow(unique_pathways) > 0) {
          display_df <- data.frame(
            `Pathway ID` = unique_pathways$ID,
            `Pathway Name` = unique_pathways$Description,
            `Database` = unique_pathways$enrichment_type,
            stringsAsFactors = FALSE
          )
          return(display_df)
        } else {
          return(data.frame(Message = "Unable to format pathway data"))
        }
        
      }, error = function(e) {
        cat("[DETAILS ERROR] extract_shared_pathways:", e$message, "\n")
        return(data.frame(Error = paste("Failed to extract pathways:", e$message)))
      })
    }
    
    # Helper function for null coalescing
    `%||%` <- function(a, b) if (is.null(a)) b else a
    
    # Helper function for string manipulation
    str_to_title <- function(x) {
      gsub("(^|[[:space:]])([[:alpha:]])", "\\1\\U\\2", x, perl = TRUE)
    }
    
    # Return reactive values for potential use by other modules
    return(list(
      analysis_results = reactive({ values$analysis_results }),
      analysis_running = reactive({ values$analysis_running })
    ))
  })
}