# iSCORE-PDecipher: Comprehensive Perturb-seq Analysis App
# Full-featured app for FDR-corrected pooled MixScale data
# Version 0.5.0 - Created October 25, 2025

# Load required files
source("global.R")
source("R/data_manager.R")

# Load utility modules
if (file.exists("modules/umap_cache_integration.R")) {
  source("modules/umap_cache_integration.R")
}

# Load analysis modules
source("modules/mod_de_results.R")          # DE Results with volcano plots
source("modules/mod_heatmap.R")             # Interactive heatmaps
source("modules/mod_enrichment_gene_display_v2.R")  # Enrichment gene display

# UI Definition
ui <- fluidPage(
  title = "iSCORE-PDecipher: Perturb-seq Analysis",

  # Include shinyjs for dynamic UI
  useShinyjs(),

  # Custom CSS
  tags$head(
    tags$style(HTML("
      .dataset-banner {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        color: white;
        padding: 15px;
        margin-bottom: 20px;
        border-radius: 5px;
        box-shadow: 0 2px 4px rgba(0,0,0,0.1);
      }

      .dataset-banner h3 {
        margin: 0 0 10px 0;
        color: white;
      }

      .dataset-info {
        display: flex;
        justify-content: space-between;
        align-items: center;
        flex-wrap: wrap;
        gap: 10px;
      }

      .dataset-stats {
        display: flex;
        gap: 20px;
        flex-wrap: wrap;
      }

      .dataset-stat {
        display: flex;
        align-items: center;
        gap: 5px;
      }

      .pval-badge {
        background: white;
        color: #667eea;
        padding: 5px 10px;
        border-radius: 3px;
        font-weight: bold;
        font-size: 0.9em;
      }

      .alert-perturb {
        background-color: #e7f3ff;
        border-left: 4px solid #2196F3;
        padding: 10px;
        margin: 10px 0;
      }
    "))
  ),

  # Dataset banner
  div(class = "dataset-banner",
    div(class = "dataset-info",
      div(
        h3(textOutput("app_title", inline = TRUE)),
        div(class = "dataset-stats",
          div(class = "dataset-stat",
            icon("database"),
            textOutput("dataset_type_display", inline = TRUE)
          ),
          div(class = "dataset-stat",
            icon("dna"),
            textOutput("perturbation_count", inline = TRUE)
          ),
          div(class = "dataset-stat",
            icon("project-diagram"),
            textOutput("cluster_count", inline = TRUE)
          ),
          div(class = "pval-badge",
            icon("calculator"),
            textOutput("pval_correction_display", inline = TRUE)
          )
        )
      ),
      div(
        actionButton("reload_data_btn", "Reload Data", icon = icon("refresh"), class = "btn-light"),
        actionButton("change_settings_btn", "Settings", icon = icon("cog"), class = "btn-light")
      )
    )
  ),

  # Main navigation tabs
  navbarPage(
    "Analysis",
    id = "main_tabs",

    # Overview Tab
    tabPanel(
      "Overview",
      icon = icon("home"),
      br(),
      fluidRow(
        column(12,
          div(class = "alert-perturb",
            h4(icon("info-circle"), " About This Dataset"),
            p("This app analyzes FDR-corrected pooled MixScale results from Perturb-seq experiments."),
            p(strong("Key Features:")),
            tags$ul(
              tags$li("Pooled analysis (NOT experiment-split)"),
              tags$li("Three p-value correction methods: Benjamini-Hochberg (BH), Uncorrected, Bonferroni"),
              tags$li("Consolidated enrichment results (GO, KEGG, Reactome, WikiPathways, STRING, GSEA)"),
              tags$li("Interactive visualizations and downloadable results")
            )
          )
        )
      ),
      fluidRow(
        column(6,
          wellPanel(
            h4("Dataset Configuration"),
            uiOutput("overview_config")
          )
        ),
        column(6,
          wellPanel(
            h4("Data Summary"),
            tableOutput("data_summary_table")
          )
        )
      ),
      fluidRow(
        column(12,
          wellPanel(
            h4("Available Perturbations"),
            DT::dataTableOutput("perturbation_overview_table")
          )
        )
      )
    ),

    # DE Results Tab
    tabPanel(
      "DE Results",
      icon = icon("chart-line"),
      br(),
      fluidRow(
        column(3,
          wellPanel(
            h4("Select Data"),
            selectInput("de_perturbation", "Perturbation:", choices = NULL),
            selectInput("de_cluster", "Cluster:", choices = NULL),
            hr(),
            h5("Volcano Plot Settings"),
            sliderInput("de_pval_threshold", "P-value threshold:",
                       min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            sliderInput("de_lfc_threshold", "Log2FC threshold:",
                       min = 0, max = 2, value = 0.25, step = 0.05),
            hr(),
            checkboxInput("de_show_labels", "Show gene labels", value = TRUE),
            sliderInput("de_top_genes", "Top genes to label:",
                       min = 5, max = 50, value = 10, step = 5),
            hr(),
            downloadButton("de_download_plot", "Download Plot"),
            downloadButton("de_download_data", "Download Data")
          )
        ),
        column(9,
          tabsetPanel(
            tabPanel("Volcano Plot",
                    br(),
                    plotOutput("de_volcano_plot", height = "600px")),
            tabPanel("Data Table",
                    br(),
                    DT::dataTableOutput("de_results_table")),
            tabPanel("Summary",
                    br(),
                    plotOutput("de_summary_plot", height = "400px"))
          )
        )
      )
    ),

    # Enrichment Tab
    tabPanel(
      "Enrichment",
      icon = icon("sitemap"),
      br(),
      fluidRow(
        column(3,
          wellPanel(
            h4("Filter Enrichment"),
            selectInput("enrich_perturbation", "Perturbation:", choices = NULL),
            selectInput("enrich_cluster", "Cluster:", choices = NULL),
            selectInput("enrich_type", "Enrichment Type:",
                       choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG",
                                 "Reactome", "WikiPathways", "STRING", "GSEA"),
                       selected = "GO_BP"),
            selectInput("enrich_direction", "Gene Direction:",
                       choices = c("All" = "ALL", "Up-regulated" = "UP",
                                 "Down-regulated" = "DOWN"),
                       selected = "ALL"),
            sliderInput("enrich_pval_max", "Max Adjusted P-value:",
                       min = 0.001, max = 0.05, value = 0.05, step = 0.001),
            sliderInput("enrich_max_terms", "Max Terms to Display:",
                       min = 10, max = 100, value = 30, step = 10),
            hr(),
            textInput("enrich_search", "Search Terms:",
                     placeholder = "e.g., mitochondr, synap"),
            hr(),
            downloadButton("enrich_download", "Download Results")
          )
        ),
        column(9,
          tabsetPanel(
            tabPanel("Dot Plot",
                    br(),
                    plotOutput("enrich_dotplot", height = "700px")),
            tabPanel("Bar Plot",
                    br(),
                    plotOutput("enrich_barplot", height = "700px")),
            tabPanel("Data Table",
                    br(),
                    DT::dataTableOutput("enrich_table")),
            tabPanel("Gene Display",
                    br(),
                    mod_enrichment_gene_display_ui("enrichment_genes"))
          )
        )
      )
    ),

    # Heatmaps Tab
    tabPanel(
      "Heatmaps",
      icon = icon("th"),
      br(),
      fluidRow(
        column(3,
          wellPanel(
            h4("Heatmap Settings"),
            selectInput("hm_cluster", "Cluster:", choices = NULL),
            selectInput("hm_type", "Heatmap Type:",
                       choices = c("P-values (-log10)" = "pvalue",
                                 "Gene Ratio" = "generatio",
                                 "Gene Count" = "count"),
                       selected = "pvalue"),
            checkboxGroupInput("hm_enrich_types", "Enrichment Types:",
                             choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG",
                                       "Reactome", "WikiPathways", "STRING", "GSEA"),
                             selected = c("GO_BP", "KEGG", "Reactome")),
            sliderInput("hm_max_terms", "Max Terms:",
                       min = 10, max = 100, value = 30, step = 5),
            selectInput("hm_direction", "Gene Direction:",
                       choices = c("All" = "ALL_DIR", "Up" = "UP",
                                 "Down" = "DOWN"),
                       selected = "ALL_DIR"),
            hr(),
            h5("Visual Settings"),
            checkboxInput("hm_cluster_rows", "Cluster rows", value = TRUE),
            checkboxInput("hm_cluster_cols", "Cluster columns", value = TRUE),
            selectInput("hm_color_scheme", "Color Scheme:",
                       choices = c("Blue-Red" = "RdBu", "Viridis" = "viridis",
                                 "Plasma" = "plasma", "Inferno" = "inferno"),
                       selected = "RdBu"),
            hr(),
            actionButton("hm_generate", "Generate Heatmap",
                        icon = icon("chart-bar"), class = "btn-primary"),
            br(), br(),
            downloadButton("hm_download", "Download Heatmap")
          )
        ),
        column(9,
          wellPanel(
            plotOutput("heatmap_plot", height = "800px")
          )
        )
      )
    ),

    # P-value Correction Comparison Tab
    tabPanel(
      "Correction Comparison",
      icon = icon("balance-scale"),
      br(),
      div(class = "alert-perturb",
        h4(icon("info-circle"), " Compare P-value Correction Methods"),
        p("This tab allows you to compare results across three FDR correction methods:"),
        tags$ul(
          tags$li(strong("Benjamini-Hochberg (BH):"), " Recommended for most analyses. Balances sensitivity and specificity."),
          tags$li(strong("Uncorrected:"), " Maximum sensitivity but higher false discovery rate."),
          tags$li(strong("Bonferroni:"), " Very conservative. Minimizes false positives but may miss true positives.")
        ),
        p("Use this comparison to understand how FDR correction affects your results.")
      ),
      fluidRow(
        column(3,
          wellPanel(
            h4("Select Data"),
            selectInput("comp_perturbation", "Perturbation:", choices = NULL),
            selectInput("comp_cluster", "Cluster:", choices = NULL),
            sliderInput("comp_pval_threshold", "P-value threshold:",
                       min = 0.001, max = 0.1, value = 0.05, step = 0.001),
            sliderInput("comp_lfc_threshold", "Log2FC threshold:",
                       min = 0, max = 2, value = 0.25, step = 0.05)
          )
        ),
        column(9,
          tabsetPanel(
            tabPanel("DEG Counts",
                    br(),
                    plotOutput("comp_deg_barplot", height = "400px"),
                    br(),
                    tableOutput("comp_deg_table")),
            tabPanel("Venn Diagram",
                    br(),
                    plotOutput("comp_venn", height = "600px")),
            tabPanel("Enrichment Comparison",
                    br(),
                    p("Compare enriched terms across correction methods"),
                    plotOutput("comp_enrich_plot", height = "600px"))
          )
        )
      )
    ),

    # Settings Tab
    tabPanel(
      "Settings",
      icon = icon("cog"),
      br(),
      fluidRow(
        column(6,
          wellPanel(
            h4("Dataset Configuration"),
            selectInput("settings_dataset", "Dataset Type:",
                       choices = c("FPD" = "FPD", "CRISPRi" = "CRISPRi"),
                       selected = "FPD"),
            selectInput("settings_pval", "P-value Correction:",
                       choices = c("Benjamini-Hochberg (BH)" = "p_weight_BH",
                                 "Uncorrected" = "p_weight",
                                 "Bonferroni" = "p_weight_bonferroni"),
                       selected = "p_weight_BH"),
            actionButton("apply_settings", "Apply Settings",
                        icon = icon("check"), class = "btn-primary")
          ),
          wellPanel(
            h4("Data Paths"),
            p(strong("MixScale Directory:")),
            verbatimTextOutput("settings_mixscale_dir"),
            p(strong("Enrichment Directory:")),
            verbatimTextOutput("settings_enrich_dir")
          )
        ),
        column(6,
          wellPanel(
            h4("Performance Settings"),
            checkboxInput("settings_cache", "Enable plot caching", value = TRUE),
            checkboxInput("settings_preview", "Preview mode for large datasets", value = FALSE),
            sliderInput("settings_max_cells", "Max cells for interactive plots:",
                       min = 10000, max = 500000, value = 100000, step = 10000)
          ),
          wellPanel(
            h4("Actions"),
            actionButton("clear_cache_btn", "Clear Cache",
                        icon = icon("trash"), class = "btn-warning"),
            br(), br(),
            actionButton("reset_app_btn", "Reset App",
                        icon = icon("redo"), class = "btn-danger")
          )
        )
      )
    )
  )
)

# Server Definition
server <- function(input, output, session) {

  # Reactive values
  rv <- reactiveValues(
    data = NULL,
    enrichment = NULL,
    dataset_type = "FPD",
    pval_correction = "p_weight_BH",
    perturbations = NULL,
    clusters = NULL,
    data_loaded = FALSE
  )

  # Initialize app with default data
  observeEvent(session$clientData, {
    if (!rv$data_loaded) {
      load_initial_data()
    }
  }, once = TRUE)

  # Load initial data function
  load_initial_data <- function() {
    withProgress(message = "Loading Perturb-seq data...", value = 0, {

      incProgress(0.3, detail = "Loading MixScale results...")

      # Determine directory
      base_dir <- "E:/THESIS/scRNASeq/mixscale"
      if (rv$dataset_type == "FPD") {
        mixscale_dir <- file.path(base_dir, "CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit")
      } else {
        mixscale_dir <- file.path(base_dir, "CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit")
      }

      # Load MixScale data
      rv$data <- get_pooled_mixscale_data(
        mixscale_dir = mixscale_dir,
        pval_column = rv$pval_correction,
        dataset_type = rv$dataset_type,
        force_reload = TRUE
      )

      incProgress(0.3, detail = "Loading enrichment results...")

      # Load enrichment data
      pval_type <- switch(rv$pval_correction,
        "p_weight" = "none",
        "p_weight_BH" = "BH",
        "p_weight_bonferroni" = "bonferroni"
      )

      rv$enrichment <- import_enrichment_with_correction(
        dataset = rv$dataset_type,
        pval_correction = pval_type
      )

      incProgress(0.3, detail = "Preparing UI...")

      # Extract perturbations and clusters
      if (!is.null(rv$data)) {
        rv$perturbations <- names(rv$data)
        rv$clusters <- unique(unlist(lapply(rv$data, names)))
        rv$data_loaded <- TRUE
      }

      incProgress(0.1, detail = "Complete!")
    })
  }

  # Update UI elements when data loads
  observe({
    req(rv$data_loaded, rv$perturbations, rv$clusters)

    # Update all perturbation selectors
    updateSelectInput(session, "de_perturbation", choices = rv$perturbations)
    updateSelectInput(session, "enrich_perturbation", choices = rv$perturbations)
    updateSelectInput(session, "comp_perturbation", choices = rv$perturbations)

    # Update all cluster selectors
    cluster_choices <- setNames(rv$clusters, paste("Cluster", rv$clusters))
    updateSelectInput(session, "de_cluster", choices = cluster_choices)
    updateSelectInput(session, "enrich_cluster", choices = cluster_choices)
    updateSelectInput(session, "comp_cluster", choices = cluster_choices)
    updateSelectInput(session, "hm_cluster", choices = c("All" = "all", cluster_choices))
  })

  # Banner outputs
  output$app_title <- renderText({
    "iSCORE-PDecipher: Perturb-seq Analysis (v0.5.0)"
  })

  output$dataset_type_display <- renderText({
    req(rv$dataset_type)
    rv$dataset_type
  })

  output$perturbation_count <- renderText({
    req(rv$perturbations)
    paste(length(rv$perturbations), "perturbations")
  })

  output$cluster_count <- renderText({
    req(rv$clusters)
    paste(length(rv$clusters), "clusters")
  })

  output$pval_correction_display <- renderText({
    req(rv$pval_correction)
    switch(rv$pval_correction,
      "p_weight" = "Uncorrected",
      "p_weight_BH" = "BH-corrected",
      "p_weight_bonferroni" = "Bonferroni"
    )
  })

  # Overview tab outputs
  output$overview_config <- renderUI({
    tagList(
      p(strong("Dataset:"), rv$dataset_type),
      p(strong("P-value Correction:"), rv$pval_correction),
      p(strong("Data Directory:"), "E:/THESIS/scRNASeq/mixscale")
    )
  })

  output$data_summary_table <- renderTable({
    req(rv$data_loaded)
    data.frame(
      Metric = c("Total Perturbations", "Total Clusters", "Loaded At"),
      Value = c(length(rv$perturbations), length(rv$clusters), as.character(Sys.time()))
    )
  })

  output$perturbation_overview_table <- DT::renderDataTable({
    req(rv$data)

    pert_summary <- data.frame(
      Perturbation = rv$perturbations,
      Clusters = sapply(rv$perturbations, function(p) length(rv$data[[p]]))
    )

    DT::datatable(pert_summary, options = list(pageLength = 20), filter = "top")
  })

  # Settings tab outputs
  output$settings_mixscale_dir <- renderText({
    paste0("E:/THESIS/scRNASeq/mixscale/CRISPRi_PerturbSeq_Reports_all_",
           rv$dataset_type, "_no_multiplets_noExptSplit")
  })

  output$settings_enrich_dir <- renderText({
    pval_suffix <- switch(rv$pval_correction,
      "p_weight" = "_p_weight",
      "p_weight_BH" = "_p_weight_BH",
      "p_weight_bonferroni" = "_p_weight_bonferroni"
    )
    paste0("E:/THESIS/scRNASeq/mixscale/enrichment_results_",
           rv$dataset_type, pval_suffix)
  })

  # Placeholder outputs for tabs (to be implemented)
  output$de_volcano_plot <- renderPlot({
    plot(1, 1, type = "n", main = "Volcano plot - Coming soon", xlab = "", ylab = "")
    text(1, 1, "Select perturbation and cluster to view volcano plot")
  })

  output$de_results_table <- DT::renderDataTable({
    data.frame(
      Message = "Select perturbation and cluster above",
      stringsAsFactors = FALSE
    )
  })

  output$enrich_dotplot <- renderPlot({
    plot(1, 1, type = "n", main = "Enrichment dot plot - Coming soon", xlab = "", ylab = "")
  })

  output$heatmap_plot <- renderPlot({
    plot(1, 1, type = "n", main = "Click 'Generate Heatmap' to create visualization", xlab = "", ylab = "")
  })

  # Apply settings button
  observeEvent(input$apply_settings, {
    rv$dataset_type <- input$settings_dataset
    rv$pval_correction <- input$settings_pval
    rv$data_loaded <- FALSE
    load_initial_data()
  })

  # Reload data button
  observeEvent(input$reload_data_btn, {
    rv$data_loaded <- FALSE
    load_initial_data()
  })

}

# Run the app
shinyApp(ui = ui, server = server)
