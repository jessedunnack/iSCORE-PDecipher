#' Perturb-seq Only Analysis Module
#'
#' Dedicated interface for analyzing Perturb-seq datasets WITHOUT mutation data
#' Focuses on perturbation comparisons and p-value correction comparisons
#'
#' Created: October 24, 2025
#' @author Claude

#' Perturb-seq Only UI
#' @param id Module namespace ID
#' @export
mod_perturbseq_only_ui <- function(id) {
  ns <- NS(id)

  tagList(
    fluidRow(
      column(12,
        h3("Perturb-seq Analysis (Pooled Data)"),
        p("Analyze FDR-corrected pooled MixScale results from Perturb-seq experiments")
      )
    ),

    fluidRow(
      # Left panel: Data configuration
      column(3,
        wellPanel(
          h4("Dataset Configuration"),

          # Dataset type selector
          selectInput(
            ns("dataset_type"),
            "Dataset:",
            choices = c("FPD" = "FPD", "CRISPRi" = "CRISPRi"),
            selected = "FPD"
          ),

          # P-value correction selector
          selectInput(
            ns("pval_correction"),
            "P-value Correction:",
            choices = c(
              "Benjamini-Hochberg (Recommended)" = "p_weight_BH",
              "Original (Uncorrected)" = "p_weight",
              "Bonferroni (Very Conservative)" = "p_weight_bonferroni"
            ),
            selected = "p_weight_BH"
          ),

          hr(),

          # Data loading
          actionButton(
            ns("load_data"),
            "Load Data",
            class = "btn-primary",
            icon = icon("database")
          ),

          br(), br(),

          # Data status
          uiOutput(ns("data_status"))
        ),

        # Help panel
        wellPanel(
          h5(icon("info-circle"), " About P-value Corrections"),
          tags$small(
            tags$ul(
              tags$li(tags$strong("BH (Recommended):"), " Balances sensitivity and false discovery rate"),
              tags$li(tags$strong("Uncorrected:"), " Maximum sensitivity, higher false positives"),
              tags$li(tags$strong("Bonferroni:"), " Most conservative, fewer false positives")
            )
          )
        )
      ),

      # Right panel: Analysis and visualization
      column(9,
        tabsetPanel(
          id = ns("analysis_tabs"),

          # Tab 1: Data Overview
          tabPanel(
            "Overview",
            br(),
            fluidRow(
              column(6,
                h4("Dataset Summary"),
                tableOutput(ns("data_summary"))
              ),
              column(6,
                h4("P-value Distribution"),
                plotOutput(ns("pval_distribution"))
              )
            ),
            br(),
            h4("Available Perturbations"),
            DT::dataTableOutput(ns("perturbation_table"))
          ),

          # Tab 2: Perturbation Comparison
          tabPanel(
            "Perturbation Analysis",
            br(),
            fluidRow(
              column(4,
                selectInput(ns("selected_pert"), "Select Perturbation:", choices = NULL),
                selectInput(ns("selected_cluster"), "Select Cluster:", choices = NULL)
              ),
              column(8,
                h4("Differential Expression Results"),
                DT::dataTableOutput(ns("de_results_table"))
              )
            )
          ),

          # Tab 3: P-value Correction Comparison
          tabPanel(
            "Correction Comparison",
            br(),
            p("Compare the same perturbation across different p-value corrections"),
            fluidRow(
              column(12,
                h4("Significant Gene Counts by Correction Method"),
                plotOutput(ns("correction_comparison_plot"))
              )
            ),
            br(),
            fluidRow(
              column(12,
                h4("Venn Diagram: Overlap Between Corrections"),
                plotOutput(ns("correction_venn"))
              )
            )
          )
        )
      )
    )
  )
}

#' Perturb-seq Only Server
#' @param id Module namespace ID
#' @export
mod_perturbseq_only_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Reactive values
    values <- reactiveValues(
      data = NULL,
      perturbations = NULL,
      clusters = NULL,
      loading = FALSE
    )

    # Load data when button clicked
    observeEvent(input$load_data, {
      values$loading <- TRUE

      withProgress(message = "Loading pooled MixScale data...", value = 0, {

        tryCatch({
          # Determine directory based on dataset type
          # NEW SHORT PATH (under Windows 260 char limit)
          # Use Windows format (E:/) which also works in WSL
          base_dir <- "E:/THESIS/scRNASeq/mixscale"

          if (input$dataset_type == "FPD") {
            mixscale_dir <- file.path(base_dir, "CRISPRi_PerturbSeq_Reports_all_FPD_no_multiplets_noExptSplit")
          } else {
            mixscale_dir <- file.path(base_dir, "CRISPRi_PerturbSeq_Reports_all_CRISPRi_no_multiplets_noExptSplit")
          }

          incProgress(0.3, detail = "Loading data files...")

          # Load data using data manager function
          data <- get_pooled_mixscale_data(
            mixscale_dir = mixscale_dir,
            pval_column = input$pval_correction,
            dataset_type = input$dataset_type,
            force_reload = TRUE
          )

          incProgress(0.5, detail = "Processing...")

          values$data <- data
          values$perturbations <- names(data)
          values$clusters <- unique(unlist(lapply(data, names)))

          incProgress(0.2, detail = "Complete!")

          showNotification("Data loaded successfully!", type = "message")

        }, error = function(e) {
          showNotification(paste("Error loading data:", e$message), type = "error")
        })

      })

      values$loading <- FALSE
    })

    # Update perturbation selector
    observe({
      if (!is.null(values$perturbations)) {
        updateSelectInput(session, "selected_pert", choices = values$perturbations)
      }
    })

    # Update cluster selector
    observe({
      if (!is.null(values$clusters)) {
        updateSelectInput(session, "selected_cluster", choices = values$clusters)
      }
    })

    # Data status output
    output$data_status <- renderUI({
      if (is.null(values$data)) {
        tags$div(
          class = "alert alert-info",
          icon("info-circle"),
          " No data loaded. Click 'Load Data' to begin."
        )
      } else {
        tags$div(
          class = "alert alert-success",
          icon("check-circle"),
          " Data loaded: ",
          tags$strong(length(values$perturbations)), " perturbations, ",
          tags$strong(length(values$clusters)), " clusters"
        )
      }
    })

    # Data summary table
    output$data_summary <- renderTable({
      req(values$data)

      data.frame(
        Metric = c(
          "Dataset Type",
          "P-value Correction",
          "Total Perturbations",
          "Total Clusters",
          "Loaded At"
        ),
        Value = c(
          input$dataset_type,
          input$pval_correction,
          length(values$perturbations),
          length(values$clusters),
          as.character(Sys.time())
        )
      )
    })

    # Perturbation table
    output$perturbation_table <- DT::renderDataTable({
      req(values$data)

      # Create summary table of perturbations
      pert_summary <- data.frame(
        Perturbation = values$perturbations,
        Clusters = sapply(values$perturbations, function(p) {
          length(values$data[[p]])
        })
      )

      DT::datatable(
        pert_summary,
        options = list(
          pageLength = 20,
          scrollX = TRUE
        ),
        filter = "top"
      )
    })

    # DE results table for selected perturbation/cluster
    output$de_results_table <- DT::renderDataTable({
      req(values$data, input$selected_pert, input$selected_cluster)

      # Get data for selected combination
      pert_data <- values$data[[input$selected_pert]]
      if (is.null(pert_data)) return(NULL)

      cluster_data <- pert_data[[input$selected_cluster]]
      if (is.null(cluster_data)) return(NULL)

      results <- cluster_data$results

      # Filter to significant genes (p < 0.05)
      pval_col <- input$pval_correction
      if (pval_col %in% colnames(results)) {
        sig_results <- results[results[[pval_col]] < 0.05, ]

        # Select key columns
        display_cols <- c("gene_ID", "log2FC", pval_col)
        display_results <- sig_results[, display_cols, drop = FALSE]

        DT::datatable(
          display_results,
          options = list(
            pageLength = 25,
            scrollX = TRUE,
            order = list(list(3, 'asc'))  # Sort by p-value
          ),
          rownames = FALSE
        ) %>%
          DT::formatRound(c("log2FC", pval_col), digits = 4)
      } else {
        NULL
      }
    })

    # P-value distribution plot
    output$pval_distribution <- renderPlot({
      req(values$data, input$selected_pert, input$selected_cluster)

      pert_data <- values$data[[input$selected_pert]]
      if (is.null(pert_data)) return(NULL)

      cluster_data <- pert_data[[input$selected_cluster]]
      if (is.null(cluster_data)) return(NULL)

      results <- cluster_data$results
      pval_col <- input$pval_correction

      if (pval_col %in% colnames(results)) {
        hist(
          results[[pval_col]],
          breaks = 50,
          main = paste("P-value Distribution -", input$pval_correction),
          xlab = "P-value",
          ylab = "Frequency",
          col = "steelblue",
          border = "white"
        )
        abline(v = 0.05, col = "red", lwd = 2, lty = 2)
      }
    })

    # Placeholder for comparison plots
    output$correction_comparison_plot <- renderPlot({
      plot(1, 1, type = "n", main = "Coming soon: Correction comparison")
    })

    output$correction_venn <- renderPlot({
      plot(1, 1, type = "n", main = "Coming soon: Venn diagram")
    })

  })
}
