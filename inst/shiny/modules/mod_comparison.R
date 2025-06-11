# Module: Method Comparison
# Compares MAST and MixScale results

mod_comparison_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Comparison Settings"),
      
      selectInput(ns("gene_select"),
                  "Select Gene:",
                  choices = character(0),
                  selected = NULL),
      
      selectInput(ns("cluster_select"),
                  "Select Cluster:",
                  choices = character(0),
                  selected = NULL),
      
      selectInput(ns("enrichment_type"),
                  "Enrichment Type:",
                  choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "WikiPathways", "STRING"),
                  selected = "GO_BP"),
      
      radioButtons(ns("direction"),
                   "Direction:",
                   choices = c("ALL", "UP", "DOWN"),
                   selected = "ALL"),
      
      radioButtons(ns("comparison_method"),
                   "Comparison Method:",
                   choices = list(
                     "Intersection (terms in both)" = "intersection",
                     "Union (terms in either)" = "union",
                     "MAST-specific" = "mast_only",
                     "MixScale-specific" = "mixscale_only"
                   ),
                   selected = "intersection"),
      
      br(),
      actionButton(ns("compare_btn"), "Compare Methods", 
                   class = "btn-primary", width = "100%")
    ),
    
    mainPanel(
      width = 9,
      fluidRow(
        column(4, uiOutput(ns("mast_terms_box"))),
        column(4, uiOutput(ns("mixscale_terms_box"))),
        column(4, uiOutput(ns("shared_terms_box")))
      ),
      br(),
      tabsetPanel(
        tabPanel("Overview",
                 br(),
                 plotOutput(ns("venn_plot"), height = "400px"),
                 br(),
                 DT::dataTableOutput(ns("comparison_table"))),
        
        tabPanel("Convergent Pathways",
                 br(),
                 plotOutput(ns("convergent_plot"), height = "600px")),
        
        tabPanel("Method-Specific",
                 br(),
                 fluidRow(
                   column(6, 
                          h4("MAST-Specific Terms"),
                          DT::dataTableOutput(ns("mast_specific_table"))),
                   column(6,
                          h4("MixScale-Specific Terms"),
                          DT::dataTableOutput(ns("mixscale_specific_table")))
                 ))
      )
    )
  )
}

mod_comparison_server <- function(id, app_data, pval_threshold) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Initialize reactive values
    comparison_data <- reactiveValues(
      mast_results = NULL,
      mixscale_results = NULL,
      comparison_results = NULL
    )
    
    # Update gene choices
    observe({
      req(app_data$available_genes)
      updateSelectInput(session, "gene_select",
                        choices = app_data$available_genes,
                        selected = app_data$available_genes[1])
    })
    
    # Update cluster choices based on gene
    observe({
      req(input$gene_select)
      clusters <- get_clusters_for_gene(input$gene_select)
      updateSelectInput(session, "cluster_select",
                        choices = clusters,
                        selected = clusters[1])
    })
    
    # Helper function to get clusters for a gene
    get_clusters_for_gene <- function(gene) {
      base_path <- APP_CONFIG$enrichment_results_path
      mast_path <- file.path(base_path, "MAST", gene)
      mixscale_path <- file.path(base_path, "MixScale", gene)
      
      clusters <- character()
      if (dir.exists(mast_path)) {
        clusters <- c(clusters, list.dirs(mast_path, full.names = FALSE, recursive = FALSE))
      }
      if (dir.exists(mixscale_path)) {
        clusters <- c(clusters, list.dirs(mixscale_path, full.names = FALSE, recursive = FALSE))
      }
      
      unique(clusters[clusters != ""])
    }
    
    # Helper function to extract significant terms
    extract_significant_terms <- function(result, enrichment_type, threshold = 0.05) {
      if (is.null(result)) return(data.frame())
      
      # Handle different enrichment result formats
      if (enrichment_type == "STRING") {
        if (is.list(result) && "enrichment" %in% names(result)) {
          df <- result$enrichment
          if (!is.null(df) && nrow(df) > 0) {
            df <- df[df$fdr < threshold, ]
            return(data.frame(
              ID = df$term,
              Description = df$description,
              pvalue = df$p_value,
              p.adjust = df$fdr,
              Count = df$number_of_genes,
              stringsAsFactors = FALSE
            ))
          }
        }
      } else if (enrichment_type == "GSEA") {
        if (is.data.frame(result)) {
          df <- result[result$padj < threshold, ]
          return(data.frame(
            ID = df$pathway,
            Description = df$pathway,
            pvalue = df$pval,
            p.adjust = df$padj,
            NES = df$NES,
            stringsAsFactors = FALSE
          ))
        }
      } else {
        # Standard enrichResult S4 object
        if (inherits(result, "enrichResult")) {
          df <- result@result
          df <- df[df$p.adjust < threshold, ]
          return(df[, c("ID", "Description", "pvalue", "p.adjust", "Count")])
        }
      }
      
      return(data.frame())
    }
    
    # Load and compare results
    observeEvent(input$compare_btn, {
      req(input$gene_select, input$cluster_select, input$enrichment_type, input$direction)
      
      showNotification("Loading comparison data...", type = "message", duration = NULL, id = "loading")
      
      base_path <- APP_CONFIG$enrichment_results_path
      
      # Load MAST results
      mast_file <- file.path(base_path, "MAST", input$gene_select, input$cluster_select,
                             "default", input$enrichment_type,
                             paste0(input$enrichment_type, "_", input$direction, ".rds"))
      
      if (file.exists(mast_file)) {
        comparison_data$mast_results <- readRDS(mast_file)
      } else {
        comparison_data$mast_results <- NULL
      }
      
      # Load MixScale results
      mixscale_base <- file.path(base_path, "MixScale", input$gene_select, input$cluster_select)
      mixscale_results <- list()
      
      if (dir.exists(mixscale_base)) {
        experiments <- list.dirs(mixscale_base, full.names = FALSE, recursive = FALSE)
        for (exp in experiments) {
          mixscale_file <- file.path(mixscale_base, exp, input$enrichment_type,
                                     paste0(input$enrichment_type, "_", input$direction, ".rds"))
          if (file.exists(mixscale_file)) {
            mixscale_results[[exp]] <- readRDS(mixscale_file)
          }
        }
      }
      
      # Combine MixScale results
      if (length(mixscale_results) > 0) {
        # For simplicity, take the first experiment or combine them
        comparison_data$mixscale_results <- mixscale_results[[1]]
      } else {
        comparison_data$mixscale_results <- NULL
      }
      
      # Extract terms
      mast_terms <- extract_significant_terms(comparison_data$mast_results, 
                                               input$enrichment_type, 
                                               pval_threshold())
      mixscale_terms <- extract_significant_terms(comparison_data$mixscale_results,
                                                   input$enrichment_type,
                                                   pval_threshold())
      
      # Store processed results
      comparison_data$comparison_results <- list(
        mast_terms = mast_terms,
        mixscale_terms = mixscale_terms,
        mast_ids = if(nrow(mast_terms) > 0) mast_terms$ID else character(),
        mixscale_ids = if(nrow(mixscale_terms) > 0) mixscale_terms$ID else character()
      )
      
      removeNotification("loading")
      showNotification("Comparison complete!", type = "message", duration = 2)
    })
    
    # Value boxes
    output$mast_terms_box <- renderUI({
      value <- ifelse(is.null(comparison_data$comparison_results),
                     0,
                     length(comparison_data$comparison_results$mast_ids))
      
      div(class = "small-box bg-blue",
        div(class = "inner",
          h3(value),
          p("MAST Terms")
        ),
        div(class = "icon",
          icon("dna")
        )
      )
    })
    
    output$mixscale_terms_box <- renderUI({
      value <- ifelse(is.null(comparison_data$comparison_results),
                     0,
                     length(comparison_data$comparison_results$mixscale_ids))
      
      div(class = "small-box bg-green",
        div(class = "inner",
          h3(value),
          p("MixScale Terms")
        ),
        div(class = "icon",
          icon("microscope")
        )
      )
    })
    
    output$shared_terms_box <- renderUI({
      shared_count <- 0
      if (!is.null(comparison_data$comparison_results)) {
        shared_count <- length(intersect(comparison_data$comparison_results$mast_ids,
                                         comparison_data$comparison_results$mixscale_ids))
      }
      
      div(class = "small-box bg-yellow",
        div(class = "inner",
          h3(shared_count),
          p("Shared Terms")
        ),
        div(class = "icon",
          icon("link")
        )
      )
    })
    
    # Venn diagram
    output$venn_plot <- renderPlot({
      req(comparison_data$comparison_results)
      
      # Get term IDs for each method
      mast_ids <- comparison_data$comparison_results$mast_ids
      mixscale_ids <- comparison_data$comparison_results$mixscale_ids
      
      if (length(mast_ids) > 0 || length(mixscale_ids) > 0) {
        # Calculate overlaps
        only_mast <- length(setdiff(mast_ids, mixscale_ids))
        only_mixscale <- length(setdiff(mixscale_ids, mast_ids))
        both <- length(intersect(mast_ids, mixscale_ids))
        
        # Create custom Venn diagram
        par(mar = c(2, 2, 3, 2))
        plot.new()
        
        # Draw circles
        symbols(c(0.35, 0.65), c(0.5, 0.5), circles = c(0.25, 0.25),
                inches = FALSE, add = TRUE,
                fg = c("#e74c3c", "#3498db"), 
                bg = adjustcolor(c("#e74c3c", "#3498db"), alpha = 0.3),
                lwd = 3)
        
        # Add labels
        text(0.2, 0.5, only_mast, cex = 2, font = 2)
        text(0.5, 0.5, both, cex = 2, font = 2)
        text(0.8, 0.5, only_mixscale, cex = 2, font = 2)
        
        # Add method names
        text(0.35, 0.85, "MAST", cex = 1.5, font = 2, col = "#e74c3c")
        text(0.65, 0.85, "MixScale", cex = 1.5, font = 2, col = "#3498db")
        
        # Add title
        title(main = "Overlap of Enriched Terms", cex.main = 1.4)
        
        # Add summary text
        text(0.5, 0.15, paste("Total unique terms:", 
                             length(unique(c(mast_ids, mixscale_ids)))),
             cex = 1.2)
        
      } else {
        plot.new()
        text(0.5, 0.5, "No enriched terms found", 
             cex = 2, col = "gray50")
      }
    })
    
    # Comparison table
    output$comparison_table <- DT::renderDataTable({
      req(comparison_data$comparison_results)
      
      # Create comparison based on method
      if (input$comparison_method == "intersection") {
        shared_ids <- intersect(comparison_data$comparison_results$mast_ids,
                                comparison_data$comparison_results$mixscale_ids)
        mast_subset <- comparison_data$comparison_results$mast_terms[
          comparison_data$comparison_results$mast_terms$ID %in% shared_ids, ]
        mixscale_subset <- comparison_data$comparison_results$mixscale_terms[
          comparison_data$comparison_results$mixscale_terms$ID %in% shared_ids, ]
        
        if (nrow(mast_subset) > 0) {
          comparison_df <- data.frame(
            ID = mast_subset$ID,
            Description = mast_subset$Description,
            MAST_pvalue = mast_subset$pvalue,
            MixScale_pvalue = mixscale_subset$pvalue[match(mast_subset$ID, mixscale_subset$ID)],
            stringsAsFactors = FALSE
          )
        } else {
          comparison_df <- data.frame()
        }
        
      } else if (input$comparison_method == "union") {
        all_ids <- union(comparison_data$comparison_results$mast_ids,
                         comparison_data$comparison_results$mixscale_ids)
        comparison_df <- data.frame(ID = all_ids, stringsAsFactors = FALSE)
        
      } else if (input$comparison_method == "mast_only") {
        comparison_df <- comparison_data$comparison_results$mast_terms
        
      } else if (input$comparison_method == "mixscale_only") {
        comparison_df <- comparison_data$comparison_results$mixscale_terms
      }
      
      if (nrow(comparison_df) > 0) {
        DT::datatable(comparison_df,
                      options = list(pageLength = 10, scrollX = TRUE),
                      rownames = FALSE)
      } else {
        DT::datatable(data.frame(Message = "No terms found for selected comparison"),
                      options = list(dom = 't'),
                      rownames = FALSE)
      }
    })
    
    # Convergent pathways plot
    output$convergent_plot <- renderPlot({
      req(comparison_data$comparison_results)
      
      # Get shared terms
      shared_ids <- comparison_data$comparison_results$shared_ids
      
      if (length(shared_ids) > 0) {
        # Get terms that appear in both MAST and MixScale
        mast_terms <- comparison_data$comparison_results$mast_terms
        mixscale_terms <- comparison_data$comparison_results$mixscale_terms
        
        # Filter to shared terms
        mast_shared <- mast_terms[mast_terms$ID %in% shared_ids, ]
        mixscale_shared <- mixscale_terms[mixscale_terms$ID %in% shared_ids, ]
        
        # Combine data for plotting
        plot_data <- data.frame(
          ID = shared_ids,
          Description = mast_shared$Description[match(shared_ids, mast_shared$ID)],
          MAST_pvalue = -log10(mast_shared$p.adjust[match(shared_ids, mast_shared$ID)]),
          MixScale_pvalue = -log10(mixscale_shared$p.adjust[match(shared_ids, mixscale_shared$ID)]),
          stringsAsFactors = FALSE
        )
        
        # Take top 20 terms by average significance
        plot_data$avg_significance <- (plot_data$MAST_pvalue + plot_data$MixScale_pvalue) / 2
        plot_data <- plot_data[order(plot_data$avg_significance, decreasing = TRUE), ]
        plot_data <- head(plot_data, 20)
        
        # Create scatter plot
        par(mar = c(5, 5, 4, 2))
        plot(plot_data$MAST_pvalue, plot_data$MixScale_pvalue,
             pch = 19, col = adjustcolor("#3c8dbc", alpha = 0.7),
             cex = 1.5,
             xlim = range(c(plot_data$MAST_pvalue, plot_data$MixScale_pvalue)),
             ylim = range(c(plot_data$MAST_pvalue, plot_data$MixScale_pvalue)),
             xlab = expression("MAST -log"[10]*"(p.adjust)"),
             ylab = expression("MixScale -log"[10]*"(p.adjust)"),
             main = "Convergent Pathways: MAST vs MixScale",
             cex.lab = 1.2, cex.main = 1.4)
        
        # Add diagonal line
        abline(0, 1, col = "gray50", lty = 2)
        
        # Add grid
        grid(col = "gray90")
        
        # Label top 5 points
        top_points <- head(plot_data, 5)
        text(top_points$MAST_pvalue, top_points$MixScale_pvalue,
             labels = substr(top_points$Description, 1, 30),
             pos = 3, cex = 0.8, col = "black")
        
        # Add legend
        legend("bottomright", 
               legend = c(paste("Shared terms:", length(shared_ids)),
                         "Equal significance"),
               pch = c(19, NA),
               lty = c(NA, 2),
               col = c("#3c8dbc", "gray50"),
               bty = "n")
        
      } else {
        plot.new()
        text(0.5, 0.5, "No convergent pathways found\nbetween MAST and MixScale", 
             cex = 2, col = "gray50")
      }
    })
    
    # Method-specific tables
    output$mast_specific_table <- DT::renderDataTable({
      req(comparison_data$comparison_results)
      
      mast_only <- setdiff(comparison_data$comparison_results$mast_ids,
                           comparison_data$comparison_results$mixscale_ids)
      
      mast_specific <- comparison_data$comparison_results$mast_terms[
        comparison_data$comparison_results$mast_terms$ID %in% mast_only, ]
      
      DT::datatable(mast_specific,
                    options = list(pageLength = 5, scrollX = TRUE),
                    rownames = FALSE)
    })
    
    output$mixscale_specific_table <- DT::renderDataTable({
      req(comparison_data$comparison_results)
      
      mixscale_only <- setdiff(comparison_data$comparison_results$mixscale_ids,
                               comparison_data$comparison_results$mast_ids)
      
      mixscale_specific <- comparison_data$comparison_results$mixscale_terms[
        comparison_data$comparison_results$mixscale_terms$ID %in% mixscale_only, ]
      
      DT::datatable(mixscale_specific,
                    options = list(pageLength = 5, scrollX = TRUE),
                    rownames = FALSE)
    })
    
    return(comparison_data)
  })
}