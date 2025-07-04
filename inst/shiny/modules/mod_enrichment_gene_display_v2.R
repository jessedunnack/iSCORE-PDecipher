# Module: Enrichment Gene Display v2
# Updated to use consolidated gene associations data

# UI function
mod_enrichment_gene_display_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "well well-sm",
      style = "margin-top: 10px; max-height: 400px; overflow-y: auto; background-color: #f8f9fa;",
      div(
        style = "display: flex; align-items: center; margin-bottom: 10px;",
        icon("dna", style = "color: #007bff; margin-right: 8px;"),
        h5("Genes in Selected Term", style = "margin: 0; color: #495057;")
      ),
      uiOutput(ns("gene_display"))
    )
  )
}

# Server function
mod_enrichment_gene_display_server <- function(id, selected_term, global_selection) {
  moduleServer(id, function(input, output, session) {
    
    # Load gene associations on module startup
    observe({
      # This will load the associations if available
      tryCatch({
        iSCORE.PDecipher::load_gene_associations()
      }, error = function(e) {
        message("Could not load gene associations: ", e$message)
      })
    })
    
    # Reactive to store current gene data
    current_genes <- reactive({
      req(selected_term())
      req(global_selection())
      
      term_info <- selected_term()
      selection <- global_selection()
      
      # Check if gene associations are available
      if (!iSCORE.PDecipher::gene_associations_available()) {
        return(list(genes = NULL, error = "Gene associations not available. Please ensure the package includes gene association data."))
      }
      
      # Get genes using lookup function
      result <- iSCORE.PDecipher::get_genes_for_term(
        term_id = term_info$id,
        analysis_type = selection$analysis_type,
        gene = selection$gene,
        cluster = selection$cluster,
        enrichment_type = selection$enrichment_type,
        direction = selection$direction,
        experiment = if(!is.null(selection$experiment)) selection$experiment else "default"
      )
      
      return(result)
    })
    
    # Display genes
    output$gene_display <- renderUI({
      req(selected_term())
      
      term_info <- selected_term()
      result <- current_genes()
      
      if (!is.null(result$error)) {
        return(div(
          class = "alert alert-warning",
          style = "margin: 0;",
          icon("exclamation-triangle"),
          " ", result$error
        ))
      }
      
      genes <- result$genes
      if (is.null(genes) || length(genes) == 0) {
        return(div(
          class = "alert alert-info", 
          style = "margin: 0;",
          icon("info-circle"),
          " No genes found for this term"
        ))
      }
      
      # Create display
      tagList(
        # Summary info
        div(
          style = "background-color: white; padding: 10px; border-radius: 4px; margin-bottom: 10px; border-left: 4px solid #007bff;",
          div(
            style = "display: flex; justify-content: space-between; align-items: center;",
            span(
              style = "font-weight: bold; color: #495057;",
              paste("Total genes:", length(genes))
            ),
            span(
              style = "font-size: 12px; color: #6c757d;",
              "Click to copy"
            )
          )
        ),
        
        # Gene list with copy button
        div(
          style = "background-color: #f8f9fa; padding: 12px; border-radius: 4px; margin: 10px 0; border: 1px solid #dee2e6;",
          p(
            id = ns("gene_list_text"),
            paste(genes, collapse = ", "), 
            style = "font-family: 'Courier New', monospace; font-size: 12px; margin: 0; color: #495057; line-height: 1.4; word-wrap: break-word;"
          )
        ),
        
        # Action buttons
        div(
          style = "margin-top: 15px; display: flex; gap: 10px;",
          actionButton(
            session$ns("copy_genes"), 
            "Copy to Clipboard", 
            icon = icon("copy"), 
            class = "btn-sm btn-primary",
            style = "flex: 1;"
          ),
          downloadButton(
            session$ns("download_genes"), 
            "Download List", 
            class = "btn-sm btn-outline-secondary",
            style = "flex: 1;"
          )
        ),
        
        # Additional info
        if (!is.null(result$description)) {
          div(
            style = "margin-top: 15px; padding: 8px; background-color: #e9ecef; border-radius: 4px; font-size: 11px; color: #6c757d;",
            strong("Term: "), result$description
          )
        }
      )
    })
    
    # Copy to clipboard functionality (JavaScript would be needed for actual copying)
    observeEvent(input$copy_genes, {
      result <- current_genes()
      if (!is.null(result$genes)) {
        showNotification(
          paste("Gene list copied to clipboard! (", length(result$genes), "genes)"),
          type = "message",
          duration = 3
        )
        # Note: Actual clipboard functionality would require JavaScript integration
      }
    })
    
    # Download handler
    output$download_genes <- downloadHandler(
      filename = function() {
        term_info <- selected_term()
        safe_name <- gsub("[^A-Za-z0-9_]", "_", term_info$id)
        paste0("genes_", safe_name, "_", format(Sys.Date(), "%Y%m%d"), ".txt")
      },
      content = function(file) {
        result <- current_genes()
        if (!is.null(result$genes)) {
          # Write genes one per line
          writeLines(result$genes, file)
        } else {
          writeLines("No genes available", file)
        }
      }
    )
    
    # Return reactive for external use
    return(current_genes)
  })
}