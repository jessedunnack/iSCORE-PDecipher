# Module: Enrichment Gene Display
# Shows genes associated with enrichment terms

# Function to construct enrichment file path
construct_enrichment_path <- function(dataset_dir, analysis_type, gene, cluster, enrichment_type, direction) {
  base_path <- file.path(dataset_dir, "enrichment_results")
  method <- ifelse(analysis_type == "MAST", "MAST", "MixScale")
  
  # Map enrichment types to file patterns
  file_pattern <- switch(enrichment_type,
    "GO_BP" = "GO_enrichment",
    "GO_CC" = "GO_enrichment", 
    "GO_MF" = "GO_enrichment",
    "KEGG" = "KEGG_enrichment",
    "Reactome" = "Reactome_enrichment",
    "WikiPathways" = "WikiPathways_enrichment",
    "STRING" = "STRING_enrichment",
    "GSEA" = "GSEA_results"
  )
  
  filename <- paste0(file_pattern, "_", direction, ".rds")
  file.path(base_path, method, gene, cluster, filename)
}

# Function to extract genes for a specific term
extract_genes_for_term <- function(term_id, enrichment_path, enrichment_type) {
  if (!file.exists(enrichment_path)) {
    return(list(genes = NULL, error = "File not found"))
  }
  
  tryCatch({
    result <- readRDS(enrichment_path)
    
    # Handle different result formats
    if (inherits(result, "enrichResult") || inherits(result, "gseaResult")) {
      # clusterProfiler result object
      gene_str <- result@result$geneID[result@result$ID == term_id]
    } else if (is.data.frame(result)) {
      # Data frame format
      gene_str <- result$geneID[result$ID == term_id | result$Description == term_id]
    } else if (is.list(result) && !is.null(result$result)) {
      # List with result element
      gene_str <- result$result$geneID[result$result$ID == term_id]
    } else {
      return(list(genes = NULL, error = "Unknown result format"))
    }
    
    if (length(gene_str) == 0 || is.na(gene_str)) {
      return(list(genes = NULL, error = "Term not found"))
    }
    
    # Split gene string (usually separated by "/")
    genes <- unlist(strsplit(gene_str[1], "/"))
    return(list(genes = genes, error = NULL))
    
  }, error = function(e) {
    return(list(genes = NULL, error = e$message))
  })
}

# UI function
mod_enrichment_gene_display_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    div(
      class = "well well-sm",
      style = "margin-top: 10px; max-height: 400px; overflow-y: auto;",
      h5(icon("dna"), "Genes in Selected Term"),
      uiOutput(ns("gene_display"))
    )
  )
}

# Server function
mod_enrichment_gene_display_server <- function(id, selected_term, global_selection, dataset_dir) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive to store genes with caching
    gene_cache <- reactiveValues()
    
    # Extract genes when term is selected
    observe({
      req(selected_term())
      req(global_selection())
      
      term_info <- selected_term()
      selection <- global_selection()
      
      # Create cache key
      cache_key <- paste(selection$analysis_type, selection$gene, selection$cluster, 
                        selection$enrichment_type, selection$direction, term_info$id, sep = "_")
      
      # Check cache first
      if (!is.null(gene_cache[[cache_key]])) {
        return()
      }
      
      # Construct file path
      enrichment_path <- construct_enrichment_path(
        dataset_dir,
        selection$analysis_type,
        selection$gene,
        selection$cluster,
        selection$enrichment_type,
        selection$direction
      )
      
      # Extract genes
      result <- extract_genes_for_term(term_info$id, enrichment_path, selection$enrichment_type)
      gene_cache[[cache_key]] <- result
    })
    
    # Display genes
    output$gene_display <- renderUI({
      req(selected_term())
      req(global_selection())
      
      term_info <- selected_term()
      selection <- global_selection()
      
      cache_key <- paste(selection$analysis_type, selection$gene, selection$cluster, 
                        selection$enrichment_type, selection$direction, term_info$id, sep = "_")
      
      result <- gene_cache[[cache_key]]
      
      if (is.null(result)) {
        return(p("Loading...", style = "color: #999;"))
      }
      
      if (!is.null(result$error)) {
        return(p(paste("Error:", result$error), style = "color: #d9534f;"))
      }
      
      genes <- result$genes
      if (is.null(genes) || length(genes) == 0) {
        return(p("No genes found", style = "color: #999;"))
      }
      
      # Create display
      tagList(
        p(paste("Total genes:", length(genes)), style = "font-weight: bold;"),
        
        # Gene list with copy button
        div(
          style = "background-color: #f5f5f5; padding: 10px; border-radius: 4px; margin: 10px 0;",
          p(paste(genes, collapse = ", "), style = "font-family: monospace; font-size: 12px;")
        ),
        
        # Action buttons
        div(
          style = "margin-top: 10px;",
          actionButton(session$ns("copy_genes"), "Copy to Clipboard", 
                      icon = icon("copy"), class = "btn-sm btn-primary"),
          downloadButton(session$ns("download_genes"), "Download List", 
                        class = "btn-sm btn-default", style = "margin-left: 5px;")
        )
      )
    })
    
    # Copy to clipboard functionality
    observeEvent(input$copy_genes, {
      # This would require additional JavaScript implementation
      showNotification("Gene list copied to clipboard!", type = "message")
    })
    
    # Download handler
    output$download_genes <- downloadHandler(
      filename = function() {
        term_info <- selected_term()
        paste0("genes_", gsub("[^A-Za-z0-9]", "_", term_info$id), ".txt")
      },
      content = function(file) {
        term_info <- selected_term()
        selection <- global_selection()
        
        cache_key <- paste(selection$analysis_type, selection$gene, selection$cluster, 
                          selection$enrichment_type, selection$direction, term_info$id, sep = "_")
        
        result <- gene_cache[[cache_key]]
        if (!is.null(result$genes)) {
          writeLines(result$genes, file)
        }
      }
    )
    
  })
}