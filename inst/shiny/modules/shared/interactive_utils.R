# Interactive Utilities
# Cross-module communication and shared reactive features

#' Create shared reactive values for cross-module communication
#' 
#' @return ReactiveValues object with shared selections
create_shared_selections <- function() {
  reactiveValues(
    selected_genes = NULL,
    selected_pathways = NULL,
    current_comparison = NULL,
    highlighted_conditions = NULL,
    last_updated_module = NULL,
    update_timestamp = Sys.time()
  )
}

#' Update shared gene selection
#' 
#' @param shared_selections ReactiveValues object
#' @param genes Vector of selected genes
#' @param module_name Name of module making the update
update_shared_genes <- function(shared_selections, genes, module_name) {
  shared_selections$selected_genes <- genes
  shared_selections$last_updated_module <- module_name
  shared_selections$update_timestamp <- Sys.time()
}

#' Update shared pathway selection
#' 
#' @param shared_selections ReactiveValues object  
#' @param pathways Vector of selected pathways
#' @param module_name Name of module making the update
update_shared_pathways <- function(shared_selections, pathways, module_name) {
  shared_selections$selected_pathways <- pathways
  shared_selections$last_updated_module <- module_name
  shared_selections$update_timestamp <- Sys.time()
}

#' Get genes associated with a pathway
#' 
#' @param pathway_name Name of the pathway
#' @param enrichment_data Enrichment data frame
#' @return Vector of gene names
get_pathway_genes <- function(pathway_name, enrichment_data) {
  
  if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
    return(character())
  }
  
  # Find the pathway in the data
  pathway_rows <- enrichment_data[enrichment_data$Description == pathway_name, ]
  
  if (nrow(pathway_rows) == 0) {
    return(character())
  }
  
  # Extract genes from geneID column if available
  if ("geneID" %in% names(pathway_rows)) {
    gene_string <- pathway_rows$geneID[1]
    if (!is.na(gene_string) && gene_string != "") {
      # Split by "/" separator commonly used in enrichment results
      genes <- strsplit(gene_string, "/")[[1]]
      return(trimws(genes))
    }
  }
  
  # Fallback: try other common gene columns
  gene_columns <- c("core_enrichment", "leading_edge", "genes")
  for (col in gene_columns) {
    if (col %in% names(pathway_rows)) {
      gene_string <- pathway_rows[[col]][1]
      if (!is.na(gene_string) && gene_string != "") {
        genes <- strsplit(gene_string, "[/,;]")[[1]]
        return(trimws(genes))
      }
    }
  }
  
  return(character())
}

#' Get pathways associated with a gene
#' 
#' @param gene_name Name of the gene
#' @param enrichment_data Enrichment data frame
#' @return Vector of pathway names
get_gene_pathways <- function(gene_name, enrichment_data) {
  
  if (is.null(enrichment_data) || nrow(enrichment_data) == 0) {
    return(character())
  }
  
  pathways <- character()
  
  # Search in geneID column
  if ("geneID" %in% names(enrichment_data)) {
    matching_rows <- enrichment_data[grepl(gene_name, enrichment_data$geneID, fixed = TRUE), ]
    pathways <- c(pathways, matching_rows$Description)
  }
  
  # Search in other gene columns
  gene_columns <- c("core_enrichment", "leading_edge", "genes")
  for (col in gene_columns) {
    if (col %in% names(enrichment_data)) {
      matching_rows <- enrichment_data[grepl(gene_name, enrichment_data[[col]], fixed = TRUE), ]
      pathways <- c(pathways, matching_rows$Description)
    }
  }
  
  return(unique(pathways))
}

#' Create click event handler for gene-pathway cross-referencing
#' 
#' @param input Shiny input object
#' @param output Shiny output object
#' @param session Shiny session object
#' @param shared_selections Shared reactive values
#' @param enrichment_data Enrichment data reactive
#' @param event_source Click event source (plotly_click, etc.)
create_gene_pathway_click_handler <- function(input, output, session, shared_selections, 
                                             enrichment_data, event_source = "plotly_click") {
  
  observeEvent(input[[event_source]], {
    click_data <- input[[event_source]]
    
    if (is.null(click_data)) return()
    
    # Extract clicked item
    clicked_item <- NULL
    
    if ("y" %in% names(click_data)) {
      # Likely a pathway name from a dot plot or bar chart
      clicked_item <- click_data$y
      
      # Find associated genes
      associated_genes <- get_pathway_genes(clicked_item, enrichment_data())
      
      if (length(associated_genes) > 0) {
        update_shared_genes(shared_selections, associated_genes, session$ns.prefix)
        
        showNotification(
          paste("Selected pathway:", clicked_item, "with", length(associated_genes), "genes"),
          type = "message",
          duration = 3
        )
      }
      
    } else if ("x" %in% names(click_data)) {
      # Might be a gene name
      clicked_item <- click_data$x
      
      # Find associated pathways
      associated_pathways <- get_gene_pathways(clicked_item, enrichment_data())
      
      if (length(associated_pathways) > 0) {
        update_shared_pathways(shared_selections, associated_pathways, session$ns.prefix)
        
        showNotification(
          paste("Selected gene:", clicked_item, "found in", length(associated_pathways), "pathways"),
          type = "message",
          duration = 3
        )
      }
    }
  })
}

#' Create hover information for enhanced interactivity
#' 
#' @param data Data frame for the plot
#' @param x_col X-axis column name
#' @param y_col Y-axis column name
#' @param additional_info Additional columns to include in hover
#' @return Character vector of hover text
create_enhanced_hover_text <- function(data, x_col, y_col, additional_info = c()) {
  
  if (nrow(data) == 0) {
    return(character())
  }
  
  # Base hover text
  hover_text <- paste0(
    y_col, ": ", data[[y_col]], "<br>",
    x_col, ": ", data[[x_col]]
  )
  
  # Add additional information
  for (info_col in additional_info) {
    if (info_col %in% names(data)) {
      hover_text <- paste0(
        hover_text, "<br>",
        info_col, ": ", data[[info_col]]
      )
    }
  }
  
  return(hover_text)
}

#' Create brush selection handler for multi-gene selection
#' 
#' @param input Shiny input object
#' @param shared_selections Shared reactive values
#' @param data_reactive Reactive data frame
#' @param y_column Column containing gene/pathway names
create_brush_selection_handler <- function(input, shared_selections, data_reactive, 
                                         y_column = "Description", event_source = "plot_brush") {
  
  observeEvent(input[[event_source]], {
    brush_data <- input[[event_source]]
    
    if (is.null(brush_data)) return()
    
    # Get current data
    data <- data_reactive()
    if (is.null(data) || nrow(data) == 0) return()
    
    # Find items within brush selection
    if ("y" %in% names(brush_data)) {
      y_min <- brush_data$ymin
      y_max <- brush_data$ymax
      
      # For discrete y-axis (categorical), get items in range
      if (y_column %in% names(data)) {
        # Convert factor levels to numeric if needed
        if (is.factor(data[[y_column]])) {
          y_levels <- levels(data[[y_column]])
          y_indices <- as.numeric(data[[y_column]])
          
          # Find items in brush range
          selected_indices <- which(y_indices >= y_min & y_indices <= y_max)
          selected_items <- data[[y_column]][selected_indices]
          
        } else {
          # For continuous y-axis
          selected_items <- data[[y_column]][data[[y_column]] >= y_min & data[[y_column]] <= y_max]
        }
        
        if (length(selected_items) > 0) {
          # Determine if these are genes or pathways
          if (grepl("gene|ENSG|^[A-Z0-9]+$", selected_items[1])) {
            update_shared_genes(shared_selections, unique(selected_items), "brush_selection")
          } else {
            update_shared_pathways(shared_selections, unique(selected_items), "brush_selection")
          }
          
          showNotification(
            paste("Selected", length(unique(selected_items)), "items via brush"),
            type = "message",
            duration = 2
          )
        }
      }
    }
  })
}

#' Apply highlighting based on shared selections
#' 
#' @param plot_data Data frame for plotting
#' @param shared_selections Shared reactive values
#' @param highlight_column Column to check for highlighting
#' @return Data frame with highlight_flag column added
apply_shared_highlighting <- function(plot_data, shared_selections, highlight_column) {
  
  if (is.null(plot_data) || nrow(plot_data) == 0) {
    return(plot_data)
  }
  
  # Initialize highlight flag
  plot_data$highlight_flag <- FALSE
  
  # Check for shared gene selections
  if (!is.null(shared_selections$selected_genes) && 
      length(shared_selections$selected_genes) > 0 &&
      highlight_column %in% names(plot_data)) {
    
    plot_data$highlight_flag <- plot_data[[highlight_column]] %in% shared_selections$selected_genes
  }
  
  # Check for shared pathway selections
  if (!is.null(shared_selections$selected_pathways) && 
      length(shared_selections$selected_pathways) > 0 &&
      highlight_column %in% names(plot_data)) {
    
    plot_data$highlight_flag <- plot_data$highlight_flag | 
                               plot_data[[highlight_column]] %in% shared_selections$selected_pathways
  }
  
  return(plot_data)
}

#' Create summary statistics from metadata
#' 
#' @param condition Condition identifier (gene_cluster_source)
#' @param metadata Cell metadata data frame
#' @return List with summary statistics
generate_condition_summary <- function(condition, metadata = NULL) {
  
  # Parse condition
  parts <- strsplit(condition, "_")[[1]]
  if (length(parts) < 3) {
    return(list(
      condition = condition,
      cell_count = NA,
      mean_nCount_RNA = NA,
      mean_percent_mt = NA,
      error = "Invalid condition format"
    ))
  }
  
  gene <- parts[1]
  cluster <- paste(parts[2:(length(parts)-1)], collapse = "_")
  source <- parts[length(parts)]
  
  # If no metadata provided, return basic structure
  if (is.null(metadata)) {
    return(list(
      condition = condition,
      gene = gene,
      cluster = cluster,
      source = source,
      cell_count = NA,
      mean_nCount_RNA = NA,
      mean_percent_mt = NA,
      note = "No metadata provided"
    ))
  }
  
  # Filter metadata for this condition
  condition_cells <- metadata
  
  # Apply filters based on available columns
  if ("gene" %in% names(metadata)) {
    condition_cells <- condition_cells[condition_cells$gene == gene, ]
  }
  
  if ("seurat_clusters" %in% names(metadata)) {
    condition_cells <- condition_cells[condition_cells$seurat_clusters == cluster, ]
  }
  
  if ("dataset" %in% names(metadata)) {
    if (source == "iSCORE-PD") {
      condition_cells <- condition_cells[condition_cells$dataset == "iSCORE-PD", ]
    } else {
      condition_cells <- condition_cells[condition_cells$dataset == "PerturbSeq", ]
    }
  }
  
  # Calculate summary statistics
  summary_stats <- list(
    condition = condition,
    gene = gene,
    cluster = cluster,
    source = source,
    cell_count = nrow(condition_cells)
  )
  
  # Add quality metrics if available
  if ("nCount_RNA" %in% names(condition_cells)) {
    summary_stats$mean_nCount_RNA <- mean(condition_cells$nCount_RNA, na.rm = TRUE)
  }
  
  if ("nFeature_RNA" %in% names(condition_cells)) {
    summary_stats$mean_nFeature_RNA <- mean(condition_cells$nFeature_RNA, na.rm = TRUE)
  }
  
  if ("percent.mt" %in% names(condition_cells)) {
    summary_stats$mean_percent_mt <- mean(condition_cells$percent.mt, na.rm = TRUE)
  }
  
  return(summary_stats)
}