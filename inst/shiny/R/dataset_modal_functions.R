# Dataset Modal Functions for PerturbSeq Analysis App
# Helper functions for dataset selection modals

#' Create dataset selector modal
#' @param available_datasets List of available datasets
#' @param current_dataset Currently selected dataset
#' @return Modal dialog UI
create_dataset_selector_modal <- function(available_datasets, current_dataset = NULL) {
  modalDialog(
    title = "Select Dataset",
    size = "l",
    
    div(style = "margin-bottom: 20px;",
      h4("Available Datasets"),
      p("Choose a dataset to analyze. Each dataset contains differential expression and enrichment analysis results.")
    ),
    
    # Dataset cards
    div(class = "row",
      lapply(names(available_datasets), function(dataset_key) {
        dataset <- available_datasets[[dataset_key]]
        is_selected <- !is.null(current_dataset) && dataset_key == current_dataset
        
        div(class = "col-md-6", style = "margin-bottom: 15px;",
          div(class = paste0("card", if(is_selected) " border-primary" else ""),
              style = "cursor: pointer;",
              onclick = paste0("Shiny.setInputValue('selected_dataset', '", dataset_key, "', {priority: 'event'})"),
              
              div(class = "card-body",
                h5(class = "card-title", 
                   dataset$name,
                   if(is_selected) tags$span(class = "badge badge-primary float-right", "Current")
                ),
                
                # Dataset description
                p(class = "card-text", dataset$description),
                
                # Dataset stats
                if (!is.null(dataset$stats)) {
                  tags$ul(class = "list-unstyled",
                    if (!is.null(dataset$stats$n_cells)) {
                      tags$li(icon("cell"), sprintf(" %s cells", format(dataset$stats$n_cells, big.mark = ",")))
                    },
                    if (!is.null(dataset$stats$n_genes)) {
                      tags$li(icon("dna"), sprintf(" %s genes", format(dataset$stats$n_genes, big.mark = ",")))
                    },
                    if (!is.null(dataset$stats$n_clusters)) {
                      tags$li(icon("project-diagram"), sprintf(" %d clusters", dataset$stats$n_clusters))
                    },
                    if (!is.null(dataset$stats$has_crispr)) {
                      tags$li(icon("cut"), " CRISPRi/CRISPRa data")
                    }
                  )
                },
                
                # Data availability indicators
                div(style = "margin-top: 10px;",
                  if (!is.null(dataset$has_de) && dataset$has_de) {
                    tags$span(class = "badge badge-success", "DE Results")
                  },
                  " ",
                  if (!is.null(dataset$has_enrichment) && dataset$has_enrichment) {
                    tags$span(class = "badge badge-info", "Enrichment")
                  },
                  " ",
                  if (!is.null(dataset$has_seurat) && dataset$has_seurat) {
                    tags$span(class = "badge badge-secondary", "Seurat Object")
                  }
                )
              )
          )
        )
      })
    ),
    
    # Custom upload option
    hr(),
    div(style = "text-align: center;",
      p("Have a custom dataset?"),
      actionButton("show_upload_option", 
                  "Upload Custom Dataset", 
                  class = "btn-outline-secondary",
                  icon = icon("upload"))
    ),
    
    footer = tagList(
      modalButton("Cancel"),
      actionButton("load_dataset", 
                  "Load Selected Dataset", 
                  class = "btn-primary",
                  icon = icon("database"))
    )
  )
}

#' Create custom upload modal
#' @return Modal dialog UI for custom uploads
create_custom_upload_modal <- function() {
  modalDialog(
    title = "Upload Custom Dataset",
    size = "l",
    
    div(style = "margin-bottom: 20px;",
      h4("Upload Your Data"),
      p("Upload preprocessed data files for analysis. The app expects specific file formats.")
    ),
    
    # File upload sections
    wellPanel(
      h5("Required Files"),
      
      # DE results file
      div(style = "margin-bottom: 15px;",
        fileInput("custom_de_file",
                 label = "Differential Expression Results (RDS)",
                 accept = c(".rds", ".RDS"),
                 placeholder = "full_DE_results.rds"),
        helpText("RDS file containing DE analysis results from MAST and/or MixScale")
      ),
      
      # Enrichment file
      div(style = "margin-bottom: 15px;",
        fileInput("custom_enrichment_file",
                 label = "Enrichment Results (RDS)",
                 accept = c(".rds", ".RDS"),
                 placeholder = "all_enrichment_padj005_complete_with_direction.rds"),
        helpText("RDS file containing functional enrichment analysis results")
      )
    ),
    
    wellPanel(
      h5("Optional Files"),
      
      # Seurat object
      div(style = "margin-bottom: 15px;",
        fileInput("custom_seurat_file",
                 label = "Seurat Object (RDS)",
                 accept = c(".rds", ".RDS"),
                 placeholder = "seurat_object.rds"),
        helpText("Seurat object with UMAP coordinates and clustering")
      ),
      
      # Metadata
      div(style = "margin-bottom: 15px;",
        fileInput("custom_metadata_file",
                 label = "Cell Metadata (RDS/CSV)",
                 accept = c(".rds", ".RDS", ".csv", ".tsv"),
                 placeholder = "metadata.rds"),
        helpText("Cell metadata including clusters and conditions")
      )
    ),
    
    # Validation messages
    div(id = "upload_validation_message"),
    
    footer = tagList(
      actionButton("back_to_dataset_selector", 
                  "Back", 
                  class = "btn-secondary",
                  icon = icon("arrow-left")),
      modalButton("Cancel"),
      actionButton("load_custom_dataset", 
                  "Load Custom Dataset", 
                  class = "btn-primary",
                  icon = icon("check"),
                  disabled = TRUE)
    )
  )
}

#' Validate uploaded dataset
#' @param file_path Path to uploaded file
#' @return List with validation result
validate_dataset <- function(file_path) {
  result <- list(valid = FALSE, message = "")
  
  tryCatch({
    # Check file exists
    if (!file.exists(file_path)) {
      result$message <- "File not found"
      return(result)
    }
    
    # Check file size
    file_size <- file.info(file_path)$size / 1024^2  # MB
    if (file_size > 500) {
      result$message <- sprintf("File too large (%.1f MB). Maximum size is 500 MB.", file_size)
      return(result)
    }
    
    # Try to load the file
    data <- readRDS(file_path)
    
    # Basic validation based on data structure
    if (is.data.frame(data) || is.list(data)) {
      result$valid <- TRUE
      result$message <- "File validated successfully"
      
      # Additional checks for specific data types
      if ("mutation_perturbation" %in% names(data) || "gene" %in% names(data)) {
        result$data_type <- "enrichment"
      } else if ("avg_log2FC" %in% names(data) || "p_val_adj" %in% names(data)) {
        result$data_type <- "de_results"
      }
      
    } else {
      result$message <- "Invalid data format. Expected data frame or list."
    }
    
  }, error = function(e) {
    result$message <- paste("Error reading file:", e$message)
  })
  
  return(result)
}

#' Get default datasets configuration
#' @return List of available datasets with metadata
get_default_datasets <- function() {
  base_dir <- "/mnt/e/ASAP/scRNASeq/PerturbSeq/final/"
  
  datasets <- list(
    iscore_pd = list(
      name = "iSCORE-PD (Mutations Only)",
      description = "iSCORE-PD dataset with 13 PD mutations analyzed across ~210K cells",
      path = file.path(base_dir, "iSCORE-PD"),
      stats = list(
        n_cells = 210000,
        n_genes = 36601,
        n_clusters = 16,
        has_crispr = FALSE
      ),
      has_de = file.exists(file.path(base_dir, "iSCORE-PD", "full_DE_results.rds")),
      has_enrichment = file.exists(file.path(base_dir, "iSCORE-PD", "all_enrichment_padj005_complete_with_direction.rds")),
      has_seurat = file.exists(file.path(base_dir, "iSCORE-PD", "iSCORE-PD_final.rds"))
    ),
    
    iscore_pd_plus_crispri = list(
      name = "iSCORE-PD + CRISPRi",
      description = "Combined dataset with mutations and CRISPRi perturbations (~230K cells)",
      path = file.path(base_dir, "iSCORE-PD_plus_CRISPRi"),
      stats = list(
        n_cells = 231874,
        n_genes = 36601,
        n_clusters = 16,
        has_crispr = TRUE
      ),
      has_de = file.exists(file.path(base_dir, "iSCORE-PD_plus_CRISPRi", "full_DE_results.rds")),
      has_enrichment = file.exists(file.path(base_dir, "iSCORE-PD_plus_CRISPRi", "all_enrichment_padj005_complete_with_direction.rds")),
      has_seurat = file.exists(file.path(base_dir, "iSCORE-PD_plus_CRISPRi", "iSCORE-PD_plus_CRISPRi_final.rds"))
    )
  )
  
  return(datasets)
}

cat("Dataset modal functions loaded successfully\n")