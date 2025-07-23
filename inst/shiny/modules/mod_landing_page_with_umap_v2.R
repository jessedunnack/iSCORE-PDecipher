# Landing Page Module with UMAP Visualization (Version 2)
# Compact layout with UMAP on left, summary boxes on right
# UMAP dataset determined by loaded data, not user selection

# Ensure required packages are loaded
if (!requireNamespace("dplyr", quietly = TRUE)) {
  stop("dplyr package is required")
}
library(dplyr)

# Source the UMAP viewer module
source("modules/mod_umap_viewer.R")

#' Landing Page UI with UMAP (Version 2)
#' 
#' @param id Module namespace
landingPageWithUmapUI <- function(id) {
  ns <- NS(id)
  
  tagList(
    # Add styling for markers section and welcome sticky note
    tags$style(HTML(paste0("
      #", ns(""), " .form-group {
        margin-bottom: 8px !important;
      }
      #", ns(""), " .control-label {
        margin-bottom: 3px !important;
        font-size: 12px !important;
      }
      #", ns("markers_table"), " .dataTables_wrapper {
        padding: 0px !important;
      }
      #", ns("markers_table"), " table {
        font-size: 11px !important;
      }
      
      /* Welcome Sticky Note Styling */
      .welcome-sticky-note {
        position: fixed;
        left: 20px;
        top: 80px;
        width: 280px;
        background: linear-gradient(135deg, #fff9c4 0%, #fff3a0 100%);
        border: 1px solid #e6d73a;
        border-radius: 0 15px 15px 0;
        box-shadow: 0 4px 8px rgba(0,0,0,0.15), 0 6px 20px rgba(0,0,0,0.1);
        padding: 20px;
        font-family: 'Kalam', 'Comic Sans MS', cursive;
        transform: rotate(-1deg);
        z-index: 1050;
        transition: all 0.5s ease-in-out;
        max-height: 85vh;
        overflow-y: auto;
      }
      
      .welcome-sticky-note::before {
        content: '';
        position: absolute;
        top: -10px;
        left: 40px;
        width: 30px;
        height: 30px;
        background: radial-gradient(circle, #ff6b6b 30%, transparent 31%);
        border-radius: 50%;
        box-shadow: 0 2px 4px rgba(0,0,0,0.2);
      }
      
      .welcome-sticky-note h4 {
        color: #d63031;
        margin-top: 0;
        margin-bottom: 15px;
        font-size: 18px;
        font-weight: bold;
        text-align: center;
        text-shadow: 1px 1px 2px rgba(0,0,0,0.1);
      }
      
      .welcome-sticky-note p {
        color: #2d3436;
        font-size: 14px;
        line-height: 1.5;
        margin-bottom: 12px;
      }
      
      .welcome-sticky-note ul {
        margin: 10px 0;
        padding-left: 20px;
      }
      
      .welcome-sticky-note li {
        color: #2d3436;
        font-size: 13px;
        margin-bottom: 8px;
        line-height: 1.4;
      }
      
      .welcome-sticky-note .btn-got-it {
        background: linear-gradient(135deg, #00b894 0%, #00a085 100%);
        color: white;
        border: none;
        border-radius: 20px;
        padding: 8px 20px;
        font-size: 14px;
        font-weight: bold;
        cursor: pointer;
        display: block;
        margin: 15px auto 0;
        transition: transform 0.2s ease;
        box-shadow: 0 2px 4px rgba(0,0,0,0.2);
      }
      
      .welcome-sticky-note .btn-got-it:hover {
        transform: scale(1.05);
        box-shadow: 0 4px 8px rgba(0,0,0,0.3);
      }
      
      .welcome-sticky-note .highlight {
        background: rgba(255, 235, 59, 0.6);
        padding: 2px 4px;
        border-radius: 3px;
        font-weight: bold;
      }
      
      /* Main content adjustment when sticky note is visible */
      .main-content-with-note {
        margin-left: 320px;
        transition: margin-left 0.5s ease-in-out;
      }
      
      .main-content-full {
        margin-left: 0;
        transition: margin-left 0.5s ease-in-out;
      }
      
      /* Hide sticky note when dismissed */
      .welcome-sticky-note.hidden {
        transform: translateX(-100%) rotate(-1deg);
        opacity: 0;
      }
    "))),
    
    # Welcome Sticky Note
    div(id = ns("welcome_sticky_note"), class = "welcome-sticky-note",
      h4(icon("lightbulb"), " Welcome to iSCORE-PDecipher!"),
      
      p("This app evaluates the ", strong("effects of gene mutations and perturbations"), " in Parkinson's disease by analyzing changes at both the ", 
        span(class = "highlight", "differential gene expression (DEG)"), " level and the ", 
        span(class = "highlight", "functional enrichment"), " level to identify convergent biological signatures across methods."),
      
      div(
        p(strong("🎯 Main Analysis Sections:")),
        tags$ul(
          tags$li(span(class = "highlight", "DE Genes"), " - Compare gene expression changes between mutations/knockdowns and controls"),
          tags$li(span(class = "highlight", "Functional Enrichment"), " - Discover affected pathways and biological processes")
        )
      ),
      
      div(
        p(strong("⚙️ Global Settings Panel (Left Sidebar):")),
        p(style = "font-size: 13px; margin-bottom: 10px;", 
          em("The sidebar is ", span(class = "highlight", "collapsible"), " - click the arrow button (◀▶) at the top-right of the sidebar to expand/collapse!")),
        tags$ul(
          tags$li("Select your ", span(class = "highlight", "gene/mutation"), " of interest"),
          tags$li("Choose a ", span(class = "highlight", "cluster"), " to focus on"),
          tags$li("Pick ", span(class = "highlight", "enrichment database"), " (GO, KEGG, etc.)"),
          tags$li("Filter by ", span(class = "highlight", "direction"), " (UP/DOWN regulated genes)")
        )
      ),
      
      p(style = "font-size: 12px; font-style: italic; margin-top: 15px;", 
        "💡 These settings affect ALL analysis tabs and update plots in real-time!"),
      
      actionButton(ns("dismiss_welcome"), "Got it! Take me to the app", 
                  class = "btn-got-it")
    ),
    
    # Main content area with dynamic margin adjustment
    div(id = ns("main_content"), class = "main-content-with-note",
      # Main content area with two columns - OPTIMIZED LAYOUT
      fluidRow(style = "min-height: 700px;", # Use explicit minimum height
      # Left column - UMAP visualization (expanded width)
      column(8,  # Increased from 7 to 8 for more width
        div(class = "box box-primary", style = "margin-top: 0;",
          div(class = "box-header with-border",
            fluidRow(
              column(9,
                h3(class = "box-title", style = "margin: 0;",
                   icon("chart-line"),
                   "Dataset UMAP Visualization")
              ),
              column(3, style = "padding-top: 8px; text-align: right;",
                selectInput(ns("pc_selection"), 
                           label = NULL,
                           choices = c("30 PCs" = "30", "50 PCs" = "50", "100 PCs" = "100"),
                           selected = "30",  # Changed default to 30 PCs
                           width = "100px")
              )
            )
          ),
          div(class = "box-body", style = "padding: 5px; text-align: center;",
            withSpinner(
              plotOutput(ns("umap_plot"), 
                        height = "720px",     # Increased to align with right column
                        width = "950px"),     # Explicit width for 1.36:1 ratio (950/720)
              type = 4,
              color = "#3c8dbc"
            )
          )
        )
      ),
      
      # Right column - Summary statistics + Cluster Markers (reduced width)
      column(4,  # Reduced from 5 to 4 to give UMAP more space
        # Summary statistics cards in a grid (top half) - Compact layout
        div(class = "value-box-compact",
          fluidRow(
            column(6,
              uiOutput(ns("total_cells_box"))
            ),
            column(6,
              uiOutput(ns("total_clusters_box"))
            )
          ),
          fluidRow(
            column(6,
              uiOutput(ns("total_results_box"))
            ),
            column(6,
              uiOutput(ns("total_genes_box"))
            )
          ),
          fluidRow(
            column(6,
              uiOutput(ns("total_experiments_box"))
            ),
            column(6,
              uiOutput(ns("enrichment_types_box"))
            )
          )
        ),
        
        # Cluster Markers section (explicit height)
        div(class = "box box-info", style = "margin-top: 20px;",
          div(class = "box-header with-border",
            h3(class = "box-title", 
               icon("dna"),
               "Cluster Marker Genes")
          ),
          div(class = "box-body", style = "padding: 10px;",
            # Controls and methodology info
            fluidRow(
              column(6,
                selectInput(ns("selected_cluster"),
                           "Select Cluster:",
                           choices = NULL,
                           width = "100%")
              ),
              column(6,
                div(style = "margin-top: 5px; font-size: 12px; color: #555; line-height: 1.4;",
                    strong("MAST test:"), " LFC≥0.5, min.pct=0.25, min.diff.pct=0.2", br(),
                    strong("Filters:"), " padj<0.05, top 25 markers, both pos/neg")
              )
            ),
            # Markers table - explicit height
            div(id = ns("markers_table_container"), 
                style = "height: 400px; overflow-y: auto; margin-top: 10px; border: 1px solid #ddd; padding: 5px;",
              p("Loading marker genes...", id = ns("markers_loading_text"), style = "color: #888; text-align: center;"),
              withSpinner(
                DT::dataTableOutput(ns("markers_table"), height = "385px"),
                type = 1,
                color = "#3c8dbc"
              )
            ),
            # Debug info
            tags$script(HTML(sprintf("
              $(document).ready(function() {
                console.log('Markers table container exists:', $('#%s').length > 0);
                console.log('Markers table output exists:', $('#%s').length > 0);
                
                // Monitor for table updates
                Shiny.addCustomMessageHandler('%s', function(message) {
                  console.log('Markers table update:', message);
                });
              });
            ", ns("markers_table_container"), ns("markers_table"), ns("markers_debug"))))
          )
        )
      )
    ),
    
    # Detailed breakdown - full width below
    fluidRow(
      # By Analysis Type
      column(4,
        div(class = "box box-primary",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Analysis Type")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("analysis_type_plot"), height = "350px"))
          )
        )
      ),
      
      # By Enrichment Type
      column(4,
        div(class = "box box-success",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Enrichment Database")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("enrichment_type_plot"), height = "350px"))
          )
        )
      ),
      
      # By Direction
      column(4,
        div(class = "box box-info",
          div(class = "box-header with-border",
            h3(class = "box-title", "Results by Direction")
          ),
          div(class = "box-body", style = "height: 400px;",
            withSpinner(plotlyOutput(ns("direction_plot"), height = "350px"))
          )
        )
      )
    ),
    
    # Detailed tables
    fluidRow(
      column(12,
        div(class = "nav-tabs-custom",
          h3("Detailed Result Counts", style = "margin-top: 20px; margin-bottom: 20px;"),
          tabsetPanel(
            id = ns("detail_tabs"),
            
            # By Gene/Mutation
            tabPanel(
              "By Gene/Mutation",
              div(style = "margin-top: 15px;",
                DT::dataTableOutput(ns("gene_table"))
              )
            ),
            
            # By Cluster
            tabPanel(
              "By Cluster",
              div(style = "margin-top: 15px;",
                DT::dataTableOutput(ns("cluster_table"))
              )
            ),
            
            # Complete Matrix
            tabPanel(
              "Complete Matrix",
              div(style = "margin-top: 15px;",
                p("This table shows the number of significant terms for each combination of parameters."),
                p("Use the filters to explore specific combinations."),
                DT::dataTableOutput(ns("matrix_table"))
              )
            ),
            
            # Top Terms
            tabPanel(
              "Top Enriched Terms",
              div(style = "margin-top: 15px;",
                fluidRow(
                  column(4,
                    selectInput(ns("top_terms_method"),
                               "Filter by Method:",
                               choices = c("All Methods" = "all",
                                         "MAST" = "MAST",
                                         "MixScale" = "MixScale"),
                               selected = "all")
                  ),
                  column(4,
                    selectInput(ns("top_terms_gene"),
                               "Filter by Gene:",
                               choices = c("All Genes" = "all"),
                               selected = "all")
                  ),
                  column(4,
                    numericInput(ns("top_terms_n"),
                                "Number of top terms:",
                                value = 20,
                                min = 10,
                                max = 100,
                                step = 10)
                  )
                ),
                DT::dataTableOutput(ns("top_terms_table"))
              )
            )
          )
        )
      )
    )
    ) # Close main content div
  )
}

#' Landing Page Server with UMAP (Version 2)
#' 
#' @param id Module namespace
#' @param data Reactive data object from app
landingPageWithUmapServer <- function(id, data, selected_dataset = NULL) {
  moduleServer(id, function(input, output, session) {
    
    # Reactive values for UMAP data and markers
    umap_data <- reactiveValues(
      sce = NULL,
      dataset_name = NULL,
      loaded = FALSE,
      markers = NULL
    )
    
    # Observer to monitor markers loading
    observe({
      if (!is.null(umap_data$markers)) {
        cat("\n[MARKERS OBSERVER] Markers have been loaded!\n")
        cat("[MARKERS OBSERVER] Total rows:", nrow(umap_data$markers), "\n")
        cat("[MARKERS OBSERVER] Columns:", paste(colnames(umap_data$markers), collapse=", "), "\n")
        cat("[MARKERS OBSERVER] Unique clusters:", paste(unique(umap_data$markers$cluster), collapse=", "), "\n\n")
      } else {
        cat("[MARKERS OBSERVER] Markers are NULL\n")
      }
    })
    
    # Welcome sticky note dismissal
    observeEvent(input$dismiss_welcome, {
      # Hide the sticky note with animation
      shinyjs::addClass(id = "welcome_sticky_note", class = "hidden")
      
      # After animation completes, adjust main content margins
      shinyjs::delay(500, {
        shinyjs::removeClass(id = "main_content", class = "main-content-with-note")
        shinyjs::addClass(id = "main_content", class = "main-content-full")
      })
      
      # Optional: Store dismissal state to prevent reappearing
      session$userData$welcome_dismissed <- TRUE
    })
    
    # Determine which UMAP dataset to load based on user selection or app data
    observe({
      req(data$data_loaded)
      
      # Use provided dataset name if available, otherwise fallback to auto-detection
      if (!is.null(selected_dataset)) {
        # Map user-friendly names to UMAP file names
        dataset_mapping <- list(
          "iSCORE-PD only" = "iSCORE_PD",
          "iSCORE-PD + CRISPRi" = "iSCORE_PD_CRISPRi", 
          "iSCORE-PD + CRISPRi + CRISPRa" = "Full_Dataset"
        )
        
        dataset_to_load <- dataset_mapping[[selected_dataset]]
        if (is.null(dataset_to_load)) {
          # If mapping fails, extract from directory name
          if (grepl("iSCORE-PD_plus_CRISPRi_and_CRISPRa", selected_dataset)) {
            dataset_to_load <- "Full_Dataset"
          } else if (grepl("iSCORE-PD_plus_CRISPRi", selected_dataset)) {
            dataset_to_load <- "iSCORE_PD_CRISPRi"
          } else {
            dataset_to_load <- "iSCORE_PD"
          }
        }
        cat("[DATASET SELECTION] Using selected dataset:", selected_dataset, "-> UMAP:", dataset_to_load, "\n")
      } else {
        # Fallback to auto-detection for backwards compatibility
        has_crispri <- any(grepl("MixScale", data$consolidated_data$method))
        has_mutations <- any(grepl("MAST", data$consolidated_data$method))
        has_crispa <- any(grepl("CRISPRa", data$consolidated_data$method))
        
        cat("[DATASET DETECTION] MAST:", has_mutations, "CRISPRi:", has_crispri, "CRISPRa:", has_crispa, "\n")
        
        # More conservative detection: prefer smaller datasets unless clearly all 3 modalities
        # Check if data actually contains all 3 modalities (not just presence)
        total_rows <- nrow(data$consolidated_data)
        mast_count <- sum(grepl("MAST", data$consolidated_data$method))
        mixscale_count <- sum(grepl("MixScale", data$consolidated_data$method))
        crispa_count <- sum(grepl("CRISPRa", data$consolidated_data$method))
        
        cat("[DATASET DETECTION] Row counts - Total:", total_rows, "MAST:", mast_count, "MixScale:", mixscale_count, "CRISPRa:", crispa_count, "\n")
        
        # Only use Full_Dataset if CRISPRa is substantial (>5% of data)
        if (has_crispa && crispa_count > (total_rows * 0.05)) {
          dataset_to_load <- "Full_Dataset"
          cat("[DATASET DETECTION] Loading Full_Dataset (substantial CRISPRa data found)\n")
        } else if (has_crispri && has_mutations) {
          dataset_to_load <- "iSCORE_PD_CRISPRi"
          cat("[DATASET DETECTION] Loading iSCORE_PD_CRISPRi (MAST + CRISPRi)\n")
        } else if (has_crispri) {
          dataset_to_load <- "iSCORE_PD_CRISPRi"
          cat("[DATASET DETECTION] Loading iSCORE_PD_CRISPRi (CRISPRi only)\n")
        } else {
          dataset_to_load <- "iSCORE_PD"
          cat("[DATASET DETECTION] Loading iSCORE_PD (MAST only)\n")
        }
      }
      
      umap_data$dataset_name <- dataset_to_load
      
      # Load default 30 PC UMAP first (changed from 100)
      load_umap_data(dataset_to_load, "30")
    })
    
    # Function to load UMAP data for specific PC count
    load_umap_data <- function(dataset_name, pc_count) {
      # Construct filename with PC suffix
      filename <- paste0(dataset_name, "_umap_data_", pc_count, "pc.rds")
      
      # Try multiple possible paths
      possible_paths <- c(
        system.file("extdata", "umap_data", filename, package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", filename),
        paste0("../../inst/extdata/umap_data/", filename)
      )
      
      # Also check for legacy filename (without pc suffix) for any PC count
      # This ensures compatibility with older UMAP files that don't have PC suffixes
      legacy_filename <- paste0(dataset_name, "_umap_data.rds")
      possible_paths <- c(possible_paths,
        system.file("extdata", "umap_data", legacy_filename, package = "iSCORE.PDecipher"),
        file.path(getwd(), "inst", "extdata", "umap_data", legacy_filename),
        paste0("../../inst/extdata/umap_data/", legacy_filename)
      )
      
      for (path in possible_paths) {
        if (file.exists(path)) {
          tryCatch({
            cat("  [UMAP DEBUG] Attempting to load from:", path, "\n")
            umap_data$sce <- readRDS(path)
            
            # Debug: Check cluster count
            if (!is.null(umap_data$sce$cluster)) {
              n_clusters <- length(unique(umap_data$sce$cluster))
              cat("  [UMAP DEBUG] Loaded data has", n_clusters, "clusters\n")
              cat("  [UMAP DEBUG] Cluster names:", paste(head(sort(unique(umap_data$sce$cluster)), 10), collapse=", "), "...\n")
            }
            
            umap_data$loaded <- TRUE
            message("Loaded UMAP (", pc_count, " PCs) for ", dataset_name, " from: ", path)
            
            # Load markers if this is the first load
            if (is.null(umap_data$markers)) {
              markers_path <- file.path(dirname(path), paste0(dataset_name, "_cluster_markers.rds"))
              cat("[MARKERS LOAD DEBUG] Looking for markers at:", markers_path, "\n")
              if (file.exists(markers_path)) {
                umap_data$markers <- readRDS(markers_path)
                message("Loaded markers for ", dataset_name)
                cat("[MARKERS LOAD DEBUG] Loaded", nrow(umap_data$markers), "total markers\n")
                cat("[MARKERS LOAD DEBUG] Unique clusters:", paste(unique(umap_data$markers$cluster), collapse=", "), "\n")
              } else {
                cat("[MARKERS LOAD DEBUG] Markers file not found!\n")
              }
            }
            
            return(TRUE)
          }, error = function(e) {
            message("Failed to load UMAP from ", path, ": ", e$message)
          })
        }
      }
      
      # If we couldn't load the requested PC count
      showNotification(
        paste0("UMAP with ", pc_count, " PCs not available. Please run generate_multipc_umaps.R"),
        type = "warning",
        duration = 5
      )
      return(FALSE)
    }
    
    # React to PC selection changes
    observeEvent(input$pc_selection, {
      req(umap_data$dataset_name)
      
      # Load new UMAP data
      if (load_umap_data(umap_data$dataset_name, input$pc_selection)) {
        # Trigger plot redraw
        umap_data$loaded <- isolate(!umap_data$loaded)
        umap_data$loaded <- isolate(!umap_data$loaded)
      }
    })
    
    # Render UMAP plot
    output$umap_plot <- renderPlot({
      # Add dependency on PC selection to trigger re-render
      req(input$pc_selection)
      
      if (!umap_data$loaded || is.null(umap_data$sce)) {
        # Placeholder when no data
        plot.new()
        text(0.5, 0.5, "UMAP data not available\nPlease run extract_umap_data.R", 
             cex = 1.2, col = "gray60")
        return()
      }
      
      # Check if dittoSeq is available
      if (!requireNamespace("dittoSeq", quietly = TRUE)) {
        plot.new()
        text(0.5, 0.5, "dittoSeq package required\nInstall with: BiocManager::install('dittoSeq')", 
             cex = 1.2, col = "red")
        return()
      }
      
      library(dittoSeq)
      
      # Create UMAP plot colored by clusters - OPTIMIZED FOR LARGER DISPLAY
      tryCatch({
        # Fix cluster ordering: ensure numeric order instead of alphabetical
        sce_copy <- umap_data$sce
        cluster_levels <- natural_sort_clusters(unique(sce_copy$seurat_clusters))
        sce_copy$seurat_clusters <- factor(sce_copy$seurat_clusters, levels = cluster_levels)
        
        p <- dittoDimPlot(
          sce_copy,
          var = "seurat_clusters",
          reduction.use = "UMAP",
          size = 0.7,  # Increased point size for better visibility
          do.label = TRUE,
          labels.size = 6,  # DOUBLED label size as requested
          legend.show = TRUE,
          main = ""  # No title to save space
        )
        
        # Maximize plot to fill the entire container space
        p <- p + 
          theme(
            # MAXIMIZE plot area - remove all margins
            plot.margin = margin(0, 0, 0, 0, "pt"),  # Zero margins
            # Larger legend text positioned close to plot
            legend.text = element_text(size = 16),
            legend.title = element_text(size = 18),
            legend.position = "right",
            legend.margin = margin(0, 0, 0, 5, "pt"),
            legend.key.size = unit(1.2, "lines"),
            # Larger axis text for the bigger plot
            axis.text = element_text(size = 14),
            axis.title = element_text(size = 16),
            # Clean backgrounds
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            # Minimal grid
            panel.grid.major = element_line(color = "grey96", size = 0.3),
            panel.grid.minor = element_blank(),
            # Remove axis ticks to save space
            axis.ticks = element_blank()
          ) +
          # Let plot fill the explicit container dimensions without restrictions
          coord_cartesian(expand = FALSE)
        
        return(p)
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, paste("Error creating UMAP:\n", e$message), 
             cex = 1.2, col = "red")
      })
    })
    
    # Update cluster choices when UMAP data is loaded
    observe({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        clusters <- unique(SummarizedExperiment::colData(umap_data$sce)$seurat_clusters)
        
        # Sort clusters numerically using standardized function
        clusters_sorted <- natural_sort_clusters(as.character(clusters))
        
        # Clean up cluster names: remove "cluster_" prefix and create clean display names
        cluster_labels <- sapply(clusters_sorted, function(x) {
          if (grepl("^cluster_", x)) {
            # Extract number after "cluster_" and format as "Cluster X"
            cluster_num <- gsub("^cluster_", "", x)
            paste("Cluster", cluster_num)
          } else {
            # If not in "cluster_X" format, just add "Cluster" prefix
            paste("Cluster", x)
          }
        })
        
        cluster_choices <- setNames(clusters_sorted, cluster_labels)
        
        updateSelectInput(session, "selected_cluster",
                         choices = cluster_choices,
                         selected = cluster_choices[1])
      }
    })
    
    # Test: Simple text output first
    output$markers_test <- renderText({
      paste("Test: Cluster", input$selected_cluster, "selected. Markers loaded:", !is.null(umap_data$markers))
    })
    
    # Render markers table
    output$markers_table <- DT::renderDataTable({
      cat("\n[MARKERS RENDER] Starting renderDataTable function\n")
      
      # Check if cluster is selected
      if (is.null(input$selected_cluster)) {
        cat("[MARKERS RENDER] No cluster selected\n")
        return(NULL)
      }
      cat("[MARKERS RENDER] Selected cluster:", input$selected_cluster, "\n")
      
      # Check if markers are loaded
      if (is.null(umap_data$markers)) {
        cat("[MARKERS RENDER] No markers data loaded\n")
        return(NULL)
      }
      cat("[MARKERS RENDER] Markers loaded with", nrow(umap_data$markers), "rows\n")
      
      # Debug logging
      cat("[MARKERS DEBUG] Available clusters in markers:", 
          paste(unique(umap_data$markers$cluster), collapse=", "), "\n")
      
      # Handle cluster name format mismatch
      # If selected_cluster is "0" but markers have "cluster_0", adjust accordingly
      cluster_to_find <- input$selected_cluster
      if (!cluster_to_find %in% umap_data$markers$cluster) {
        # Try with "cluster_" prefix
        cluster_with_prefix <- paste0("cluster_", cluster_to_find)
        if (cluster_with_prefix %in% umap_data$markers$cluster) {
          cluster_to_find <- cluster_with_prefix
          cat("[MARKERS DEBUG] Using prefixed cluster name:", cluster_to_find, "\n")
        } else {
          cat("[MARKERS DEBUG] ERROR: Cluster not found with or without prefix\n")
          return(NULL)
        }
      }
      
      # Filter markers for selected cluster
      cat("[MARKERS DEBUG] Filtering for cluster:", cluster_to_find, "\n")
      cluster_markers <- tryCatch({
        umap_data$markers %>%
          filter(cluster == cluster_to_find) %>%
          arrange(desc(avg_log2FC)) %>%
          head(25) %>%  # Fixed to top 25 markers
          select(gene, avg_log2FC, p_val_adj, pct.1, pct.2) %>%
          mutate(
            avg_log2FC = round(avg_log2FC, 3),
            p_val_adj = formatC(p_val_adj, format = "e", digits = 2),
            pct.1 = round(pct.1, 3),
            pct.2 = round(pct.2, 3)
          )
      }, error = function(e) {
        cat("[MARKERS DEBUG] ERROR in data processing:", e$message, "\n")
        return(NULL)
      })
      
      # Check if we have any markers for this cluster
      if (is.null(cluster_markers) || nrow(cluster_markers) == 0) {
        cat("[MARKERS DEBUG] No markers found for cluster:", cluster_to_find, "\n")
        return(NULL)
      }
      
      cat("[MARKERS DEBUG] Found", nrow(cluster_markers), "markers for display\n")
      cat("[MARKERS DEBUG] First 3 genes:", paste(head(cluster_markers$gene, 3), collapse=", "), "\n")
      
      # Try to create the table
      cat("[MARKERS DEBUG] Creating DT::datatable...\n")
      dt_result <- tryCatch({
        DT::datatable(
          cluster_markers,
          options = list(
          pageLength = 12,  # Show fewer rows to fit in compact space
          scrollY = "360px",  # Increased to match new table height
          scrollCollapse = TRUE,
          dom = 't',  # Only show table (no search/pagination)
          autoWidth = TRUE,  # Enable automatic width calculation
          columnDefs = list(
            list(width = '25%', targets = 0),  # Gene column (use percentage)
            list(width = '20%', targets = 1),  # Log2FC
            list(width = '20%', targets = 2),  # P-val
            list(width = '17.5%', targets = 3),  # % in cluster
            list(width = '17.5%', targets = 4),  # % in other
            list(className = 'dt-center', targets = 1:4)  # Center align numeric columns
          ),
          # Add initialization callback to fix column widths after render
          initComplete = DT::JS(
            "function(settings, json) {",
            "  setTimeout(function() {",
            "    $(this.api().table().container()).find('.dataTables_scrollBody table').css('width', '100%');",
            "    this.api().columns.adjust();",
            "  }, 100);",
            "}"
          )
        ),
        rownames = FALSE,
        colnames = c('Gene', 'Log2FC', 'P-adj', '% this clust', '% in others')  # Clear, descriptive names
      ) %>%
        DT::formatStyle(
          'avg_log2FC',
          background = DT::styleColorBar(cluster_markers$avg_log2FC, 'lightblue'),
          backgroundSize = '100% 90%',
          backgroundRepeat = 'no-repeat',
          backgroundPosition = 'center'
        )
      }, error = function(e) {
        cat("[MARKERS DEBUG] ERROR creating DT::datatable:", e$message, "\n")
        return(NULL)
      })
      
      cat("[MARKERS DEBUG] DT::datatable creation complete\n")
      return(dt_result)
    })
    
    # Cell count box
    output$total_cells_box <- renderUI({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        n_cells <- ncol(umap_data$sce)
        valueBox(
          value = format(n_cells, big.mark = ","),
          subtitle = "Total Cells",
          icon = icon("microscope"),
          color = "aqua"
        )
      } else {
        valueBox(
          value = "N/A",
          subtitle = "Total Cells",
          icon = icon("microscope"),
          color = "gray"
        )
      }
    })
    
    # Cluster count box
    output$total_clusters_box <- renderUI({
      if (umap_data$loaded && !is.null(umap_data$sce)) {
        n_clusters <- length(unique(SummarizedExperiment::colData(umap_data$sce)$seurat_clusters))
        valueBox(
          value = n_clusters,
          subtitle = "Cell Clusters",
          icon = icon("object-group"),
          color = "yellow"
        )
      } else {
        n_clusters <- length(unique(data$consolidated_data$cluster))
        valueBox(
          value = n_clusters,
          subtitle = "Cell Clusters",
          icon = icon("object-group"),
          color = "yellow"
        )
      }
    })
    
    # Results count box
    output$total_results_box <- renderUI({
      valueBox(
        value = format(nrow(data$consolidated_data), big.mark = ","),
        subtitle = "Enrichment Results",
        icon = icon("chart-bar"),
        color = "blue"
      )
    })
    
    # Genes count box
    output$total_genes_box <- renderUI({
      n_genes <- length(unique(data$consolidated_data$gene))
      valueBox(
        value = n_genes,
        subtitle = "Genes Analyzed",
        icon = icon("dna"),
        color = "green"
      )
    })
    
    # Experiments count box
    output$total_experiments_box <- renderUI({
      n_exp <- length(unique(data$consolidated_data$method))
      valueBox(
        value = n_exp,
        subtitle = "Analysis Methods",
        icon = icon("flask"),
        color = "purple"
      )
    })
    
    # Enrichment types box
    output$enrichment_types_box <- renderUI({
      n_types <- length(unique(data$consolidated_data$enrichment_type))
      valueBox(
        value = n_types,
        subtitle = "Enrichment Types",
        icon = icon("database"),
        color = "orange"
      )
    })
    
    # Dataset info panel
    output$dataset_info <- renderUI({
      if (umap_data$loaded) {
        dataset_label <- switch(umap_data$dataset_name,
          "iSCORE_PD" = "iSCORE-PD (Mutations Only)",
          "iSCORE_PD_CRISPRi" = "iSCORE-PD + CRISPRi",
          "Full_Dataset" = "Full Dataset (CRISPRi + CRISPRa)",
          umap_data$dataset_name
        )
        
        tagList(
          p(strong("Active Dataset:"), dataset_label),
          p(strong("Analysis Types:"), paste(unique(data$consolidated_data$method), collapse = ", ")),
          p(strong("Data Loaded:"), format(Sys.time(), "%Y-%m-%d %H:%M"))
        )
      } else {
        p("Loading dataset information...", style = "color: gray;")
      }
    })
    
    # Analysis type plot
    output$analysis_type_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        group_by(method) %>%
        summarise(count = n(), .groups = 'drop')
      
      plot_ly(summary_data, 
              x = ~method, 
              y = ~count, 
              type = 'bar',
              width = 0.8,
              marker = list(color = c('#374E55', '#DF8F44'))) %>%
        layout(title = NULL,
               xaxis = list(title = "", automargin = TRUE),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE,
               autosize = TRUE,
               margin = list(l = 50, r = 20, t = 20, b = 40),
               bargap = 0.3)
    })
    
    # Enrichment type plot
    output$enrichment_type_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        group_by(enrichment_type) %>%
        summarise(count = n(), .groups = 'drop') %>%
        arrange(desc(count))
      
      # Define colors for each enrichment type
      colors <- c(
        "GO_BP" = "#8dd3c7",
        "GO_CC" = "#ffffb3", 
        "GO_MF" = "#bebada",
        "KEGG" = "#fb8072",
        "Reactome" = "#80b1d3",
        "WikiPathways" = "#fdb462",
        "STRING" = "#b3de69",
        "GSEA" = "#fccde5"
      )
      
      plot_ly(summary_data,
              x = ~enrichment_type,
              y = ~count,
              type = 'bar',
              width = 0.8,
              marker = list(color = unname(colors[summary_data$enrichment_type]))) %>%
        layout(title = NULL,
               xaxis = list(title = "", tickangle = -45, automargin = TRUE),
               yaxis = list(title = "Number of Results"),
               showlegend = FALSE,
               autosize = TRUE,
               margin = list(l = 50, r = 20, t = 20, b = 80),
               bargap = 0.2)
    })
    
    # Direction plot
    output$direction_plot <- renderPlotly({
      summary_data <- data$consolidated_data %>%
        filter(direction != "RANKED") %>%
        group_by(direction) %>%
        summarise(count = n(), .groups = 'drop')
      
      colors <- c("UP" = "#DC3220", "DOWN" = "#005AB5", "ALL" = "#7C4DFF")
      
      plot_ly(summary_data,
              labels = ~direction,
              values = ~count,
              type = 'pie',
              hole = 0.3,
              marker = list(colors = colors[summary_data$direction])) %>%
        layout(title = NULL,
               showlegend = TRUE,
               autosize = TRUE,
               margin = list(l = 20, r = 20, t = 20, b = 20))
    })
    
    # Gene table
    output$gene_table <- DT::renderDataTable({
      gene_summary <- data$consolidated_data %>%
        group_by(gene, method) %>%
        summarise(
          total_terms = n(),
          enrichment_types = n_distinct(enrichment_type),
          clusters = n_distinct(cluster),
          .groups = 'drop'
        ) %>%
        arrange(desc(total_terms))
      
      DT::datatable(gene_summary,
                    options = list(pageLength = 15),
                    rownames = FALSE) %>%
        formatStyle('total_terms',
                   background = styleColorBar(gene_summary$total_terms, 'lightblue'),
                   backgroundSize = '100% 90%',
                   backgroundRepeat = 'no-repeat',
                   backgroundPosition = 'center')
    })
    
    # Cluster table
    output$cluster_table <- DT::renderDataTable({
      cluster_summary <- data$consolidated_data %>%
        group_by(cluster) %>%
        summarise(
          total_terms = n(),
          genes = n_distinct(gene),
          methods = paste(unique(method), collapse = ", "),
          .groups = 'drop'
        ) %>%
        arrange(cluster)
      
      DT::datatable(cluster_summary,
                    options = list(pageLength = 15),
                    rownames = FALSE)
    })
    
    # Complete matrix table
    output$matrix_table <- DT::renderDataTable({
      matrix_data <- data$consolidated_data %>%
        group_by(gene, cluster, method, enrichment_type, direction) %>%
        summarise(n_terms = n(), .groups = 'drop') %>%
        arrange(gene, cluster, method)
      
      DT::datatable(matrix_data,
                    filter = 'top',
                    options = list(
                      pageLength = 20,
                      scrollX = TRUE
                    ),
                    rownames = FALSE)
    })
    
    # Update gene choices for top terms
    observe({
      gene_choices <- c("All Genes" = "all")
      unique_genes <- sort(unique(data$consolidated_data$gene))
      gene_choices <- c(gene_choices, setNames(unique_genes, unique_genes))
      
      updateSelectInput(session, "top_terms_gene",
                       choices = gene_choices)
    })
    
    # Top terms table
    output$top_terms_table <- DT::renderDataTable({
      top_data <- data$consolidated_data
      
      # Apply filters
      if (input$top_terms_method != "all") {
        top_data <- top_data %>%
          filter(method == input$top_terms_method)
      }
      
      if (!is.null(input$top_terms_gene) && input$top_terms_gene != "all") {
        top_data <- top_data %>%
          filter(gene == input$top_terms_gene)
      }
      
      # Get top terms by p.adjust
      top_terms <- top_data %>%
        arrange(p.adjust) %>%
        head(input$top_terms_n) %>%
        select(gene, cluster, method, enrichment_type, direction, 
               Description, p.adjust, Count) %>%
        mutate(p.adjust = format(p.adjust, scientific = TRUE, digits = 3))
      
      DT::datatable(top_terms,
                    options = list(
                      pageLength = 20,
                      scrollX = TRUE
                    ),
                    rownames = FALSE)
    })
    
  })
}

# Helper function for value boxes (if not already defined)
valueBox <- function(value, subtitle, icon = NULL, color = "blue", width = 12) {
  div(class = paste0("small-box bg-", color), style = "margin-bottom: 15px;",
    div(class = "inner",
      h3(value, style = "margin-bottom: 5px;"),
      p(subtitle, style = "margin: 0;")
    ),
    if (!is.null(icon)) {
      div(class = "icon",
        icon
      )
    }
  )
}