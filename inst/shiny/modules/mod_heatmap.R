# Module: Advanced Heatmaps
# Provides heatmap visualizations using the unified heatmap system

mod_heatmap_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Heatmap Settings"),
      
      selectInput(ns("heatmap_type"),
                  "Heatmap Type:",
                  choices = c("P-value" = "pvalue",
                              "Fold Enrichment" = "foldenrichment", 
                              "Z-score" = "zscore",
                              "Gene Count" = "count",
                              "GSEA NES" = "nes"),
                  selected = "pvalue"),
      
      selectInput(ns("method_filter"),
                  "Analysis Method:",
                  choices = c("All" = "all",
                              "MAST Only" = "MAST",
                              "MixScale Only" = "MixScale", 
                              "Intersection" = "intersection",
                              "Union" = "union"),
                  selected = "all"),
      
      checkboxGroupInput(ns("enrichment_types"),
                         "Enrichment Types:",
                         choices = c("GO_BP", "GO_CC", "GO_MF", "KEGG", "Reactome", "STRING", "GSEA"),
                         selected = c("GO_BP", "KEGG")),
      
      # GSEA-specific options
      conditionalPanel(
        condition = paste0("input['", ns("heatmap_type"), "'] == 'nes' || input['", ns("enrichment_types"), "'].includes('GSEA')"),
        wellPanel(
          style = "background-color: #f0f0f0; margin-top: 10px;",
          h5("GSEA Options"),
          checkboxInput(ns("gsea_only"), "Show GSEA results only", value = FALSE),
          sliderInput(ns("nes_threshold"), 
                      "Min |NES| threshold:",
                      min = 0, max = 3, value = 1, step = 0.1)
        )
      ),
      
      selectInput(ns("cluster_select"),
                  "Cluster:",
                  choices = c("All Clusters" = "all", paste0("cluster_", 0:9)),
                  selected = "all"),
      
      sliderInput(ns("max_terms"),
                  "Maximum Terms:",
                  min = 10,
                  max = 100,
                  value = 50,
                  step = 10),
      
      selectInput(ns("direction_filter"),
                  "Gene Direction:",
                  choices = c("All Directions" = "ALL_DIR",
                              "Up-regulated only" = "UP",
                              "Down-regulated only" = "DOWN",
                              "Both (non-ALL)" = "BOTH"),
                  selected = "ALL_DIR"),
      
      selectInput(ns("color_scale"),
                  "Color Scale:",
                  choices = c("Red (0 to max)" = "red",
                              "Blue (0 to max)" = "blue",
                              "Red-Blue (diverging)" = "RdBu",
                              "Viridis" = "viridis",
                              "Yellow-Orange-Red" = "YlOrRd"),
                  selected = "red"),
      
      selectInput(ns("scale_method"),
                  "Color Scaling Method:",
                  choices = c("Linear" = "linear",
                              "Quantile (adaptive)" = "quantile",
                              "Fixed breaks" = "fixed"),
                  selected = "linear"),
      
      checkboxInput(ns("show_annotations"),
                    "Show Row Annotations",
                    value = TRUE),
      
      checkboxInput(ns("cluster_rows"),
                    "Cluster Rows",
                    value = TRUE),
      
      checkboxInput(ns("cluster_cols"),
                    "Cluster Columns", 
                    value = TRUE),
      
      br(),
      actionButton(ns("generate_heatmap"),
                   "Generate Heatmap",
                   class = "btn-primary", 
                   width = "100%"),
      
      br(),
      br(),
      downloadButton(ns("download_heatmap"),
                     "Download Heatmap",
                     class = "btn-success",
                     style = "width: 100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        tabPanel("Heatmap",
                 br(),
                 withSpinner(
                   plotlyOutput(ns("heatmap_plot"), height = "800px", width = "100%"),
                   type = 4,
                   color = "#3c8dbc"
                 )),
        
        tabPanel("Data View",
                 br(),
                 DT::dataTableOutput(ns("heatmap_data"))),
        
        tabPanel("PDF Export",
                 br(),
                 h4("Generate Publication-Quality PDF"),
                 p("This will create a static PDF version using ComplexHeatmap with the same settings as the interactive heatmap."),
                 br(),
                 fluidRow(
                   column(6,
                     numericInput(ns("pdf_width"), "Width (inches):", value = 12, min = 6, max = 24),
                     numericInput(ns("pdf_height"), "Height (inches):", value = 10, min = 6, max = 36)
                   ),
                   column(6,
                     numericInput(ns("pdf_fontsize"), "Base Font Size:", value = 10, min = 6, max = 16),
                     numericInput(ns("pdf_dpi"), "DPI:", value = 300, min = 150, max = 600)
                   )
                 ),
                 br(),
                 actionButton(ns("generate_pdf"), "Generate PDF", class = "btn-success"),
                 br(),
                 br(),
                 uiOutput(ns("pdf_status"))
        ),
        
        tabPanel("Settings Info",
                 br(),
                 verbatimTextOutput(ns("settings_info")))
      )
    )
  )
}

mod_heatmap_server <- function(id, app_data, pval_threshold) {
  moduleServer(id, function(input, output, session) {
    
    # Helper function to validate and create color vectors
    create_validated_colors <- function() {
      # Define complete color vectors with validation
      type_colors <- c(
        "GO_BP" = "#8DD3C7", 
        "GO_CC" = "#FFFFB3", 
        "GO_MF" = "#BEBADA",
        "KEGG" = "#FB8072", 
        "Reactome" = "#80B1D3", 
        "WikiPathways" = "#FDB462", 
        "STRING" = "#B3DE69",
        "GSEA" = "#FCCDE5"
      )
      
      dir_colors <- c(
        "UP" = "#FF6B6B", 
        "DOWN" = "#4ECDC4", 
        "ALL" = "#95E1D3", 
        "RANKED" = "#45B7D1"
      )
      
      # Validate that all colors are proper hex codes
      validate_hex <- function(color_vec) {
        valid_colors <- sapply(color_vec, function(x) {
          grepl("^#[0-9A-Fa-f]{6}$", x)
        })
        if (!all(valid_colors)) {
          warning("Invalid hex colors found, using defaults")
          # Replace invalid colors with gray
          color_vec[!valid_colors] <- "#808080"
        }
        return(color_vec)
      }
      
      return(list(
        type = validate_hex(type_colors),
        direction = validate_hex(dir_colors)
      ))
    }
    
    # Helper function to create safe row annotations
    create_safe_row_annotations <- function(mat, annotation_data, show_annotations) {
      if (!show_annotations || is.null(annotation_data) || nrow(annotation_data) == 0) {
        return(list(colors = NULL, palette = NULL))
      }
      
      # Get validated color vectors
      color_vectors <- create_validated_colors()
      
      # Ensure annotation data matches matrix rows
      if (nrow(annotation_data) != nrow(mat)) {
        message("Row annotation data doesn't match matrix dimensions")
        return(list(colors = NULL, palette = NULL))
      }
      
      # Check required columns exist
      if (!all(c("enrichment_type", "direction") %in% names(annotation_data))) {
        message("Required annotation columns missing")
        return(list(colors = NULL, palette = NULL))
      }
      
      # Create annotation data frame with proper row names
      row_colors <- data.frame(
        Type = as.character(annotation_data$enrichment_type),
        Direction = as.character(annotation_data$direction),
        stringsAsFactors = FALSE
      )
      rownames(row_colors) <- rownames(mat)
      
      # Ensure all types and directions have colors
      missing_types <- setdiff(unique(row_colors$Type), names(color_vectors$type))
      missing_dirs <- setdiff(unique(row_colors$Direction), names(color_vectors$direction))
      
      if (length(missing_types) > 0) {
        message("Adding default colors for missing types: ", paste(missing_types, collapse = ", "))
        new_colors <- rainbow(length(missing_types))
        names(new_colors) <- missing_types
        color_vectors$type <- c(color_vectors$type, new_colors)
      }
      
      if (length(missing_dirs) > 0) {
        message("Adding default colors for missing directions: ", paste(missing_dirs, collapse = ", "))
        new_colors <- rainbow(length(missing_dirs))
        names(new_colors) <- missing_dirs
        color_vectors$direction <- c(color_vectors$direction, new_colors)
      }
      
      # For heatmaply, ensure palette is a simple named vector, not a list
      # Extract just the colors that are actually used
      used_types <- unique(row_colors$Type)
      used_dirs <- unique(row_colors$Direction)
      
      type_palette <- color_vectors$type[used_types]
      dir_palette <- color_vectors$direction[used_dirs]
      
      # Ensure no NA values
      type_palette <- type_palette[!is.na(names(type_palette))]
      dir_palette <- dir_palette[!is.na(names(dir_palette))]
      
      return(list(
        colors = row_colors,
        palette = list(
          Type = type_palette,
          Direction = dir_palette
        )
      ))
    }
    
    # Reactive values
    heatmap_data <- reactiveValues(
      data = NULL,
      plot = NULL
    )
    
    # Load and prepare data for heatmap
    observeEvent(input$generate_heatmap, {
      req(input$enrichment_types)
      
      showNotification("Loading data for heatmap...", type = "default", duration = NULL, id = "loading")
      
      tryCatch({
        # Use the consolidated data from app_data
        if (!is.null(app_data$consolidated_data) && nrow(app_data$consolidated_data) > 0) {
          # Use the data already loaded in the app
          all_data <- app_data$consolidated_data
          message("Using consolidated data with ", nrow(all_data), " rows")
          message("Columns available: ", paste(names(all_data), collapse = ", "))
          
          # Filter by enrichment types
          filtered_data <- all_data %>%
            dplyr::filter(enrichment_type %in% input$enrichment_types)
          message("After enrichment type filter: ", nrow(filtered_data), " rows")
          
          # Filter by cluster if not "all"
          if (input$cluster_select != "all") {
            filtered_data <- filtered_data %>%
              dplyr::filter(cluster == input$cluster_select)
            message("After cluster filter: ", nrow(filtered_data), " rows")
          }
          
          # Filter by method
          if (input$method_filter != "all") {
            if (input$method_filter %in% c("MAST", "MixScale")) {
              # Check which column name to use
              if ("method" %in% names(filtered_data)) {
                filtered_data <- filtered_data %>%
                  dplyr::filter(method == input$method_filter)
              } else if ("analysis_type" %in% names(filtered_data)) {
                filtered_data <- filtered_data %>%
                  dplyr::filter(analysis_type == input$method_filter)
              }
              message("After method filter: ", nrow(filtered_data), " rows")
            }
          }
          
          # GSEA-specific filtering
          if (input$gsea_only) {
            filtered_data <- filtered_data %>%
              dplyr::filter(enrichment_type == "GSEA")
            message("After GSEA-only filter: ", nrow(filtered_data), " rows")
          }
          
          # Filter by NES threshold if GSEA data
          if (input$heatmap_type == "nes" && "NES" %in% names(filtered_data)) {
            filtered_data <- filtered_data %>%
              dplyr::filter(abs(NES) >= input$nes_threshold)
            message("After NES threshold filter: ", nrow(filtered_data), " rows")
          }
          
          # Filter by direction
          if (input$direction_filter != "ALL_DIR") {
            if (input$direction_filter == "BOTH") {
              # Filter out "ALL" direction, keep only UP and DOWN
              filtered_data <- filtered_data %>%
                dplyr::filter(direction %in% c("UP", "DOWN"))
            } else {
              # Filter for specific direction
              filtered_data <- filtered_data %>%
                dplyr::filter(direction == input$direction_filter)
            }
            message("After direction filter: ", nrow(filtered_data), " rows")
          }
          
          # Filter by significance
          if ("p.adjust" %in% names(filtered_data)) {
            filtered_data <- filtered_data %>%
              dplyr::filter(p.adjust <= 0.05)
            message("After significance filter: ", nrow(filtered_data), " rows")
          }
          
          if (nrow(filtered_data) == 0) {
            removeNotification("loading")
            showNotification("No data matches the selected criteria", type = "warning", duration = 5)
            return()
          }
          
          heatmap_data$data <- filtered_data
          
          removeNotification("loading")
          showNotification("Data loaded successfully! Heatmap will be generated automatically.", 
                           type = "default", duration = 3)
          
        } else {
          removeNotification("loading")
          showNotification("No consolidated data available. Please ensure data is loaded properly.", 
                           type = "warning", duration = 5)
        }
        
      }, error = function(e) {
        removeNotification("loading")
        message("Error loading heatmap data: ", e$message)
        showNotification(paste("Error loading data:", e$message), type = "error", duration = 10)
      })
    })
    
    # Generate heatmap plot
    output$heatmap_plot <- renderPlotly({
      req(heatmap_data$data)
      
      # Check if heatmaply is available
      if (!requireNamespace("heatmaply", quietly = TRUE)) {
        showNotification("heatmaply package is required for interactive clustered heatmaps with dendrograms. Install with: install.packages('heatmaply')", 
                         type = "warning", duration = 10)
        use_heatmaply <- FALSE
      } else {
        use_heatmaply <- TRUE
      }
      
      # Create an interactive heatmap
      tryCatch({
        df <- heatmap_data$data
        message("Creating heatmap with ", nrow(df), " rows")
        message("Available columns: ", paste(names(df), collapse = ", "))
        
        # Create a frequency matrix of terms across genes/conditions
        # Use Description column for term names, and create columns for different conditions
        
        # Determine what to use as columns (x-axis)
        if ("gene" %in% names(df)) {
          x_var <- "gene"
        } else if ("mutation" %in% names(df)) {
          x_var <- "mutation"
        } else if ("cluster" %in% names(df)) {
          x_var <- "cluster"
        } else {
          # Create a composite identifier
          if ("method" %in% names(df)) {
            df$condition <- paste(df$method, df$cluster, sep = "_")
          } else if ("analysis_type" %in% names(df)) {
            df$condition <- paste(df$analysis_type, df$cluster, sep = "_")
          } else {
            df$condition <- df$cluster
          }
          x_var <- "condition"
        }
        
        # Use Description or term_name for y-axis
        if ("Description" %in% names(df)) {
          y_var <- "Description"
        } else if ("term_name" %in% names(df)) {
          y_var <- "term_name"
        } else {
          # Create generic term names
          df$term_id <- paste("Term", 1:nrow(df))
          y_var <- "term_id"
        }
        
        # Determine value column based on heatmap type
        if (input$heatmap_type == "pvalue" && "p.adjust" %in% names(df)) {
          df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
          value_var <- "heatmap_value"
          legend_title <- "-log10(p-value)"
        } else if (input$heatmap_type == "foldenrichment") {
          # Check for fold enrichment columns
          if ("FoldEnrichment" %in% names(df)) {
            value_var <- "FoldEnrichment"
            legend_title <- "Fold Enrichment"
          } else if ("fold_enrichment" %in% names(df)) {
            value_var <- "fold_enrichment"
            legend_title <- "Fold Enrichment"
          } else if ("RichFactor" %in% names(df)) {
            value_var <- "RichFactor"
            legend_title <- "Rich Factor"
          } else {
            # Calculate fold enrichment if we have the necessary columns
            if ("Count" %in% names(df) && "GeneRatio" %in% names(df)) {
              # Parse GeneRatio (format: "X/Y")
              df$heatmap_value <- sapply(df$GeneRatio, function(x) {
                parts <- as.numeric(strsplit(x, "/")[[1]])
                if (length(parts) == 2 && parts[2] > 0) parts[1]/parts[2] else 1
              })
              value_var <- "heatmap_value"
              legend_title <- "Gene Ratio"
            } else {
              df$heatmap_value <- 1
              value_var <- "heatmap_value"
              legend_title <- "Fold Enrichment (not available)"
            }
          }
        } else if (input$heatmap_type == "zscore") {
          if ("zScore" %in% names(df)) {
            value_var <- "zScore"
            legend_title <- "Z-score"
          } else if ("z_score" %in% names(df)) {
            value_var <- "z_score"
            legend_title <- "Z-score"
          } else {
            # Calculate z-score from p-values if available
            if ("p.adjust" %in% names(df)) {
              df$heatmap_value <- qnorm(1 - df$p.adjust/2) * sign(rnorm(nrow(df)))
              value_var <- "heatmap_value"
              legend_title <- "Z-score (approximated)"
            } else {
              df$heatmap_value <- 0
              value_var <- "heatmap_value"
              legend_title <- "Z-score (not available)"
            }
          }
        } else if (input$heatmap_type == "count" && "Count" %in% names(df)) {
          value_var <- "Count"
          legend_title <- "Gene Count"
        } else if (input$heatmap_type == "nes") {
          # Handle GSEA NES values
          if ("NES" %in% names(df)) {
            value_var <- "NES"
            legend_title <- "Normalized Enrichment Score"
          } else if ("nes" %in% names(df)) {
            value_var <- "nes"
            legend_title <- "Normalized Enrichment Score"
          } else if ("normalizedEnrichmentScore" %in% names(df)) {
            value_var <- "normalizedEnrichmentScore"
            legend_title <- "Normalized Enrichment Score"
          } else {
            # If no NES column, show warning and use p-values
            showNotification("NES values not found. Using p-values instead.", type = "warning", duration = 5)
            df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
            value_var <- "heatmap_value"
            legend_title <- "-log10(p-value)"
          }
        } else {
          # Default to p-value if available
          if ("p.adjust" %in% names(df)) {
            df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
            value_var <- "heatmap_value"
            legend_title <- "-log10(p-value)"
          } else {
            df$heatmap_value <- 1
            value_var <- "heatmap_value"
            legend_title <- "Presence"
          }
        }
        
        # Create matrix
        if (length(unique(df[[x_var]])) < 2) {
          # If only one condition, create a simple bar plot
          top_terms <- head(df[order(df[[value_var]], decreasing = TRUE), ], input$max_terms)
          
          p <- plot_ly(
            x = top_terms[[value_var]],
            y = substr(top_terms[[y_var]], 1, 50),
            type = "bar",
            orientation = "h",
            text = top_terms[[y_var]],
            hovertemplate = "%{text}<br>%{x:.2f}<extra></extra>",
            marker = list(color = top_terms[[value_var]], colorscale = "Viridis")
          ) %>%
            layout(
              title = paste("Top Terms -", legend_title),
              xaxis = list(title = legend_title),
              yaxis = list(title = "", tickfont = list(size = 10)),
              margin = list(l = 200)
            )
          
          return(p)
          
        } else {
          # Create proper heatmap matrix
          data_wide <- df %>%
            dplyr::select(all_of(c(y_var, x_var, value_var))) %>%
            tidyr::pivot_wider(names_from = all_of(x_var), 
                              values_from = all_of(value_var), 
                              values_fill = 0,
                              values_fn = mean)  # Handle duplicates by averaging
          
          # Convert to matrix
          mat <- as.matrix(data_wide[,-1])
          rownames(mat) <- substr(data_wide[[y_var]], 1, 80)  # Truncate long names
          
          # Store original row names and full descriptions
          original_rownames <- rownames(mat)
          full_descriptions_all <- data_wide[[y_var]]
          
          # Limit to max terms
          if (nrow(mat) > input$max_terms) {
            # Keep most significant terms (highest values)
            row_means <- rowMeans(mat, na.rm = TRUE)
            top_indices <- order(row_means, decreasing = TRUE)[1:input$max_terms]
            mat <- mat[top_indices, , drop = FALSE]
            full_descriptions <- full_descriptions_all[top_indices]
          } else {
            full_descriptions <- full_descriptions_all
          }
          
          # Use heatmaply if available, otherwise use plotly
          if (use_heatmaply) {
            # Prepare custom hover text for heatmaply
            custom_text <- matrix(
              paste0("Term: ", rep(full_descriptions, ncol(mat)), "<br>",
                     "Condition: ", rep(colnames(mat), each = nrow(mat)), "<br>",
                     "Value: "),
              nrow = nrow(mat),
              ncol = ncol(mat),
              byrow = FALSE
            )
            
            # Set color scale based on user selection
            if (input$color_scale == "red") {
              colors <- colorRampPalette(c("white", "red"))(256)
            } else if (input$color_scale == "blue") {
              colors <- colorRampPalette(c("white", "darkblue"))(256)
            } else if (input$color_scale == "RdBu") {
              colors <- colorRampPalette(c("blue", "white", "red"))(256)
            } else if (input$color_scale == "viridis") {
              colors <- viridis::viridis(256)
            } else if (input$color_scale == "YlOrRd") {
              colors <- colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(256)
            } else {
              colors <- colorRampPalette(c("white", "red"))(256)
            }
            
            # Apply scaling method if specified
            if (input$scale_method == "quantile") {
              # Use quantile breaks for adaptive scaling
              mat_values <- as.vector(mat)
              mat_values <- mat_values[!is.na(mat_values)]
              if (length(mat_values) > 0) {
                quantile_breaks <- quantile(mat_values, probs = seq(0, 1, length.out = 10))
                mat <- matrix(
                  cut(as.vector(mat), breaks = quantile_breaks, labels = FALSE, include.lowest = TRUE),
                  nrow = nrow(mat),
                  ncol = ncol(mat)
                )
              }
            } else if (input$scale_method == "fixed" && input$heatmap_type == "pvalue") {
              # Use fixed breaks for p-values
              mat[mat > 20] <- 20  # Cap at -log10(1e-20)
            }
            
            # Determine dendrogram settings
            dendrogram <- "none"
            if (input$cluster_rows && input$cluster_cols) {
              dendrogram <- "both"
            } else if (input$cluster_rows) {
              dendrogram <- "row"
            } else if (input$cluster_cols) {
              dendrogram <- "column"
            }
            
            # Prepare row annotations if requested
            annotation_result <- NULL
            row_side_colors <- NULL
            row_side_palette <- NULL
            
            if (input$show_annotations) {
              # Get annotation data for each row
              row_annotations <- df[!duplicated(df$Description), ] %>%
                dplyr::filter(Description %in% full_descriptions) %>%
                dplyr::arrange(match(Description, full_descriptions))
              
              # Use the safe annotation creation function
              annotation_result <- create_safe_row_annotations(mat, row_annotations, input$show_annotations)
              
              if (!is.null(annotation_result)) {
                row_side_colors <- annotation_result$colors
                row_side_palette <- annotation_result$palette
                
                # Debug the structure
                message("Annotation result structure:")
                message("  row_side_colors class: ", class(row_side_colors))
                message("  row_side_palette class: ", class(row_side_palette))
                if (is.list(row_side_palette)) {
                  message("  row_side_palette names: ", paste(names(row_side_palette), collapse = ", "))
                }
              }
            }
            
            # Create heatmaply heatmap with comprehensive error handling
            p <- tryCatch({
              if (!is.null(row_side_colors) && !is.null(row_side_palette) && input$show_annotations) {
                message("Creating heatmaply with annotations...")
                message("Row side colors dimensions: ", nrow(row_side_colors), " x ", ncol(row_side_colors))
                message("Row side palette Type colors: ", length(row_side_palette$Type))
                message("Row side palette Direction colors: ", length(row_side_palette$Direction))
                
                # Try with annotations first
                tryCatch({
                  # Debug palette structure
                  message("Attempting heatmaply with row_side_palette structure:")
                  message("  Type palette: ", paste(names(row_side_palette$Type), "=", row_side_palette$Type, collapse=", "))
                  message("  Direction palette: ", paste(names(row_side_palette$Direction), "=", row_side_palette$Direction, collapse=", "))
                  
                  heatmaply::heatmaply(
                    mat,
                    dendrogram = dendrogram,
                    colors = colors,
                    xlab = "",
                    ylab = "",
                    main = paste("Interactive Enrichment Heatmap -", legend_title),
                    margins = c(150, 300, 50, 50),
                    custom_hovertext = custom_text,
                    label_names = c("Row", "Column", "Value"),
                    fontsize_row = 10,
                    fontsize_col = 10,
                    showticklabels = c(TRUE, TRUE),
                    plot_method = "plotly",
                    row_side_colors = row_side_colors,
                    row_side_palette = row_side_palette
                  )
                }, error = function(e_ann) {
                  message("Heatmaply with annotations failed: ", e_ann$message)
                  message("Trying without annotations...")
                  # Fallback to heatmaply without annotations
                  heatmaply::heatmaply(
                    mat,
                    dendrogram = dendrogram,
                    colors = colors,
                    xlab = "",
                    ylab = "",
                    main = paste("Interactive Enrichment Heatmap -", legend_title),
                    margins = c(150, 250, 50, 50),
                    custom_hovertext = custom_text,
                    label_names = c("Row", "Column", "Value"),
                    fontsize_row = 10,
                    fontsize_col = 10,
                    showticklabels = c(TRUE, TRUE),
                    plot_method = "plotly"
                  )
                })
              } else {
                message("Creating heatmaply without annotations...")
                heatmaply::heatmaply(
                  mat,
                  dendrogram = dendrogram,
                  colors = colors,
                  xlab = "",
                  ylab = "",
                  main = paste("Interactive Enrichment Heatmap -", legend_title),
                  margins = c(150, 250, 50, 50),
                  custom_hovertext = custom_text,
                  label_names = c("Row", "Column", "Value"),
                  fontsize_row = 10,
                  fontsize_col = 10,
                  showticklabels = c(TRUE, TRUE),
                  plot_method = "plotly"
                )
              }
            }, error = function(e) {
              message("Heatmaply failed: ", e$message)
              message("Falling back to simple plotly heatmap...")
              
              # Fallback to basic plotly heatmap without annotations
              tryCatch({
                # Create basic hover text
                hover_text <- matrix(
                  paste0("Term: ", rep(full_descriptions, ncol(mat)), "<br>",
                         "Condition: ", rep(colnames(mat), each = nrow(mat)), "<br>",
                         "Value: ", round(as.vector(mat), 3)),
                  nrow = nrow(mat),
                  ncol = ncol(mat),
                  byrow = FALSE
                )
                
                plot_ly(
                  x = colnames(mat),
                  y = rownames(mat),
                  z = mat,
                  type = "heatmap",
                  colorscale = "Viridis",
                  hovertext = hover_text,
                  hovertemplate = "%{hovertext}<extra></extra>",
                  colorbar = list(title = legend_title)
                ) %>%
                  layout(
                    title = paste("Enrichment Heatmap -", legend_title),
                    xaxis = list(title = "", tickangle = -45),
                    yaxis = list(title = ""),
                    margin = list(l = 250, b = 100, t = 50)
                  )
              }, error = function(e2) {
                message("Plotly fallback also failed: ", e2$message)
                # Return empty plot with error message
                plot_ly() %>%
                  add_annotations(
                    text = paste("Error creating heatmap:<br>", e$message, "<br><br>Fallback error:<br>", e2$message),
                    x = 0.5, y = 0.5,
                    xref = "paper", yref = "paper",
                    showarrow = FALSE,
                    font = list(size = 12, color = "red")
                  )
              })
            })
            
            # Store the plotly object for download
            heatmap_data$plotly_object <- p
            
            return(p)
            
          } else {
            # Fallback to regular plotly without clustering
            # Create hover text with full descriptions
            hover_text <- matrix(
              paste0("Term: ", rep(full_descriptions, ncol(mat)), "<br>",
                     "Condition: ", rep(colnames(mat), each = nrow(mat)), "<br>",
                     "Value: ", round(as.vector(mat), 3)),
              nrow = nrow(mat),
              ncol = ncol(mat),
              byrow = FALSE
            )
            
            # Create color scale based on user selection
            if (input$color_scale == "red") {
              colorscale <- list(c(0, "white"), c(1, "red"))
            } else if (input$color_scale == "blue") {
              colorscale <- list(c(0, "white"), c(1, "darkblue"))
            } else if (input$color_scale == "RdBu") {
              colorscale <- "RdBu"
            } else if (input$color_scale == "viridis") {
              colorscale <- "Viridis"
            } else if (input$color_scale == "YlOrRd") {
              colorscale <- list(c(0, "white"), c(0.5, "yellow"), c(0.75, "orange"), c(1, "red"))
            } else {
              colorscale <- list(c(0, "white"), c(1, "red"))
            }
            
            # Create basic plotly heatmap
            p <- plot_ly(
              x = colnames(mat),
              y = rownames(mat),
              z = mat,
              type = "heatmap",
              colorscale = colorscale,
              hovertext = hover_text,
              hovertemplate = "%{hovertext}<extra></extra>",
              colorbar = list(title = legend_title)
            ) %>%
              layout(
                title = paste("Interactive Enrichment Heatmap -", legend_title),
                xaxis = list(
                  title = "",
                  tickangle = -45,
                  tickfont = list(size = 12)
                ),
                yaxis = list(
                  title = "",
                  tickfont = list(size = 10)
                ),
                margin = list(l = 250, b = 100, t = 50)
              )
            
            # Store the plotly object for download
            heatmap_data$plotly_object <- p
            
            return(p)
          }
        }
        
      }, error = function(e) {
        message("Heatmap error: ", e$message)
        # Return an empty plotly object with error message
        plot_ly() %>%
          add_annotations(
            text = paste("Error creating heatmap:\n", e$message),
            x = 0.5,
            y = 0.5,
            xref = "paper",
            yref = "paper",
            showarrow = FALSE,
            font = list(size = 14)
          ) %>%
          layout(
            xaxis = list(visible = FALSE),
            yaxis = list(visible = FALSE)
          )
      })
    })
    
    # Display data table
    output$heatmap_data <- DT::renderDataTable({
      req(heatmap_data$data)
      
      # Select available columns
      available_cols <- intersect(
        c("gene", "mutation_perturbation", "cluster", "enrichment_type", "Description", 
          "term_name", "p.adjust", "Count", "FoldEnrichment", "RichFactor"),
        names(heatmap_data$data)
      )
      
      DT::datatable(
        heatmap_data$data %>%
          dplyr::select(all_of(available_cols)) %>%
          dplyr::arrange(p.adjust),
        options = list(pageLength = 25, scrollX = TRUE),
        filter = 'top'
      ) %>%
        DT::formatSignif(columns = "p.adjust", digits = 3)
    })
    
    # Settings info
    output$settings_info <- renderPrint({
      cat("Current Settings:\n")
      cat("================\n")
      cat("Heatmap Type:", input$heatmap_type, "\n")
      cat("Method Filter:", input$method_filter, "\n")
      cat("Enrichment Types:", paste(input$enrichment_types, collapse = ", "), "\n")
      cat("Cluster:", input$cluster_select, "\n")
      cat("Max Terms:", input$max_terms, "\n")
      cat("Direction Filter:", input$direction_filter, "\n")
      cat("Color Scale:", input$color_scale, "\n")
      cat("Scaling Method:", input$scale_method, "\n")
      cat("Show Annotations:", input$show_annotations, "\n")
      cat("Cluster Rows:", input$cluster_rows, "\n")
      cat("Cluster Columns:", input$cluster_cols, "\n")
      
      # GSEA-specific settings
      if (input$heatmap_type == "nes" || "GSEA" %in% input$enrichment_types) {
        cat("\nGSEA Settings:\n")
        cat("GSEA Only:", input$gsea_only, "\n")
        cat("Min |NES| Threshold:", input$nes_threshold, "\n")
      }
      
      if (!is.null(heatmap_data$data)) {
        cat("\nData Summary:\n")
        cat("Total rows:", nrow(heatmap_data$data), "\n")
        cat("Unique genes:", length(unique(heatmap_data$data$gene)), "\n")
        cat("Unique terms:", length(unique(heatmap_data$data$term_name)), "\n")
      }
    })
    
    # PDF Export functionality
    observeEvent(input$generate_pdf, {
      req(heatmap_data$data)
      
      showNotification("Generating PDF heatmap...", type = "default", duration = NULL, id = "pdf_loading")
      
      # Create a temporary file
      temp_file <- tempfile(fileext = ".pdf")
      
      tryCatch({
        # Get the current data
        df <- heatmap_data$data
        
        # Same data preparation as for interactive heatmap
        if ("gene" %in% names(df)) {
          x_var <- "gene"
        } else if ("mutation" %in% names(df)) {
          x_var <- "mutation"
        } else if ("cluster" %in% names(df)) {
          x_var <- "cluster"
        } else {
          if ("method" %in% names(df)) {
            df$condition <- paste(df$method, df$cluster, sep = "_")
          } else {
            df$condition <- df$cluster
          }
          x_var <- "condition"
        }
        
        y_var <- if ("Description" %in% names(df)) "Description" else "term_id"
        
        # Get value variable based on heatmap type
        value_var <- switch(input$heatmap_type,
          "pvalue" = {
            df$heatmap_value <- -log10(pmax(df$p.adjust, 1e-300))
            "heatmap_value"
          },
          "foldenrichment" = {
            if ("FoldEnrichment" %in% names(df)) {
              "FoldEnrichment"
            } else if ("fold_enrichment" %in% names(df)) {
              "fold_enrichment"
            } else if ("RichFactor" %in% names(df)) {
              "RichFactor"
            } else {
              df$heatmap_value <- 1
              "heatmap_value"
            }
          },
          "zscore" = {
            if ("zScore" %in% names(df)) {
              "zScore"
            } else if ("z_score" %in% names(df)) {
              "z_score"
            } else {
              df$heatmap_value <- 0
              "heatmap_value"
            }
          },
          "nes" = if("NES" %in% names(df)) "NES" else "heatmap_value",
          "count" = if("Count" %in% names(df)) "Count" else "heatmap_value",
          "heatmap_value"
        )
        
        # Create matrix
        data_wide <- df %>%
          dplyr::select(all_of(c(y_var, x_var, value_var))) %>%
          tidyr::pivot_wider(names_from = all_of(x_var), 
                            values_from = all_of(value_var), 
                            values_fill = 0,
                            values_fn = mean)
        
        mat <- as.matrix(data_wide[,-1])
        rownames(mat) <- substr(data_wide[[y_var]], 1, 100)
        
        # Limit rows if needed
        if (nrow(mat) > input$max_terms) {
          row_means <- rowMeans(mat, na.rm = TRUE)
          mat <- mat[order(abs(row_means), decreasing = TRUE)[1:input$max_terms], , drop = FALSE]
        }
        
        # Create color function based on settings
        if (input$color_scale == "YlOrRd") {
          col_fun <- circlize::colorRamp2(
            breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 9),
            colors = RColorBrewer::brewer.pal(9, "YlOrRd")
          )
        } else if (input$color_scale == "RdBu" || input$heatmap_type == "nes") {
          col_fun <- circlize::colorRamp2(
            breaks = c(min(mat, na.rm = TRUE), 0, max(mat, na.rm = TRUE)), 
            colors = c("blue", "white", "red")
          )
        } else if (input$color_scale == "viridis") {
          col_fun <- circlize::colorRamp2(
            breaks = seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = 100),
            colors = viridis::viridis(100)
          )
        } else {
          col_fun <- circlize::colorRamp2(
            breaks = seq(0, max(mat, na.rm = TRUE), length.out = 100),
            colors = colorRampPalette(c("white", ifelse(input$color_scale == "blue", "darkblue", "red")))(100)
          )
        }
        
        # Create row annotations if requested
        row_ha <- NULL
        if (input$show_annotations && nrow(df) > 0) {
          tryCatch({
            # Get unique rows for annotations
            row_ann_data <- df[match(rownames(mat), substr(df[[y_var]], 1, 100)), ]
            
            if (nrow(row_ann_data) > 0 && all(c("enrichment_type", "direction") %in% names(row_ann_data))) {
              # Use validated colors
              color_vectors <- create_validated_colors()
              
              # Ensure annotation data is complete
              row_ann_data$enrichment_type[is.na(row_ann_data$enrichment_type)] <- "Unknown"
              row_ann_data$direction[is.na(row_ann_data$direction)] <- "Unknown"
              
              # Add colors for unknown categories if needed
              if (!"Unknown" %in% names(color_vectors$type)) {
                color_vectors$type["Unknown"] <- "#808080"
              }
              if (!"Unknown" %in% names(color_vectors$direction)) {
                color_vectors$direction["Unknown"] <- "#808080"
              }
              
              row_ha <- ComplexHeatmap::rowAnnotation(
                Type = row_ann_data$enrichment_type,
                Direction = row_ann_data$direction,
                col = list(Type = color_vectors$type, Direction = color_vectors$direction),
                annotation_name_gp = grid::gpar(fontsize = 8),
                simple_anno_size = grid::unit(0.3, "cm")
              )
            }
          }, error = function(e) {
            message("Error creating ComplexHeatmap annotations: ", e$message)
            row_ha <<- NULL
          })
        }
        
        # Create ComplexHeatmap
        ht <- ComplexHeatmap::Heatmap(
          mat,
          name = switch(input$heatmap_type,
            "pvalue" = "-log10(p.adj)",
            "nes" = "NES",
            "count" = "Count",
            "foldenrichment" = "Fold Enrich",
            "zscore" = "Z-score",
            "Value"
          ),
          col = col_fun,
          cluster_rows = input$cluster_rows && nrow(mat) > 1,
          cluster_columns = input$cluster_cols && ncol(mat) > 1,
          show_row_names = TRUE,
          show_column_names = TRUE,
          row_names_gp = grid::gpar(fontsize = max(4, input$pdf_fontsize - nrow(mat)/10)),
          column_names_gp = grid::gpar(fontsize = input$pdf_fontsize),
          column_names_rot = 45,
          column_title = paste("Enrichment Heatmap -", input$heatmap_type),
          column_title_gp = grid::gpar(fontsize = input$pdf_fontsize + 2),
          left_annotation = row_ha,
          use_raster = FALSE
        )
        
        # Save to PDF
        pdf(temp_file, width = input$pdf_width, height = input$pdf_height)
        ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
        dev.off()
        
        removeNotification("pdf_loading")
        
        # Update status
        output$pdf_status <- renderUI({
          tags$div(
            class = "alert alert-success",
            icon("check-circle"),
            "PDF generated successfully!",
            br(),
            downloadButton(session$ns("download_pdf"), "Download PDF", class = "btn-primary")
          )
        })
        
        # Store the temp file path
        heatmap_data$pdf_file <- temp_file
        
      }, error = function(e) {
        removeNotification("pdf_loading")
        output$pdf_status <- renderUI({
          tags$div(
            class = "alert alert-danger",
            icon("exclamation-triangle"),
            "Error generating PDF: ", e$message
          )
        })
      })
    })
    
    # Download handler for PDF
    output$download_pdf <- downloadHandler(
      filename = function() {
        paste0("heatmap_", format(Sys.Date(), "%Y%m%d"), "_", 
               input$heatmap_type, "_", 
               gsub("[^A-Za-z0-9]", "", paste(input$enrichment_types, collapse = "_")), 
               ".pdf")
      },
      content = function(file) {
        if (!is.null(heatmap_data$pdf_file) && file.exists(heatmap_data$pdf_file)) {
          file.copy(heatmap_data$pdf_file, file)
        } else {
          showNotification("No PDF file available. Please generate one first.", type = "error")
        }
      }
    )
    
    # Download handler for interactive heatmap
    output$download_heatmap <- downloadHandler(
      filename = function() {
        paste0("heatmap_interactive_", Sys.Date(), ".html")
      },
      content = function(file) {
        if (!is.null(heatmap_data$plotly_object)) {
          htmlwidgets::saveWidget(
            heatmap_data$plotly_object,
            file,
            selfcontained = TRUE
          )
        } else {
          showNotification("No heatmap available. Please generate one first.", type = "error")
        }
      }
    )
    
    return(heatmap_data)
  })
}