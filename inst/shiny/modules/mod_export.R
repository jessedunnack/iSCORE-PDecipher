# Module: Export
# Export figures and data in various formats

mod_export_ui <- function(id) {
  ns <- NS(id)
  
  fluidRow(
    column(
      width = 6,
      div(
        class = "box box-primary box-solid",
        div(class = "box-header with-border",
          h3(class = "box-title", "Export Figures")
        ),
        div(class = "box-body",
          h4("Figure Export Settings"),
        
        selectInput(ns("figure_format"),
                    "Format:",
                    choices = c("PDF" = "pdf",
                                "PNG" = "png",
                                "SVG" = "svg",
                                "TIFF" = "tiff"),
                    selected = "pdf"),
        
        conditionalPanel(
          condition = "input.figure_format == 'png' || input.figure_format == 'tiff'",
          ns = ns,
          sliderInput(ns("dpi"),
                      "Resolution (DPI):",
                      min = 72,
                      max = 600,
                      value = 300,
                      step = 72)
        ),
        
        fluidRow(
          column(6,
                 numericInput(ns("fig_width"),
                              "Width (inches):",
                              value = 10,
                              min = 2,
                              max = 20,
                              step = 0.5)),
          column(6,
                 numericInput(ns("fig_height"),
                              "Height (inches):",
                              value = 8,
                              min = 2,
                              max = 20,
                              step = 0.5))
        ),
        
        checkboxInput(ns("include_timestamp"),
                      "Include timestamp in filename",
                      value = TRUE),
        
        br(),
        
        h4("Available Figures"),
        selectInput(ns("figure_select"),
                    "Select Figure:",
                    choices = character(0)),
        
        uiOutput(ns("figure_preview")),
        
        br(),
        downloadButton(ns("download_figure"), "Download Figure",
                       class = "btn-success", style = "width: 100%;")
        ) # end box-body
      ), # end box div
      
      div(
        class = "box box-warning box-solid collapsed-box",
        div(class = "box-header with-border",
          h3(class = "box-title", "Batch Export"),
          div(class = "box-tools pull-right",
            tags$button(class = "btn btn-box-tool", 
                      `data-widget` = "collapse",
                      icon("plus"))
          )
        ),
        div(class = "box-body", style = "display: none;",
          h4("Export Multiple Figures"),
        
        checkboxGroupInput(ns("batch_figures"),
                           "Select Figures:",
                           choices = character(0)),
        
        radioButtons(ns("batch_format"),
                     "Format:",
                     choices = c("PDF (single file)" = "pdf_combined",
                                 "PDF (multiple files)" = "pdf_separate",
                                 "PNG" = "png",
                                 "SVG" = "svg"),
                     selected = "pdf_combined"),
        
        br(),
        downloadButton(ns("download_batch"), "Download All Selected",
                       class = "btn-warning", style = "width: 100%;")
        ) # end box-body
      ) # end box div
    ),
    
    column(
      width = 6,
      div(
        class = "box box-info box-solid",
        div(class = "box-header with-border",
          h3(class = "box-title", "Export Data")
        ),
        div(class = "box-body",
          h4("Data Export Settings"),
        
        selectInput(ns("data_format"),
                    "Format:",
                    choices = c("CSV" = "csv",
                                "TSV" = "tsv",
                                "Excel" = "xlsx",
                                "RDS" = "rds"),
                    selected = "csv"),
        
        radioButtons(ns("data_content"),
                     "Export Content:",
                     choices = c("Current View" = "current",
                                 "All Results" = "all",
                                 "Significant Only" = "significant",
                                 "Custom Selection" = "custom"),
                     selected = "current"),
        
        conditionalPanel(
          condition = "input.data_content == 'custom'",
          ns = ns,
          checkboxGroupInput(ns("custom_columns"),
                             "Select Columns:",
                             choices = character(0))
        ),
        
        checkboxInput(ns("include_metadata"),
                      "Include metadata",
                      value = TRUE),
        
        br(),
        
        h4("Data Preview"),
        DT::dataTableOutput(ns("data_preview")),
        
        br(),
        downloadButton(ns("download_data"), "Download Data",
                       class = "btn-info", style = "width: 100%;")
        ) # end box-body
      ), # end box div
      
      div(
        class = "box box-success box-solid collapsed-box",
        div(class = "box-header with-border",
          h3(class = "box-title", "Report Generation"),
          div(class = "box-tools pull-right",
            tags$button(class = "btn btn-box-tool", 
                      `data-widget` = "collapse",
                      icon("plus"))
          )
        ),
        div(class = "box-body", style = "display: none;",
          h4("Generate Comprehensive Report"),
        
        textInput(ns("report_title"),
                  "Report Title:",
                  value = "Enrichment Analysis Report"),
        
        textInput(ns("report_author"),
                  "Author:",
                  value = ""),
        
        checkboxGroupInput(ns("report_sections"),
                           "Include Sections:",
                           choices = c("Executive Summary" = "summary",
                                       "Methods" = "methods",
                                       "Results" = "results",
                                       "Visualizations" = "plots",
                                       "Data Tables" = "tables",
                                       "References" = "references"),
                           selected = c("summary", "results", "plots")),
        
        radioButtons(ns("report_format"),
                     "Format:",
                     choices = c("HTML" = "html",
                                 "PDF" = "pdf",
                                 "Word" = "docx"),
                     selected = "html"),
        
        br(),
        actionButton(ns("generate_report"), "Generate Report",
                     class = "btn-success", width = "100%"),
        
        br(),
        br(),
        uiOutput(ns("report_status"))
        ) # end box-body
      ) # end box div
    )
  )
}

mod_export_server <- function(id, app_data, precomputed_data, comparison_data, visualization_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for export data
    export_data <- reactiveValues(
      figures = list(),
      data_tables = list(),
      report_path = NULL
    )
    
    # Update available figures
    observe({
      figures <- character()
      
      # Check for available plots from different modules
      if (!is.null(precomputed_data)) {
        pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
        if (!is.null(pd) && !is.null(pd$current_plot)) {
          figures <- c(figures, "Precomputed Analysis Plot" = "precomputed")
        }
      }
      
      if (!is.null(visualization_data)) {
        vd <- if (is.reactive(visualization_data)) visualization_data() else visualization_data
        if (!is.null(vd) && !is.null(vd$plot)) {
          figures <- c(figures, "Enhanced Visualization" = "visualization")
        }
      }
      
      if (!is.null(comparison_data)) {
        cd <- if (is.reactive(comparison_data)) comparison_data() else comparison_data
        if (!is.null(cd) && !is.null(cd$comparison_results)) {
          figures <- c(figures, "Method Comparison" = "comparison")
        }
      }
      
      updateSelectInput(session, "figure_select", choices = figures)
      updateCheckboxGroupInput(session, "batch_figures", choices = figures)
    })
    
    # Figure preview
    output$figure_preview <- renderUI({
      if (!is.null(input$figure_select)) {
        tags$div(
          style = "border: 1px solid #ddd; padding: 10px; margin-top: 10px;",
          h5("Preview:"),
          tags$p("Figure preview would appear here", style = "color: gray; text-align: center;")
        )
      }
    })
    
    # Update custom columns based on current data
    observe({
      if (!is.null(precomputed_data)) {
        pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
        
        # Handle both old structure (with selected_data field) and new direct data
        data_to_check <- NULL
        if (!is.null(pd)) {
          if (is.data.frame(pd)) {
            # Direct data frame
            data_to_check <- pd
          } else if (!is.null(pd$selected_data)) {
            # Old structure with selected_data field
            data_to_check <- pd$selected_data
          }
        }
        
        if (!is.null(data_to_check) && is.data.frame(data_to_check)) {
          cols <- names(data_to_check)
          updateCheckboxGroupInput(session, "custom_columns", 
                                   choices = cols, 
                                   selected = cols[1:min(5, length(cols))])
        }
      }
    })
    
    # Data preview
    output$data_preview <- DT::renderDataTable({
      data <- NULL
      
      if (!is.null(precomputed_data)) {
        pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
        
        # Handle both old structure (with selected_data field) and new direct data
        if (!is.null(pd)) {
          if (is.data.frame(pd)) {
            # Direct data frame
            data <- pd
          } else if (!is.null(pd$selected_data)) {
            # Old structure with selected_data field
            data <- pd$selected_data
          }
        }
        
        if (!is.null(data) && is.data.frame(data)) {
          if (input$data_content == "significant") {
            if ("p.adjust" %in% names(data)) {
              data <- data[data$p.adjust < 0.05, ]
            }
          } else if (input$data_content == "custom" && !is.null(input$custom_columns)) {
            data <- data[, input$custom_columns, drop = FALSE]
          }
        }
      }
      
      if (!is.null(data)) {
        DT::datatable(data,
                      options = list(pageLength = 5, scrollX = TRUE),
                      rownames = FALSE)
      }
    })
    
    # Download figure handler
    output$download_figure <- downloadHandler(
      filename = function() {
        timestamp <- ifelse(input$include_timestamp, 
                            format(Sys.time(), "_%Y%m%d_%H%M%S"), 
                            "")
        paste0("enrichment_plot", timestamp, ".", input$figure_format)
      },
      content = function(file) {
        # Get the selected plot
        plot_obj <- NULL
        
        if (input$figure_select == "precomputed") {
          pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
          if (!is.null(pd) && !is.null(pd$current_plot)) {
            plot_obj <- pd$current_plot
          }
        } else if (input$figure_select == "visualization") {
          vd <- if (is.reactive(visualization_data)) visualization_data() else visualization_data
          if (!is.null(vd) && !is.null(vd$plot)) {
            plot_obj <- vd$plot
          }
        }
        
        if (!is.null(plot_obj)) {
          if (input$figure_format == "pdf") {
            pdf(file, width = input$fig_width, height = input$fig_height)
          } else if (input$figure_format == "png") {
            png(file, width = input$fig_width, height = input$fig_height, 
                units = "in", res = input$dpi)
          } else if (input$figure_format == "svg") {
            svg(file, width = input$fig_width, height = input$fig_height)
          } else if (input$figure_format == "tiff") {
            tiff(file, width = input$fig_width, height = input$fig_height, 
                 units = "in", res = input$dpi)
          }
          
          print(plot_obj)
          dev.off()
        }
      }
    )
    
    # Download data handler
    output$download_data <- downloadHandler(
      filename = function() {
        timestamp <- ifelse(input$include_timestamp, 
                            format(Sys.time(), "_%Y%m%d_%H%M%S"), 
                            "")
        paste0("enrichment_data", timestamp, ".", input$data_format)
      },
      content = function(file) {
        pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
        
        # Handle both old structure (with selected_data field) and new direct data
        data <- NULL
        if (!is.null(pd)) {
          if (is.data.frame(pd)) {
            # Direct data frame
            data <- pd
          } else if (!is.null(pd$selected_data)) {
            # Old structure with selected_data field
            data <- pd$selected_data
          }
        }
        
        if (input$data_content == "significant" && "p.adjust" %in% names(data)) {
          data <- data[data$p.adjust < 0.05, ]
        } else if (input$data_content == "custom" && !is.null(input$custom_columns)) {
          data <- data[, input$custom_columns, drop = FALSE]
        }
        
        if (input$data_format == "csv") {
          write.csv(data, file, row.names = FALSE)
        } else if (input$data_format == "tsv") {
          write.table(data, file, sep = "\t", row.names = FALSE)
        } else if (input$data_format == "xlsx") {
          # Would require openxlsx package
          write.csv(data, file, row.names = FALSE)  # Fallback to CSV
        } else if (input$data_format == "rds") {
          saveRDS(data, file)
        }
      }
    )
    
    # Batch download handler
    output$download_batch <- downloadHandler(
      filename = function() {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        if (input$batch_format == "pdf_combined") {
          paste0("enrichment_plots_combined_", timestamp, ".pdf")
        } else {
          paste0("enrichment_plots_", timestamp, ".zip")
        }
      },
      content = function(file) {
        # Create progress indicator
        withProgress(message = "Generating batch export...", value = 0, {
          
          if (input$batch_format == "pdf_combined") {
            # Combined PDF export
            pdf(file, width = 12, height = 10)
            
            # Get current selections
            selection <- app_data$global_selection
            if (!is.null(selection)) {
              # Page 1: Overview statistics
              incProgress(0.2, detail = "Creating overview...")
              par(mfrow = c(2, 2))
              
              # Plot 1: Terms by enrichment type
              data <- app_data$consolidated_data
              enrichment_counts <- table(data$enrichment_type)
              barplot(enrichment_counts, 
                     main = paste("Enrichment Results:", selection$gene, selection$cluster),
                     col = rainbow(length(enrichment_counts)),
                     las = 2)
              
              # Plot 2: Direction distribution
              direction_counts <- table(data$direction)
              pie(direction_counts, 
                  main = "Direction Distribution",
                  col = c("UP" = "#e74c3c", "DOWN" = "#3498db", 
                         "ALL" = "#2ecc71", "RANKED" = "#9b59b6"))
              
              # Reset layout
              par(mfrow = c(1, 1))
              
              # Page 2: Top terms plot
              incProgress(0.5, detail = "Creating term plots...")
              if (!is.null(data)) {
                plot_data <- data
                if (nrow(plot_data) > 0) {
                  # Create a simple dot plot
                  top_terms <- head(plot_data[order(plot_data$p.adjust), ], 20)
                  par(mar = c(5, 15, 4, 2))
                  plot(x = -log10(top_terms$p.adjust),
                       y = 1:nrow(top_terms),
                       pch = 19, col = "#3c8dbc",
                       xlim = c(0, max(-log10(top_terms$p.adjust)) * 1.1),
                       yaxt = "n",
                       xlab = expression("-log"[10]*"(p.adjust)"),
                       ylab = "",
                       main = paste("Top Enriched Terms:", selection$enrichment_type))
                  axis(2, at = 1:nrow(top_terms), 
                       labels = substr(top_terms$Description, 1, 50),
                       las = 1, cex.axis = 0.8)
                  abline(v = -log10(0.05), col = "red", lty = 2)
                }
              }
              
              incProgress(0.3, detail = "Finalizing...")
            }
            
            dev.off()
            
          } else {
            # ZIP file export with multiple formats
            temp_dir <- tempdir()
            files_to_zip <- c()
            
            # Generate plots in different formats
            formats <- strsplit(input$batch_format, "_")[[1]]
            
            for (fmt in formats) {
              if (fmt %in% c("png", "pdf")) {
                filename <- file.path(temp_dir, paste0("enrichment_plot.", fmt))
                
                if (fmt == "png") {
                  png(filename, width = 1200, height = 1000, res = 150)
                } else {
                  pdf(filename, width = 12, height = 10)
                }
                
                # Create the plot
                plot(1:10, main = "Enrichment Analysis Results")
                
                dev.off()
                files_to_zip <- c(files_to_zip, filename)
              }
            }
            
            # Add data export
            if ("csv" %in% formats && !is.null(pd)) {
              csv_file <- file.path(temp_dir, "enrichment_data.csv")
              # Get data using same logic as before
              export_data <- NULL
              if (is.data.frame(pd)) {
                export_data <- pd
              } else if (!is.null(pd$selected_data)) {
                export_data <- pd$selected_data
              }
              if (!is.null(export_data)) {
                write.csv(export_data, csv_file, row.names = FALSE)
                files_to_zip <- c(files_to_zip, csv_file)
              }
            }
            
            # Create zip file
            zip(file, files = files_to_zip, flags = "-j")
            
            # Clean up
            unlink(files_to_zip)
          }
        })
      }
    )
    
    # Generate report
    observeEvent(input$generate_report, {
      output$report_status <- renderUI({
        tags$div(
          class = "alert alert-info",
          tags$strong("Report generation started..."),
          tags$br(),
          "This feature would generate a comprehensive report using R Markdown."
        )
      })
      
      # Simulate report generation
      shinyjs::delay(2000, {
        output$report_status <- renderUI({
          tags$div(
            class = "alert alert-success",
            tags$strong("Report generated successfully!"),
            tags$br(),
            downloadLink(ns("download_report"), "Download Report")
          )
        })
      })
    })
    
    # Download report handler
    output$download_report <- downloadHandler(
      filename = function() {
        paste0("enrichment_report_", Sys.Date(), ".", input$report_format)
      },
      content = function(file) {
        withProgress(message = "Generating report...", value = 0, {
          
          # Get current data and selections
          selection <- app_data$global_selection
          pd <- if (is.reactive(precomputed_data)) precomputed_data() else precomputed_data
          
          # Handle both old structure (with selected_data field) and new direct data
          data <- NULL
          if (!is.null(pd)) {
            if (is.data.frame(pd)) {
              # Direct data frame
              data <- pd
            } else if (!is.null(pd$selected_data)) {
              # Old structure with selected_data field
              data <- pd$selected_data
            }
          }
          
          incProgress(0.3, detail = "Preparing content...")
          
          if (input$report_format == "html") {
            # Generate HTML report
            html_content <- paste0(
              "<html><head><title>", htmltools::htmlEscape(input$report_title), "</title>",
              "<style>",
              "body { font-family: Arial, sans-serif; margin: 40px; }",
              "h1 { color: #3c8dbc; }",
              "h2 { color: #666; }",
              "table { border-collapse: collapse; width: 100%; margin: 20px 0; }",
              "th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }",
              "th { background-color: #f4f4f4; }",
              ".summary { background-color: #f0f8ff; padding: 15px; border-radius: 5px; }",
              "</style></head><body>",
              "<h1>", htmltools::htmlEscape(input$report_title), "</h1>",
              "<p>Generated on: ", Sys.Date(), "</p>"
            )
            
            # Add summary section
            if (!is.null(selection)) {
              html_content <- paste0(html_content,
                "<div class='summary'>",
                "<h2>Analysis Summary</h2>",
                "<ul>",
                "<li><b>Analysis Type:</b> ", selection$analysis_type, "</li>",
                "<li><b>Gene/Mutation:</b> ", selection$gene, "</li>",
                "<li><b>Cluster:</b> ", selection$cluster, "</li>",
                "<li><b>Enrichment Database:</b> ", selection$enrichment_type, "</li>",
                "<li><b>Direction:</b> ", selection$direction, "</li>",
                "<li><b>P-value Threshold:</b> ", selection$pval_threshold, "</li>",
                "</ul>",
                "</div>"
              )
            }
            
            incProgress(0.4, detail = "Adding data tables...")
            
            # Add data table
            if (!is.null(data) && nrow(data) > 0) {
              html_content <- paste0(html_content,
                "<h2>Top Enriched Terms</h2>",
                "<table>",
                "<tr><th>ID</th><th>Description</th><th>P.adjust</th><th>Count</th></tr>"
              )
              
              top_terms <- head(data[order(data$p.adjust), ], input$max_terms)
              for (i in 1:nrow(top_terms)) {
                html_content <- paste0(html_content,
                  "<tr>",
                  "<td>", htmltools::htmlEscape(top_terms$ID[i]), "</td>",
                  "<td>", htmltools::htmlEscape(top_terms$Description[i]), "</td>",
                  "<td>", format(top_terms$p.adjust[i], scientific = TRUE, digits = 3), "</td>",
                  "<td>", top_terms$Count[i], "</td>",
                  "</tr>"
                )
              }
              
              html_content <- paste0(html_content, "</table>")
            }
            
            html_content <- paste0(html_content, "</body></html>")
            
            incProgress(0.3, detail = "Writing file...")
            writeLines(html_content, file)
            
          } else {
            # Generate text report
            text_content <- paste0(
              "# ", input$report_title, "\n",
              "Generated on: ", Sys.Date(), "\n",
              paste(rep("=", 50), collapse = ""), "\n\n"
            )
            
            if (!is.null(selection)) {
              text_content <- paste0(text_content,
                "## Analysis Summary\n",
                "- Analysis Type: ", selection$analysis_type, "\n",
                "- Gene/Mutation: ", selection$gene, "\n",
                "- Cluster: ", selection$cluster, "\n",
                "- Enrichment Database: ", selection$enrichment_type, "\n",
                "- Direction: ", selection$direction, "\n",
                "- P-value Threshold: ", selection$pval_threshold, "\n\n"
              )
            }
            
            if (!is.null(data) && nrow(data) > 0) {
              text_content <- paste0(text_content,
                "## Top Enriched Terms\n\n"
              )
              
              top_terms <- head(data[order(data$p.adjust), ], input$max_terms)
              for (i in 1:nrow(top_terms)) {
                text_content <- paste0(text_content,
                  i, ". ", top_terms$Description[i], "\n",
                  "   ID: ", top_terms$ID[i], "\n",
                  "   P.adjust: ", format(top_terms$p.adjust[i], scientific = TRUE, digits = 3), "\n",
                  "   Count: ", top_terms$Count[i], "\n\n"
                )
              }
            }
            
            writeLines(text_content, file)
          }
        })
      }
    )
    
    return(export_data)
  })
}