# Module: Pathview KEGG Pathway Visualization
# Provides interactive KEGG pathway diagrams with DE gene highlighting

# UI function
mod_pathview_ui <- function(id) {
  ns <- NS(id)
  
  tagList(
    fluidRow(
      # Control panel
      column(3,
        div(
          class = "box box-primary box-solid",
          div(class = "box-header with-border",
            h3(class = "box-title", "Pathway Settings")
          ),
          div(class = "box-body",
        
        # KEGG pathway selection
        selectInput(
          ns("pathway_id"),
          "KEGG Pathway:",
          choices = character(0),
          selected = NULL
        ),
        
        # Gene data type
        radioButtons(
          ns("gene_data_type"),
          "Color Genes By:",
          choices = c(
            "Log2 Fold Change" = "logfc",
            "P-value Significance" = "pvalue",
            "Binary (DE vs Non-DE)" = "binary"
          ),
          selected = "logfc"
        ),
        
        
        # Log2FC threshold for binary mode
        conditionalPanel(
          condition = paste0("input['", ns("gene_data_type"), "'] == 'binary'"),
          numericInput(
            ns("logfc_threshold"),
            "Log2FC Threshold:",
            value = 1.0,
            min = 0.1,
            max = 5.0,
            step = 0.1
          )
        ),
        
        # Species selection
        selectInput(
          ns("species"),
          "Species:",
          choices = c("Human" = "hsa", "Mouse" = "mmu", "Rat" = "rno"),
          selected = "hsa"
        ),
        
        # Output format
        selectInput(
          ns("output_format"),
          "Output Format:",
          choices = c("PNG" = "png", "PDF" = "pdf"),
          selected = "png"
        ),
        
        hr(),
        
        # Action buttons
        actionButton(
          ns("generate_pathway"),
          "Generate Pathway",
          class = "btn-primary btn-block",
          icon = icon("image")
        ),
        
        br(),
        
        actionButton(
          ns("test_button"),
          "Test Module",
          class = "btn-info btn-block",
          icon = icon("bug")
        ),
        
        br(),
        
        downloadButton(
          ns("download_pathway"),
          "Download Pathway",
          class = "btn-success",
          style = "width: 100%;"
        )
        ) # end box-body
      ) # end box div
      ), # end column
      
      # Visualization panel
      column(9,
        div(
          class = "box box-primary",
          div(class = "box-header with-border",
            h3(class = "box-title", "KEGG Pathway Diagram")
          ),
          div(class = "box-body",
            tabsetPanel(
              id = ns("pathway_tabs"),
          
          # Pathway diagram tab
          tabPanel(
            "Pathway Diagram",
            value = "diagram",
            br(),
            withSpinner(
              uiOutput(ns("pathway_image")),
              type = 4,
              color = "#3c8dbc"
            )
          ),
          
          # Gene mapping tab
          tabPanel(
            "Gene Mapping",
            value = "mapping",
            br(),
            DT::dataTableOutput(ns("gene_mapping_table"))
          ),
          
          # Pathway info tab
          tabPanel(
            "Pathway Info",
            value = "info",
            br(),
            verbatimTextOutput(ns("pathway_info"))
          )
        )
        ) # end box-body
      ) # end box div
      ) # end column
    ),
    
    # Additional info row
    fluidRow(
      column(12,
        div(
          class = "box box-default collapsed-box",
          div(class = "box-header with-border",
            h3(class = "box-title", "Instructions"),
            div(class = "box-tools pull-right",
              tags$button(class = "btn btn-box-tool", 
                        `data-widget` = "collapse",
                        icon("plus"))
            )
          ),
          div(class = "box-body", style = "display: none;",
        
        h4("How to Use Pathview:"),
        tags$ol(
          tags$li("Select your data using the sidebar controls (Analysis Type, Gene, Cluster, etc.)"),
          tags$li("Choose a KEGG pathway from the dropdown"),
          tags$li("Select how to color genes (fold change, p-value, or binary)"),
          tags$li("Click 'Generate Pathway' to create the diagram"),
          tags$li("View the annotated pathway with your differentially expressed genes highlighted")
        ),
        
        h4("Data Source:"),
        tags$ul(
          tags$li(strong("Real DE Data:"), " Uses actual differential expression results from your MAST/MixScale analyses with real log2 fold changes and statistical significance values.")
        ),
        
        h4("Color Schemes:"),
        tags$ul(
          tags$li(strong("Log2 Fold Change:"), " Red = higher enrichment/upregulation, Blue = lower"),
          tags$li(strong("P-value:"), " Intensity based on significance"),
          tags$li(strong("Binary:"), " Green = significant, Gray = not significant")
        ),
        
        h4("Supported File Formats:"),
        p("Upload CSV, TSV, or Excel files with columns: gene_symbol, log2FoldChange, padj")
        ) # end box-body
      ) # end box div
      ) # end column
    )
  )
}

# Server function
mod_pathview_server <- function(id, app_data, selected_enrichment_data, global_selection) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    message("=== Pathview Module Server Started ===")
    
    # Reactive values for this module
    module_data <- reactiveValues(
      de_data = NULL,
      pathway_data = NULL,
      current_pathway_file = NULL,
      available_pathways = NULL
    )
    
    # Load pathview library with fallback
    check_pathview <- reactive({
      if (!requireNamespace("pathview", quietly = TRUE)) {
        showNotification(
          "Pathview package not installed. Please install with: BiocManager::install('pathview')",
          type = "error", duration = 10
        )
        return(FALSE)
      }
      return(TRUE)
    })
    
    # Load KEGG pathway list
    observe({
      if (check_pathview()) {
        tryCatch({
          # Load KEGG pathway information using KEGGREST if available
          if (requireNamespace("KEGGREST", quietly = TRUE)) {
            # Get pathway list from KEGG API
            pathway_list <- KEGGREST::keggList("pathway", input$species)
            pathway_choices <- names(pathway_list)
            names(pathway_choices) <- paste0(pathway_choices, ": ", substr(pathway_list, 1, 60))
          } else {
            # Fallback to manual pathway list for common pathways
            common_pathways <- c(
              "hsa04010" = "MAPK signaling pathway",
              "hsa04020" = "Calcium signaling pathway", 
              "hsa04060" = "Cytokine-cytokine receptor interaction",
              "hsa04110" = "Cell cycle",
              "hsa04115" = "p53 signaling pathway",
              "hsa04120" = "Ubiquitin mediated proteolysis",
              "hsa04141" = "Protein processing in endoplasmic reticulum",
              "hsa04144" = "Endocytosis",
              "hsa04150" = "mTOR signaling pathway",
              "hsa04210" = "Apoptosis",
              "hsa04310" = "Wnt signaling pathway",
              "hsa04330" = "Notch signaling pathway",
              "hsa04340" = "Hedgehog signaling pathway",
              "hsa04350" = "TGF-beta signaling pathway",
              "hsa04360" = "Axon guidance",
              "hsa04370" = "VEGF signaling pathway",
              "hsa04510" = "Focal adhesion",
              "hsa04520" = "Adherens junction",
              "hsa04530" = "Tight junction",
              "hsa04540" = "Gap junction",
              "hsa04550" = "Signaling pathways regulating pluripotency of stem cells",
              "hsa04610" = "Complement and coagulation cascades",
              "hsa04620" = "Toll-like receptor signaling pathway",
              "hsa04630" = "Jak-STAT signaling pathway",
              "hsa04640" = "Hematopoietic cell lineage",
              "hsa04650" = "Natural killer cell mediated cytotoxicity",
              "hsa04660" = "T cell receptor signaling pathway",
              "hsa04670" = "Leukocyte transendothelial migration",
              "hsa04710" = "Circadian rhythm",
              "hsa04720" = "Long-term potentiation",
              "hsa04730" = "Long-term depression",
              "hsa04910" = "Insulin signaling pathway",
              "hsa04920" = "Adipocytokine signaling pathway",
              "hsa04930" = "Type II diabetes mellitus",
              "hsa05010" = "Alzheimer disease",
              "hsa05012" = "Parkinson disease",
              "hsa05014" = "Amyotrophic lateral sclerosis (ALS)",
              "hsa05016" = "Huntington disease",
              "hsa05020" = "Prion diseases",
              "hsa05030" = "Cocaine addiction",
              "hsa05031" = "Amphetamine addiction",
              "hsa05032" = "Morphine addiction",
              "hsa05033" = "Nicotine addiction",
              "hsa05034" = "Alcoholism",
              "hsa05200" = "Pathways in cancer",
              "hsa05210" = "Colorectal cancer",
              "hsa05212" = "Pancreatic cancer",
              "hsa05213" = "Endometrial cancer",
              "hsa05214" = "Glioma",
              "hsa05215" = "Prostate cancer",
              "hsa05216" = "Thyroid cancer",
              "hsa05217" = "Basal cell carcinoma",
              "hsa05218" = "Melanoma",
              "hsa05219" = "Bladder cancer",
              "hsa05220" = "Chronic myeloid leukemia",
              "hsa05221" = "Acute myeloid leukemia",
              "hsa05222" = "Small cell lung cancer",
              "hsa05223" = "Non-small cell lung cancer"
            )
            
            # Adjust for species
            if (input$species == "mmu") {
              pathway_choices <- gsub("hsa", "mmu", names(common_pathways))
              names(pathway_choices) <- common_pathways
            } else if (input$species == "rno") {
              pathway_choices <- gsub("hsa", "rno", names(common_pathways))
              names(pathway_choices) <- common_pathways
            } else {
              pathway_choices <- names(common_pathways)
              names(pathway_choices) <- common_pathways
            }
          }
          
          module_data$available_pathways <- pathway_choices
          updateSelectInput(session, "pathway_id", choices = pathway_choices)
          
        }, error = function(e) {
          showNotification(paste("Error loading KEGG pathways:", e$message), 
                          type = "error")
        })
      }
    })
    
    # Test button for debugging
    observeEvent(input$test_button, {
      message("=== TEST BUTTON CLICKED ===")
      showNotification("Pathview module is responding to button clicks!", type = "message")
    })
    
    # Prepare DE data from global selection
    observe({
      req(global_selection())
      selection <- global_selection()
      
      # Check that we have valid selection values
      req(selection$analysis_type)
      req(selection$gene)
      req(selection$cluster)
      
      message("=== Loading DE data from global selection ===")
      message("Analysis type: ", selection$analysis_type)
      message("Gene: ", selection$gene)
      message("Cluster: ", selection$cluster)
      
      # Load the actual DE results from full_DE_results.rds
      tryCatch({
        # Look for the DE results file
      de_env_path <- Sys.getenv("ISCORE_DE_FILE", "")
      possible_paths <- c(
        de_env_path,
        "/Users/hockemeyer/Desktop/Functional Enrichment/full_DE_results.rds",
        "/Users/hockemeyer/Desktop/Functional Enrichment/shiny_vis/full_DE_results.rds",
        file.path(dirname(Sys.getenv("ISCORE_ENRICHMENT_DIR", "")), "full_DE_results.rds")
      )
      # Remove empty paths
      possible_paths <- possible_paths[possible_paths != ""]
      
      de_file_path <- NULL
      for (path in possible_paths) {
        if (file.exists(path)) {
          de_file_path <- path
          break
        }
      }
      
      if (!is.null(de_file_path)) {
        message("Loading DE results from: ", de_file_path)
        full_de <- readRDS(de_file_path)
        
        # Navigate to the correct data based on global selection
        if (selection$analysis_type == "MAST") {
          de_results <- full_de$iSCORE_PD_MAST[[selection$gene]][[selection$cluster]]$results
        } else if (selection$analysis_type == "MixScale") {
          de_results <- full_de$CRISPRi_Mixscale[[selection$gene]][[selection$cluster]]$results
        }
        
        if (!is.null(de_results) && is.data.frame(de_results)) {
          message("Found DE results: ", nrow(de_results), " genes")
          
          # Create proper DE data frame based on analysis type
          if (selection$analysis_type == "MAST") {
            # MAST has avg_log2FC and p_val_adj columns
            de_data <- data.frame(
              gene_symbol = rownames(de_results),
              log2FoldChange = de_results$avg_log2FC,
              padj = de_results$p_val_adj,
              stringsAsFactors = FALSE
            )
          } else if (selection$analysis_type == "MixScale") {
            # MixScale has different structure - use first available log2FC and p-value columns
            log2fc_cols <- grep("^log2FC_", names(de_results), value = TRUE)
            p_cols <- grep("^p_", names(de_results), value = TRUE)
            
            if (length(log2fc_cols) > 0 && length(p_cols) > 0) {
              # Use the first log2FC column and first p-value column
              de_data <- data.frame(
                gene_symbol = if("gene_ID" %in% names(de_results)) de_results$gene_ID else rownames(de_results),
                log2FoldChange = de_results[[log2fc_cols[1]]],
                padj = de_results[[p_cols[1]]],
                stringsAsFactors = FALSE
              )
            } else {
              stop("Could not find log2FC or p-value columns in MixScale data")
            }
          }
          
          message("DE data created with ", nrow(de_data), " genes")
          message("Log2FC range: ", round(min(de_data$log2FoldChange, na.rm = TRUE), 2), 
                 " to ", round(max(de_data$log2FoldChange, na.rm = TRUE), 2))
          
          sig_genes <- sum(de_data$padj <= 0.05, na.rm = TRUE)
          message("Significant genes (padj <= 0.05): ", sig_genes)
          
          module_data$de_data <- de_data
          showNotification(paste("Loaded real DE data:", nrow(de_data), "genes,", sig_genes, "significant"), type = "message")
          
        } else {
          message("Could not find DE results for this condition")
          showNotification("No DE results found for this analysis condition", type = "error")
        }
        
      } else {
        message("Could not find full_DE_results.rds file")
        showNotification("Could not locate DE results file", type = "warning")
      }
      
      }, error = function(e) {
        message("Error loading DE data: ", e$message)
        showNotification(paste("Error loading DE data:", e$message), type = "error")
      })
    })
    
    # Generate pathway diagram
    observeEvent(input$generate_pathway, {
      message("=== Generate Pathway Button Clicked ===")
      message("Pathway ID: ", input$pathway_id)
      message("Using global data selection")
      message("DE data available: ", !is.null(module_data$de_data))
      message("Pathview available: ", check_pathview())
      
      # Check requirements with detailed messages
      if (is.null(input$pathway_id) || input$pathway_id == "") {
        showNotification("Please select a KEGG pathway first", type = "error")
        return()
      }
      
      if (is.null(module_data$de_data)) {
        showNotification("No differential expression data available. Please load enrichment results first.", type = "error")
        return()
      }
      
      if (!check_pathview()) {
        showNotification("Pathview package not available", type = "error")
        return()
      }
      
      message("All requirements met, proceeding with pathway generation...")
      
      withProgress(message = 'Generating pathway diagram...', value = 0, {
        
        tryCatch({
          incProgress(0.2, detail = "Preparing gene data...")
          
          # Prepare gene data based on selected type
          de_data <- module_data$de_data
          
          # Filter genes based on global significance threshold
          selection <- global_selection()
          sig_genes <- de_data[de_data$padj <= selection$pval_threshold, ]
          
          message("Total DE genes: ", nrow(de_data))
          message("Significant genes (p <= ", selection$pval_threshold, "): ", nrow(sig_genes))
          
          if (nrow(sig_genes) == 0) {
            showNotification(paste("No significant genes found with threshold", selection$pval_threshold), 
                            type = "warning")
            return()
          }
          
          # Prepare gene data vector for pathview
          if (input$gene_data_type == "logfc") {
            gene_data <- sig_genes$log2FoldChange
            names(gene_data) <- sig_genes$gene_symbol
          } else if (input$gene_data_type == "pvalue") {
            gene_data <- -log10(sig_genes$padj)
            names(gene_data) <- sig_genes$gene_symbol
          } else if (input$gene_data_type == "binary") {
            gene_data <- ifelse(abs(sig_genes$log2FoldChange) >= input$logfc_threshold, 1, 0)
            names(gene_data) <- sig_genes$gene_symbol
          }
          
          message("Gene data prepared for ", length(gene_data), " genes")
          message("Gene symbols: ", paste(head(names(gene_data), 10), collapse = ", "))
          
          incProgress(0.4, detail = "Converting gene symbols...")
          
          # Enhanced gene symbol to ENTREZ ID conversion
          if (requireNamespace("clusterProfiler", quietly = TRUE) && 
              requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
            
            # Clean gene symbols - remove any special characters and ensure proper format
            gene_symbols <- names(gene_data)
            gene_symbols_clean <- toupper(trimws(gene_symbols))
            
            message("Starting gene mapping for ", length(gene_symbols_clean), " genes")
            message("Sample gene symbols: ", paste(head(gene_symbols_clean, 5), collapse = ", "))
            
            # Primary mapping attempt
            gene_mapping <- tryCatch({
              clusterProfiler::bitr(
                gene_symbols_clean, 
                fromType = "SYMBOL", 
                toType = "ENTREZID", 
                OrgDb = org.Hs.eg.db::org.Hs.eg.db
              )
            }, error = function(e) {
              message("Primary mapping failed: ", e$message)
              data.frame(SYMBOL = character(0), ENTREZID = character(0))
            })
            
            message("Primary mapping result: ", nrow(gene_mapping), " genes mapped")
            
            # Try alternative mapping for unmapped genes
            if (nrow(gene_mapping) < length(gene_symbols_clean) * 0.5) {
              message("Low mapping success, trying alternative approaches...")
              
              # Try ALIAS mapping for unmapped genes
              unmapped_symbols <- setdiff(gene_symbols_clean, gene_mapping$SYMBOL)
              
              if (length(unmapped_symbols) > 0) {
                alias_mapping <- tryCatch({
                  clusterProfiler::bitr(
                    unmapped_symbols, 
                    fromType = "ALIAS", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db::org.Hs.eg.db
                  )
                }, error = function(e) {
                  message("Alias mapping failed: ", e$message)
                  data.frame(ALIAS = character(0), ENTREZID = character(0))
                })
                
                if (nrow(alias_mapping) > 0) {
                  # Rename to match primary mapping
                  names(alias_mapping)[1] <- "SYMBOL"
                  gene_mapping <- rbind(gene_mapping, alias_mapping)
                  message("Additional genes mapped via ALIAS: ", nrow(alias_mapping))
                }
              }
              
              # Try ENSEMBL mapping if available
              if (length(unmapped_symbols) > 0) {
                ensemble_mapping <- tryCatch({
                  clusterProfiler::bitr(
                    unmapped_symbols, 
                    fromType = "ENSEMBL", 
                    toType = "ENTREZID", 
                    OrgDb = org.Hs.eg.db::org.Hs.eg.db
                  )
                }, error = function(e) {
                  message("ENSEMBL mapping failed: ", e$message)
                  data.frame(ENSEMBL = character(0), ENTREZID = character(0))
                })
                
                if (nrow(ensemble_mapping) > 0) {
                  names(ensemble_mapping)[1] <- "SYMBOL"
                  gene_mapping <- rbind(gene_mapping, ensemble_mapping)
                  message("Additional genes mapped via ENSEMBL: ", nrow(ensemble_mapping))
                }
              }
            }
            
            # Remove duplicates and ensure unique mapping
            gene_mapping <- gene_mapping[!duplicated(gene_mapping$SYMBOL), ]
            gene_mapping <- gene_mapping[!duplicated(gene_mapping$ENTREZID), ]
            
            message("Final mapping result: ", nrow(gene_mapping), " genes mapped")
            message("Mapping success rate: ", round(nrow(gene_mapping)/length(gene_symbols_clean)*100, 1), "%")
            
            if (nrow(gene_mapping) == 0) {
              showNotification("No genes could be mapped to ENTREZ IDs", type = "error")
              return()
            }
            
            # Create mapping between original symbols and clean symbols
            symbol_map <- data.frame(
              original = gene_symbols,
              clean = gene_symbols_clean,
              stringsAsFactors = FALSE
            )
            
            # Map gene data to ENTREZ IDs using clean symbols
            mapped_gene_data <- merge(
              data.frame(clean = names(gene_data), value = gene_data, stringsAsFactors = FALSE),
              symbol_map,
              by = "clean"
            )
            
            # Now map to ENTREZ IDs
            entrez_mapping <- merge(
              mapped_gene_data,
              gene_mapping,
              by.x = "clean",
              by.y = "SYMBOL"
            )
            
            if (nrow(entrez_mapping) == 0) {
              showNotification("No genes matched between data and mapping", type = "error")
              return()
            }
            
            # Create final ENTREZ data vector
            entrez_data <- entrez_mapping$value
            names(entrez_data) <- entrez_mapping$ENTREZID
            
            message("Final ENTREZ data: ", length(entrez_data), " genes")
            message("Sample ENTREZ IDs: ", paste(head(names(entrez_data), 5), collapse = ", "))
            
          } else {
            showNotification("Gene annotation packages not available", type = "error")
            return()
          }
          
          incProgress(0.6, detail = "Generating pathway diagram...")
          
          # Validate entrez_data before pathview
          if (length(entrez_data) == 0) {
            showNotification("No genes successfully mapped to ENTREZ IDs", type = "error")
            return()
          }
          
          # Remove any NA or infinite values
          entrez_data <- entrez_data[!is.na(entrez_data) & is.finite(entrez_data)]
          
          if (length(entrez_data) == 0) {
            showNotification("No valid gene values after cleaning", type = "error")
            return()
          }
          
          message("Final gene data for pathview:")
          message("- Number of genes: ", length(entrez_data))
          message("- Value range: ", min(entrez_data, na.rm = TRUE), " to ", max(entrez_data, na.rm = TRUE))
          message("- Sample ENTREZ IDs: ", paste(head(names(entrez_data), 3), collapse = ", "))
          message("- Sample values: ", paste(head(entrez_data, 3), collapse = ", "))
          
          # Extract pathway ID
          pathway_id <- gsub("^[a-z]{3}([0-9]+).*", "\\1", input$pathway_id)
          message("Using pathway ID: ", pathway_id)
          
          # Generate pathway with pathview
          pathway_result <- pathview::pathview(
            gene.data = entrez_data,
            pathway.id = pathway_id,
            species = input$species,
            out.suffix = paste0("_", format(Sys.time(), "%Y%m%d_%H%M%S")),
            kegg.native = TRUE,
            same.layer = FALSE,
            limit = list(gene = c(min(entrez_data, na.rm = TRUE), max(entrez_data, na.rm = TRUE))),
            low = c("blue"),
            mid = c("white"), 
            high = c("red")
          )
          
          message("Pathview completed")
          message("Pathway result object: ", class(pathway_result))
          
          incProgress(0.8, detail = "Processing results...")
          
          # Find generated pathway file
          pathway_files <- list.files(pattern = paste0("^", input$species, pathway_id, ".*\\.", input$output_format, "$"))
          
          if (length(pathway_files) > 0) {
            module_data$current_pathway_file <- pathway_files[1]
            module_data$pathway_data <- list(
              pathway_id = pathway_id,
              gene_data = gene_data,
              entrez_data = entrez_data,
              gene_mapping = gene_mapping,
              entrez_mapping = entrez_mapping,
              mapping_stats = list(
                total_input = length(gene_symbols_clean),
                mapped_genes = nrow(gene_mapping),
                final_entrez = length(entrez_data),
                success_rate = round(nrow(gene_mapping)/length(gene_symbols_clean)*100, 1)
              )
            )
            
            showNotification("Pathway diagram generated successfully!", type = "message")
          } else {
            showNotification("No pathway file generated", type = "error")
          }
          
          incProgress(1.0, detail = "Complete!")
          
        }, error = function(e) {
          showNotification(paste("Error generating pathway:", e$message), type = "error")
        })
      })
    })
    
    # Display pathway image
    output$pathway_image <- renderUI({
      req(module_data$current_pathway_file)
      
      if (file.exists(module_data$current_pathway_file)) {
        if (input$output_format == "png") {
          tags$img(
            src = base64enc::dataURI(file = module_data$current_pathway_file, mime = "image/png"),
            style = "max-width: 100%; height: auto;"
          )
        } else {
          tags$p("PDF generated. Use download button to view.")
        }
      } else {
        tags$p("No pathway image available", style = "text-align: center; color: gray;")
      }
    })
    
    # Gene mapping table
    output$gene_mapping_table <- DT::renderDataTable({
      req(module_data$pathway_data)
      
      # Use enhanced mapping data if available
      if (!is.null(module_data$pathway_data$entrez_mapping)) {
        result_table <- module_data$pathway_data$entrez_mapping[, c("original", "clean", "ENTREZID", "value")]
        colnames(result_table) <- c("Original Symbol", "Clean Symbol", "ENTREZ ID", paste("Value (", input$gene_data_type, ")"))
      } else {
        # Fallback to old method
        mapping_data <- module_data$pathway_data$gene_mapping
        gene_data <- module_data$pathway_data$gene_data
        
        result_table <- merge(
          mapping_data,
          data.frame(
            SYMBOL = names(gene_data),
            Value = gene_data,
            stringsAsFactors = FALSE
          ),
          by = "SYMBOL"
        )
        colnames(result_table) <- c("Gene Symbol", "ENTREZ ID", paste("Value (", input$gene_data_type, ")"))
      }
      
      # Add mapping statistics at the top
      stats <- module_data$pathway_data$mapping_stats
      if (!is.null(stats)) {
        stats_row <- data.frame(
          V1 = paste("MAPPING STATS - Success Rate:", stats$success_rate, "%"),
          V2 = paste("Input:", stats$total_input),
          V3 = paste("Mapped:", stats$mapped_genes),
          V4 = paste("Final:", stats$final_entrez),
          stringsAsFactors = FALSE
        )
        names(stats_row) <- names(result_table)
        result_table <- rbind(stats_row, result_table)
      }
      
      DT::datatable(
        result_table,
        options = list(pageLength = 25, scrollX = TRUE),
        filter = 'top',
        caption = "Gene Symbol to ENTREZ ID Mapping Results"
      ) %>%
        DT::formatRound(columns = ncol(result_table), digits = 3)
    })
    
    # Pathway info
    output$pathway_info <- renderPrint({
      cat("=== Pathview Module Status ===\n")
      cat("Module loaded:", TRUE, "\n")
      cat("Data source:", input$data_source, "\n")
      cat("Pathview available:", check_pathview(), "\n")
      cat("DE data loaded:", !is.null(module_data$de_data), "\n")
      
      if (!is.null(input$pathway_id) && input$pathway_id != "") {
        cat("\nSelected Pathway:\n")
        cat("================\n")
        cat("ID:", input$pathway_id, "\n")
        cat("Species:", input$species, "\n")
        cat("Gene Data Type:", input$gene_data_type, "\n")
        cat("P-value Threshold:", input$pval_threshold, "\n")
        
        if (!is.null(module_data$de_data)) {
          sig_genes <- sum(module_data$de_data$padj <= input$pval_threshold, na.rm = TRUE)
          cat("Significant Genes:", sig_genes, "\n")
          cat("Total Genes:", nrow(module_data$de_data), "\n")
        } else {
          cat("No DE data available\n")
        }
        
        if (!is.null(module_data$pathway_data)) {
          cat("Mapped Genes:", length(module_data$pathway_data$entrez_data), "\n")
          
          if (!is.null(module_data$pathway_data$mapping_stats)) {
            stats <- module_data$pathway_data$mapping_stats
            cat("\nGene Mapping Statistics:\n")
            cat("======================\n")
            cat("Input genes:", stats$total_input, "\n")
            cat("Successfully mapped:", stats$mapped_genes, "\n")
            cat("Final ENTREZ genes:", stats$final_entrez, "\n")
            cat("Success rate:", stats$success_rate, "%\n")
          }
        }
      } else {
        cat("\nNo pathway selected\n")
      }
    })
    
    # Download handler
    output$download_pathway <- downloadHandler(
      filename = function() {
        if (!is.null(module_data$current_pathway_file)) {
          basename(module_data$current_pathway_file)
        } else {
          paste0("pathway_", Sys.Date(), ".", input$output_format)
        }
      },
      content = function(file) {
        if (!is.null(module_data$current_pathway_file) && 
            file.exists(module_data$current_pathway_file)) {
          file.copy(module_data$current_pathway_file, file)
        }
      }
    )
    
    return(module_data)
  })
}