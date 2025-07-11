# Script to integrate gene display into the Plot Details data table
# This will add gene lists directly to the enrichment terms table

library(stringr)

cat("=== Integrating Gene Display into Plot Details Table ===\n")

# Read the visualization module
viz_file <- "inst/shiny/modules/mod_visualization_enhanced.R"
viz_content <- readLines(viz_file)

# Find the renderDataTable section (around line 663)
render_start <- which(str_detect(viz_content, "output\\$plot_data <- DT::renderDataTable"))
if (length(render_start) == 0) {
  stop("Could not find renderDataTable section")
}

cat("Found renderDataTable at line", render_start, "\n")

# Create the updated renderDataTable code with gene integration
new_render_code <- c(
  '    output$plot_data <- DT::renderDataTable({',
  '      result <- processed_data()',
  '      req(result$data)',
  '      ',
  '      # Load gene associations if available',
  '      gene_data_available <- FALSE',
  '      tryCatch({',
  '        if (!exists(".gene_associations")) {',
  '          source("R/gene_association_lookup.R")',
  '          load_gene_associations()',
  '        }',
  '        gene_data_available <- gene_associations_available()',
  '      }, error = function(e) {',
  '        message("Gene associations not available: ", e$message)',
  '      })',
  '      ',
  '      if (result$is_gsea) {',
  '        # GSEA-specific columns',
  '        display_data <- result$data %>%',
  '          select(ID, Description, NES, p.adjust, setSize, any_of(c("enrichmentScore", "rank"))) %>%',
  '          arrange(desc(abs(NES)))',
  '      } else {',
  '        # Standard enrichment columns',
  '        display_data <- result$data %>%',
  '          select(ID, Description, p.adjust, Count, any_of(c("FoldEnrichment", "GeneRatio"))) %>%',
  '          arrange(p.adjust)',
  '      }',
  '      ',
  '      # Add gene lists if available',
  '      if (gene_data_available && !is.null(global_selection())) {',
  '        selection <- global_selection()',
  '        ',
  '        # Function to get genes for a term',
  '        get_gene_list <- function(term_id) {',
  '          result <- get_genes_for_term(',
  '            term_id = term_id,',
  '            analysis_type = selection$analysis_type,',
  '            gene = selection$gene,',
  '            cluster = selection$cluster,',
  '            enrichment_type = selection$enrichment_type,',
  '            direction = selection$direction,',
  '            experiment = "default"',
  '          )',
  '          ',
  '          if (!is.null(result$genes) && length(result$genes) > 0) {',
  '            # Limit display to first 20 genes for space',
  '            genes_display <- head(result$genes, 20)',
  '            if (length(result$genes) > 20) {',
  '              genes_display <- c(genes_display, paste("...", length(result$genes) - 20, "more"))',
  '            }',
  '            return(paste(genes_display, collapse = ", "))',
  '          } else {',
  '            return("No genes found")',
  '          }',
  '        }',
  '        ',
  '        # Add gene column',
  '        display_data$Associated_Genes <- sapply(display_data$ID, get_gene_list)',
  '      }',
  '      ',
  '      # Create datatable with enhanced features',
  '      dt <- DT::datatable(',
  '        display_data,',
  '        options = list(',
  '          pageLength = 20,',
  '          scrollX = TRUE,',
  '          columnDefs = list(',
  '            list(width = "200px", targets = which(names(display_data) == "Description") - 1),',
  '            list(width = "300px", targets = which(names(display_data) == "Associated_Genes") - 1)',
  '          )',
  '        ),',
  '        rownames = FALSE,',
  '        selection = "single"  # Allow row selection',
  '      ) %>%',
  '        DT::formatSignif(columns = "p.adjust", digits = 3)',
  '      ',
  '      # Format gene column if it exists',
  '      if ("Associated_Genes" %in% names(display_data)) {',
  '        dt <- dt %>%',
  '          DT::formatStyle("Associated_Genes",',
  '                         fontSize = "90%",',
  '                         color = "#007bff")',
  '      }',
  '      ',
  '      return(dt)',
  '    })'
)

# Find the end of the current renderDataTable block
current_end <- render_start
brace_count <- 0
for (i in render_start:length(viz_content)) {
  line <- viz_content[i]
  brace_count <- brace_count + str_count(line, "\\{") - str_count(line, "\\}")
  if (brace_count == 0 && i > render_start) {
    current_end <- i
    break
  }
}

cat("Current renderDataTable ends at line", current_end, "\n")

# Replace the section
viz_content <- c(
  viz_content[1:(render_start - 1)],
  new_render_code,
  viz_content[(current_end + 1):length(viz_content)]
)

# Add reactive for selected term (for additional gene display panel if desired)
# Find where to insert it (before the return statement)
return_line <- which(str_detect(viz_content, "return\\(processed_data\\)"))
if (length(return_line) > 0) {
  insert_at <- return_line[1] - 1
  
  selected_term_code <- c(
    '    # Reactive for selected term from data table',
    '    selected_term <- reactive({',
    '      req(input$plot_data_rows_selected)',
    '      result <- processed_data()',
    '      req(result$data)',
    '      ',
    '      selected_row <- input$plot_data_rows_selected',
    '      if (length(selected_row) > 0) {',
    '        term_data <- result$data[selected_row, ]',
    '        return(list(',
    '          id = term_data$ID,',
    '          description = term_data$Description',
    '        ))',
    '      }',
    '      return(NULL)',
    '    })',
    '    '
  )
  
  viz_content <- c(
    viz_content[1:insert_at],
    selected_term_code,
    viz_content[(insert_at + 1):length(viz_content)]
  )
}

# Write the updated file
writeLines(viz_content, viz_file)

cat("✅ Updated visualization module with integrated gene display\n")

# Also update the UI to show instructions
ui_update_line <- which(str_detect(viz_content, 'tabPanel\\("Plot Details"'))
if (length(ui_update_line) > 0) {
  cat("Consider adding instructions to the Plot Details tab\n")
}

cat("\n=== Integration Complete ===\n")
cat("The Plot Details table now includes an 'Associated_Genes' column\n")
cat("Users will see gene lists directly in the data table\n")
cat("Table is scrollable horizontally to accommodate gene lists\n")

# Create a test script
test_script <- '# Test Gene Display in Data Table
# Run this after restarting the Shiny app

cat("=== Testing Gene Display in Plot Details ===\\n")
cat("1. Launch the Shiny app\\n")
cat("2. Go to Functional Enrichment Results\\n")
cat("3. Click on the Plot Details tab\\n")
cat("4. Check for Associated_Genes column in the table\\n")
cat("5. Verify genes are displayed for each term\\n")
cat("\\nExpected behavior:\\n")
cat("- Table should have an Associated_Genes column\\n")
cat("- Each term shows its DE genes (up to 20 + count of remaining)\\n")
cat("- Gene text is blue and slightly smaller\\n")
cat("- Table is horizontally scrollable\\n")
'

writeLines(test_script, "test_gene_display_table.R")
cat("\nCreated test script: test_gene_display_table.R\n")