# Module: Enhanced Visualization with Gene Display
# Advanced visualizations with integrated gene lists in both tables and tooltips

mod_visualization_ui <- function(id) {
  ns <- NS(id)
  
  sidebarLayout(
    sidebarPanel(
      width = 3,
      h4("Visualization Settings"),
      
      # Conditional UI based on enrichment type
      conditionalPanel(
        condition = "output.is_gsea == false",
        ns = ns,
        
        selectInput(ns("plot_type"),
                    "Visualization Type:",
                    choices = list(
                      "Dot Plot" = "dotplot",
                      "Bar Plot" = "barplot",
                      "Lollipop Plot" = "lollipop"
                    ),
                    selected = "dotplot"),
        
        numericInput(ns("top_terms"),
                     "Number of Top Terms:",
                     value = 20,
                     min = 5,
                     max = 50,
                     step = 5),
        
        br(),
        
        selectInput(ns("x_axis"),
                    "X-axis:",
                    choices = c("-log10(p.adjust)" = "neg_log10_pval",
                              "Fold Enrichment" = "FoldEnrichment", 
                              "Gene Count" = "Count",
                              "Gene Ratio" = "GeneRatio",
                              "Rich Factor" = "RichFactor"),
                    selected = "neg_log10_pval"),
        
        selectInput(ns("color_by"),
                    "Color By:",
                    choices = c("-log10(adjusted p-value)" = "p-value", 
                               "Fold Enrichment" = "Fold Enrichment"),
                    selected = "p-value"),
        
        checkboxInput(ns("show_labels"),
                      "Show Labels",
                      value = TRUE),
        
        checkboxInput(ns("show_genes_in_hover"),
                      "Show Genes in Hover",
                      value = TRUE)
      ),
      
      # GSEA-specific settings
      conditionalPanel(
        condition = "output.is_gsea == true",
        ns = ns,
        
        selectInput(ns("gsea_plot_type"),
                    "GSEA Plot Type:",
                    choices = list(
                      "Enrichment Plot" = "gseaplot",
                      "Dot Plot" = "dotplot",
                      "Ridge Plot" = "ridgeplot",
                      "GSEA Table" = "table"
                    ),
                    selected = "gseaplot"),
        
        conditionalPanel(
          condition = "input.gsea_plot_type == 'gseaplot'",
          ns = ns,
          
          selectInput(ns("gsea_term_select"),
                      "Select Gene Set:",
                      choices = NULL),
          
          checkboxInput(ns("gsea_show_pval"),
                        "Show p-value on plot",
                        value = TRUE),
          
          checkboxInput(ns("gsea_show_running"),
                        "Show running score",
                        value = TRUE),
          
          sliderInput(ns("gsea_base_size"),
                      "Base font size:",
                      min = 8,
                      max = 16,
                      value = 11,
                      step = 1)
        ),
        
        conditionalPanel(
          condition = "input.gsea_plot_type != 'gseaplot'",
          ns = ns,
          
          numericInput(ns("gsea_top_terms"),
                       "Number of Top Terms:",
                       value = 20,
                       min = 5,
                       max = 50,
                       step = 5)
        )
      ),
      
      br(),
      helpText("Plot updates automatically when you change global settings"),
      helpText("Gene lists are shown in Plot Details table and hover tooltips (if enabled)"),
      
      br(),
      downloadButton(ns("download_plot"), "Download Plot",
                     class = "btn-success", style = "width: 100%;")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        id = ns("main_tabs"),
        selected = "primary",
        
        tabPanel("Primary Plot",
                 value = "primary",
                 br(),
                 conditionalPanel(
                   condition = "output.is_gsea == false",
                   ns = ns,
                   plotlyOutput(ns("interactive_plot"), height = "700px") %>%
                     withSpinner(color = "#0073b7")
                 ),
                 conditionalPanel(
                   condition = "output.is_gsea == true",
                   ns = ns,
                   plotOutput(ns("gsea_plot"), height = "700px") %>%
                     withSpinner(color = "#0073b7")
                 )),
        
        tabPanel("Plot Details",
                 value = "details",
                 br(),
                 div(
                   style = "background-color: #e8f4fd; padding: 10px; border-radius: 5px; margin-bottom: 15px;",
                   icon("info-circle", style = "color: #0073b7;"),
                   strong(" Gene Lists: ", style = "color: #0073b7;"),
                   "Each enrichment term shows its associated DE genes in the table below."
                 ),
                 verbatimTextOutput(ns("plot_info")),
                 br(),
                 DT::dataTableOutput(ns("plot_data")))
      )
    )
  )
}
