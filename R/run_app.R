#' Run iSCORE-PDecipher Shiny Application
#'
#' @description
#' Launches the iSCORE-PDecipher Shiny application for interactive analysis
#' of Parkinson's Disease single-cell RNA-seq data including MAST and MixScale
#' differential expression results.
#'
#' @param app_type Character string specifying which app to launch:
#'   - "main" (default): Full app with dataset selector (recommended)
#'   - "basic": Basic app without dataset selector
#'   - "perturbseq": Perturb-seq only module (v0.5.0+)
#' @param launch.browser Logical, whether to launch in browser (default TRUE)
#' @param port Integer, port number for the app (default NULL for automatic)
#' @param host Character, host IP address (default "127.0.0.1")
#' @param ... Additional arguments passed to shiny::runApp()
#'
#' @return Starts the Shiny application
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch the main app
#' run_app()
#'
#' # Launch the Perturb-seq only module
#' run_app(app_type = "perturbseq")
#'
#' # Launch on specific port
#' run_app(port = 8080)
#' }
run_app <- function(app_type = "main",
                    launch.browser = TRUE,
                    port = NULL,
                    host = "127.0.0.1",
                    ...) {

  # Validate app_type
  valid_types <- c("main", "basic", "perturbseq")
  if (!app_type %in% valid_types) {
    stop("app_type must be one of: ", paste(valid_types, collapse = ", "))
  }

  # Determine which app file to use
  app_file <- switch(app_type,
    "main" = "app_with_dataset_selector.R",
    "basic" = "app.R",
    "perturbseq" = NULL  # Special case handled below
  )

  # Handle Perturb-seq module separately
  if (app_type == "perturbseq") {
    message("Launching Perturb-seq Only Module (v0.5.0)...")

    # Get package installation directory
    app_dir <- system.file("shiny", package = "iSCORE.PDecipher")

    if (app_dir == "") {
      stop("Could not find iSCORE.PDecipher installation directory")
    }

    # Source required files
    module_file <- file.path(app_dir, "modules/mod_perturbseq_only.R")
    data_manager_file <- file.path(app_dir, "R/data_manager.R")

    if (!file.exists(module_file)) {
      stop("Perturb-seq module not found at: ", module_file)
    }

    source(module_file)
    source(data_manager_file)

    # Create standalone app for the module
    ui <- shiny::fluidPage(
      title = "iSCORE-PDecipher: Perturb-seq Analysis",
      mod_perturbseq_only_ui("perturbseq")
    )

    server <- function(input, output, session) {
      mod_perturbseq_only_server("perturbseq")
    }

    # Launch the app
    app <- shiny::shinyApp(ui, server)
    shiny::runApp(app,
                  launch.browser = launch.browser,
                  port = port,
                  host = host,
                  ...)

  } else {
    # Launch main or basic app
    message("Launching iSCORE-PDecipher Shiny App (", app_type, ")...")

    # Get the app directory
    app_dir <- system.file("shiny", package = "iSCORE.PDecipher")

    if (app_dir == "") {
      stop("Could not find iSCORE.PDecipher installation directory")
    }

    # Full path to app file
    app_path <- file.path(app_dir, app_file)

    if (!file.exists(app_path)) {
      stop("App file not found at: ", app_path)
    }

    # Set working directory to app directory for proper sourcing
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(app_dir)

    # Launch the app
    shiny::runApp(app_file,
                  launch.browser = launch.browser,
                  port = port,
                  host = host,
                  ...)
  }
}
