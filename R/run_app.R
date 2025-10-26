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

  # Handle Perturb-seq full app separately
  if (app_type == "perturbseq") {
    message("Launching Perturb-seq Full-Featured App (v0.5.0)...")

    # Load required packages
    if (!requireNamespace("shiny", quietly = TRUE)) {
      stop("Package 'shiny' is required but not installed.")
    }

    # Get package installation directory
    app_dir <- system.file("shiny", package = "iSCORE.PDecipher")

    if (app_dir == "") {
      stop("Could not find iSCORE.PDecipher installation directory")
    }

    # Check for comprehensive app file
    app_file <- file.path(app_dir, "app_perturbseq_full.R")

    if (!file.exists(app_file)) {
      stop("Perturb-seq app not found at: ", app_file)
    }

    # Set working directory to app directory (required for sourcing files)
    old_wd <- getwd()
    on.exit(setwd(old_wd))
    setwd(app_dir)

    # Run the comprehensive app
    shiny::runApp(app_file,
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
