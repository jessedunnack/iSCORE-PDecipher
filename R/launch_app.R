#' Launch iSCORE-PDecipher Shiny Application
#'
#' This function launches the interactive Shiny application for exploring
#' Parkinson's disease mutation and perturbation enrichment results.
#'
#' @param data_file Path to the RDS file containing consolidated enrichment results.
#'   If NULL, the app will prompt for file upload.
#' @param port Port number for the Shiny app (default: auto-select)
#' @param launch.browser Whether to launch the app in browser (default: TRUE)
#' @param ... Additional arguments passed to shiny::runApp
#'
#' @return Launches the Shiny application
#' @export
#'
#' @examples
#' \dontrun{
#' # Launch with data file
#' launch_iscore_app("path/to/enrichment_results.rds")
#' 
#' # Launch with file upload interface
#' launch_iscore_app()
#' }
launch_iscore_app <- function(data_file = NULL, port = getOption("shiny.port"), 
                              launch.browser = getOption("shiny.launch.browser", TRUE), ...) {
  
  # Get the path to the Shiny app
  app_dir <- system.file("shiny", package = "iSCORE.PDecipher")
  
  if (app_dir == "") {
    stop("Could not find Shiny app. Try re-installing `iSCORE.PDecipher`.", call. = FALSE)
  }
  
  # Set environment variable for data file if provided
  if (!is.null(data_file)) {
    if (!file.exists(data_file)) {
      stop("Data file not found: ", data_file, call. = FALSE)
    }
    Sys.setenv(ISCORE_DATA_FILE = data_file)
  }
  
  # Launch the app
  shiny::runApp(appDir = app_dir, 
                port = port, 
                launch.browser = launch.browser, 
                ...)
}

#' @rdname launch_iscore_app
#' @export
launch_app <- launch_iscore_app