\name{launch_iscore_app}
\alias{launch_iscore_app}
\alias{launch_app}
\title{Launch iSCORE-PDecipher Shiny Application}
\usage{
launch_iscore_app(
  data_file = NULL,
  port = getOption("shiny.port"),
  launch.browser = getOption("shiny.launch.browser", TRUE),
  ...
)

launch_app(
  data_file = NULL,
  port = getOption("shiny.port"),
  launch.browser = getOption("shiny.launch.browser", TRUE),
  ...
)
}
\arguments{
\item{data_file}{Path to the RDS file containing consolidated enrichment results.
If NULL, the app will prompt for file upload.}

\item{port}{Port number for the Shiny app (default: auto-select)}

\item{launch.browser}{Whether to launch the app in browser (default: TRUE)}

\item{...}{Additional arguments passed to shiny::runApp}
}
\value{
Launches the Shiny application
}
\description{
This function launches the interactive Shiny application for exploring
Parkinson's disease mutation and perturbation enrichment results.
}
\examples{
\dontrun{
# Launch with data file
launch_iscore_app("path/to/enrichment_results.rds")

# Launch with file upload interface
launch_iscore_app()
}
}