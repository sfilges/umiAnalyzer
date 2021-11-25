#' Function to run the umiVisualizer shiny app
#' @export
#' 
#' @import shinyFiles
#' @import shinyWidgets
#' @import shinydashboard
#' @import DT
#' 
#' @importFrom shiny runApp
#' 
#' @return Opens the umiVisualizer app
#' 
#' @examples 
#' \dontrun{
#' library(umiAnalyzer)
#' 
#' runUmiVisualizer()
#' }
#' 
runUmiVisualizer <- function() {
  appDir <- system.file("shiny", package = "umiAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umiAnalyzer`.", call. = FALSE)
  }

  # run app
  shiny::runApp(appDir, display.mode = "normal")
}
