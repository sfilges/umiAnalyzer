#' Function to run the umiVisualiser shiny app
#' @export
#' @importFrom shiny runApp
runUmiVisualiser <- function() {
  appDir <- system.file("shiny", package = "umiAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umiAnalyzer`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
