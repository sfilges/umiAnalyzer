#' Function to run the umiVisualiser shiny app
#' @export
#' @importFrom shiny runApp
#' @param docker Boolean. Are you running in a docker container?
runUmiVisualiser <- function(docker=FALSE) {
  appDir <- system.file("shiny", package = "umiAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umiAnalyzer`.", call. = FALSE)
  }

  if(docker == FALSE){
    shiny::runApp(appDir, display.mode = "normal")
  } else {
    shiny::runApp(appDir, launch.browser = FALSE, port = 8080, host = "0.0.0.0")
  }

}
