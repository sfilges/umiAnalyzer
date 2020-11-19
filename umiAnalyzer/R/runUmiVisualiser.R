#' Function to run the umiVisualiser shiny app
#' @export
#' @importFrom shiny runApp
#' @param docker Boolean. Are you running in a docker container?
#' @param path Path to directory containing umierrorocorrect output.
runUmiVisualiser <- function(docker=FALSE, path=NULL) {
  appDir <- system.file("shiny", package = "umiAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umiAnalyzer`.", call. = FALSE)
  }

  if(docker == FALSE){
    # Set path as global variable
    .GlobalEnv$path_to_umierrorcorrect_data <- path
    on.exit(rm(path_to_umierrorcorrect_data, envir=.GlobalEnv))

    print(path_to_umierrorcorrect_data)
    if(!is.null(path)) {
      print(paste("Loading data from: ", path_to_umierrorcorrect_data, sep =""))
    }

    # run app
    shiny::runApp(appDir, display.mode = "normal")
  } else {
    # run shiny app
    shiny::runApp(appDir, launch.browser = FALSE, port = 8080, host = "0.0.0.0")
  }
}
