#' Function to run the umiVisualizer shiny app
#' @export
#' 
#' @importFrom shiny runApp
#' @importFrom utils globalVariables
#' 
#' @param docker Boolean. Are you running in a docker container?
#' @param path Path to directory containing UMIErrorCorrect output.
runUmiVisualiser <- function(docker=FALSE, path=NULL) {
  appDir <- system.file("shiny", package = "umiAnalyzer")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `umiAnalyzer`.", call. = FALSE)
  }
  
  #utils::globalVariables(names("path_to_umierrorcorrect_data"), add = FALSE)

  if(docker == FALSE){
    # Set path as global variable
    
    .GlobalEnv$path_to_umierrorcorrect_data <- path
    on.exit(rm(.GlobalEnv$path_to_umierrorcorrect_data, envir=.GlobalEnv))

    if(!is.null(path)) {
      print(paste("Loading data from: ", .GlobalEnv$path_to_umierrorcorrect_data, sep =""))
    } else {
      print("Start app without loading data.")
    }

    # run app
    shiny::runApp(appDir, display.mode = "normal")
  } else {
    # run shiny app
    shiny::runApp(appDir, launch.browser = FALSE, port = 8080, host = "0.0.0.0")
  }
}
