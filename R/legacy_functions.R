#--------------------------// Functionality for analyzing Debarcer output //------------------------

#' Define sample class
#' 
#' @noRd
#' 
DebarcerSample <- setClass(
  "DebarcerSample",
  slots = list(
    name = "character",
    cons.data = "tbl_df"
  )
)

#' Define experiment class
#' 
#' @import tibble
#' 
#' @noRd
#' 
DebarcerExperiment <- setClass(
  "DebarcerExperiment",

  # Define the slots
  slots = list(
    name = "character",
    cons.data = "tbl_df"
  ),

  # Set the default values for the slots. (optional)
  prototype = list(name = NULL,
                   cons.data = NULL)
)



#' Method for creating a UMI sample
#' @export
#' @importFrom readr read_delim
#' @importFrom methods new
#' @importFrom utils read.csv
#' @importFrom dplyr rename
#' @param sample.name UMI sample object name
#' @param sample.dir Path to UMI sample
#' @param cons Consensus depth. Needs to be string; default is 10.
createUMIsample_Debarcer <- function(sample.name,sample.dir,cons = "10"){

  cons.file <- list.files(path = file.path(sample.dir,"tables"),pattern = paste("\\.cons",
                                                                                cons,".txt$",
                                                                                sep = ""))

  cons.table <- readr::read_delim(file = file.path(sample.dir,"tables",cons.file),
                                  delim = "\t",
                                  col_types = cols(
                                    `#AmpliconChromStart` = col_character(),
                                    Alias = col_character(),
                                    Position = col_double(),
                                    ProbableRef = col_character(),
                                    rawA = col_double(),
                                    rawC = col_double(),
                                    rawD = col_double(),
                                    rawG = col_double(),
                                    rawI = col_double(),
                                    rawN = col_double(),
                                    rawT = col_double(),
                                    rawDepth = col_double(),
                                    consA = col_double(),
                                    consC = col_double(),
                                    consD = col_double(),
                                    consG = col_double(),
                                    consI = col_double(),
                                    consN = col_double(),
                                    consT = col_double(),
                                    consDepth = col_double()
                                  ))


  UMI.sample <- DebarcerSample(name = sample.name,
                               cons.data = cons.table)

  return(UMI.sample)
}

#' Method for creating a UMI experiment object
#' @export
#' @import tibble
#' @importFrom utils read.csv
#' @importFrom methods new
#' @importFrom dplyr bind_rows
#' @param experiment.name Name of the experiment
#' @param main.dir Main experiment directory
#' @param dir.names List of sample names
createUMIexperiment_Debarcer <- function(experiment.name,main.dir,dir.names){
  main = main.dir
  cons.data.merged = tibble()

  for(i in 1:length(dir.names)){

    sample <- createUMIsample_Debarcer(dir.names[i], file.path(main,dir.names[i]))

    cons <- sample@cons.data
    cons$sample <- dir.names[i]
    cons.data.merged <- dplyr::bind_rows(cons.data.merged,cons)
  }

  UMIexperiment <- DebarcerExperiment(name = experiment.name,
                                      cons.data = cons.data.merged)
  return(UMIexperiment)
}
