#' createUmiSample
#'
#' Method for creating a UMI sample from UMIErrorCorrect output.
#'
#' @param sampleName UMI sample object name
#' @param sampleDir Path to UMI sample folders. Must be a folder generated by UMIErrorCorrect
#' @param importBam Logical. Should BAM files be imported at object initialization? Default is False.
#'
#' @export
#'
#' @import tibble
#' @importFrom readr read_delim cols col_character col_double
#' @importFrom methods new
#' @importFrom utils read.csv
#' @importFrom dplyr rename
#' 
#' @return An object of class UMIsample
#' 
#' @examples 
#' library(umiAnalyzer)
#' 
#' main = system.file('extdata', package = 'umiAnalyzer')
#' samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' s1 <- createUmiSample('s1',sampleDir = paste(main,"/",samples[1],sep=""))
#'
createUmiSample <- function(
  sampleName,
  sampleDir,
  importBam = FALSE
) {
  
  if(!dir.exists(sampleDir)){
    warning("You need to provide a valid path.")
    return(NULL)
  } else if(!is.logical(importBam)){
    warning("importBam needs to be of type boolean. Using defaults instead.")
    importbam <- FALSE
  }
  
  consFile <- list.files(path = sampleDir, pattern = "\\.cons$")
  summaryFile <- list.files(path = sampleDir, pattern = "\\_summary_statistics.txt$")
  
  if(length(consFile) == 0 | length(summaryFile) == 0){
    warning(
      paste(
        "No consensus or summary file found for sample: ",
        sampleName,". Did you supply a correct data folder?"
      )
    )
    return(NULL)
  }
  
  consTable <- readr::read_delim(
    file = file.path(sampleDir, consFile),
    delim = "\t",
    col_types = readr::cols(
      `Sample Name` = readr::col_character(),
      Contig = readr::col_character(),
      Position = readr::col_double(),
      Name = readr::col_character(),
      Reference = readr::col_character(),
      A = readr::col_double(),
      C = readr::col_double(),
      G = readr::col_double(),
      T = readr::col_double(),
      I = readr::col_double(),
      D = readr::col_double(),
      N = readr::col_double(),
      Coverage = readr::col_double(),
      `Consensus group size` = readr::col_double(),
      `Max Non-ref Allele Count` = readr::col_double(),
      `Max Non-ref Allele Frequency` = readr::col_double(),
      `Max Non-ref Allele` = readr::col_character()
    )
  )
  
  summaryTable <- readr::read_delim(
    file = file.path(sampleDir, summaryFile),
    delim = "\t",
    col_names = FALSE,
    col_types = readr::cols(
      X1 = readr::col_character(),
      X2 = readr::col_character(),
      X3 = readr::col_character(),
      X4 = readr::col_double(),
      X5 = readr::col_double(),
      X6 = readr::col_double(),
      X7 = readr::col_double()
    )
  ) %>%
    dplyr::rename(
      ID = .data$X1,
      region = .data$X2,
      assay = .data$X3,
      depth = .data$X4,
      fraction = .data$X5,
      totalCount = .data$X6,
      UMIcount = .data$X7
    )
  
  if (importBam) {
    readsTable <- readBamFile(sampleDir = sampleDir)
    
    UMIsample <- UMIsample(
      name = sampleName,
      cons.data = consTable,
      summary.data = summaryTable,
      reads = readsTable
    )
    
    return(UMIsample)
  }
  else {
    UMIsample <- UMIsample(
      name = sampleName,
      cons.data = consTable,
      summary.data = summaryTable,
      reads = tibble()
    )
    
    return(UMIsample)
  }
}

#' Method for creating a UMI experiment object
#'
#'
#' @param experimentName Name of the experiment
#' @param mainDir Main experiment directory
#' @param sampleNames List of sample names. Can be either NULL or list. If NULL all subdirectories of mainDir will be searched.
#' @param importBam Logical. Should bam files be imported on creation? Default is False.
#' @param as.shiny Set to TRUE if run within a shiny::withProgress environment
#'
#' @export
#'
#' @import tibble
#' @importFrom stringr str_remove
#' @importFrom utils read.csv
#' @importFrom methods new
#' @importFrom dplyr bind_rows
#' @importFrom shiny incProgress
#'
#' @return An object of class UMIexperiment 
#'
#' @examples
#' library(umiAnalyzer)
#' 
#' main = system.file('extdata', package = 'umiAnalyzer')
#' 
#' samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' 
#' exp1 <- createUmiExperiment(experimentName = 'exp1',mainDir = main,sampleNames = samples)
#'
createUmiExperiment <- function(
  mainDir,
  experimentName = NULL,
  sampleNames = NULL,
  importBam = FALSE,
  as.shiny = FALSE
){
  
  if (!dir.exists(mainDir)) {
    warning("Must provide a valid directory.")
    return(NULL)
  } else if(!is.null(experimentName)) {
    if( !is.character(experimentName) && length(experimentName) == 1 ) {
      warning("experimentName needs to be a string or NULL. Using default.")
      experimentName <- NULL
    }
  } else if(!is.logical(importBam)){
    warning("importBam needs to be of type boolean. Using default.")
    importBam <- FALSE
  } else if (!is.null(sampleNames) && !is.list(sampleNames)){
    warning("sampleNames must be NULL or a list. Using default.")
    sampleNames <- NULL
  }
  
  # Get sample names
  if (is.null(sampleNames)){
    sampleNames <- list.dirs(
      path = mainDir,
      full.names = FALSE,
      recursive = FALSE
    )
  }
  
  # Initialise dataframes
  cons.data.merged <- tibble::tibble()
  summary.data.merged <- tibble::tibble()
  reads.merged <- tibble::tibble()
  
  for (i in 1:length(sampleNames)) {
    
    if(as.shiny){
      n <- length(sampleNames)
      shiny::incProgress(1/n, detail = paste("Loading sample", i, " of ",n))
    }
    
    if(!dir.exists(file.path(mainDir, sampleNames[i]))){
      warning("Sample directory not found. Did you provide a top level directory
           containing you sample folders?")
    }
    
    # Find .cons file
    consFile <- list.files(
      path = file.path(mainDir, sampleNames[i]),
      pattern = "\\.cons$")
    
    if( identical(consFile, character(0))  ){
      next
    }
    
    # Remove file ending
    consFile <- stringr::str_remove(consFile, '.cons')
    
    # Create UMI sample
    sample <- createUmiSample(
      sampleName = consFile,
      sampleDir= file.path(mainDir, sampleNames[i]),
      importBam = importBam
    )
    
    # Sample returns NULL if no consensus file is present, skip that sample
    if(is.null(sample)){
      next
    }
    
    cons <- sample@cons.data
    cons$sample <- consFile
    cons.data.merged <- dplyr::bind_rows(cons.data.merged, cons)
    
    summary <- sample@summary.data
    summary$sample <- consFile
    summary.data.merged <- dplyr::bind_rows(summary.data.merged, summary)
    
    if(importBam) {
      reads <- sample@reads
      reads$sample <- consFile
      reads.merged <- dplyr::bind_rows(reads.merged, reads)
    }
  }
  
  if(identical(cons.data.merged, tibble())){
    return(NULL)
  } else {
    # Save experiment object
    UMIexperiment <- UMIexperiment(
      name = experimentName,
      cons.data = cons.data.merged,
      summary.data = summary.data.merged,
      reads = reads.merged
    )
    
    return(UMIexperiment) 
  }
}

#' Add UMI sample to an existing experiment object
#'
#' @param object UMIexperiment object
#' @param sampleName Name of new sample
#' @param sampleDir Directory to new sample
#' @param clearData Should other data in UMIexperiment be cleared
#'
#' @importFrom dplyr bind_rows
#'
#' @export
#' 
#' @return A UMIexperiment object
#'
addUmiSample <- function(
  object,
  sampleName,
  sampleDir,
  clearData = FALSE) {
  
  if(missing(x = object)) {
    stop("No UMIexperiment object supplied.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if (!dir.exists(sampleDir)){
    stop("No valid path provided.")
  } else if(!is.logical(clearData)) {
    warning("clearData needs to be of type boolean. Using defaults instead.")
    clearData <- FALSE
  }
  
  newSample <- createUmiSample(
    sampleName = sampleName,
    sampleDir = sampleDir,
    importBam = FALSE
  )
  
  newConsData <- newSample@cons.data
  newConsData$sample <- sampleName
  
  object@cons.data <- dplyr::bind_rows(object@cons.data, newConsData)
}