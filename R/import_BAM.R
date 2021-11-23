#' readBamFile
#'
#' Method for reading bam files using scanBam from the Rsamtools package.
#'
#' @param sampleDir Path to UMI sample
#' @param consDepth Only retain consensus reads of at least cons.depth. Default is 0.
#'
#' @import tibble
#' @import magrittr
#' 
#' @importFrom Rsamtools scanBam
#' @importFrom tidyr separate unite
#' @importFrom dplyr filter
#' @importFrom BiocManager install
#' 
#' @noRd
#' 
#' @return A tibble containing reads extracted from BAM file
#'
readBamFile <- function(
  sampleDir,
  consDepth = 0
) {
  
  if (!requireNamespace("Rsamtools", quietly = TRUE)){
    BiocManager::install("Rsamtools")
  } 
  
  if(!dir.exists(sampleDir)){
    stop("You must supply a valid path.")
  }
  
  bam.file <- list.files(
    path = sampleDir,
    pattern = "\\_consensus_reads.bam$"
  )
  
  # If no bam files are located print an error message and return NULL
  if(length(bam.file) == 0){
    warning(paste("The directory you supplied does not seem to contain bam files: "),
            sampleDir, sep = '')
    return(NULL)
  } else {
    # Load bam file
    bam <- Rsamtools::scanBam(
      file = file.path(sampleDir, bam.file[1])
    )
    
    # Extract sequence information from bam object
    sequences <- tibble(
      qname = bam[[1]]$qname,
      chrom = bam[[1]]$rname,
      pos = bam[[1]]$pos,
      seq = as.data.frame(bam[[1]]$seq)$x
    )
    
    sequences <- tidyr::separate(
      sequences,
      col = .data$qname,
      into = c(NA, NA, NA, 'barcode', 'count'),
      sep = '_',
      remove = TRUE
    ) %>%
      tidyr::separate(
        col = .data$count,
        sep = '=',
        into = c(NA, 'count')
      ) %>%
      tidyr::unite(
        col = 'position',
        .data$chrom,
        .data$pos,
        sep = ':'
      )
    
    sequences$count %<>% as.integer
    
    sequences <- dplyr::filter(sequences, count >= consDepth)
    
    return(sequences)
  }
}

#' Function to parse bam files
#' @export
#' 
#' @import tibble
#' @importFrom dplyr bind_rows progress_estimated
#' @importFrom graphics plot
#' 
#' @param mainDir Directory containing UMIErrorCorrect output folders.
#' @param sampleNames A list of sample names.
#' @param consDepth Only retain consensus reads of at least cons.depth. Default is 0.
#' @param as.shiny Set to TRUE if run within a shiny::withProgress environment
#' 
#' @return A data table
#'
parseBamFiles <- function(
  mainDir,
  sampleNames = NULL,
  consDepth = 0,
  as.shiny = FALSE) {
  
  if (missing(x = mainDir)) {
    stop("Must provide a working directory and sample names.")
  } else if(!is.numeric(consDepth) | consDepth < 0){
    stop("consDepth needs to be a positive integer.")
  } else if(!dir.exists(mainDir)) {
    stop("You must supply a valid directory.")
  } else if(!is.null(sampleNames) && !is.list(sampleNames)){
    warning("sampleNames must be NULL or list. Using defaults.")
    sampleNames = NULL
  }
  
  # Get sample names
  if (is.null(sampleNames)){
    dir.names <- list.dirs(
      path = mainDir,
      full.names = FALSE,
      recursive = FALSE
    )
  }
  
  seq.Data <- tibble()
  
  if(length(dir.names) == 0){
    stop("No sample folders found. Have you supplied a valid top level
         directory containing umierrorcorrect output folders?")
  }
  
  for (i in 1:length(dir.names)) {
    
    seq.Table <- readBamFile(
      sampleDir = file.path(mainDir, dir.names[i]),
      consDepth = consDepth
    )
    
    # If NULL is returned by readBamFile no bam file was found
    if(is.null(seq.Table)){
      print(
        paste(
          "No bam file was found for sample:",dir.names[i],
          "in directory:",mainDir,
          sep=" "
        )
      )
    }
    
    # If running in shiny app, displat a loading bar
    if(as.shiny){
      n <- length(dir.names)
      shiny::incProgress(1/n, detail = paste("Parsing reads", i, " of ", n))
    }
    
    seq.Table$sample <- dir.names[i]
    
    seq.Data <- dplyr::bind_rows(seq.Data, seq.Table)
  }
  
  return(seq.Data)
}