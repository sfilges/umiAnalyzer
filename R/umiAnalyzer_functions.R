#' Save consensus data
#' If save is set to TRUE data will be written to a csv file otherwise consensus data will
#' be returned as a tibble.
#' @export
#' @importFrom readr write_excel_csv write_delim
#' @param object UMIexperiment object
#' @param outDir output directory, defaults to working directory
#' @param save Logical. Should data be saved to file? Default is FALSE.
#' @param delim Single character string, either ';' or ',' or tab
#' @param fileName String. Name of the file to be saved. Default is 'consensus_data.csv'
#'
#' @return A data table
#'
saveConsData <- function(
  object,
  save = FALSE,
  fileName = 'consensus_data.csv',
  outDir = getwd(),
  delim = ';'
  ){

  if(missing(x = object)) {
    stop("No UMIexperiment object supplied.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!is.logical(save)){
    stop("Save needs to be of type boolean.")
  } else if(!(is.character(fileName) && length(fileName)==1)){
    stop("Invalid file name.")
  } else if(!dir.exists(outDir)) {
    stop("Output directory does not exist.")
  } else if(! delim %in% c(';',',','\t')){
    stop("Invalid delimeter, needs to be comma, semicolon or tab.")
  }

  consData <- object@cons.data

  if (save) {
    path <- file.path(outDir,fileName)
    if (delim == ';') {
      readr::write_excel_csv(consData, path, delim = ';')
    } else if (delim == ',') {
      readr::write_excel_csv(consData, path)
    } else if (delim == '\t') {
      readr::write_delim(consData, path, delim = delim)
    }
  } else {
    return(consData)
  }
}

#' Find consensus reads
#' A function to analyze consensus read tables generated with parseBamFiles or
#' a UMIexperiment object containing reads.
#'
#' @export
#'
#' @import tibble
#'
#' @param object Either a tibble generated with parseBamFiles or a UMIexperiment object
#' @param pattern Regular expression
#' @param consDepth Minimum consensus depth to keep. Default is 0.
#' @param groupBy Should data be grouped by position, sample, both or not at all.
#'
#' @return A data table
#'
findConsensusReads <- function(
  object,
  consDepth = 0,
  groupBy = c("none", "sample", "position", "both"),
  pattern = NULL
  ){

  if(missing(x = object)){
    stop("No object supplied")
  } else if(!tibble::is_tibble(object) || !class(object) == "UMIexperiment") {
    stop("Need to supply either a UMIexperiment object tibble generated
         with parseBamFiles.")
  }

  if (class(object)[1] == "UMIexperiment") {
    readsTable <- object@reads
  } else {
    readsTable <- object
  }
}

#' Method for filtering UMIexperiment and sample objects
#' @export
#'
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom utils data
#' @importFrom dplyr filter
#' @importFrom tibble as_tibble
#'
#' @param object Requires a UMI sample or UMI experiment object.
#' @param name String. Name of the filter. Default is "default".
#' @param minDepth Consensus depth to analyze. Default is 3.
#' @param minCoverage Minimum coverage required for amplicons. Default is 1.
#' @param minFreq Minimum variant allele frequency to keep. Default is 0.
#' @param minCount Minimum variant allele count to keep. Default is 3.
#'
#' @return A UMI sample or UMI experiment object.
#'
#'
filterUmiObject <- function(
  object,
  name = "default",
  minDepth = 3,
  minCoverage = 100,
  minFreq = 0,
  minCount = 0) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object and filter name.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(minDepth < 3){
    warning("You set minDepth to a value below 3. This will severely impact
            error correction.")
  } else if(minCoverage < 100){
    warning("Minimum coverage is below 50 consensus reads. Data with so few
            reads may be very unreliable.")
  } else if( tibble::is_tibble(object@filters[name][[1]]) ){
    warning("Filter ", name, " already exists. Will be overwritten.")
  }

  cons.table <- object@cons.data

  raw.error <- cons.table %>%
    dplyr::filter(
      .data$`Consensus group size` == 0,
      .data$Coverage >= minCoverage,
      .data$Name != '',
      .data$`Max Non-ref Allele Frequency` >= minFreq,
      .data$`Max Non-ref Allele Count` >= minCount
    )

  cons.table <- cons.table %>%
    dplyr::filter(
      .data$`Consensus group size` == minDepth,
      .data$Coverage >= minCoverage,
      .data$Name != '',
      .data$`Max Non-ref Allele Frequency` >= minFreq,
      .data$`Max Non-ref Allele Count` >= minCount
    )

  object@filters[[name]] <- cons.table
  object@raw.error <- raw.error

  return(object)
}

#' Method for retrieving filtered data
#' @export
#' @importFrom readr write_excel_csv write_delim
#' @param object Requires a UMI sample or UMI experiment object.
#' @param name String. Name of the filter. Default is "default".
#' @param save Logical, should data be saved as csv file? Default is FALSE.
#' @param outDir Output directory
#' @param fileName Filename to be used, default is the same as 'name'
#' @param delim Character string denoting delimiter to be used, default is ';'.
#' @return A filtered consensus table, as a tibble.
#'
getFilteredData <- function(
  object,
  name = 'default',
  save = FALSE,
  outDir = getwd(),
  fileName = NULL,
  delim = ';'
  ) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object and filter name.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!is.logical(save)){
    warning("save needs to be of type boolean. Using defaults instead.")
    save = FALSE
  } else if(!dir.exists(outDir)){
    warning("outDir needs to be a valid path. Using working directory.")
    outDir = getwd()
  } else if(!is.character(fileName) && !is.null(fileName)){
    stop("fileName needs to be a string or NULL")
  } else if(! delim %in% c(';',',','\t')){
    warning("Invalid delimeter, needs to be comma, semicolon or tab.
            Using comma instead.")
    delim = ','
  } else if(is.null(object@filters[name][[1]])) {
    if(!is.null(object@filters$default)){
      warning("Requested filter not found, using default.")
      name = 'default'
    } else {
      stop("Filter not found. Have you run filterUmiObject?")
    }
  }

  filter <- object@filters[name][[1]]

  if (is.null(fileName)) {
    outFile <- paste(name, ".csv", sep = '')
  } else {
    outFile <- paste(fileName, ".csv", sep = '')
  }

  if (save) {
    path <- file.path(outDir, outFile)
    if (delim == ';') {
      readr::write_excel_csv(filter, path, delim = ';')
    } else if (delim == ',') {
      readr::write_excel_csv(filter, path)
    } else {
      readr::write_delim(filter, path, delim = delim)
    }
  } else {
    return(filter)
  }
}


#' Filter variants based on p values or depth
#'
#' You can filter variants called with the the "callVariants" function based
#' on adjusted p-value, minimum variant allele count and supply a list
#' of assays and samples to plot.
#'
#' @export
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr select filter between
#' @importFrom tibble as_tibble
#' @param object A UMIexperiment object
#' @param p.adjust Numeric. Adjusted p value (FDR). Default is 0.2.
#' @param minVarCount Integer. Minimum variant allele count. Default is 5.
#' @param amplicons NULL or list of assays to plot. NULL uses all.
#' @param samples NULL or list of samples to plot. NULL uses all.
#' @seealso \code{\link{callVariants}} on how to call variants.
#' @return A UMIexperiment object with filtered variants. Can be used to
#'   generate VCF files.
#'
filterVariants <- function(
  object,
  p.adjust = 0.2,
  minVarCount = 5,
  amplicons = NULL,
  samples = NULL
  ) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!dplyr::between(p.adjust, 0, 1)) {
    warning("Adjusted p-value cutoff needs to be between 0 and 1, using defaults.")
    p.adjust = 0.2
  } else if(minVarCount < 0) {
    warning("minVarCount must be a positive integer. Using defaults instead.")
    minVarCount = 5
  }

  # TODO update the check for presence of the variant data to checking object@variants instead of attributes

  if ("varCalls" %in% names(attributes(object))) {
    # Load the consensus data from object
    vars.to.print <- object@variants

    # Filter based on p-value and minimum variant allele depth and select important columns
    # using .data also prevents R CMD check from giving a NOTE about undefined global variables
    # (provided that you have also imported rlang::.data with @importFrom rlang .data).

    vars.to.print <- filterConsensusTable(
      consensus.data = vars.to.print,
      amplicons =  amplicons,
      samples = samples
    )

    vars.to.print <- vars.to.print %>%
      dplyr::filter(
        .data$`Max Non-ref Allele Count` >= minVarCount,
        .data$p.adjust <= p.adjust
      ) %>%
      dplyr::select(
        .data$`Sample Name`,
        .data$Contig,
        .data$Position,
        .data$Name,
        .data$Reference,
        .data$`Max Non-ref Allele`,
        .data$p.adjust,
        .data$Coverage,
        .data$`Max Non-ref Allele Count`,
        .data$`Max Non-ref Allele Frequency`,
        .data$sample
      )

    print(vars.to.print)
    object@variants <- vars.to.print

    return(object)
  }
  else {
    stop("You need to run callVariants before running filterVariants.")
  }
}

#' Import experimental design meta data such as replicates, treatments, categorical variables.
#' @export
#' 
#' @importFrom utils read.table
#' 
#' @param object UMI.experiment to which to add metadata
#' @param file File containing meta data
#' @param delim Column separator. Default is NULL (automatically determine delimiter)
#'
#' @return A UMIexperiment object
#'
importDesign <- function(
  object,
  file,
  delim = NULL
  ){

  # Error Handling
  if (missing(x = object) || missing(x = file)) {
    stop("Must provide a umiExperiment object and file name.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!is.character(file)) {
    stop("File must be a valid name.")
  }

  if(is.null(delim)) {
    # Automatically determine file type if delim = NULL
    # Import data using all three delimiters and then check the number of
    # columns. If the delimiter is wrong, there will only be one column.

    # TODO this works for these delimiters but doesn't handle other
    # delimiters well.
    comma <- read.table(file = file, sep = ',', header = TRUE)
    semicolon <- read.table(file = file, sep = ';', header = TRUE)
    tab <- read.table(file = file, sep = '\t', header = TRUE)

    if(ncol(comma) > 1){
      metaData <- comma
      print("Uploading comma separated meta data file.")
    } else if(ncol(semicolon) > 1) {
      metaData <- semicolon
      print("Uploading semicolon separated meta data file.")
    } else if(ncol(tab) > 1) {
      metaData <- tab
      print("Uploading tab separated meta data file.")
    } else {
      warning('Automatic delimiter selection failed: It seems like your metadata
              file is not delimited by either comma, semicolon or tab.')
    }
  } else if (!delim %in% c(',', ';', '\t')) {
    # If delim is not NULL and not comma, simicolon or tab, throw exception
    stop("Delimiter needs to be one of: c(',', ';', '\t')")
  } else {
    # Import file using user defined delimiter
    metaData <- read.table(
      file = file,
      sep = delim,
      header = TRUE
    )
  }

  # Add imported table to the object meta.data slot
  object@meta.data <- metaData
  object <- addMetaData(object = object, attributeName = 'design', metaData)

  return(object)
}

#' mergeTechnicalReplicates
#'
#' A function to merge replicates in UMIexperiment object. This will result in a merged data set
#' accessible from the UMIexperiment object using merged.data. This is meant to provide statistical
#' information across multiple replicates. If you want to merge multiple sequencing runs of the
#' sample into a single sample using the collapseReplicates function instead.
#'
#' @param object UMI.experiment to which to add metadata
#' @param filter.name Name of the filter to use. Defaults to "default".
#' @param do.plot Should normalization plot be shown. Default is TRUE.
#' @param group.by Variable used to group data. If NULL sample names will be used.
#' @param amplicons List of amplicons to use
#' @param samples List of samples to use
#' @param normalise.by.sample If TRUE, normalizes reads depth by both samples and assays. Otherwise only assays are used.
#' @param remove.singletons Remove variants only found in one replicate.
#' @param zero.counts Number between 0 and 1. What values should negative counts get?
#' @param option Color scale for plotting.
#' @param direction Direction of color scale if using viridis package.
#'
#' @export
#'
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats sd
#'
#' @return A UMIexperiment object
mergeTechnicalReplicates <- function(
  object,
  filter.name = 'default',
  do.plot = TRUE,
  group.by = NULL,
  amplicons = NULL,
  samples = NULL,
  normalise.by.sample = FALSE,
  remove.singletons = TRUE,
  zero.counts = 0.5,
  option = c('viridis', 'magma', 'plasma', 'inferno'),
  direction = c(1, -1)
  ) {

  # Error handling
  if (missing(x = object)) {
    stop("Must provide a umiExperiment object.")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(is.null(object@filters[filter.name][[1]])) {
    if(!is.null(object@filters$default)){
      warning("Requested filter ", filter.name, " not found, using default.")
      filter.name = 'default'
    } else {
      stop("Filter not found. Have you run filterUmiObject?")
    }
  } else if(!is.logical(do.plot)) {
    warning("do.plot needs to  be of type boolean. Using default instead.")
    do.plot = TRUE
  }

  # Import filtered data
  consData <- getFilteredData(
    object = object,
    name = filter.name
  )

  consData$Position %<>% as.factor

  # Select amplicons and samples
  consData <- filterConsensusTable(
    consData,
    amplicons = amplicons,
    samples = samples
  )

  # Get metadata and
  metaData <- as_tibble(object@meta.data)
  if (is.null(group.by)){
    metaData$group.by <- dplyr::pull(metaData, 1)
    metaData$group.by %<>% as.factor
  } else {
    metaData$group.by <- dplyr::pull(metaData, group.by)
    metaData$group.by %<>% as.factor
  }

  # Join meta data and consData replicate ID column to consData
  # Change to left_join?
  # First column of meta data needs to contain sample names
  consData <- dplyr::inner_join(consData, metaData, by = c(`Sample Name` = colnames(metaData[,1]) ))

  # Calculate normalization factor
  consData <- consData %>% group_by(.data$Name) %>%
    mutate(normFac= mean(.data$Coverage)/.data$Coverage) %>% # Normalization factor
    mutate(normCoverage = .data$Coverage*.data$normFac)  %>% # Normalized Coverage
    ungroup()

  # Plot coverage before and after normalization
  plot.norm <- vizNormalization(consData)

  # Summarize normalized data for output
  consData <- consData %>%
    # Group data by factors
    dplyr::group_by(
      .data$Name,
      .data$Contig,
      .data$Position,
      .data$Reference,
      .data$group.by
    ) %>%
    dplyr::summarise(
      avg.A = mean(.data$A*.data$normFac),
      avg.T = mean(.data$T*.data$normFac),
      avg.C = mean(.data$C*.data$normFac),
      avg.G = mean(.data$G*.data$normFac),
      avg.N = mean(.data$N*.data$normFac),
      avg.I = mean(.data$I*.data$normFac),
      avg.D = mean(.data$D*.data$normFac),
      avg.Depth = mean(.data$Coverage*.data$normFac),
      std.Depth = sd(.data$Coverage*.data$normFac),
      avg.Max.AF = mean(.data$`Max Non-ref Allele Frequency`),
      std.MaxAF = sd(.data$`Max Non-ref Allele Frequency`),
      avg.MaxAC = mean(.data$`Max Non-ref Allele Count`),
      std.MaxAC = sd(.data$`Max Non-ref Allele Count`)
    )

  # Plot normalised counts stacked by variant allele
  stacked.counts <- vizStackedCounts(
    consData,
    option = option,
    direction = direction
  )

  if(do.plot){
    # Return object
    object@plots$stacked_counts <- stacked.counts
    object@plots$norm_plot <- plot.norm
    object@merged.data <- consData

    print(object@plots$stacked_counts)

    return(object)

  } else {
    # Return object
    object@plots$stacked_counts <- stacked.counts
    object@plots$norm_plot <- plot.norm
    object@merged.data <- consData
    return(object)
  }
}
