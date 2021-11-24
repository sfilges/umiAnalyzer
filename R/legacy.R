#-------// Legacy functions //---------


#' Analyze time-course data
#'
#'
#' @param object UMIexperiment object containing meta data
#' @param filter.name Name of the filter to use.
#' @param time.var String. Name of the time variable. Default is "time"
#' @param use.variants Logical. Should pre computed variants be used? Default is FALSE.
#' @param group.by String. Variable for grouping data, e.g. replicates. Default is NULL.
#' @param do.plot Should plot be shown?
#'
#' @noRd
#'
#' @importFrom magrittr "%>%" "%<>%"
#' @import dplyr
#' @importFrom stats sd
#'
#' @return A UMIexperiment object
#' 
#' 
analyzeTimeSeries <- function(
  object,
  filter.name = "default",
  time.var = "time",
  use.variants = FALSE,
  group.by = NULL,
  do.plot = TRUE
) {
  
  if (missing(x = object)) {
    stop("Must provide a umiExperiment object and filter names")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(is.null(object@filters$default)) {
    stop("No data filter found.")
  }
  
  # Check if variant caller has been run on object
  if (use.variants == FALSE) {
    consData <- getFilteredData(
      object = object,
      name = filter.name
    )
    consData$Position %<>% as.factor
    consData$Variants <- ifelse(consData$`Max Non-ref Allele Count` >= 5, "Variant", "Background")
  }
  else {
    consData <- object@variants
    consData$Variants <- ifelse(consData$p.adjust <= 0.05, "Variant", "Background")
  }
  
  metaData <- as_tibble(object@meta.data)
  
  
  # If no group.by info is provided use sample name instead, else use the
  # group.by column from meta.data.
  if (is.null(group.by)){
    metaData$group.by <- dplyr::pull(metaData, 1)
    metaData$group.by %<>% as.factor
  } else {
    metaData$group.by <- dplyr::pull(metaData, group.by)
    metaData$group.by %<>% as.factor
  }
  
  
  # Join meta data and consData replicate ID column to consData
  # Change to left_join?
  summaryData <- dplyr::inner_join(
    consData,
    metaData,
    by = c(`Sample Name` ="Sample_Name")
  )
  
  print(summaryData)
  
  summaryData <- summaryData %>%
    dplyr::filter(.data$Variants == "Variant",
                  !.data$`Max Non-ref Allele` %in% c("I", "D")) %>%
    tidyr::unite("Position", .data$Contig, .data$Position, sep = ":") %>%
    tidyr::unite("Change", .data$Reference, .data$`Max Non-ref Allele`, sep = ">") %>%
    tidyr::unite("Change", .data$Change, .data$Position, sep = "*") %>%
    tidyr::unite("Change", .data$Change, .data$Name, sep = "*") %>%
    tidyr::unite("Change", .data$Change, .data$group.by, sep = "*")
  
  summaryData <- summaryData[duplicated(summaryData$Change),]
  
  summaryData$time_var <- dplyr::pull(summaryData, time.var)
  summaryData$time_var %<>% as.factor
  
  summaryData <- summaryData %>% dplyr::group_by(.data$Change, .data$time_var) %>%
    dplyr::summarise(VAF = 100 * mean(.data$`Max Non-ref Allele Frequency`)) %>%
    dplyr::ungroup() %>%
    tidyr::separate(.data$Change, c("Change", "Position", "Name", "Sample"), sep = "\\*")
  
  
  time_course <- ggplot(
    summaryData, aes_(
      x = ~time_var,
      y = ~VAF,
      group = ~Change,
      shape = ~Name)) +
    theme_bw() +
    xlab("Time") +
    ylab("VAF [%]") +
    geom_point() +
    geom_line(aes(color = Position))
  
  if(do.plot){
    print(time_course)
    
    object@plots$time_course <- time_course
    return(object)
  } else {
    object@plots$time_course <- time_course
    return(object)
  }
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
#' @noRd
#'
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats sd
#'
#' @return A UMIexperiment object
#' 
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


#-------// Functions analyzing Debarcer output //---------

#' Define sample class
#' 
#' @noRd
#' 
DebarcerSample <- setClass(
  'DebarcerSample',
  slots = list(
    name = 'character',
    cons.data = 'tbl_df'
  )
)

#' Define experiment class
#' 
#' @import tibble
#' 
#' @noRd
#' 
DebarcerExperiment <- setClass(
  'DebarcerExperiment',

  # Define the slots
  slots = list(
    name = 'character',
    cons.data = 'tbl_df'
  ),

  # Set the default values for the slots. (optional)
  prototype = list(
    name = NULL,
    cons.data = NULL
  )
)



#' Method for creating a UMIsample object
#' @export
#' 
#' @importFrom readr read_delim
#' @importFrom methods new
#' @importFrom utils read.csv
#' @importFrom dplyr rename
#' 
#' @param sample.name UMI sample object name
#' @param sample.dir Path to UMI sample
#' @param cons Consensus depth. Needs to be string; default is 10.
#' 
#' @return A UMIsample object
#' 
createUMIsample_Debarcer <- function(sample.name,sample.dir,cons = '10'){

  cons.file <- list.files(
    path = file.path(sample.dir,'tables'),
    pattern = paste('\\.cons',cons,'.txt$',sep = '')
  )

  cons.table <- readr::read_delim(
    file = file.path(sample.dir,'tables',cons.file),
    delim = '\t',
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


  UMI.sample <- DebarcerSample(
    name = sample.name,
    cons.data = cons.table
  )

  return(UMI.sample)
}

#' Method for creating a UMI experiment object
#' @export
#' 
#' @import tibble
#' @importFrom utils read.csv
#' @importFrom methods new
#' @importFrom dplyr bind_rows
#' 
#' @param experiment.name Name of the experiment
#' @param main.dir Main experiment directory
#' @param dir.names List of sample names
#' 
#' @return A UMIexperiment object
#' 
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
