#' Save consensus data
#' 
#' If save is set to TRUE data will be written to a csv file otherwise consensus 
#' data will be returned as a tibble.
#' 
#' @export
#' 
#' @importFrom readr write_excel_csv write_delim
#' 
#' @param object UMIexperiment object
#' @param outDir output directory, defaults to working directory
#' @param save Logical. Should data be saved to file? Default is FALSE.
#' @param delim Single character string, either ';' or ',' or tab
#' @param fileName String. Name of the file to be saved. Default 
#' is 'consensus_data.csv'
#'
#' @return A data table
#' 
#' @examples
#' library(umiAnalyzer)
#' 
#' main = system.file('extdata', package = 'umiAnalyzer')
#' 
#' samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' 
#' example <- createUmiExperiment(experimentName = 'example',mainDir = main,sampleNames = samples)
#' 
#' consensus_data <- saveConsData(object = example)
#' consensus_data
#' 
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
#' @examples
#' library(umiAnalyzer)
#' 
#' main = system.file('extdata', package = 'umiAnalyzer')
#' 
#' samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' 
#' simsen <- createUmiExperiment(experimentName = 'simsen',mainDir = main,sampleNames = samples)
#' 
#' simsen <- filterUmiObject(simsen)
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
#' 
#' @importFrom readr write_excel_csv write_delim
#' 
#' @param object Requires a UMI sample or UMI experiment object.
#' @param name String. Name of the filter. Default is "default".
#' @param save Logical, should data be saved as csv file? Default is FALSE.
#' @param outDir Output directory
#' @param fileName Filename to be used, default is the same as 'name'
#' @param delim Character string denoting delimiter to be used, default is ';'.
#' 
#' @return A filtered consensus table, as a tibble.
#' 
#' @examples
#' library(umiAnalyzer)
#' 
#' main = system.file('extdata', package = 'umiAnalyzer')
#' 
#' samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' 
#' simsen <- createUmiExperiment(experimentName = 'simsen',mainDir = main,sampleNames = samples)
#' simsen <- filterUmiObject(simsen)
#' 
#' myfilter <- getFilteredData(simsen)
#' myfilter
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
    outDir <- getwd()
  } else if(!is.character(fileName) && !is.null(fileName)){
    stop("fileName needs to be a string or NULL")
  } else if(! delim %in% c(';',',','\t')){
    warning("Invalid delimeter, needs to be comma, semicolon or tab.
            Using comma instead.")
    delim <- ','
  } else if(is.null(object@filters[name][[1]])) {
    if(!is.null(object@filters$default)){
      warning("Requested filter not found, using default.")
      name <- 'default'
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
#' @examples
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#' 
#' metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#'
#' simsen <- importDesign(object = simsen,file = metaData)
#' 
#' # Retrieve meta data
#' design <- getMetaData(object = simsen, attributeName = "design")
#' design
#' 
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


#' Beta binomial model
#' 
#' Code was obtained from VGAM package function VGAM::rbetabinom.ab. The VGAM
#' package is available under the GPL-3 license and maintained by
#' Thomas Yee <t.yee at auckland.ac.nz>. Source code of the function is identical
#' to rbetabinom.ab, but the function name was changed to beta_binom.
#' 
#' @references Yee TW (2015). Vector Generalized Linear and Additive Models: With an Implementation in R. Springer, New York, USA.
#' 
#' @export
#' 
#' @importFrom stats rbinom rbeta
#' 
#' @param n n
#' @param size size
#' @param shape1 alpha
#' @param shape2 beta
#' @param limit.prob 0.5
#' @param .dontuse.prob NULL
#' 
#' @return Numeric 
#' 
#' @examples
#' beta_binom(10,5, 0.5, 1)
#' beta_binom(10,2, 0.5, 1)
#' 
beta_binom <- function(
  n,
  size,
  shape1,
  shape2,
  limit.prob = 0.5,
  .dontuse.prob = NULL
){
  use.n <- if((length.n <- length(n)) > 1){
    length.n
  } else if(!is_Numeric(n, integer.valued = TRUE, length.arg = 1, positive = TRUE)) {
    stop("bad input for argument 'n'")
  } else {
    n
  }
  
  if (length(size) != use.n) {
    size <- rep_len(size, use.n)
  } 
    
  if (length(shape1) != use.n){
    shape1 <- rep_len(shape1, use.n)
  }
    
  if (length(shape2) != use.n){
    shape2 <- rep_len(shape2, use.n)
  }
  
  if (length(limit.prob) != use.n) {
    limit.prob <- rep_len(limit.prob, use.n)
  }
    
  ans <- rep_len(NA_real_, use.n)
  ind3 <- !is.na(shape1) & !is.na(shape2) & ((is.infinite(shape1) & 
                                                is.infinite(shape2)))
  if (sum.ind3 <- sum(ind3)) {
        ans[ind3] <- stats::rbinom(
          n = sum.ind3,
          size = size[ind3], 
          prob = limit.prob[ind3]
        )
  }

  if (ssum.ind3 <- sum(!ind3)) {
        ans[!ind3] <- stats::rbinom(
          n = ssum.ind3,
          size = size[!ind3],
          prob = stats::rbeta(
            n = ssum.ind3, 
            shape1 = shape1[!ind3],
            shape2 = shape2[!ind3]
          )
        )
  }

  ans[is.na(shape1) | shape1 < 0] <- NaN
  ans[is.na(shape2) | shape2 < 0] <- NaN
  
  ans
}

#' Is numeric
#' 
#' VGAM package function VGAM:::is.Numeric. The VGAM
#' package is available under the GPL-3 license and maintained by
#' Thomas Yee <t.yee at auckland.ac.nz>. Source code of the function is identical
#' to is.Numeric., but the function name was changed to is_Numeric.
#' 
#' @references Yee TW (2015). Vector Generalized Linear and Additive Models: With an Implementation in R. Springer, New York, USA.
#' 
#' @noRd
#' 
#' @importFrom stats rbinom rbeta
#' 
#' @param x x
#' @param length.arg Inf
#' @param integer.valued FALSE
#' @param positive FALSE
#' 
#' @return Boolean
#' 
is_Numeric <- function(
  x,
  length.arg = Inf,
  integer.valued = FALSE,
  positive = FALSE
  ){
    if (all(is.numeric(x)) && 
        all(is.finite(x)) && 
        (if (is.finite(length.arg)) length(x) == length.arg else TRUE) && 
        (if (integer.valued) all(x == round(x)) else TRUE) && 
        (if (positive) all(x > 0) else TRUE)){
      TRUE
    } else {
      FALSE
    } 

}
  


#' Download meta data template
#'
#' Function for downloading a template file containing metadata.
#'
#' @param object A UMIexperiment object
#'
#' @importFrom tibble enframe
#' @importFrom dplyr rename
#'
#' @export
#'
#' @return A tibble containing a metadata template
#' 
#' @examples
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#' 
#' download_template(simsen)
#'
download_template <- function(object){
  data <- object@cons.data

  samples <- tibble::enframe(unique(data$`Sample Name`), name = NULL)
  samples <- dplyr::rename(samples, Sample_Name = .data$value)

  return(samples)
}

#' Theme selection
#'
#' Function to select plotting theme based on user choice.
#'
#' @param theme User supplied theme selection
#'
#' @import ggplot2
#' 
#' @noRd
#'
#' @return A ggplot theme.
#'
select_theme <- function(theme){
  if(theme == 'classic'){
    use_theme <- ggplot2::theme_classic()
  } else if(theme == 'umiVisualiser') {

    use_theme <- ggplot2::theme_bw() +
      theme(
        # remove the plot background
        plot.background = ggplot2::element_blank(),
        panel.background = ggplot2::element_blank(),

        # make the legend and strip background transparent
        legend.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        ),
        legend.key = ggplot2::element_rect(
          fill = "transparent",
          colour = NA),
        strip.background = ggplot2::element_rect(
          fill = "transparent",
          colour = NA
        ),

        panel.grid.major = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),

    # remove the axis tick marks and hide axis lines
    axis.ticks = ggplot2::element_line(colour = "#454545", size = 0.3),
    axis.line = ggplot2::element_line(color = "#454545", size = 0.3),

    # modify the bottom margins of the title and subtitle
    plot.title = ggplot2::element_text(
      size = 18, colour = "#454545",
      hjust = 0.5,
      margin = ggplot2::margin(b = 10)
    ),
    plot.subtitle = ggplot2::element_text(
      size = 12, colour = "#454545",
      hjust = 0.5,
      margin = ggplot2::margin(b = 10)
    ),

    # add padding to the caption
    plot.caption = ggplot2::element_text(
      size = 10, colour = "#454545",
      hjust = 1,
      margin = ggplot2::margin(t = 15)
      ),

    # Adjust text size and axis title position
    axis.title = ggplot2::element_text(
      size = 13,
      colour = "#454545",
      hjust = 0.95
    ),
    axis.text = ggplot2::element_text(
      size = 10,
      colour = "#212121"
    ),
    legend.title = ggplot2::element_text(
      size = 12,
      colour = "#454545"
    ),
    legend.text = ggplot2::element_text(
      size = 10,
      colour = "#454545"
    ),
    strip.text = ggplot2::element_text(
      size = 12, colour = "#454545",
      margin = ggplot2::margin(10, 10, 10, 10, "pt")
      )
    )

  } else if(theme == 'bw'){
    use_theme <- ggplot2::theme_bw()
  } else if(theme == 'gray'){
    use_theme <- ggplot2::theme_gray()
  } else if(theme == 'minimal'){
    use_theme <- ggplot2::theme_minimal()
  } else if(theme == 'light'){
    use_theme <- ggplot2::theme_light()
  } else{
    warning('Invalid theme chosen, using classic theme.')
    use_theme <- ggplot2::theme_classic()
  }

  return(use_theme)
}

#' Filter samples and amplicons from a consensus table
#'
#' @importFrom dplyr filter
#' @param consensus.data A consensus fdata table.
#' @param amplicons Null or list of amplicons to use.
#' @param samples Null or a list of samples to use.
#' @param positions Null or a list of positions to use.
#'
#' @noRd
#'
#' @return A consensus table.
filterConsensusTable <- function(
  consensus.data,
  amplicons = NULL,
  samples = NULL,
  positions = NULL
  ) {

  #if (!is.null(amplicons)) {
  #  consensus.data <- consensus.data %>%
  #  dplyr::filter(.data$Name %in% amplicons)
  #}

  consensus.data.bind <- consensus.data

  if (!is.null(amplicons)) {
    consensus.data.bind <- tibble()

    for (name in amplicons ){
      bind <- consensus.data[stringr::str_detect(string = consensus.data$Name, pattern = name),]
      consensus.data.bind <- dplyr::bind_rows(consensus.data.bind,bind)
    }
  }

  if (!is.null(samples)) {
    consensus.data.bind <- consensus.data.bind %>%
      dplyr::filter(.data$`Sample Name` %in% samples)
  }

  if (!is.null(positions)) {
    consensus.data.bind <- consensus.data.bind %>%
      dplyr::filter(.data$Position %in% positions)
  }


  return(consensus.data.bind)
}

#' Merge assays
#'
#' Merge assays together by name. Requires a name of the new assay and
#' a list of assays that will be merged.
#'
#' @param object A UMIexperiment object
#' @param name Name of the new assay
#' @param assay.list List of assays to merge
#'
#' @export
#'
#' @importFrom dplyr mutate
#'
#' @return merged consensus data
#' 
#' @examples
#' \donttest{
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#'
#' simsen <- mergeAssays(object = simsen,name = "new",assay.list = c("PIK3CA_123", "PIK3CA_234"))
#' }
#'
mergeAssays <- function(object, name, assay.list){

  # Update consensus data
  data <- object@cons.data
  data <- data %>%
    dplyr::mutate(Name = as.character(.data$Name)) %>%
    dplyr::mutate(Name = replace(.data$Name, .data$Name %in% assay.list, name))

  # Update summary data
  summary.data <- object@summary.data
  summary.data <- summary.data %>%
    dplyr::mutate(assay = as.character(.data$assay)) %>%
    dplyr::mutate(assay = replace(.data$assay, .data$assay %in% assay.list, name))

  # Update umiExperiment object
  object@cons.data <- data
  object@summary.data <- summary.data

  return(object)
}



#' Add metaData
#'
#' @param object R object to which meta data should be added
#' @param attributeName Name of the meta data attribute.
#' @param attributeValue Meta data to be saved.
#'
#' @export
#' 
#' @return A UMIexperiment object
#'
#' @examples
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#' 
#' metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#'
#' simsen <- addMetaData(simsen,'metaData',metaData)
#'
addMetaData <- function(object,attributeName,attributeValue){
  attr(x = object, attributeName) <- attributeValue
  return(object)
}

#' Retrieve meta data by name.
#' @export
#' 
#' @param object R object from which to get meta data.
#' @param attributeName Name of the meta data attribute.
#' 
#' @return Metadata
#' 
#' @examples
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#' 
#' metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#'
#' simsen <- importDesign(object = simsen,file = metaData)
#' design <- getMetaData(object = simsen, attributeName = "design")
#'
getMetaData <- function(object,attributeName){
  if(attributeName %in% names(attributes(object))){
    value <- attributes(object)[names(attributes(object)) == attributeName][[1]]
    return(value)
  }
  else{
    warning("Attribute not found in object.")
  }
}

#' Generate VCF file from UMI sample or UMI experiment object
#'
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param outDir String. Output directory, defaults to working directory.
#' @param outFile String. Name of the output file
#' @param printAll Logical. Should all or only trusted variant be printed?
#'
#' @export
#' 
#' @return A VCF file
#' 
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' main <- system.file("extdata", package = "umiAnalyzer")
#'
#' simsen <- createUmiExperiment(main)
#' 
#' simsen <- filterUmiObject(simsen)
#'
#' generateVCF(simsen,'simsen.vcf', printAll = FALSE, save = FALSE)
#' }
#'
generateVCF <- function(object, outDir = getwd(), outFile, printAll = FALSE) {
  cons.table <- object@cons.table
  cons.table$Variants <- ifelse(cons.table$`Max Non-ref Allele Count` >= 5, "Variant", "Background")

  header <- c(
    "##fileformat=VCFv4.3",
    paste("##fileDate=", Sys.Date(), sep = ""),
    "##source=umiAnalyzerv0.3.0",
    "##INFO=<ID=DP,Number=1,Type=Integer,Description='Total Depth'>",
    "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele Frequency'>",
    "##INFO=<ID=PADJ,Number=A,Type=Float,Description='FDR-adjusted p-value'>",
    # paste("##INFO=<ID=SAMPLE,Number=A,Type=Float,Description=",
    #      sample.name,">"),
    paste("#",
          paste("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                collapse = "\t"
          ),
          sep = ""
    )
  )
  lines <- c()
  lines <- append(lines, header)

  if (printAll == FALSE) {
    cons.table <- cons.table %>% dplyr::filter(.data$Variants == "Variant")
  }

  for (i in 1:nrow(cons.table)) {
    row <- cons.table[i, ]

    if (row$Variants == "Variant") {
      vcfRow <- paste(row$Contig, row$Position, ".", row$Reference, row$`Max Non-ref Allele`,
                      "PASS", paste("DP=", row$Coverage, ";",
                                    "AF=", row$`Max Non-ref Allele Frequency`,
                                    sep = ""
                      ),
                      collapse = "\t"
      )
      lines <- append(lines, vcfRow)
    }
    else {
      if (is.na(row$`Max Non-ref Allele`)) {
        vcfRow <- paste(row$Contig, row$Position, ".", row$Reference, ".",
                        "FAIL", paste("DP=", row$Coverage, ";",
                                      "AF=", row$`Max Non-ref Allele Frequency`,
                                      sep = ""
                        ),
                        collapse = "\t"
        )
        lines <- append(lines, vcfRow)
      }
      else {
        vcfRow <- paste(row$Contig, row$Position, ".", row$Reference, row$`Max Non-ref Allele`,
                        "FAIL", paste("DP=", row$Coverage, ";",
                                      "AF=", row$`Max Non-ref Allele Frequency`,
                                      sep = ""
                        ),
                        collapse = "\t"
        )
        lines <- append(lines, vcfRow)
      }
    }
  }

  fileConn <- file(file.path(outDir, paste(outFile, ".vcf", sep = "")))
  writeLines(lines, fileConn)
  close(fileConn)
}


#' Import bed file
#'
#' @param path path to bed file
#'
#'
#' @import readr
#' @importFrom dplyr rename
#'
#' @export
#' 
#' @return A table containing genome positions
#' 
#' @examples 
#' library(umiAnalzyer)
#' 
#' bed_dir <- system.file("extdata", "simple.bed", package = "umiAnalyzer")
#' bed <- importBedFile(path = bed_dir)
#'
importBedFile <- function(path){

  bed <- readr::read_delim(
    file = path,
    delim = '\t',
    col_names = FALSE,
    col_types = readr::cols(
      X1 = readr::col_character(),
      X2 = readr::col_integer(),
      X3 = readr::col_integer(),
      X4 = readr::col_character()
    )
  ) %>%
    dplyr::rename(
      chrom = .data$X1,
      chromStart = .data$X2,
      chromEnd = .data$X3,
      Variant = .data$X4
    )

  l <- NULL

  for(i in 1:nrow(bed)){

    l <- append(l, bed$chromStart[i]:bed$chromEnd[i])

  }

  return(l)
}
