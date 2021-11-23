#' Beta binomial model
#' 
#' VGAM package function VGAM::rbetabinom.ab
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

#' Beta binomial model
#' 
#' VGAM package function VGAM:::is.Numeric
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
#' @export
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

  consensus.data.bind = consensus.data

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
#' @export
#'
#' @importFrom magrittr "%>%" "%<>%"
#' @import dplyr
#' @importFrom stats sd
#'
#' @return A UMIexperiment object
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
