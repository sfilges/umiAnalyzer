#' Download meta data template
#'
#' Function for downloading a template file containing metadata.
#'
#' @param object A UmiExperiment object
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
#' @return A ggplot theme.
#'
select_theme <- function(theme){
  if(theme == 'classic'){
    use_theme <- ggplot2::theme_classic()
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

#' Filter samples and amplicons from a consensus tabel
#'
#' @importFrom dplyr filter
#' @param consensus.data A consensus fdata table.
#' @param amplicons Null or list of amplicons to use.
#' @param samples Null or a list of samples to use.
#' @return A consensus table.
filterConsensusTable <- function(
  consensus.data,
  amplicons = NULL,
  samples = NULL
  ) {

  if (!is.null(amplicons)) {
    consensus.data <- consensus.data %>%
      dplyr::filter(.data$Name %in% amplicons)
  }

  if (!is.null(samples)) {
    consensus.data <- consensus.data %>%
      dplyr::filter(.data$`Sample Name` %in% samples)
  }

  return(consensus.data)
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
#' @param time.var String. Name of thethe time variable. Default is "time"
#' @param use.variants Logical. Should pre computed variants be used? Default is FALSE.
#' @param group.by String. Variable for grouping data, e.g. replicates. Default is NULL.
#' @param do.plot Should plot be shown?
#'
#' @export
#'
#' @importFrom magrittr "%>%" "%<>%"
#' @import dplyr
#' @importFrom rlang .data
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
    by = c(`Sample Name` =  colnames(metaData[,1]))
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
      group = ~Position,
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
addMetaData <- function(object,attributeName,attributeValue){
  attr(x = object, attributeName) <- attributeValue
  return(object)
}

#' Retrieve meta data by name.
#' @export
#' @param object R object from which to get meta data.
#' @param attributeName Name of the meta data attribute.
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
#' @param outDir String. Output directory, defaults to wokring directory.
#' @param outFile String. Name of the output file
#' @param printAll Logical. Should all or only trusted variant be printed?
#'
#' @export
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
