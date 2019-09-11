#' Define sample class
UMIsample <- setClass("UMIsample",
                      slots = list(name = "character",
                                   cons.data = "tbl_df",
                                   summary.data = "tbl_df",
                                   reads = "tbl_df")
)

setOldClass(c("tbl_df", "tbl", "data.frame"))

#' Define experiment class
#' @import tibble
UMIexperiment <- setClass(
  "UMIexperiment",

  # Define the slots
  slots = list(name = "character",
               cons.data = "tbl_df",
               summary.data = "tbl_df",
               reads = "tbl_df",
               meta.data = "data.frame",
               filters = "list",
               variants = "tbl_df",
               merged.data = "tbl_df"),

  # Set the default values for the slots. (optional)
  prototype = list(name = NULL,
               cons.data = NULL,
               summary.data = NULL,
               reads = NULL,
               meta.data = data.frame(),
               filters = list(),
               variants = tibble(),
               merged.data = tibble())
)

# Add function to append UMIexperiment

# Add slot for plots

#' Method for creating a UMI sample
#' @export
#' @import readr
#' @importFrom methods new
#' @importFrom utils read.csv
#' @importFrom dplyr rename
#' @param sample.name UMI sample object name
#' @param sample.dir Path to UMI sample
createUMIsample <- function(sample.name,sample.dir){
  cons.file <- list.files(path = sample.dir,pattern = "\\.cons$")

  cons.table <- readr::read_delim(file = file.path(sample.dir,cons.file),
                                  delim = "\t",
                                  col_types = cols(
                                    `Sample Name` = col_character(),
                                    Contig = col_character(),
                                    Position = col_double(),
                                    Name = col_character(),
                                    Reference = col_character(),
                                    A = col_double(),
                                    C = col_double(),
                                    G = col_double(),
                                    T = col_double(),
                                    I = col_double(),
                                    D = col_double(),
                                    N = col_double(),
                                    Coverage = col_double(),
                                    `Consensus group size` = col_double(),
                                    `Max Non-ref Allele Count` = col_double(),
                                    `Max Non-ref Allele Frequency` = col_double(),
                                    `Max Non-ref Allele` = col_character()
                                  ))

  summary.file <- list.files(path = sample.dir,pattern = "\\_summary_statistics.txt$")

  summary.table <- readr::read_delim(file = file.path(sample.dir,summary.file),
                                     delim = "\t",
                                     col_names = FALSE,
                                     col_types = cols(
                                       X1 = col_character(),
                                       X2 = col_character(),
                                       X3 = col_character(),
                                       X4 = col_double(),
                                       X5 = col_double(),
                                       X6 = col_double(),
                                       X7 = col_double()
                                       )) %>%
                                       dplyr::rename(ID = .data$X1,
                                                     region = .data$X2,
                                                     assay = .data$X3,
                                                     depth = .data$X4,
                                                     fraction = .data$X5,
                                                     totalCount = .data$X6,
                                                     UMIcount = .data$X7)

  reads.table = readBamFile(sample.dir = sample.dir)

  UMI.sample <- UMIsample(name = sample.name,
                          cons.data = cons.table,
                          summary.data = summary.table,
                          reads = reads.table)
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
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' main = system.file("extdata", package = "umiAnalyzer")
#' sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#'
#' exp1 <- create.UMIexperiment(experiment.name = "exp1", main.dir = main, dir.names = sample.names)
#' }
createUMIexperiment <- function(experiment.name,main.dir,dir.names){
  main = main.dir
  cons.data.merged = tibble()
  summary.data.merged = tibble()
  reads.merged = tibble()

  for(i in 1:length(dir.names)){

    sample <- createUMIsample(dir.names[i], file.path(main,dir.names[i]))

    cons <- sample@cons.data
    cons$sample <- dir.names[i]
    cons.data.merged <- dplyr::bind_rows(cons.data.merged,cons)

    summary <- sample@summary.data
    summary$sample <- dir.names[i]
    summary.data.merged <- dplyr::bind_rows(summary.data.merged,summary)

    reads = sample@reads
    reads$sample <- dir.names[i]
    reads.merged <- dplyr::bind_rows(reads.merged,reads)
  }

  UMIexperiment <- UMIexperiment(name = experiment.name,
                                 cons.data = cons.data.merged,
                                 summary.data = summary.data.merged,
                                 reads = reads.merged)
  return(UMIexperiment)
}

#' Method for reading bam files
#' @import tibble
#' @import magrittr
#' @import Rsamtools
#' @importFrom tidyr separate unite
#' @param sample.dir Path to UMI sample
readBamFile <- function(sample.dir){
  bam.file <- list.files(path = sample.dir,pattern = "\\_consensus_reads.bam$")

  bam <- scanBam(file.path(sample.dir, bam.file))

  sequences <- tibble(qname = bam[[1]]$qname,
                      chrom = bam[[1]]$rname,
                      pos = bam[[1]]$pos,
                      seq = as.data.frame(bam[[1]]$seq)$x)

  sequences <- tidyr::separate(sequences,
                               col = qname,
                               into = c(NA, NA, NA, "barcode", "count"),
                               sep = "_",
                               remove = TRUE) %>%
    tidyr::separate(col=.data$count,sep="=",into=c(NA,"count"))%>%
    tidyr::unite(col="position",.data$chrom,.data$pos,sep=":")

  sequences$count %<>% as.integer

  return(sequences)
}

#' Method for filtering UMIexperiment and sample objects
#' @export
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom utils data
#' @importFrom dplyr filter
#' @param object Requires a UMI sample or UMI experiment object.
#' @param name String. Name of the filter.
#' @param minDepth Consensus depth to analyze. Default is 3.
#' @param minCoverage Mininum coverage required for amplicons. Default is 1.
#' @param minFreq Minimum variant allele frequency to keep. Default is 0.
#' @param minCount Minimum variant allele count to keep. Default is 3.
#' @return A UMI sample or UMI experiment object.
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' data <- simsen
#' data <- filterUMIobject(data)
#' }
filterUMIobject <- function(object, name, minDepth=3, minCoverage=50, minFreq=0, minCount=0){

  cons.table <- object@cons.data

  cons.table <- cons.table %>%
    dplyr::filter(.data$`Consensus group size` == minDepth,
                  .data$Coverage >= minCoverage,
                  .data$Name != "",
                  .data$`Max Non-ref Allele Frequency` >= minFreq,
                  .data$`Max Non-ref Allele Count` >= minCount)

  object@filters[[name]] <- cons.table

  return(object)
}

#' Method for retrieving filtered data
#' @export
#' @param object Requires a UMI sample or UMI experiment object.
#' @param name String. Name of the filter.
#' @return A filtered consensus table, as a tibble.
getFilter <- function(object, name){
  filter <- object@filters[name]
  return(filter)
}


#' Calculates the negative log likelihood for the beta distribution.
#' @importFrom stats dbeta
#' @param params Non-negative parameters of the Beta distribution.
#' @param data consensus.data table of a UMisample or UMIexperiment object.
#' @return Negative log-likelihood for beta distribution.
betaNLL <- function(params,data){
  a<-params[1]
  b<-params[2]

  print(paste("a= ",a," b= ",b,sep=""))

  # negative log likelihood for beta
  return(-sum(dbeta(data,shape1=a, shape2=b, log=TRUE)))
}

#' Calculate variant p-values using permutation-based testing. A prior is fitted to model the background
#' error using maximum likelihood estimation of a beta distribution. The maximum likelihood estimate
#' of the beta distribution is then used to define the shape of a beta-binomial distribution used
#' to estimate variant P-Values. This can be interpreted as a probability for a variant to not have
#' arisen by chance
#' @export
#' @importFrom dplyr mutate progress_estimated
#' @importFrom tibble as_tibble
#' @importFrom VGAM rbetabinom.ab
#' @importFrom stats nlm var p.adjust
#' @importFrom utils install.packages
#' @importFrom graphics plot
#' @param object A UMierrorcorrect object.
#' @return Object containing raw and FDR-adjusted P-Values
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' data <- simsen
#' data <- callVariants(data)
#' }
callVariants <- function(object){

  if("filter" %in% names(attributes(object))){
    warning("It appears the UMI experiment object has been filtered. Consider
            running callVariants on an unfiltered object instead.")
    print(attributes(object)$filter)
  }

  if(!requireNamespace("VGAM", quietly = TRUE)){
    install.packages("VGAM")
  }

  object <- filterUMIobject(object = object, name = "varCalls",
                            minDepth = 3, # Require consensus 3
                            minCoverage = 100, # Require at least 100 cons reads
                            minFreq = 0, # no minimum allele freq
                            minCount = 0) # no minimum variant allele count

  cons.table <- object@filters["varCalls"][[1]]

  a1 <- cons.table$Coverage*cons.table$`Max Non-ref Allele Frequency` # No. of variant alleles
  b1 <- cons.table$Coverage # Total coverage

  m <- mean(a1/b1) # average background count
  v <- var(a1/b1) # variance of background counts

  # Calculate initial values
  a0 <- m*(m * (1-m) / v-1 )
  b0 <- (1-m)*(m * (1-m) / v-1 )
  params0 = c(a0,b0)

  fit <- nlm(betaNLL,params0,a0/b0)
  a<-fit$estimate[1]
  b<-fit$estimate[2]
  pval<-NULL

  pbar <- dplyr::progress_estimated(length(a1))
  print("Estimating p-values using permutation test. This might take a while.")
  for (i in 1:length(a1)){ # for each named amplicon position:
    pbar$tick()$print()
    r1 <- rbetabinom.ab(10000,b1[i],shape1=a,shape2=b) # Calculate probability of success
    pval[i] = sum(r1>a1[i])/10000                      # Estimate p value of variant

    graphics::plot(r1)
  }

  padj <- p.adjust(pval,method="fdr")

  cons.table <- dplyr::mutate(cons.table, pval = pval)
  cons.table <- dplyr::mutate(cons.table, p.adjust = padj)

  #object@cons.data <- cons.table
  object@variants <- cons.table

  object <- addMetaData(object = object, attributeName = "varCalls", "varCalls")

  return(object)
}


#' Filter variants based on p values or depth
#' @export
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr select filter
#' @importFrom tibble as_tibble
#' @importFrom rlang .data
#' @param object A UMIexperiment object
#' @param p.adjust Numeric. Adjusted p value (FDR). Default is 0.2.
#' @param minVarCount Integer. Minimum variant allele count. Default is 5.
#' @return A UMIexperiment object with filtered variants. Can be used to generate vcf files.
filterVariants <- function(object, p.adjust = 0.2, minVarCount = 5){
  if("varCalls" %in% names(attributes(object))){
    # Load the consensus data from object
    vars.to.print <- object@variants

    # Filter based on p-value and minimum variant allele depth and select important columns
    # using .data also prevents R CMD check from giving a NOTE about undefined global variables
    # (provided that you???ve also imported rlang::.data with @importFrom rlang .data).

    vars.to.print <- vars.to.print %>%
      dplyr::filter(.data$`Max Non-ref Allele Count` >= minVarCount,
                    .data$p.adjust <= p.adjust) %>%
      dplyr::select(.data$`Sample Name`,
                    .data$Contig,
                    .data$Position,
                    .data$Name,
                    .data$Reference,
                    .data$`Max Non-ref Allele`,
                    .data$Coverage,
                    .data$`Max Non-ref Allele Count`,
                    .data$`Max Non-ref Allele Frequency`)

    print(vars.to.print)
    object@variants <- vars.to.print

    return(object)
  }
  else{
    message <- simpleError("You need to run callVariants before running filterVariants.")
    print(message)
  }
}

#' Import experimental design meta data such as replicates, treatments, categorical variables.
#' @export
#' @importFrom utils read.table
#' @param object UMI.experiment to which to add metadata
#' @param file File containing meta data
#' @param sep Column separator. Default is tab.
importDesign <- function(object,file,sep="\t"){
  mData <- read.table(file = file, sep = sep, header = TRUE)

  object@meta.data <- mData
  object <- addMetaData(object = object, attributeName = "design", mData)

  return(object)
}

#' A function to merge replicates in UMIexperiment object. This will result in a merged data set
#' accessible from the UMIexperiment object using merged.data. This is meant to provide statistical
#' information across multiple replicates. If you want to merge multiple sequencing runs of the
#' sample into a single sample using the collapseReplicates function instead.
#' @export
#' @importFrom magrittr "%>%" "%<>%"
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats sd
#' @param object UMI.experiment to which to add metadata
#' @param filter.name Name of the filter to use.
mergeTechnicalReplicates <- function(object, filter.name) {

  consData <- getFilter(object = object, name = filter.name)
  consData <- consData[[1]]
  consData$Position %<>% as.factor

  mData <- as_tibble(object@meta.data)
  mData$replicate %<>% as.factor

  # Join meta data and consData replicate ID column to consData
  # Change to left_join?
  consData <- dplyr::inner_join(consData, mData, by = c(`Sample Name` = "sample"))

  # Calculate normalization factor
  consData <- consData %>% group_by(.data$Name) %>%
    mutate(normFac= mean(.data$Coverage)/.data$Coverage) %>% # Normalization factor
    mutate(normCoverage = .data$Coverage*.data$normFac)  %>% # Normalized Coverage
    ungroup()

  # Plot coverage before and after normalization
  plot.norm <- viz_Normalization(consData)
  print(plot.norm)

  # Summarize normalized data for output
  consData <- consData %>%
    # Group data by factors
    dplyr::group_by(.data$Name,
                    .data$Contig,
                    .data$Position,
                    .data$Reference,
                    .data$replicate) %>%
    dplyr::summarise(avg.A = mean(.data$A*.data$normFac),
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
                     std.MaxAC = sd(.data$`Max Non-ref Allele Count`))

  # Plot normalised counts stacked by variant allele
  stacked.counts <- viz_stacked_counts(consData)
  print(stacked.counts)

  # Return object
  object@merged.data <- consData
  return(object)
}

#' Analyze time-course data
#' @export
#' @importFrom magrittr "%>%" "%<>%"
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom stats sd
#' @param object UMI.experiment containing meta data
#' @param filter.name Name of the filter to use.
#' @param time.var String. Name of thethe time variable. Default is "time"
#' @param use.variants Logical. Should pre computed variants be used? Default is FALSE.
#' @param group.by String. Variable for grouping data, eg.g. replicates. Default is NULL.
analyzeTimeSeries <- function(object,
                              filter.name,
                              time.var = "time",
                              use.variants = FALSE,
                              group.by = NULL){


  #data <- filterUMIobject(object = data, name = "myfilter", minDepth = 3,
  #                        minCoverage = 100, minFreq = 0, minCount = 0)
  #myfilter <- getFilter(object = data, name = "myfilter")
  #metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
  #data <- importDesign(object = data, file = metaData)

  # Check if variant caller has been run on object
  if( use.variants == FALSE ) {
    consData <- getFilter(object = object, name = filter.name)
    consData <- consData[[1]]
    consData$Position %<>% as.factor

    consData$Variants <- ifelse(consData$`Max Non-ref Allele Count` >= 5, "Variant","Background")
  }
  else{
    consData <- object@variants
    consData$Variants <- ifelse(consData$p.adjust <= 0.05, "Variant","Background")
  }

  metaData <- as_tibble(object@meta.data)
  metaData$replicate %<>% as.factor
  metaData$sample %<>% as.character

  # Join meta data and consData replicate ID column to consData
  # Change to left_join?
  summaryData <- dplyr::inner_join(consData, metaData, by = c(`Sample Name` = "sample"))
  summaryData <- summaryData %>%
    dplyr::filter(.data$Variants == "Variant") %>%
    tidyr::unite(.data$Position, .data$Contig, .data$Position, sep = ":") %>%
    dplyr::group_by(.data$`Sample Name`, .data$Position, .data$Name) %>%
    dplyr::summarise()

}


#' Add metaData
#' @export
#' @param object R object to which meta data should be added
#' @param attributeName Name of the meta data attribute.
#' @param attributeValue Meta data to be saved.
addMetaData <- function(object,attributeName,attributeValue){
  attr(x = object, attributeName) <- attributeValue
  return(object)
}

#' Retrieve meta data by name.
#' @export
#' @param object R object from which to get meta data.
#' @param attributeName Name of the meta data attribute.
getMetaData <- function(object,attributeName){
  if(attributeName %in% names(attributes(object))){
    value <- attributes(object)[names(attributes(object)) == attributeName][[1]]
    return(value)
  }
  else{
    message <- simpleError("Attribute not found in object.")
    print(message)
  }
}

#' Generate VCF file from UMI sample or UMI experiment object
#' @export
#' @param object Requires a UMI sample or UMI experiment object
#' @param outDir String. Output directory, defaults to wokring directory.
#' @param outFile String. Name of the output file
#' @param printAll Logical. Should all or only trusted variant be printed?
generateVCF <- function(object, outDir = getwd(), outFile, printAll = FALSE){

  cons.table <- object@cons.table
  cons.table$Variants <- ifelse(cons.table$`Max Non-ref Allele Count` >= 5, "Variant","Background")

  header <- c("##fileformat=VCFv4.3",
              paste("##fileDate=", Sys.Date(),sep=""),
              "##source=umiAnalyzerv0.3.0",
              "##INFO=<ID=DP,Number=1,Type=Integer,Description='Total Depth'>",
              "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele Frequency'>",
              "##INFO=<ID=PADJ,Number=A,Type=Float,Description='FDR-adjusted p-value'>",
              #paste("##INFO=<ID=SAMPLE,Number=A,Type=Float,Description=",
              #      sample.name,">"),
              paste("#",
                    paste("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO",
                          collapse = "\t"), sep=""))
  lines = c()
  lines = append(lines, header)

  if(printAll == FALSE){

    cons.table = cons.table %>% dplyr::filter(.data$Variants == "Variant")

  }

  for(i in 1:nrow(cons.table)){
    row = cons.table[i,]

    if(row$Variants == "Variant"){
      vcfRow = paste(row$Contig, row$Position , ".", row$Reference, row$`Max Non-ref Allele`,
                     "PASS", paste("DP=",row$Coverage,";",
                                   "AF=",row$`Max Non-ref Allele Frequency`,
                                   sep=""),
                     collapse = "\t")
      lines = append(lines, vcfRow)
    }
    else {

      if(is.na(row$`Max Non-ref Allele`)){
        vcfRow = paste(row$Contig, row$Position , ".", row$Reference, ".",
                       "FAIL", paste("DP=",row$Coverage,";",
                                     "AF=",row$`Max Non-ref Allele Frequency`,
                                     sep=""),
                       collapse = "\t")
        lines = append(lines, vcfRow)
      }
      else {
        vcfRow = paste(row$Contig, row$Position , ".", row$Reference, row$`Max Non-ref Allele`,
                       "FAIL", paste("DP=",row$Coverage,";",
                                     "AF=",row$`Max Non-ref Allele Frequency`,
                                     sep=""),
                       collapse = "\t")
        lines = append(lines, vcfRow)
      }
    }
  }

  fileConn<-file(file.path(outDir,paste(outFile,".vcf",sep="")))
  writeLines(lines, fileConn)
  close(fileConn)

  return(object)
}



