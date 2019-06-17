#' Define sample class
UMIsample <- setClass("UMIsample",
                      slots = list(name = "character",
                                   cons.data = "data.frame",
                                   summary.data = "data.frame")
)

#' Define experiment class
UMIexperiment <- setClass("UMIexperiment",
                          slots = list(name = "character",
                                       sample.list = "list",
                                       cons.data = "data.frame",
                                       summary.data = "data.frame")
)

#' Method for creating a UMI sample
#' @export
#' @importFrom methods new
#' @importFrom utils read.csv
#' @param sample.name UMI sample object name
#' @param sample.dir Path to UMI sample
create.UMIsample <- function(sample.name,sample.dir){
  cons.file <- list.files(path = sample.dir,pattern = "\\.cons$")

  cons.table <- read.csv(file = file.path(sample.dir,cons.file),
                         sep = "\t", row.names = NULL,
                         header = TRUE)

  summary.file <- list.files(path = sample.dir,pattern = "\\.txt$")
  summary.table <- read.csv(file = file.path(sample.dir,summary.file),
                            sep = "\t", row.names = NULL,
                            header = FALSE)

  UMI.sample <- UMIsample(name = sample.name,
                          cons.data = cons.table,
                          summary.data = summary.table)
}

#' Method for creating a UMI experiment object
#' @export
#' @importFrom utils read.csv
#' @importFrom methods new
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
create.UMIexperiment <- function(experiment.name,main.dir,dir.names){
  main = main.dir
  sample.list = list()
  cons.data.merged = data.frame()
  summary.data.merged = data.frame()

  for(i in 1:length(dir.names)){

    sample.list[[i]] <- create.UMIsample(dir.names[i], file.path(main,dir.names[i]))

    sample <- sample.list[[i]]

    cons <- sample@cons.data
    cons$sample <- dir.names[i]
    cons.data.merged <- rbind(cons.data.merged,cons)

    summary <- sample@summary.data
    summary$sample <- dir.names[i]
    summary.data.merged <- rbind(summary.data.merged,summary)
  }
  colnames(summary.data.merged) <- c("ID","region","assay","depth","fraction","totalCount","UMIcount","sample")

  UMIexperiment <- UMIexperiment(name = experiment.name,
                                 sample.list = sample.list,
                                 cons.data = cons.data.merged,
                                 summary.data = summary.data.merged)
  return(UMIexperiment)
}


#' Method for filtering UMIexperiment and sample objects
#' @export
#' @param object Requires a UMI sample or UMI experiment object.
#' @param minDepth Consensus depth to analyze. Default is 3.
#' @param minCoverage Mininum coverage required for amplicons. Default is 1.
#' @param minFreq Minimum variant allele frequency to keep. Default is 0.
#' @param minCount Minimum variant allele count to keep. Default is 3.
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' data <- simsen
#' data <- filterUMIobject(data)
#' }
filterUMIobject <- function(object, minDepth=3, minCoverage=50, minFreq=0, minCount=3){
  cons.table <- object@cons.data

  cons.table <- cons.table[cons.table$Consensus.group.size == minDepth,]
  cons.table <- cons.table[cons.table$Coverage >= minCoverage,]
  cons.table <- cons.table[cons.table$Name != "",]
  cons.table <- cons.table[cons.table$Max.Non.ref.Allele.Frequency >= minFreq,]
  cons.table <- cons.table[cons.table$Max.Non.ref.Allele.Count >= minCount,]

  object@cons.data <- cons.table

  filter = list(minDepth=minDepth,
                minCoverage=minCoverage,
                minFreq=minFreq,
                minCount=minCount)

  object <- addMetaData(object = object,
                        attributeName = "filter",
                        attributeValue = filter)
  return(object)
}

#' Calculates the negative log likelihood for the beta distribution.
#' @importFrom stats dbeta
#' @param params Non-negative parameters of the Beta distribution.
#' @param data consensus.data table of a UMisample or UMIexperiment object.
#' @return Negative log-likelihood for beta distribution.
betaNLL <- function(params,data){
  a<-params[1]
  b<-params[2]
  # negative log likelihood for beta
  return(-sum(dbeta(data,shape1=a, shape2=b, log=TRUE)))
}

#' Calculate variant p-values using permutation-based testing. A prior is fitted to model the background
#' error using maximum likelihood estimation of a beta distribution. The maximum likelihood estimate
#' of the beta distribution is then used to define the shape of a beta-binomial distribution used
#' to estimate variant P-Values. This can be interpreted as a probability for a variant to not have
#' arisen by chance
#' @export
#' @importFrom VGAM rbetabinom.ab
#' @importFrom stats nlm var p.adjust
#' @importFrom utils install.packages
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

  object <- filterUMIobject(object = object,
                            minDepth = 3, # Require consensus 3
                            minCoverage = 100, # Require at least 100 cons reads
                            minFreq = 0, # no minimum allele freq
                            minCount = 0) # no minimum variant allele count

  cons.table <- object@cons.data

  a1 <- cons.table$Coverage*cons.table$Max.Non.ref.Allele.Frequency # No. of variant alleles
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

  for (i in 1:length(a1)){ # for each named amplicon position:
    r1 <- rbetabinom.ab(10000,b1[i],shape1=a,shape2=b) # Calculate probability of success
    pval[i] = sum(r1>a1[i])/10000                      # Estimate p value of variant
  }

  cons.table$pval <- pval
  cons.table$p.adjust <-p.adjust(pval,method="fdr")

  object@cons.data <- cons.table

  object <- addMetaData(object = object, attributeName = "varCalls", "varCalls")

  return(object)
}


#' Filter variants based on p values or depth
#' @export
#' @param object A UMIexperiment object
#' @param p.adjust Adjusted p value (FDR)
#' @param minDepth Minimum consensus depth required
#' @return A UMIexperiment object with filtered variants
filterVariants <- function(object, p.adjust = 0.2, minDepth = 50){
  if("varCalls" %in% names(attributes(object))){
    vars <- object@cons.data
    vars <- vars[vars$p.adjust <= p.adjust,]
    vars <- vars[vars$minDepth >= minDepth,]

    object@cons.data <- vars
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

  object <- addMetaData(object = object, attributeName = "design", mData)
  return(object)
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
#' @param object Requires a UMI sample or UMI experiment object
#' @param outFile String. Name of the output file
generateVCF <- function(object, outFile){

  # This should be moved to a function annotateVCF

  #if (!requireNamespace("BiocManager", quietly = TRUE)) {
  #  install.packages("BiocManager")
  #  if (!requireNamespace("VariantAnnotation"), quietly = TRUE)){
  #     BiocManager::install("VariantAnnotation")
  #  }
  #}

  cons.table <- object@cons.table


  header <- c("##fileformat=VCFv4.3",
              paste("##fileDate=", Sys.Date(),sep=""),
              "##source=umiAnalyzerv0.3.0",
              "##INFO=<ID=DP,Number=1,Type=Integer,Description='Total Depth'>",
              "##INFO=<ID=AF,Number=A,Type=Float,Description='Allele Frequency'>",
              "##INFO=<ID=PADJ,Number=A,Type=Float,Description='FDR-adjusted p-value'>",
              #paste("##INFO=<ID=SAMPLE,Number=A,Type=Float,Description=",
              #      sample.name,">"),
              "#CHROM POS ID  REF ALT QUAL  FILTER  INFO  FORMAT")

  #fileConn<-file("output.txt")
  #writeLines(header, fileConn)
  #close(fileConn)


  return(object)
}

#' Downsample raw reads
#' @importFrom ShortRead readFastq
#' @param fastq Path to fastq containing UMIs in header
downSample <- function(fastq){

  #import fastq and extract header
  reads <- readFastq(fastq)

  #samples <- sample(1:ncol(data), N_samples, rep=TRUE, prob=data[row,])


}


