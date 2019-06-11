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
#' library(umiAnalyzer)
#' main = system.file("extdata", package = "umiAnalyzer")
#' sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#' exp1 <- create.UMIexperiment(experiment.name = "exp1", main.dir = main, dir.names = sample.names)
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
#' library(umiAnalyzer)
#' data <- UMIexperiment
#' data <- filterUMIobject(data)
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

  attr(object@cons.data, "filter") <- filter
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
#' @param object A UMierrorcorrect object.
#' @return Object containing raw and FDR-adjusted P-Values
#' @examples
#' library(umiAnalyzer)
#' data <- UMIexperiment
#' data <- callVariants(data)
callVariants <- function(object){

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
  return(object)
}


#' Generate QC plots
#' @param object Requires a UMI sample or UMI experiment object
generateQCplots <- function(object){
  cons.table <- object@cons.data
  hist.table <- object@hist.data
  summary.table <- object@summary.data

  # Consensus depth plot per assay

  # Plot consensus depth distribution

  # Downsampling plots

  return(object)
}


#' Generate VCF file from UMI sample or UMI experiment object
#' @param object Requires a UMI sample or UMI experiment object
generateVCF <- function(object){
  cons.table <- object@cons.table


  return(object)
}


