#' Define sample class
UMIsample <- setClass("UMIsample",
                      slots = list(name = "character",
                                   cons.data = "data.frame",
                                   hist.data = "data.frame",
                                   summary.data = "data.frame")
)

#' Define experiment class
UMIexperiment <- setClass("UMIexperiment",
                          slots = list(name = "character",
                                       sample.list = "list",
                                       cons.data = "data.frame",
                                       hist.data = "data.frame",
                                       summary.data = "data.frame")
)

#' Method for creating a UMI sample
#' @param sample.name UMI sample object name
#' @param sample.dir Path to UMI sample
create.UMIsample <- function(sample.name,sample.dir){
  cons.file <- list.files(path = sample.dir,pattern = "\\.cons$")

  cons.table <- read.csv(file = file.path(sample.dir,cons.file),
                         sep = "\t", row.names = NULL,
                         header = TRUE)

  hist.file <- list.files(path = sample.dir,pattern = "\\.hist$")
  hist.table <- read.csv(file = file.path(sample.dir,hist.file),
                         sep = "\t", row.names = 1,
                         header = FALSE)

  summary.file <- list.files(path = sample.dir,pattern = "\\.txt$")
  summary.table <- read.csv(file = file.path(sample.dir,summary.file),
                            sep = "\t", row.names = NULL,
                            header = FALSE)

  UMI.sample <- UMIsample(name = sample.name,
                          cons.data = cons.table,
                          hist.data = hist.table,
                          summary.data = summary.table)
}

#' Method for creating a UMI experiment object
#' @param experiment.name Name of the experiment
#' @param main.dir Main experiment directory
#' @param dir.names List of sample names
create.UMIexperiment <- function(experiment.name,main.dir,dir.names){
  main = main.dir
  sample.list = list()
  cons.data.merged = data.frame()
  hist.data.merged = data.frame()
  summary.data.merged = data.frame()

  for(i in 1:length(dir.names)){

    sample.list[[i]] <- create.UMIsample(dir.names[i], file.path(main,dir.names[i]))

    sample <- sample.list[[i]]

    cons <- sample@cons.data
    cons$sample <- dir.names[i]
    cons.data.merged <- rbind(cons.data.merged,cons)

    hist <- sample@hist.data
    hist$sample <- dir.names[i]
    hist.data.merged <- rbind(hist.data.merged,hist)

    summary <- sample@summary.data
    summary$sample <- dir.names[i]
    summary.data.merged <- rbind(summary.data.merged,summary)
  }
  colnames(hist.data.merged) <- c("coordinates","assay","consReads","singletons","sample")
  colnames(summary.data.merged) <- c("ID","region","assay","depth","fraction","counts1","counts2","sample")

  UMIexperiment <- UMIexperiment(name = experiment.name,
                                 sample.list = sample.list,
                                 cons.data = cons.data.merged,
                                 hist.data = hist.data.merged,
                                 summary.data = summary.data.merged)
  return(UMIexperiment)
}


#' Method for filtering UMIexperiment and sample objects
#' @param object Requires a UMI sample or UMI experiment object
#' @param minDepth Consensus depth to analyze. Default is 3
#' @param minCoverage Mininum coverage required for amplicons. Default is 1
#' @param minFreq Minimum variant allele frequency to keep. Default is 0
filterUMIobject <- function(object, minDepth=3, minCoverage=1, minFreq=0){
  cons.table <- object@cons.data

  cons.table <- cons.table[cons.table$Consensus.group.size == minDepth,]
  cons.table <- cons.table[cons.table$Coverage >= minCoverage,]
  cons.table <- cons.table[cons.table$Name != "",]
  cons.table <- cons.table[cons.table$Max.Non.ref.Allele.Frequency >= minFreq,]

  object@cons.data <- cons.table
  return(object)
}

#' Generate QC plots
#' @param object Requires a UMI sample or UMI experiment object
generateQCplots <- function(object){
  requireNamespace("ggplot2", quietly = TRUE)

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



