# Define sample class
UMIsample <- setClass("UMIsample",
                      slots = list(name = "character",
                                   cons.data = "data.frame",
                                   hist.data = "data.frame",
                                   summary.data = "data.frame")
)

# Define experiment class
UMIexperiment <- setClass("UMIexperiment",
                          slots = list(name = "character",
                                       sample.list = "list",
                                       cons.data = "data.frame",
                                       hist.data = "data.frame",
                                       summary.data = "data.frame")
)

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

create.UMIexperiment <- function(experiment.name,dir.names){

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
    summary.data.merged <- rbind(summary.data.merged,cons)
  }

  UMIexperiment <- UMIexperiment(name = experiment.name,
                                 sample.list = sample.list,
                                 cons.data = cons.data.merged,
                                 hist.data = hist.data.merged,
                                 summary.data = summary.data.merged)
  return(UMIexperiment)
}


# Filter data
filterUMIobject <- function(object, minDepth, minCoverage, minFreq){
  cons.table <- object@cons.data

  cons.table <- cons.table[cons.table$Consensus.group.size == minDepth,]
  cons.table <- cons.table[cons.table$Coverage >= minCoverage,]
  cons.table <- cons.table[cons.table$Name != "",]
  cons.table <- cons.table[cons.table$Max.Non.ref.Allele.Frequency >= minFreq,]

  object@cons.data <- cons.table
  return(object)
}

# Generate QC plots
generateQCplots <- function(object){
  cons.table <- object@cons.table

  return(object)
}
