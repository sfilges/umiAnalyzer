## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, 
                        comment = "#>",
                        fig.width=9, 
                        fig.height=6)

## ----example1, eval=FALSE------------------------------------------------
#  library(umiAnalyzer)
#  
#  main = system.file("extdata", package = "umiAnalyzer")
#  
#  sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#  
#  exp1 <- createUMIexperiment(experiment.name = "exp1",
#                               main.dir = main,
#                               dir.names = sample.names)
#  
#  #simsen <- exp1
#  #save(simsen, file = "data/simsen.RData")
#  
#  exp1 <- generateQCplots(exp1, do.plot = TRUE, group.by = "assay")
#  
#  exp1 <- filterUMIobject(object = exp1, name = "myfilter", minDepth = 3,
#                          minCoverage = 100, minFreq = 0, minCount = 0)
#  
#  myfilter <- getFilter(object = exp1, name = "myfilter")
#  myfilter
#  
#  exp1 <- generateAmpliconPlots(object = exp1, filter.name = "myfilter", do.plot = TRUE)
#  
#  metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#  exp1 <- importDesign(object = exp1, file = metaData)
#  
#  exp1 <- merge_technical_replicates(object = exp1, filter.name = "myfilter")
#  exp1@merged.data
#  
#  viz_Merged_data(exp1)

## ----example2, eval=FALSE------------------------------------------------
#  #data <- simsen
#  #data <- callVariants(data)
#  
#  #data <- filterVariants(object = data, p.adjust = 0.2, minDepth = 5)

## ----design, eval=FALSE--------------------------------------------------
#  data <- simsen
#  data <- callVariants(data)
#  
#  metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#  
#  data <- importDesign(object = data, file = metaData)

## ----getmetadata, eval=FALSE---------------------------------------------
#  design <- getMetaData(object = data, attributeName = "design")
#  
#  design

## ----addmetadata, eval=FALSE---------------------------------------------
#  comment <- "fix this"
#  data <- addMetaData(object = data, attributeName = "my-comment", attributeValue = comment)
#  
#  myattribute <- getMetaData(object = data, attributeName = "my-comment")
#  myattribute

