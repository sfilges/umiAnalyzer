## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, 
                        comment = "#>",
                        fig.width=9, 
                        fig.height=6)

## ----example1, eval=FALSE------------------------------------------------
#  library(umiAnalyzer)
#  
#  main <- system.file("extdata", package = "umiAnalyzer")
#  
#  sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
#  
#  simsen <- createUMIexperiment(
#    experiment.name = "test",
#    main.dir = main,
#    dir.names = sample.names,
#    importBam = FALSE
#  )
#  
#  # save(simsen, file = "data/simsen.RData")
#  
#  reads <- parseBamFiles(main.dir = main)
#  
#  simsen <- generateQCplots(simsen, do.plot = TRUE, group.by = "assay")
#  
#  simsen <- filterUMIobject(
#    object = simsen, name = "myfilter", minDepth = 3,
#    minCoverage = 100, minFreq = 0, minCount = 0
#  )
#  
#  myfilter <- getFilter(object = simsen, name = "myfilter")
#  myfilter
#  
#  simsen <- generateAmpliconPlots(object = simsen, filter.name = "myfilter", do.plot = TRUE)
#  
#  metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#  simsen <- importDesign(object = exp1, file = metaData)
#  
#  simsen <- merge_technical_replicates(object = simsen, filter.name = "myfilter")
#  simsen@merged.data
#  
#  viz_Merged_data(simsen)

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

