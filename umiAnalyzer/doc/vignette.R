## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, 
                        comment = "#>",
                        fig.width=9, 
                        fig.height=6)

## ----runApp, eval=FALSE--------------------------------------------------
#  library(umiAnalyzer)
#  
#  runUmiVisualiser()

## ----example1, eval=TRUE-------------------------------------------------
library(umiAnalyzer)

main <- system.file("extdata", package = "umiAnalyzer")

samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

simsen <- createUmiExperiment(experimentName = "simsen",
                              mainDir = main,
                              sampleNames = samples)

reads <- parseBamFiles(mainDir = main, sampleNames = samples, consDepth = 10)

plotFamilyHistogram(reads)

## ----example1continued, eval=TRUE----------------------------------------
simsen <- generateQCplots(simsen, do.plot = TRUE, group.by = "assay")

simsen <- filterUmiobject(
  object = simsen, name = "myfilter", minDepth = 3,
  minCoverage = 100, minFreq = 0, minCount = 0
)

myfilter <- getFilterdData(object = simsen, name = "myfilter")
myfilter

## ----ampliconPlots, eval=TRUE--------------------------------------------
simsen <- generateAmpliconPlots(
  object = simsen,
  filter.name = "myfilter",
  do.plot = TRUE)

simsen <- generateAmpliconPlots(
  object = simsen,
  filter.name = "myfilter",
  do.plot = TRUE, 
  amplicons = c("PIK3CA_123", "PIK3CA_234"),
  samples = "VAF-1-5ng-1-10x")

## ----replicates, eval=TRUE-----------------------------------------------
metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
simsen <- importDesign(object = simsen, file = metaData)

simsen <- mergeTechnicalReplicates(object = simsen, filter.name = "myfilter")
simsen@merged.data

vizMergedData(simsen)

## ----example2, eval=FALSE------------------------------------------------
#  data <- simsen
#  data <- callVariants(data)
#  
#  data <- filterVariants(object = data, p.adjust = 0.2, minDepth = 5)

## ----design, eval=FALSE--------------------------------------------------
#  metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#  
#  data <- importDesign(object = simsen, file = metaData)

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

## ----vcf, eval=FALSE-----------------------------------------------------
#  generateVCF(object = simsen, outFile = 'simsen.vcf', printAll = FALSE, save = FALSE)

## ----csv, eval=FALSE-----------------------------------------------------
#  consensus_data <- saveConsData(object = simsen)
#  
#  outDir <- "~/Documents/"
#  saveConsData(object = simsen, outDir = outDir, delim = ";", save = TRUE)

