## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=9,
  fig.height=6
)

## ----runApp, eval=FALSE-------------------------------------------------------
#  library(umiAnalyzer)
#  
#  runUmiVisualiser()

## ----example1, eval=TRUE------------------------------------------------------
library(umiAnalyzer)

main <- system.file("extdata", package = "umiAnalyzer")

simsen <- createUmiExperiment(main)

simsen <- mergeAssays(
  object = simsen,
  name = "new",
  assay.list = c("PIK3CA_123", "PIK3CA_234")
)

## ----bam-files, eval=TRUE-----------------------------------------------------
reads <- parseBamFiles(main, consDepth = 10)

plotFamilyHistogram(reads)

## ----example1continued, eval=TRUE---------------------------------------------
simsen <- generateQCplots(simsen, group.by = 'assay')

## ----example1continued_2, eval=TRUE-------------------------------------------
simsen <- filterUmiObject(simsen)

myfilter <- getFilteredData(simsen)
myfilter

## ----example1continued_3, eval=TRUE-------------------------------------------
simsen <- plotUmiCounts(object = simsen)

## ---- eval=TRUE---------------------------------------------------------------
simsen <- generateAmpliconPlots(
  object = simsen,
  do.plot = TRUE,
  amplicons = 'KIT_125',
  plot.ref = TRUE,
  plot.text = FALSE
)

## ---- eval=TRUE---------------------------------------------------------------
simsen <- generateAmpliconPlots(
  object = simsen,
  amplicons = c("new"),
  samples = "VAF-1-5ng-1-10x"
)

## ----metaDataImport, eval=TRUE------------------------------------------------
metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")

simsen <- importDesign(
  object = simsen,
  file = metaData
)

design <- getMetaData(
  object = simsen, 
  attributeName = "design"
)

design

## ----time_course, eval=FALSE--------------------------------------------------
#  simsen <- analyzeTimeSeries(
#    object = simsen,
#    time.var = 'time',
#    group.by = 'replicate'
#  )

## ----replicates, eval=TRUE----------------------------------------------------
simsen <- mergeTechnicalReplicates(
  object = simsen,
  group.by = 'replicate',
  option = 'magma',
  direction = -1, 
  amplicons = c('PIK3CA_234', 'TP53_1', 'TP53_7')
)

#viewNormPlot(simsen)
#simsen@plots$stacked_counts
#simsen@merged.data

## ----vizMerge, eval=TRUE------------------------------------------------------

vizMergedData(
  simsen,
  amplicons = c('PIK3CA_234', 'TP53_1', 'TP53_7')
)


## ----example2, eval=FALSE-----------------------------------------------------
#  simsen <- callVariants(simsen)
#  
#  simsen <- filterVariants(simsen)
#  
#  simsen <- generateAmpliconPlots(
#    object = simsen
#  )

## ----vcf, eval=FALSE----------------------------------------------------------
#  generateVCF(object = simsen, outFile = 'simsen.vcf', printAll = FALSE, save = FALSE)

## ----csv, eval=FALSE----------------------------------------------------------
#  consensus_data <- saveConsData(object = simsen)
#  
#  outDir <- getwd()
#  saveConsData(object = simsen, outDir = outDir, delim = ';', save = TRUE)

## ----sessionInfo--------------------------------------------------------------
sessionInfo()

