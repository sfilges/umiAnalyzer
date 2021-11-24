## ---- echo = FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=9,
  fig.height=6
)

## ----runApp, eval=FALSE-------------------------------------------------------
#  umiAnalyzer::runUmiVisualizer()

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

BarcodeFamilyHistogram(reads)

## ----example1continued, eval=TRUE---------------------------------------------
QCplot(simsen)

## ----filerobject, eval=TRUE---------------------------------------------------
simsen <- filterUmiObject(simsen)

## ----example1continued_2, eval=TRUE-------------------------------------------
# This is optional
simsen <- callVariants(simsen)

## ----getFilter, eval=TRUE-----------------------------------------------------
myfilter <- getFilteredData(simsen)
myfilter

## ----example1continued_3, eval=TRUE-------------------------------------------
UmiCountsPlot(simsen)

## ---- eval=FALSE--------------------------------------------------------------
#  AmpliconPlot(simsen)

## ---- eval=FALSE--------------------------------------------------------------
#  AmpliconPlot(
#    object = simsen,
#    amplicons = 'KIT_125',
#    samples = 'VAF-1-5ng-1-10x'
#  )

## ----metaDataImport, eval=FALSE-----------------------------------------------
#  metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")
#  
#  simsen <- importDesign(
#    object = simsen,
#    file = metaData
#  )
#  
#  design <- getMetaData(
#    object = simsen,
#    attributeName = "design"
#  )
#  
#  design

## ----time_course, eval=FALSE--------------------------------------------------
#  simsen <- umiAnalyzer::analyzeTimeSeries(
#    object = simsen,
#    time.var = 'time',
#    group.by = 'replicate'
#  )

## ----replicates, eval=FALSE---------------------------------------------------
#  simsen <- mergeTechnicalReplicates(
#    object = simsen,
#    group.by = 'replicate',
#    option = 'Set1',
#    amplicons = c('new', 'TP53_1', 'TP53_7')
#  )

## ----vizMerge, eval=FALSE-----------------------------------------------------
#  vizMergedData(
#    simsen,
#    amplicons = c('new')
#  )

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

