## ---- echo = FALSE-------------------------------------------------------
  knitr::opts_chunk$set(collapse = TRUE, 
                        comment = "#>")

## ----example1------------------------------------------------------------
library(umiAnalyzer)

main = system.file("extdata", package = "umiAnalyzer")

sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

exp1 <- create.UMIexperiment(experiment.name = "exp1",
                             main.dir = main,
                             dir.names = sample.names)

## ----example2------------------------------------------------------------
data <- simsen
data <- callVariants(data)

vars <- data@cons.data
vars <- vars[vars$p.adjust <= 0.1,]

## ----design--------------------------------------------------------------
data <- simsen
data <- callVariants(data)

metaData <- system.file("extdata", "metadata.txt", package = "umiAnalyzer")

data <- importDesign(object = data, file = metaData)

## ----getmetadata---------------------------------------------------------
design <- getMetaData(object = data, attributeName = "design")

design

## ----addmetadata---------------------------------------------------------
comment <- "fix this"
data <- addMetaData(object = data, attributeName = "my-comment", attributeValue = comment)

myattribute <- getMetaData(object = data, attributeName = "my-comment")
myattribute

## ----QC plots------------------------------------------------------------
data <- generateQCplots(data, do.plot = TRUE, group.by = "assay")



