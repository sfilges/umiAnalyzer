## ---- echo = FALSE-------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)

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

