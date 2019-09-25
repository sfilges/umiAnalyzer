# umiAnalyzer

Tools for analyzing sequencing data containing UMI generated by umierrorcorrect using comprehensive easy-to-use methods. The package allows merging of multiple samples into a single UMIexperiment object which can be easily manipulated using build-in functions to generate tabular and graphical output.

This readme serves as a basic introduction, for more complete informations see the R vignette using:
```
browseVignettes("umiAnalyzer")
```

Requirements
------------

- R (>= 3.5), which can be downloaded and installed via The Comprehensive R Archive Network [CRAN](https://cran.r-project.org/).
- Installation from R using install_github requires the devtools package
- Installation of the package will, if necessary, automatically install its dependencies: ggplot2, methods, utils, VGAM, stats

Installation from R using devtools
------------

In the R studio console, type:

```
devtools::install_github("ozimand1as/umiAnalyzer")
```

How to make your own umiExperiment object
---------------------

Define a variable containing the path to the directory with all the umierrorcorrect output folders 
belonging to your experiment. umiAnalyzer comes with raw test data generated with umierrorcorrect that 
you can import if you don't have any of your own.

Call the create.UMIexperiment to create your umiExperiment object.

```
library(umiAnalyzer)

main <- system.file("extdata", package = "umiAnalyzer")

sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

simsen <- createUMIexperiment(experimentName = "test",
                              mainDir = main,
                              sampleNames = sample.names)

# save(simsen, file = "data/simsen.RData")

reads <- parseBamFiles(mainDir = main, sampleNames = sample.names, consDepth = 10)

consDepthsPlot(reads)
```

A basic example
-----------------

umiAnalyzer comes with a build-in umiExperiment object to explore, which was generated using the code 
above, so it can be used without creating it first, if so desired.

In order to call variants using the umiAnalyzer variant caller simply load the package and test data
and use the callVariants function. You can then filter the resulting consensus data (cons.data) within
the object, e.g. for significant variants.

```
data <- simsen
data <- callVariants(data)

vars <- data@cons.data
vars <- vars[vars$p.adjust <= 0.1,]
```




