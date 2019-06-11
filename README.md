# umiAnalyzer

Tools for analyzing umierrorcorrect output

Requirements
------------

- R (>= 3.5), which can be downloaded and installed via The Comprehensive R Archive Network [CRAN](https://cran.r-project.org/).

Installation from R using devtools
------------

In the R studio console, type:

```
require(devtools)

install_github("ozimand1as/umiAnalyzer")
```

How to use umiAnalyzer
-----------

```
library(umiAnalyzer)

main = "/Path_to_Project"

sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

exp1 <- create.UMIexperiment(experiment.name = "exp1",
                             main.dir = main,
                             dir.names = sample.names)

exp1 <- filterUMIobject(exp1)

```
