# umiAnalyzer

Tools for analyzing umierrorcorrect output

Installation
------------

In the R studio console, type:

```
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
