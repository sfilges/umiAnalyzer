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
# Define main directory
main = "/Users/Project/"

# List umierrorcorrect output directories within main
dir.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

# Create a UMI experiment object
experiment <- create.UMIexperiment(experiment.name = "exp1", dir.names = dir.names)
```
