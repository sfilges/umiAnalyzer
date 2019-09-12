pkgname <- "umiAnalyzer"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
base::assign(".ExTimings", "umiAnalyzer-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('umiAnalyzer')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("callVariants")
### * callVariants

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: callVariants
### Title: Calculate variant p-values using permutation-based testing. A
###   prior is fitted to model the background error using maximum
###   likelihood estimation of a beta distribution. The maximum likelihood
###   estimate of the beta distribution is then used to define the shape of
###   a beta-binomial distribution used to estimate variant P-Values. This
###   can be interpreted as a probability for a variant to not have arisen
###   by chance
### Aliases: callVariants

### ** Examples

## Not run: 
##D library(umiAnalyzer)
##D 
##D data <- simsen
##D data <- callVariants(data)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("callVariants", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("createUMIexperiment")
### * createUMIexperiment

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: createUMIexperiment
### Title: Method for creating a UMI experiment object
### Aliases: createUMIexperiment

### ** Examples

## Not run: 
##D library(umiAnalyzer)
##D 
##D main = system.file("extdata", package = "umiAnalyzer")
##D sample.names <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)
##D 
##D exp1 <- create.UMIexperiment(experiment.name = "exp1", main.dir = main, dir.names = sample.names)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("createUMIexperiment", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("filterUMIobject")
### * filterUMIobject

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: filterUMIobject
### Title: Method for filtering UMIexperiment and sample objects
### Aliases: filterUMIobject

### ** Examples

## Not run: 
##D library(umiAnalyzer)
##D 
##D data <- simsen
##D data <- filterUMIobject(data)
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("filterUMIobject", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
