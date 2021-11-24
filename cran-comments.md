## Second resubmission

The following points were addressed:

* Please add small executable examples in your Rd-files if possible.

Added examples to exported functions. Longer running or less important
examples were enclosed in `\dontrun{}`.

* Please do not modify the .GlobalEnv. 

Removed modification of the .GlobalEnv in the function runUmiVisualizer().


## First resubmission

* Please omit "+ file LICENSE" 

Fixed.

* The Date field is over a month old.

This was changed to 2021-11-23.


## Tested environments
* Local macOS installation (Intel), Monterey 12.0.1, R 4.1.1
* Apple Silicon (M1), macOS 11.6 Big Sur, R-release
* Windows Server 2022, UCRT R-devel, 64 bit
* Debian Linux, R-devel, GCC
* Ubuntu Linux 20.04.1 LTS, R-devel, GCC

## R CMD check results
There were no ERRORs, WARNINGs or NOTEs.
