setOldClass(c("tbl_df", "tbl", "data.frame"))

#' UMIsample class
#'
#' @exportClass UMIsample
#'
#' @slot name Sample name
#' @slot cons.data Raw consensus data
#' @slot summary.data Summary data from umierrorcorrect
#' @slot reads Consensus reads imported from a bam file.
UMIsample <- setClass(
  "UMIsample",

  slots = list(
    name = "character",
    cons.data = "tbl_df",
    summary.data = "tbl_df",
    reads = "tbl_df"
  )
)

#' UMIexperiment class
#'
#' The UMIexperiment is the core data object, storing all data and relevant
#' analysis data associated with your experiment. Each object has number of
#' slots storing raw data, graphs and processed data.
#'
#' @exportClass UMIexperiment
#' @import tibble
#'
#' @slot name Optional project name for record keeping.
#' @slot cons.data The raw consensus data supplied by the user.
#' @slot summary.data Summary data from umierrorcorrect.
#' @slot reads Consensus reads imported using the parseBamFiles function.
#' @slot meta.data Sample data optionally supplied by the user.
#' @slot filters A list of filtered cons.data, which can be accessed separately.
#' @slot plots A list of generated plots.
#' @slot variants Consensus table generated with the umiAnalyzer variant caller.
#' @slot merged.data Data generated using the mergeTechnicalReplicates function.
UMIexperiment <- setClass(
  Class = "UMIexperiment",

  # Define the slots
  slots = list(
    name = "ANY",
    cons.data = "tbl_df",
    summary.data = "tbl_df",
    reads = "tbl_df",
    meta.data = "data.frame",
    filters = "list",
    plots = "list",
    variants = "tbl_df",
    merged.data = "tbl_df"
  ),

  # Set the default values for the slots. (optional)
  prototype = list(
    name = NULL,
    cons.data = NULL,
    summary.data = NULL,
    reads = tibble(),
    meta.data = data.frame(),
    filters = list(),
    plots = list(),
    variants = tibble(),
    merged.data = tibble()
  )
)
