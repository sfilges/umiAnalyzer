#' Generate QC plots
#' @export
#' @import ggplot2
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param group.by String. Which variable should be used as a factor on the x-axis. Default is assay.
generateQCplots <- function(object, do.plot = TRUE, group.by = "assay"){
  cons.table <- object@cons.data
  summary.table <- object@summary.data

  # Consensus depth plot per assay

  cdepths <- summary.table[summary.table$assay != "",]
  cdepths <- cdepths[cdepths$depth == 3,]

  # From the ggplot2 vignette:
  # https://github.com/tidyverse/ggplot2/releases
  # aes_() replaces aes_q(). It also supports formulas, so the most concise
  # SE version of aes(carat, price) is now aes_(~carat, ~price). You may
  # want to use this form in packages, as it will avoid spurious R CMD check
  # warnings about undefined global variables.

  if(group.by == "assay") {
    depth_plot <- ggplot(cdepths, aes_(x=~as.factor(assay), y=~UMIcount)) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,
                   outlier.size=4) +
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle("Consensus 3 depths by assay")
  } else if (group.by == "sample") {
    depth_plot <- ggplot(cdepths, aes_(x=~as.factor(sample), y=~UMIcount)) +
      geom_boxplot(outlier.colour="red", outlier.shape=8,
                   outlier.size=4) +
      theme(axis.text.x = element_text(angle = 90)) +
      ggtitle("Consensus 3 depths by sample")
  }

  # Plot consensus depth distribution

  if(do.plot){
    print(depth_plot)
    object <- addMetaData(object = object, attributeName = "depth_plot", depth_plot)
  }
  else{
    object <- addMetaData(object = object, attributeName = "depth_plot", depth_plot)
  }

  return(object)
}

#' Generate Amplicon plots
#' @export
#' @import ggplot2
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
generateAmpliconPlots <- function(object, do.plot = TRUE){

  cons.table <- object@variants
  cons.table$Name <- as.factor(cons.table$Name)
  cons.table$sample<- as.factor(cons.table$sample)

  cons.table$Variants <- ifelse(cons.table$p.adjust <= 0.05, "Variant","Background")

  amplicon_plot <- ggplot(cons.table, aes_(x=~as.factor(Position),
                                           y=~Max.Non.ref.Allele.Frequency,
                                           fill=~Variants)) +
    geom_bar(stat="identity") +
    theme(axis.text.x = element_text(size = 2, angle = 90)) +
    facet_grid(sample ~ Name, scales = "free_x", space = "free_x")

  if(do.plot){
    print(amplicon_plot)
    object <- addMetaData(object = object, attributeName = "amplicon_plot", amplicon_plot)
  }
  else{
    object <- addMetaData(object = object, attributeName = "amplicon_plot", amplicon_plot)
  }
  return(object)
}




