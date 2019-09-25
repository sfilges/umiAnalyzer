#' Generate QC plots
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats median
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param group.by String. Which variable should be used as a factor on the x-axis. Default is assay.
#'
generateQCplots <- function(object,
                            do.plot = TRUE,
                            group.by = "assay") {
  cons.table <- object@cons.data
  summary.table <- object@summary.data

  # Consensus depth plot per assay

  cdepths <- summary.table %>% dplyr::filter(
    .data$assay != "",
    .data$depth == 3
  )

  cdepths$assay %<>% as.factor
  cdepths$sample %<>% as.factor

  # From the ggplot2 vignette:
  # https://github.com/tidyverse/ggplot2/releases
  # aes_() replaces aes_q(). It also supports formulas, so the most concise
  # SE version of aes(carat, price) is now aes_(~carat, ~price). You may
  # want to use this form in packages, as it will avoid spurious R CMD check
  # warnings about undefined global variables.

  if (group.by == "assay") {
    depth_plot <- ggplot(cdepths, aes_(x = ~assay, y = ~UMIcount)) +
      geom_boxplot(outlier.colour = "black", outlier.shape = 10, outlier.size = 3) +
      theme(axis.text.x = element_text(angle = 90)) +
      geom_hline(yintercept = median(cdepths$UMIcount), linetype = "dashed", color = "red") +
      geom_hline(yintercept = mean(cdepths$UMIcount), linetype = "dashed", color = "blue") +
      labs(
        title = "Consensus 3 depths by assay",
        subtitle = paste(
          "Mean depth: ", round(mean(cdepths$UMIcount)),
          "Median depth: ", round(median(cdepths$UMIcount))
        ),
        caption = ""
      )
  } else if (group.by == "sample") {
    depth_plot <- ggplot(cdepths, aes_(x = ~sample, y = ~UMIcount)) +
      geom_boxplot(outlier.colour = "black", outlier.shape = 10, outlier.size = 3) +
      theme(axis.text.x = element_text(angle = 90)) +
      geom_hline(yintercept = median(cdepths$UMIcount), linetype = "dashed", color = "red") +
      geom_hline(yintercept = mean(cdepths$UMIcount), linetype = "dashed", color = "blue") +
      labs(
        title = "Consensus 3 depths by sample",
        subtitle = paste(
          "Mean depth: ", round(mean(cdepths$UMIcount)),
          "Median depth: ", round(median(cdepths$UMIcount))
        ),
        caption = ""
      )
  }

  summary.table <- as_tibble(summary.table)

  cons0.depths <- summary.table %>%
    dplyr::filter(.data$assay != "", .data$depth == 0) %>%
    dplyr::select(.data$assay, .data$region, .data$sample, .data$totalCount)

  cons3.UMIcount <- summary.table %>%
    dplyr::filter(.data$assay != "", .data$depth == 3) %>%
    dplyr::select(.data$assay, .data$region, .data$sample, .data$UMIcount)

  avg.depths <- dplyr::left_join(cons0.depths, cons3.UMIcount, c(
    "region" = "region",
    "assay" = "assay",
    "sample" = "sample"
  ))

  avg.depths <- avg.depths %>% dplyr::mutate(avg.FamDepth = .data$totalCount / .data$UMIcount)

  avg.depths_plot <- ggplot(avg.depths, aes_(x = ~assay, y = ~avg.FamDepth)) +
    geom_boxplot(
      outlier.colour = "red", outlier.shape = 8,
      outlier.size = 4
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Average family size per assay")

  # Plot consensus depth distribution

  if (do.plot) {
    print(depth_plot)
    print(avg.depths_plot)
    object <- addMetaData(object = object, attributeName = "depth_plot", depth_plot)
  }
  else {
    object <- addMetaData(object = object, attributeName = "depth_plot", depth_plot)
  }

  return(object)
}

#' Generate Amplicon plots
#' @export
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr filter
#' @param object Requires a UMI sample or UMI experiment object
#' @param filter.name Name of the filter to be plotted.
#' @param do.plot Logical. Should plots be shown.
generateAmpliconPlots <- function(object,filter.name, do.plot = TRUE){

  # Check if variant caller has been run on object
  if( identical(dim(object@variants),dim(tibble())) ) {
    cons.table <- getFilter(object = object, name = filter.name)
    cons.table <- cons.table[[1]]

    cons.table$Variants <- ifelse(cons.table$`Max Non-ref Allele Count` >= 5, "Variant","Background")
  }
  else{
    cons.table <- object@variants
    cons.table$Variants <- ifelse(cons.table$p.adjust <= 0.05, "Variant","Background")
  }

  cons.table$Position %<>% as.factor
  cons.table$Name %<>% as.factor
  cons.table$sample %<>% as.factor

  # If the plot is too big, limit number of positions plotted;
  # also output tabular output as an html table
  if(length(unique(cons.table$Name)) > 3){

    amplicon_plot <- ggplot(cons.table, aes_(x=~Name,
                                             y=~(100 * `Max Non-ref Allele Frequency`) )) +
      geom_point(aes(col=.data$Variants, size = .data$`Max Non-ref Allele Count`)) + theme_bw() +
      theme(axis.text.x = element_text(size = 8, angle = 90)) +
      ylab("Variant Allele Frequency (%)") +
      xlab("Assay") +
      labs(title = "Maximum variant allele frequency by assay",
           caption = "Each dot represenst a position. All samples are included.
           Blue dots represent positions with at least 5 variant alleles.")
  }
  else{
    amplicon_plot <- ggplot(cons.table, aes_(x=~Position,
                                             y=~`Max Non-ref Allele Frequency`,
                                             fill=~Variants)) +
      geom_bar(stat="identity") +
      theme(axis.text.x = element_text(size = 2, angle = 90)) +
      facet_grid(`Sample Name` ~ Name, scales = "free_x", space = "free_x")

  }

  # Show plot and add ggplot object to the UMIexperiment
  if(do.plot){
    print(amplicon_plot)
    object <- addMetaData(object = object, attributeName = "amplicon_plot", amplicon_plot)
  }
  else{
    object <- addMetaData(object = object, attributeName = "amplicon_plot", amplicon_plot)
  }
  return(object)
}

#' Generate Amplicon plots
#' @export
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
viz_Merged_data <- function(object, do.plot = TRUE){

  # Plotting maximum alternate alle count on merged data
  data <- object@merged.data
  data$Position %<>% as.factor

  data$Variants <- ifelse(data$avg.MaxAC > 5, "Variant","Background")

  plot <- ggplot(data, aes_(x=~Position, y=~avg.MaxAC,fill=~Variants)) +
    geom_bar(stat="identity") +
    geom_errorbar(aes_(ymin=~avg.MaxAC, ymax=~avg.MaxAC+std.MaxAC), width=.1) +
    theme(axis.text.x = element_text(size = 2, angle = 90)) +
    facet_grid(replicate ~ Name, scales = "free_x", space = "free_x")

  print(plot)
}

#' Generate consensus depths plots
#' @export
#' @import ggplot2
#' @param object Requires a UMI sample or UMI experiment object
plotFamilyHistogram <- function(object) {
  if (class(object)[1]== "UMIexperiment") {
    reads <- object@reads

    cons_depth_plot <- ggplot(reads, aes(x = count, fill = sample, color = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      xlim(0, 100) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")
      ) +
      facet_wrap(~sample)
  } else {
    cons_depth_plot <- ggplot(object, aes(x = count, fill = sample, color = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      xlim(0, 100) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")
      ) +
      facet_wrap(~sample)
  }

  plot(cons_depth_plot)
}


#' Plot coverage before and after normalization
#' @importFrom gridExtra grid.arrange
#' @importFrom magrittr "%>%" "%<>%"
#' @param cons.data Consensus data table
#' @return A ggplot object.
viz_Normalization <- function(cons.data){

  cons.data$Name %<>% as.factor

  # plot coverage before nomralization per assay and group by replicate
  p1 <- ggplot(cons.data, aes_(x=~Name, y=~Coverage)) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("Before normalization") +
    facet_grid(. ~ replicate, scales = "free_x", space = "free_x")

  # plot coverage after normalization per assay and group by replicate
  p2 <- ggplot(cons.data, aes_(x=~Name, y=~normCoverage)) +
    geom_boxplot(outlier.colour="red", outlier.shape=8,
                 outlier.size=4) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle("After normalization") +
    facet_grid(. ~ replicate, scales = "free_x", space = "free_x")

  # group replicates
  merged <- grid.arrange(p1, p2, nrow = 1)

  return(merged)
}

#' Plot coverage before and after normalization
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @param cons.data A consensus data table
#' @return A ggplot object.
# Remove counts for reference allele for plotting
viz_stacked_counts <- function(cons.data){

  out.file <- tibble()
  for(j in 1:nrow(cons.data)) {
    row <- cons.data[j,]
    if( row$Reference == "A" ) {
      row$avg.A <- 0
    }
    else if( row$Reference == "C" ) {
      row$avg.C <- 0
    }
    else if( row$Reference == "G" ) {
      row$avg.G <- 0
    }
    else if( row$Reference == "T" ) {
      row$avg.T <- 0
    }
    out.file <- dplyr::bind_rows(out.file, row)
  }

  # Plot normalised counts
  out.file <- out.file %>% dplyr::select(.data$Name, .data$Position, .data$replicate,
                                         .data$avg.Depth,.data$Reference,
                                         .data$avg.A,.data$avg.T,.data$avg.C,.data$avg.G,
                                         .data$avg.N,.data$avg.I,.data$avg.D) %>%
    dplyr::rename(">A"= .data$avg.A,
                  ">T"= .data$avg.T,
                  ">C"= .data$avg.C,
                  ">G"= .data$avg.G,
                  ">N"= .data$avg.N,
                  ">I"= .data$avg.I,
                  ">D"= .data$avg.D) %>%
    tidyr::gather("variant","count", -c(.data$Name,
                                        .data$Position,
                                        .data$replicate,
                                        .data$Reference,
                                        .data$avg.Depth))

  # Stacked
  stacked <- ggplot(out.file, aes_(fill=~variant, y=~count, x=~Position)) +
    geom_bar( stat="identity") +
    facet_grid(replicate ~ Name, scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90))

  return(stacked)
}


