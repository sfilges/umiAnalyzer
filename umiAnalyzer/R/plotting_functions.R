#' Generate QC plots
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats median
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param group.by String. Which variable should be used as a factor on the x-axis. Default is assay.
#' @param plotDepth Which consensus depth to plot
#' @param assays (Optional) user-supplied list of assays to plot. Default is all.
#' @param samples (Optional) user-supplied list of samples to plot. Default is all.
#'
generateQCplots <- function(
  object,
  do.plot = TRUE,
  group.by = 'assay',
  plotDepth = 3,
  assays = NULL,
  samples = NULL) {

  cons.table <- object@cons.data
  summary.table <- object@summary.data

  # Consensus depth plot per assay

  cdepths <- summary.table %>% dplyr::filter(
    .data$assay != '',
    .data$depth == plotDepth
  )

  if (!is.null(assays)) {
    cdepths <- cdepths %>% dplyr::filter(.data$assay %in% assays)
  }

  if (!is.null(samples)) {
    cdepths <- cdepths %>% dplyr::filter(.data$sample %in% samples)
  }

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
      geom_boxplot(
        outlier.colour = "black",
        outlier.shape = 10,
        outlier.size = 3) +
      geom_jitter(
        size = 8,
        shape = 16,
        position = position_jitter(0.2)) +
      theme_bw() +
      theme(
        axis.text.x = element_text(size = 14, angle = 90),
        axis.text.y = element_text(size = 14)) +
      geom_hline(
        yintercept = median(cdepths$UMIcount),
        linetype = "dashed", color = "red") +
      geom_hline(
        yintercept = mean(cdepths$UMIcount),
        linetype = "dashed", color = "blue") +
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
      geom_boxplot(
        outlier.colour = "black",
        outlier.shape = 10,
        outlier.size = 3) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90)) +
      geom_hline(
        yintercept = median(cdepths$UMIcount),
        linetype = "dashed",
        color = "red") +
      geom_hline(
        yintercept = mean(cdepths$UMIcount),
        linetype = "dashed",
        color = "blue") +
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
    dplyr::filter(
      .data$assay != "",
      .data$depth == 0
    ) %>%
    dplyr::select(
      .data$assay,
      .data$region,
      .data$sample,
      .data$totalCount
    )

  cons3.UMIcount <- summary.table %>%
    dplyr::filter(
      .data$assay != '',
      .data$depth == 3
    ) %>%
    dplyr::select(
      .data$assay,
      .data$region,
      .data$sample,
      .data$UMIcount
    )

  avg.depths <- dplyr::left_join(cons0.depths, cons3.UMIcount, c(
    'region' = 'region',
    'assay' = 'assay',
    'sample' = 'sample'
  ))

  avg.depths <- avg.depths %>% dplyr::mutate(avg.FamDepth = .data$totalCount / .data$UMIcount)

  avg.depths_plot <- ggplot(avg.depths, aes_(x = ~assay, y = ~avg.FamDepth)) +
    geom_boxplot(
      outlier.colour = 'red', outlier.shape = 8,
      outlier.size = 4
    ) +
    theme(axis.text.x = element_text(angle = 90)) +
    ggtitle('Average family size per assay')

  # Plot consensus depth distribution
  if (do.plot) {
    print(depth_plot)
    object <- addMetaData(
      object = object,
      attributeName = "depth_plot",
      depth_plot
    )
  }
  else {
    object <- addMetaData(
      object = object,
      attributeName = "depth_plot",
      depth_plot
    )
  }

  return(object)
}

#' Plot UMI counts
#' @export
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats median
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param amplicons (Optional) user-supplied list of assays to plot. Default is all.
#' @param samples (Optional) user-supplied list of samples to plot. Default is all.
plotUmiCounts <- function(
  object,
  do.plot = TRUE,
  amplicons = NULL,
  samples = NULL
  ) {

  # Read summary data from object
  data <- object@summary.data %>%
    dplyr::filter(!is.na(.data$assay),
                  .data$depth > 0)

  # Select amplicons
  if (!is.null(amplicons)) {
    data <- data %>% dplyr::filter(.data$Name %in% amplicons)
  }

  # Select samples
  if (!is.null(samples)) {
    data <- data %>% dplyr::filter(.data$`Sample Name` %in% samples)
  }

  # Generate ggplot object
  data$depth %<>% as.factor

  plot <- ggplot(data, aes(x=.data$depth, y=.data$UMIcount, fill=sample)) +
    theme_classic() +
    geom_col(alpha=0.6) +
    facet_grid(assay ~ sample)

  # Return plot
  if(do.plot){

    print(plot)
    return(object)

  } else {

    return(object)

  }

}

#' Generate Amplicon plots
#' Plots variant allele frequencies or alternate allele counts for chosen
#' samples and assays.
#' @export
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr filter
#' @param object Requires a UMI sample or UMI experiment object
#' @param filter.name Name of the filter to be plotted.
#' @param do.plot Logical. Should plots be shown?
#' @param cut.off How many variant reads are necessary to consider a variant above background? Default is 5 reads.
#' @param amplicons (Optional) character vector of amplicons to be plotted.
#' @param samples (Optional) character vector of samples to be plotted.
#' @param abs.count Should absolute counts be plotted instead of frequencies?
#' Default is FALSE.
#' @examples
#' \dontrun{
#' library(umiAnalyzer)
#'
#' data <- simsen
#' data <- filterUmiobject(data, "myfilter")
#'
#' data <- generateAmpliconPlots(simsen, "myfilter")
#' }
#' @return A umiExperiment object containing a ggplot object with the
#' amplicon plot.
generateAmpliconPlots <- function(
  object,
  filter.name = "default",
  do.plot = TRUE,
  cut.off = 5,
  amplicons = NULL,
  samples = NULL,
  abs.count = FALSE
  ) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object and filter names")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  } else if(!is.logical(do.plot)){
    warning("do.plot needs to be of type boolean. Using default.")
    do.plot = TRUE
  } else if(!is.logical(abs.count)){
    warning("abs.count needs to be of type boolean. Using defaults.")
    abs.count = FALSE
  } else if (!is.numeric(cut.off) || cut.off < 0) {
    warning("cut.off needs to be a positive integer. Using defaults.")
    cut.off = 5
  } else if(is.null(object@filters[filter.name][[1]])) {
    if(!is.null(object@filters$default)){
      warning("Requested filter not found, using default.")
    } else {
      stop("Filter not found. Have you run filterUmiObject?")
    }
  }

  # Check if variant caller has been run on object
  if (identical(dim(object@variants), dim(tibble()))) {
    cons.table <- getFilteredData(
      object = object,
      name = filter.name
    )

    cons.table$Variants <- ifelse(cons.table$`Max Non-ref Allele Count` >= cut.off, "Variant", "Background")
  }
  else {
    cons.table <- object@variants
    cons.table$Variants <- ifelse(cons.table$p.adjust <= 0.05, "Variant", "Background")
  }

  # Make variables factors to ensure equidistance on the x-axis
  cons.table$`Sample Name` %<>% as.factor
  cons.table$Position %<>% as.factor
  cons.table$Name %<>% as.factor
  cons.table$sample %<>% as.factor

  if (!is.null(amplicons)) {
    cons.table <- cons.table %>%
      dplyr::filter(.data$Name %in% amplicons)
  }

  if (!is.null(samples)) {
    cons.table <- cons.table %>%
      dplyr::filter(.data$`Sample Name` %in% samples)
  }

  # If the plot is too big, limit number of positions plotted;
  # also output tabular output as an html table
  if (length(unique(cons.table$`Sample Name`)) > 6) {
    amplicon_plot <- ggplot(
      cons.table, aes_(
        x = ~Name,
        y = ~ (100 * `Max Non-ref Allele Frequency`))
      ) +
      geom_point(aes(col = .data$Variants, size = .data$`Max Non-ref Allele Count`)) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 8, angle = 90)) +
      ylab("Variant Allele Frequency (%)") +
      xlab("Assay") +
      labs(
        title = "Maximum variant allele frequency by assay",
        caption = "Each dot represenst a position. All samples are included.
           Blue dots represent positions with at least 5 variant alleles."
      )
  } else {
    if(abs.count) {
      amplicon_plot <- ggplot(cons.table, aes_(
        x = ~Position,
        y = ~ (100 * `Max Non-ref Allele Count`),
        fill = ~Variants)
      ) +
        theme_bw() +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 6, angle = 90)) +
        ylab("Variant UMI count") +
        xlab("Assay") +
        facet_grid(`Sample Name` ~ Name, scales = "free_x", space = "free_x")
    } else {
      amplicon_plot <- ggplot(cons.table, aes_(
        x = ~Position,
        y = ~ (100 * `Max Non-ref Allele Frequency`),
        fill = ~Variants)
        ) +
        theme_bw() +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 6, angle = 90)) +
        ylab("Variant Allele Frequency (%)") +
        xlab("Assay") +
        facet_grid(`Sample Name` ~ Name, scales = "free_x", space = "free_x")
    }
  }

  # Show plot and add ggplot object to the UMIexperiment object
  if (do.plot) {
    print(amplicon_plot)
    object@plots$amplicon_plot <- amplicon_plot
  } else {
    object@plots$amplicon_plot <- amplicon_plot
  }
  return(object)
}

#' Generate Merged data plots
#' @export
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%" "%<>%"
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param cut.off How many variant reads are necessary to consider a variant above background? Default is 5 reads.
#' @param amplicons (Optional) character vector of amplicons to plot.
vizMergedData <- function(
  object,
  cut.off = 5,
  amplicons = NULL,
  do.plot = TRUE
  ){

  # Plotting maximum alternate alle count on merged data
  data <- object@merged.data
  data$Position %<>% as.factor

  if (!is.null(amplicons)) {
    data <- data %>% dplyr::filter(.data$Name %in% amplicons)
  }

  data$Variants <- ifelse(data$avg.MaxAC > cut.off, "Variant","Background")

  plot <- ggplot(data, aes_(x=~Position, y=~avg.MaxAC,fill=~Variants)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes_(ymin=~avg.MaxAC, ymax=~avg.MaxAC+std.MaxAC), width=.1) +
    theme(axis.text.x = element_text(size = 2, angle = 90)) +
    facet_grid(group.by ~ Name, scales = "free_x", space = "free_x")

  if(do.plot) {
    print(plot)
    object@plots$merged_amplicons <- plot
  } else {
    object@plots$merged_amplicons <- plot
  }
}

#' Generate consensus depth histograms
#' @export
#' @import ggplot2
#' @importFrom tibble is_tibble
#' @param object Requires a UMI sample or UMI experiment object
#' @param xMin Minimum consensus family size to plot, default is 0.
#' @param xMax Maximum consensus family size to plot. Default is 100.
#' @param samples List of samples to be shown.
plotFamilyHistogram <- function(
  object,
  xMin = 0,
  xMax = 100,
  samples = NULL
  ) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object.")
  } else if(!is.numeric(xMin) && !is.numeric(xMax)){
    warning("xMin and xMax needs to be numerical. Using default values instead.")
    xMin = 0
    xMax = 10
  } else if(xMin < 0 || xMax < 0){
    warning("xMin and xMax need to be positive numbers. Using default values")
    xMin = 0
    xMax = 10
  } else if(xMax < xMin){
    warning("xMax smaller than xMin, using default values instead.")
    xMin = 0
    xMax = 100
  } else if(!is.character(samples) && !is.null(samples)){
    warning("Samples need to be a character or NULL. Using default.")
    samples = NULL
  }

  # check if object is a UMIexperiment
  if (class(object)[1]== "UMIexperiment") {
    reads <- object@reads

    if (!is.null(samples)) {
      reads <- reads %>% dplyr::filter(.data$sample %in% samples)
    }

    cons_depth_plot <- ggplot(reads, aes(x = count, fill = sample, color = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      theme_classic() +
      xlim(xMin, xMax) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")
      ) +
      facet_wrap(~sample)
  } else {
    # a tibble containing read info is passed directly
    if (!is.null(samples)) {
      object <- object %>% dplyr::filter(.data$sample %in% samples)
    }

    cons_depth_plot <- ggplot(object, aes(x = count, fill = sample, color = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      theme_classic() +
      xlim(xMin, xMax) +
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
#' @return A list of ggplot objects.
vizNormalization <- function(cons.data){

  cons.data$Name %<>% as.factor
  cons.data$replicate <- cons.data$group.by

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
  merged <- list(p1,p2)

  return(merged)
}

#' View normalisation
#'
#' @export
#' @importFrom gridExtra grid.arrange
#'
#' @param object A umiExperiment object containing norm plots
#' @param do.plot should plot be shown? If false returns a grid.arrange object
#'
viewNormPlot <- function(object, do.plot = TRUE){

  plot <- grid.arrange(
    object@plots$norm_plot[[1]],
    object@plots$norm_plot[[2]],
    nrow = 1
  )

  if(do.plot){
    plot
  } else{
    return(plot)
  }
}

#' Plot counts by nucleotide change
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @param cons.data A consensus data table
#' @param do.plot Logical. Should plot be shown?
#' @return A ggplot object.
# Remove counts for reference allele for plotting
vizStackedCounts <- function(
  cons.data,
  do.plot = TRUE
  ){

  # For each row in consensus data, set the reference count to 0.
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

  # Rename variables
  out.file <- out.file %>% dplyr::select(
    .data$Name, .data$Position, .data$group.by,
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
                                        .data$group.by,
                                        .data$Reference,
                                        .data$avg.Depth))

  # Stacked count plot.
  stacked <- ggplot(out.file, aes_(fill=~variant, y=~count, x=~Position)) +
    geom_bar( stat="identity") +
    facet_grid(. ~ Name, scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90))

  if(do.plot){
    return(stacked)
  } else {
    return(NULL)
  }
}

