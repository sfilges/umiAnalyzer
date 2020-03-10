#' Generate QC plots
#'
#' Visualise the UMI count for each selected assay and sample for a given
#' consensus depth. This is useful to detect differences in coverage,
#' especially for multiplexed assays.
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param group.by String. Which variable should be used as a factor on the x-axis. Default is assay.
#' @param plotDepth Which consensus depth to plot
#' @param assays (Optional) user-supplied list of assays to plot. Default is all.
#' @param samples (Optional) user-supplied list of samples to plot. Default is all.
#' @param theme ggplot theme to use.
#' @param option Colour palette to use, etiher ggplot default or viridis colours.
#' @param direction If viridis colours are used, choose orientation of colour scale.
#'
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats median
#' @importFrom viridis scale_fill_viridis
#'
generateQCplots <- function(
  object,
  do.plot = TRUE,
  group.by = c('assay', 'sample'),
  plotDepth = 3,
  assays = NULL,
  samples = NULL,
  theme = 'classic',
  option = 'default',
  direction = 'default'
  ) {

  if (missing(x = object)) {
    stop("Must provide a umiExperiment object and filter names")
  } else if(!class(object) == "UMIexperiment"){
    stop("Object is not of class UMIexperiment.")
  }

  cons.table <- object@cons.data
  summary.table <- object@summary.data

  # Consensus depth plot per assay
  cdepths <- summary.table %>% dplyr::filter(
    .data$assay != '',
    .data$depth == plotDepth
  )

  if (!is.null(assays)) {
    cdepths <- cdepths %>%
      dplyr::filter(.data$assay %in% assays)
  }

  if (!is.null(samples)) {
    cdepths <- cdepths %>%
      dplyr::filter(.data$sample %in% samples)
  }

  # Set assay and sample to factor for better plotting
  cdepths$assay %<>% as.factor
  cdepths$sample %<>% as.factor

  # From the ggplot2 vignette:
  # https://github.com/tidyverse/ggplot2/releases
  # aes_() replaces aes_q(). It also supports formulas, so the most concise
  # SE version of aes(carat, price) is now aes_(~carat, ~price). You may
  # want to use this form in packages, as it will avoid spurious R CMD check
  # warnings about undefined global variables.

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  if (group.by == "assay") {
    depth_plot <- ggplot(cdepths, aes_(x = ~assay, y = ~UMIcount, fill=~sample)) +
      geom_bar(position = "dodge", stat = "identity") +
      use_theme +
      theme(
        axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 14)
      ) +
      geom_hline(
        yintercept = median(cdepths$UMIcount),
        linetype = "dashed", color = "red"
      ) +
      geom_hline(
        yintercept = mean(cdepths$UMIcount),
        linetype = "dashed", color = "blue"
      ) +
      labs(
        title = paste("Consensus ", plotDepth, " depths by assay", sep = ""),
        subtitle = paste(
          "Mean depth: ", round(mean(cdepths$UMIcount)),
          "Median depth: ", round(median(cdepths$UMIcount))
        ),
        caption = ""
      )
  } else if (group.by == "sample") {
    depth_plot <- ggplot(cdepths, aes_(x = ~sample, y = ~UMIcount, fill=~assay)) +
      geom_bar(position = "dodge", stat = "identity") +
      use_theme +
      theme(
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)
      ) +
      geom_hline(
        yintercept = median(cdepths$UMIcount),
        linetype = "dashed",
        color = "red") +
      geom_hline(
        yintercept = mean(cdepths$UMIcount),
        linetype = "dashed",
        color = "blue") +
      labs(
        title = paste("Consensus ", plotDepth, " depths by sample", sep = ""),
        subtitle = paste(
          "Mean depth: ", round(mean(cdepths$UMIcount)),
          "Median depth: ", round(median(cdepths$UMIcount))
        ),
        caption = ""
      )
  }

  if(option != 'default'){

    if(direction == 'default'){
      orientation = 1
    } else {
      orientation = -1
    }

    depth_plot <- depth_plot + viridis::scale_fill_viridis(
      discrete = TRUE,
      option = option,
      direction = orientation
    )
  }

  # Plot consensus depth distribution
  if (do.plot) {
    print(depth_plot)
    object@plots$qc_depth_plot <- depth_plot
  }
  else {
    object@plots$qc_depth_plot <- depth_plot
  }

  return(object)
}

#' Plot UMI counts
#'
#' Visualise the number detected UMI for each consensus depth cut-off. This may
#' may helpful in choosing the right consensus depth for your analysis, by
#' checking the number of reads still available for each assay and sample
#' for your chosen cut-off.
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param amplicons (Optional) user-supplied list of assays to plot. Default is all.
#' @param samples (Optional) user-supplied list of samples to plot. Default is all.
#' @param theme Plotting theme, default is classic
#' @param option Colour palette. Default uses ggplot standard, otherwise viridis options.
#' @param direction If using viridis colours should the scale be inverted or default?
#'
#' @export
#'
#' @import ggplot2
#' @import dplyr
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom stats median
#' @importFrom viridis scale_fill_viridis
#'
plotUmiCounts <- function(
  object,
  do.plot = TRUE,
  amplicons = NULL,
  samples = NULL,
  theme = 'classic',
  option = 'viridis',
  direction = 1
  ) {

  # Read summary data from object
  data <- object@summary.data %>%
    dplyr::filter(!is.na(.data$assay),
                  .data$depth > 0)

  # Select amplicons
  if (!is.null(amplicons)) {
    data <- data %>%
      dplyr::filter(.data$assay %in% amplicons)
  }

  # Select samples
  if (!is.null(samples)) {
    data <- data %>%
      dplyr::filter(.data$sample %in% samples)
  }

  # Generate ggplot object
  data$depth %<>% as.factor

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  # If colour option is not default use the viridis package for colour palettes.
  # https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html
  if(option != 'default') {
    plot <- ggplot(data, aes(x=.data$depth, y=.data$UMIcount, fill=sample)) +
      use_theme +
      viridis::scale_fill_viridis(
        discrete = TRUE,option = option,direction = direction) +
        geom_col(alpha=0.6) +
        facet_grid(assay ~ sample) +
        ylab("UMI count") +
        xlab("Consensus depth cut-off")
  } else {
    plot <- ggplot(data, aes(x=.data$depth, y=.data$UMIcount, fill=sample)) +
      use_theme +
      geom_col(alpha=0.6) +
      facet_grid(assay ~ sample) +
      ylab("UMI count") +
      xlab("Consensus depth cut-off")
  }

  # Return object and plot
  if(do.plot){
    print(plot)
    object@plots$umi_counts <- plot
    return(object)
  } else {
    object@plots$umi_counts <- plot
    return(object)
  }
}

#' Generate Amplicon plots
#'
#' Plots variant allele frequencies or alternate allele counts for chosen
#' samples and assays.
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param filter.name Name of the filter to be plotted.
#' @param do.plot Logical. Should plots be shown?
#' @param cut.off How many variant reads are necessary to consider a variant above background? Default is 5 reads.
#' @param amplicons (Optional) character vector of amplicons to be plotted.
#' @param samples (Optional) character vector of samples to be plotted.
#' @param abs.count Should absolute counts be plotted instead of frequencies? Default is FALSE.
#' @param theme Plotting theme to use, default is classic.
#' @param option Colour palette to use.
#' @param direction Orientation of the colour palette.
#' @param plot.text Should non-references bases be indicated above the bar?
#' @param plot.ref If true show reference base instead of position on x-axis.
#' @param stack.plot Show all variant alleles in a stacked bar plot.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom dplyr filter
#' @importFrom viridis scale_fill_viridis
#' @importFrom scales rescale_none
#'
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
  filter.name = 'default',
  do.plot = TRUE,
  cut.off = 5,
  amplicons = NULL,
  samples = NULL,
  abs.count = FALSE,
  y_min = 0,
  y_max = NULL,
  theme = 'classic',
  option = 'default',
  direction = 'default',
  plot.text = TRUE,
  plot.ref = TRUE,
  stack.plot = FALSE
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

  cons.table <- filterConsensusTable(
    cons.table,
    amplicons = amplicons,
    samples = samples
  )

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  # Set maximum y-axis limit to largest variant allele frequency in data
  # rounded up to the nearest integer.
  if(is.null(y_max)){
    y_max <- ceiling(100*max(cons.table$`Max Non-ref Allele Frequency`))
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
      use_theme +
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
        y = ~(`Max Non-ref Allele Count`),
        fill = ~Variants)
      ) +
        use_theme +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 6, angle = 90)) +
        ylab("Variant UMI count") +
        xlab("Assay") +
        facet_grid(`Sample Name` ~ Name, scales = "free_x", space = "free_x")
    } else {
      amplicon_plot <- ggplot(cons.table, aes_(
        x = ~Position,
        y = ~ (100 * `Max Non-ref Allele Frequency`),
        fill = ~Variants)) +
        use_theme +
        geom_bar(stat = "identity") +
        theme(axis.text.x = element_text(size = 6, angle = 90)) +
        ylab("Variant Allele Frequency (%)") +
        xlab("Assay") +
        scale_y_continuous(limits=c(y_min,y_max), oob = scales::rescale_none) +
        facet_grid(`Sample Name` ~ Name, scales = "free_x", space = "free_x")
    }

    if(plot.text){
      amplicon_plot <- amplicon_plot +
        geom_text(data = cons.table,
                  mapping = aes(label = `Max Non-ref Allele`),
                  position = position_dodge(width = 1),
                  size = 4
        )
    }

    if(plot.ref){
      amplicon_plot <- amplicon_plot +
        scale_x_discrete(
          breaks = cons.table$Position,
          labels = cons.table$Reference
        ) +
        theme(
          axis.text.x = element_text(size = 9, angle = 0)
        )
    }

  }

  if(direction == 'default'){
    orientation = 1
  } else {
    orientation = -1
  }

  if(option != 'default'){

    # If not using default colour scheme use either
    # (1) Colours from viridis package
    if( option %in% c('viridis','magma','plasma','inferno','cividis') ){
      amplicon_plot <- amplicon_plot + viridis::scale_fill_viridis(
        discrete = TRUE,
        option = option,
        direction = orientation
      )
    } else {
    # (2) Colours from ggplot: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
      amplicon_plot <- amplicon_plot +
        ggplot2::scale_fill_brewer(
          palette = option
        )
    }
  }

  if(stack.plot){
    amplicon_plot <- stacked_amplicon_plot(
      cons.data = cons.table,
      theme = theme,
      plot.ref = plot.ref,
      abs.count = abs.count,
      option = option,
      direction = orientation
    )
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
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param do.plot Logical. Should plots be shown.
#' @param cut.off How many variant reads are necessary to consider a variant above background? Default is 5 reads.
#' @param amplicons (Optional) character vector of amplicons to plot.
#' @param theme ggplot theme to use. Default is classic.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom dplyr filter
#' @importFrom magrittr "%>%" "%<>%"
#'
vizMergedData <- function(
  object,
  cut.off = 5,
  amplicons = NULL,
  do.plot = TRUE,
  theme = 'classic'
  ){

  # Plotting maximum alternate alle count on merged data
  data <- object@merged.data
  data$Position %<>% as.factor

  if (!is.null(amplicons)) {
    data <- data %>% dplyr::filter(.data$Name %in% amplicons)
  }

  data$Variants <- ifelse(data$avg.MaxAC > cut.off, 'Variant', 'Background')

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  plot <- ggplot(data, aes_(x=~Position, y=~avg.MaxAC,fill=~Variants)) +
    geom_bar(stat = 'identity') +
    use_theme +
    geom_errorbar(aes_(ymin=~avg.MaxAC, ymax=~avg.MaxAC+std.MaxAC), width=.1) +
    theme(axis.text.x = element_text(size = 2, angle = 90)) +
    facet_grid(group.by ~ Name, scales = 'free_x', space = 'free_x')

  if(do.plot) {
    print(plot)
    object@plots$merged_amplicons <- plot
  } else {
    object@plots$merged_amplicons <- plot
  }
}

#' Consensus depth histograms
#'
#' Generate histograms for the frequency of barcode family depths.
#'
#' @param object Requires a UMI sample or UMI experiment object
#' @param xMin Minimum consensus family size to plot, default is 0.
#' @param xMax Maximum consensus family size to plot. Default is 100.
#' @param samples List of samples to be shown.
#' @param option Colour scheme to use
#' @param direction If using viridis colours sets the orientation of colour scale.
#' @param theme ggplot theme to use. Defaults to classic.
#'
#' @export
#'
#' @import ggplot2
#' @importFrom tibble is_tibble
#' @importFrom viridis scale_fill_viridis
#'
plotFamilyHistogram <- function(
  object,
  xMin = 0,
  xMax = 100,
  samples = NULL,
  option = 'viridis',
  direction = 1,
  theme = 'classic'
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

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  # check if object is a UMIexperiment
  if (class(object)[1]== "UMIexperiment") {
    reads <- object@reads

    if (!is.null(samples)) {
      reads <- reads %>% dplyr::filter(.data$sample %in% samples)
    }

    cons_depth_plot <- ggplot(reads, aes(x = count, fill = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      use_theme +
      viridis::scale_fill_viridis(discrete = TRUE,option = option,direction = direction) +
      xlim(xMin, xMax) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")
      ) +
      facet_wrap(~sample) +
      ylab("Number of families") +
      xlab("Barcode familiy size")

    # Assign plot to UMIexperiment object
    object@plots$family_histogram <- cons_depth_plot

  } else {
    # a tibble containing read info is passed directly
    if (!is.null(samples)) {
      object <- object %>% dplyr::filter(.data$sample %in% samples)
    }

    cons_depth_plot <- ggplot(object, aes(x = count, fill = sample)) +
      geom_histogram(binwidth = 1, alpha = 0.5) +
      use_theme +
      viridis::scale_fill_viridis(discrete = TRUE,option = option,direction = direction) +
      xlim(xMin, xMax) +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 14, face = "bold")
      ) +
      facet_wrap(~sample) +
      ylab("Number of families") +
      xlab("Barcode familiy size")
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

#' Plot all variant allele bases
#'
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#'
#' @param cons.data Consensus data table
#' @param theme Plotting theme
#' @param plot.ref If true, shows reference base on x-axis
#' @param abs.count Plot absolute countsinstead of frequencies.
#'
#' @return A ggplot object.
stacked_amplicon_plot <- function(
  cons.data,
  theme = 'classic',
  plot.ref = FALSE,
  abs.count = FALSE,
  option = 'Pastel1',
  direction = 1
  ){

  out.file <- tibble()
  for(j in 1:nrow(cons.data)) {
    row <- cons.data[j,]
    if( row$Reference == "A" ) {
      row$A <- 0
    }
    else if( row$Reference == "C" ) {
      row$C <- 0
    }
    else if( row$Reference == "G" ) {
      row$G <- 0
    }
    else if( row$Reference == "T" ) {
      row$T <- 0
    }
    out.file <- dplyr::bind_rows(out.file, row)
  }

  # Rename variables
  out.file <- out.file %>% dplyr::select(
    .data$Name, .data$Position,
    .data$Coverage,.data$Reference,
    .data$A,.data$T,.data$C,.data$G,
    .data$N,.data$I,.data$D) %>%
    tidyr::gather(
      "variant",
      "count",
      -c(.data$Name,
         .data$Position,
         .data$Reference,
         .data$Coverage
        )
      )

  out.file$Name %<>% as.factor
  out.file$Position %<>% as.factor

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  if(abs.count){
    # Stacked count plot.
    stacked <- ggplot(out.file, aes_(
      fill=~variant,
      y=~count,
      x=~Position)) +
      geom_bar(position="stack", stat="identity") +
      use_theme +
      ylab("Variant Allele Frequency (%)") +
      xlab("Assay") +
      facet_grid(. ~ Name, scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_text(angle = 90))
  } else {
    # Stacked count plot.
    stacked <- ggplot(out.file, aes_(
      fill=~variant,
      y=~(100*count/Coverage),
      x=~Position)) +
      geom_bar(position="stack", stat="identity") +
      use_theme +
      ylab("Variant Allele Frequency (%)") +
      xlab("Assay") +
      facet_grid(. ~ Name, scales = "free_x", space = "free_x") +
      theme(axis.text.x = element_text(angle = 90))
  }

  if(plot.ref){
    stacked <- stacked +
      scale_x_discrete(
        breaks = out.file$Position,
        labels = out.file$Reference
      ) +
      theme(
        axis.text.x = element_text(size = 9, angle = 0)
      )
  }

  if(option == 'default'){
    option = 'Pastel1'
  }

  if(option != 'default'){
    # If not using default colour scheme use either
    # (1) Colours from viridis package
    if( option %in% c('viridis','magma','plasma','inferno','cividis') ){
      stacked <- stacked + viridis::scale_fill_viridis(
        discrete = TRUE,
        option = option,
        direction = direction
      )
    } else {
      # (2) Colours from ggplot: Accent, Dark2, Paired, Pastel1, Pastel2, Set1, Set2, Set3
      stacked <- stacked +
        ggplot2::scale_fill_brewer(
          palette = option
        )
    }
  }

  return(stacked)

}

#' View count normalisation
#'
#' @param object A umiExperiment object containing norm plots
#' @param do.plot should plot be shown? If false returns a grid.arrange object
#'
#' @export
#' @importFrom gridExtra grid.arrange
#'
viewNormPlot <- function(
  object,
  do.plot = TRUE
  ){

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
#' @param cons.data A consensus data table
#' @param do.plot Logical. Should plot be shown?
#' @param option Colour palette to use
#' @param direction If using viridis colors, choose orientation of palette
#' @param theme ggplot theme to use, default is classic.
#'
#' @import tibble
#' @import dplyr
#' @import ggplot2
#' @importFrom magrittr "%>%" "%<>%"
#' @importFrom tidyr gather
#' @importFrom rlang .data
#' @importFrom viridis scale_fill_viridis
#'
#' @return A ggplot object.
#'
vizStackedCounts <- function(
  cons.data,
  do.plot = TRUE,
  option = 'viridis',
  direction = 1,
  theme = 'classic'
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

  # Use selected plotting theme
  use_theme <- select_theme(theme = theme)

  # Stacked count plot.
  stacked <- ggplot(out.file, aes_(fill=~variant, y=~count, x=~Position)) +
    geom_bar( stat="identity") +
    use_theme +
    viridis::scale_fill_viridis(discrete = TRUE,option = option,direction = direction) +
    facet_grid(. ~ Name, scales = "free_x", space = "free_x") +
    theme(axis.text.x = element_text(angle = 90))

  if(do.plot){
    return(stacked)
  } else {
    return(NULL)
  }
}

