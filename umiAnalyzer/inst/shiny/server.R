# Define server logic
server <- function(input, output, session, plotFun) {


  #----- Download selected data csv ----

  output$downloadData <- downloadHandler(

    filename = function() {
      paste("consensus_data", ".csv", sep = "")
    },
    content = function(file) {

      filter <- umiAnalyzer::getFilteredData(
        object = filteredData()
      )

      filter <- filter %>%
        dplyr::filter(.data$Name %in% input$assays) %>%
        dplyr::filter(.data$`Sample Name` %in% input$samples) %>%
        dplyr::filter(.data$`Max Non-ref Allele Count` >= input$minCount) %>%
        dplyr::filter(.data$`Max Non-ref Allele Frequency` >= input$minFreq)

      readr::write_csv2(filter, file)
    }
  )

  #----Output_report-----

  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    file = 'report.html',
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), 'report.Rmd')
      file.copy('report.Rmd', tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list(
        data = filteredData(),
        assays = input$assays,
        samples = input$samples,
        minDepth =  input$consensus
      )

      # Knit the document, passing in the `params` list, and eval it in a
      # child of the global environment (this isolates the code in the document
      # from the code in this app).
      rmarkdown::render(
        tempReport,
        output_file = file,
        params = params,
        envir = new.env(parent = globalenv())
      )
    }
  )

  #----Download amplicon plot-----

  # Output pdf report upon button click
  output$download_plot <- downloadHandler(
    filename <- function() {
      paste('amplicon-plot-', Sys.time(),'.pdf',sep='') },
    content <- function(file) {
      pdf(file, width = 7, height = 3)
      object <- umiAnalyzer::generateAmpliconPlots(
        object = filteredData(),
        do.plot = TRUE,
        amplicons = input$assays,
        samples = input$samples,
        abs.count = input$abs_counts,
        theme = input$theme,
        option = input$colors,
        direction = input$direction
      )
      dev.off()
    }
  )

  #----Download QC plot-----

  output$download_qc_plot <- downloadHandler(
    filename <- function() {
      paste('qc-plot-', Sys.time(),'.pdf',sep='') },
    content <- function(file) {
      pdf(file, width = 7, height = 3)
      object <- umiAnalyzer::generateQCplots(
        object = experiment(),
        do.plot = TRUE,
        group.by = 'sample',
        plotDepth = input$consensus,
        assays = input$assays,
        samples = input$samples,
        theme = input$theme_qc,
        option = input$colors_qc,
        direction = input$direction_qc
      )
      dev.off()
    }
  )

  #----Download UMI plor-----

  # Output pdf report upon button click
  output$download_umi_plot <- downloadHandler(
    filename = 'umi_plot.pdf',
    content = function(file) {

      if(is.null(experiment())){
        return(NULL)
      }

      pdf(file, width = 7, height = 3)
      umiAnalyzer::plotUmiCounts(
        object = experiment(),
        do.plot = TRUE,
        amplicons = input$assays,
        samples = input$samples
      )
      dev.off()
    }
  )

  #----Shiny files setup-----

  # Define avalible volumes for shinyFiles
  volumes <- c(Home = fs::path_home(),
               'R Installation' = R.home(),
               getVolumes()())

  shinyDirChoose(
    input = input,
    id = 'dir',
    roots = volumes,
    session = session,
    restrictions = system.file(package = 'base')
  )

  shinyFileChoose(
    input = input,
    id = 'file',
    root = volumes,
    filetypes = c('.csv','.txt','.tsv')
  )

  shinyFileChoose(
    input = input,
    id = 'zipFile',
    root = volumes,
    filetypes = c('.zip')
  )

  #----Upload zipped data----

    temp_data_main <- reactive({

    zip_path <- input$zipFile$datapath

    if ( is.null(zip_path) ) {
      return(NULL)
    } else {

      temp_dir <- file.path(tempdir(),'appData')

      # Use 7zip to change windows paths to linux format
      #7z rn windows.zip $(7z l windows.zip | grep '\\' | awk '{ print $6, gensub(/\\/, "/", "g", $6); }' | paste -s)

      unzip(
        zipfile = zip_path,
        list = FALSE,
        exdir = temp_dir,
        unzip = 'internal'
      )

      return(temp_dir)
    }
  })

  #----remove this line----



  #----Meta data reactive----

  metaData <- reactive({

    if(is.null(experiment())){
      return(NULL)
    }

    #Note "file" is the name of the metadata from the inputUI

    metaData <- input$file$datapath

    if (identical(metaData, character(0))) {
      return(NULL)
    } else {

      data <- umiAnalyzer::importDesign(
        object = experiment(),
        file = metaData,
        delim = NULL # automatically select delimiter
      )

      design <- umiAnalyzer::getMetaData(
        object = data,
        attributeName = 'design'
      )

      choices <- colnames(design)

      # Updates values based on content from metadata file

      updateSelectInput(
        session = session,
        inputId = 'replicates',
        choices = choices,
        selected = head(choices,1)
      )

      updateSelectInput(
        session = session,
        inputId = 'timeVar',
        choices = choices,
        selected = head(choices,1)
      )

      return(design)
    }
  })

  # Set up user_data_main
  user_data_main <- reactive({
    # Path selected by the user
    main <- shinyFiles::parseDirPath(
      roots = volumes,
      selection = input$dir
    )

    # Create umiExperiment object
    if (identical(main, character(0))){
      return(NULL)
    } else {
      return(main)
    }
  })

  # Set up a test_data main
  test_data_main <- eventReactive(input$importTest,{
    main <- system.file('extdata', package = 'umiAnalyzer')
    return(main)
  })

  # Values is a reactive object to which a umiExperiment object is added in
  # the data slot.
  values <- reactiveValues(data=NULL, merge=FALSE)

  # Create experiment
  experiment <- reactive({
    # select directories
    if( !is.null(user_data_main()) || !is.null(temp_data_main()) ){
      if( !is.null(user_data_main())  ) {
        main <- user_data_main()
      } else {
        main <- temp_data_main()
      }
    } else {
      main <- test_data_main()
    }

    if (is.null(main)) {
      return(NULL)
    } else {

      # Check if assays have been merged. If false, initialise the umiExperiment
      # object and assing to the values object
      if( values$merge == FALSE){

        withProgress(message="Creating experiment object", value = 0, {
          values$data <- umiAnalyzer::createUmiExperiment(main, as.shiny = TRUE)
        })
      }

      #print( unique(values$data@cons.data$Name) )
      #data <- umiAnalyzer::createUmiExperiment(main)

      return(values$data)
    }
  })

  experiment_merged <- observeEvent(input$mergeAssays, {

    # Check of experiment exists and if new assay names have been defined
    if (is.null(experiment())){
      return(NULL)
    } else if(input$new_name == ""){
      return(NULL)
    }

    # Merge assays based on user input: (1) new assay name (2) list of assays to merge.
    new_data <- umiAnalyzer::mergeAssays(
      object = experiment(),
      name = input$new_name,
      assay.list = input$assay_list
    )

    # Update values object. This will trigger the reactive experiment() object
    # which will update data throughout the app with new assay info.
    values$data <- new_data
    values$merge <- TRUE

    #print( unique(values$data@cons.data$Name) )

    return(new_data)
  })

  mergedData <- observeEvent(input$mergeReplicates, {

    if (is.null(filteredData())){
      return(NULL)
    }

    if ( is.null(metaData()) ) {
      return(NULL)
    }

    if( input$replicates == "" ){
      replicates = NULL
    } else {
      replicates <- input$replicates
    }

    data <- filteredData()

    data@meta.data <- metaData()

    data <- umiAnalyzer::mergeTechnicalReplicates(
      object = data,
      do.plot = FALSE,
      group.by = input$replicates,
      amplicons = input$assays,
      samples = input$samples
    )

    out_data <- data@merged.data %>%
      dplyr::mutate_if(is.numeric, round, 1)

    output$mergedDataTable <- DT::renderDataTable({
      out_data
    })

    output$normPlot <- renderPlot({
      umiAnalyzer::viewNormPlot(data)
    })

    output$stackPlot <- renderPlot({
      data@plots$stacked_counts
    })

    output$mergePlot <- renderPlot({
      umiAnalyzer::vizMergedData(data)
    })

    return(data)
  })

  variantCalls <- observeEvent(input$runVarCaller, {

    if (is.null(filteredData())){
      return(NULL)
    }

    # Call and filter variants based on user input
    data <- filteredData()
    data <- umiAnalyzer::callVariants(data)

    data <- umiAnalyzer::filterVariants(
      object = data,
      p.adjust = input$pVal,
      minVarCount = input$minVarCount
    )

    out_data <- data@variants

    withProgress(message = 'Rendering outputs', value = 0.25, {

      # Render variant call table in app
      output$varDataTable <- DT::renderDataTable({
        out_data
      })

      # Render amplicon plot for computed variants
      output$varPlot <- renderPlot({
        umiAnalyzer::generateAmpliconPlots(
          object = data,
          amplicons = input$assays,
          samples = input$samples,
          abs.count = input$abs_counts
        )
      })

      shiny::incProgress(1, detail = paste("Rendering plot"))

    })

    return(data)
  })

  # filteredData returns an updated version of the experimen() object containing
  # a single filter called "user_filter" which is used downstream
  filteredData <- reactive({

    if (is.null(experiment())){
      return(NULL)
    }

    withProgress(message = 'Filtering object', value = 0.25, {

      data <- umiAnalyzer::filterUmiObject(
        object = experiment(),
        minDepth = input$consensus,
        minCoverage = 100,
        minFreq = input$minFreq/100,
        minCount = input$minCount
      )

      shiny::incProgress(1, detail = paste("Filtering"))

    })

    return(data)
  })

  # Update assay and sample list based on initially loaded object, meaning that
  # the lists will be visible even if filter are applied
  observe({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- umiAnalyzer::saveConsData( experiment() )

    updateSelectInput(
      session = session,
      inputId = 'assay_list',
      choices = unlist(strsplit(unique(data$Name), split = ',')),
      selected = head(unlist(strsplit(unique(data$Name), split = ',')),1)
    )

    updateSelectInput(
      session = session,
      inputId = 'assays',
      choices = unlist(strsplit(unique(data$Name), split = ',')),
      selected = head(unlist(strsplit(unique(data$Name), split = ',')),1)
    )

    updateSelectInput(
      session = session,
      inputId = 'samples',
      choices = unlist(strsplit(unique(data$`Sample Name`), split = ',')),
      selected = head(unlist(strsplit(unique(data$`Sample Name`), split = ',')),1)
    )

  })

  #--------- Output the consensus data --------
  output$dataTable <- DT::renderDataTable({

    if (is.null(filteredData())){
      return(NULL)
    }

    filter <- umiAnalyzer::getFilteredData(
      object = filteredData()
    )

    filter %>%
      dplyr::filter(.data$Name %in% input$assays) %>%
      dplyr::filter(.data$`Sample Name` %in% input$samples)

  }, options = list(
    orderClasses = TRUE,
    pageLength = 5,
    lengthMenu = c(5, 10, 50, 100)
  ))


  output$metaDataTable <- DT::renderDataTable({

    if (is.null(metaData())){
      return(NULL)
    }

    # TODO this table can now be edited by the user, but the changes are
    # not used by downstream functions.
    DT::datatable(metaData(), editable = TRUE)

  }, options = list(
    orderClasses = TRUE,
    pageLenght = 50,
    lengthMenu = c(10, 50, 100)
  ))

  # make reactive expresion of input values
  amplicon_settings <- reactive({input$assays})
  sample_settings <- reactive({input$samples})


  # delay amplicon plot until reactive stop changing
  amplicon_settings_d <- amplicon_settings %>% debounce(500)
  sample_settings_d <- sample_settings %>% debounce(500)

  # plot amplicon plot reactive value
  output$amplicon_plot <- renderPlot({

    if(is.null(filteredData())){
      return(NULL)
    }

    # TODO this generatesa new progress bar each time the plot changes. Consider
    # moving everything into an umbrella reactive object?

    withProgress(message = 'Rendering amplicon plot', value = 0.25, {

      object <- umiAnalyzer::generateAmpliconPlots(
        object = filteredData(),
        do.plot = TRUE,
        amplicons = amplicon_settings_d(),
        samples = sample_settings_d(),
        abs.count = input$abs_counts,
        theme = input$theme,
        option = input$colors,
        direction = input$direction
      )

      shiny::incProgress(1, detail = paste("Rendering complete"))

    })
    object@plots$amplicon_plot


  })


  #------ Output the QC plot -------

  output$qcPlot <- renderPlot({

    if(is.null(experiment())){
      return(NULL)
    }

    shiny::withProgress(message = 'Rendering QC plot', value = 0.25, {
      umiAnalyzer::generateQCplots(
        object = experiment(),
        do.plot = TRUE,
        group.by = 'sample',
        plotDepth = input$consensus,
        assays = input$assays,
        samples = input$samples,
        theme = input$theme_qc,
        option = input$colors_qc,
        direction = input$direction_qc
      )
      shiny::incProgress(1, detail = paste("Rendering QC plot"))
    })
  })

  observeEvent(input$timeSeries, {

    output$timeSeriesPlot <- renderPlot({

      if(is.null(filteredData())){
        return(NULL)
      }

      if(is.null(metaData())){
        return(NULL)
      }

      data <- filteredData()
      data@meta.data <- metaData()

      umiAnalyzer::analyzeTimeSeries(
        object = data,
        time.var = input$timeVar,
        group.by = input$replicates
      )
    })
  })

  output$umiCounts <- renderPlot({

    if(is.null(experiment())){
      return(NULL)
    }

    if(input$direction_umi == "default") {
      direction = 1
    } else {
      direction = -1
    }

    # Initialise progress bar
    shiny::withProgress(
      message = 'Rendering UMI plot',
      value = 0.25, {

        umiAnalyzer::plotUmiCounts(
          object = experiment(),
          do.plot = TRUE,
          amplicons = input$assays,
          samples = input$samples,
          theme = input$theme_umi,
          option = input$colors_umi,
          direction = direction
        )

        # Update progress bar
        shiny::incProgress(
          amount = 1,
          detail = paste("Rendering UMIs")
        )
      })

  })

  #----- Import BAM files ------

  # Import consensus read bam file upon button click to generate histograms
  # of barcode distribution. It is possible to import directly into the umiExperiment object
  # by setting the importBam flag, but file and this might become very large, so this
  # function is outsourced to parseBamFiles
  observeEvent(input$importBam, {

    # select between main
    if(!is.null(user_data_main())){
      main <- user_data_main()
    } else {
      main <- test_data_main()
    }

    if (identical(main, character(0))) {
      return(NULL)
    } else {
      # List sample names in main directory
      samples <- list.dirs(
        path = main,
        full.names = FALSE,
        recursive = FALSE
      )

      shiny::withProgress(
        message = 'Parsing consensus reads',
        value = 0, {

          # Parse each sample folder for .bam files containing consensus reads
          reads <- umiAnalyzer::parseBamFiles(
            mainDir = main,
            sampleNames = samples,
            consDepth = 0,
            as.shiny= TRUE
          )

        })

      # Output barcode family histogram
      output$histogram <- renderPlot({

        # TODO progress bar initialises at 0.25 and finishes at 1 when plot is
        # rendered. Implement continuous bar?
        shiny::withProgress(
          message = 'Rendering histograms',
          value = 0.25, {

            # Generate histogram plot using user defined parameters
            umiAnalyzer::plotFamilyHistogram(
              object = reads,
              xMin = input$famSize[1],
              xMax = input$famSize[2],
              samples = input$samples
            )

            # Update progress bar
            shiny::incProgress(
              amount = 1,
              detail = paste("Rendering histograms")
            )
          })

      })
    }
  })
}

# Run the application
