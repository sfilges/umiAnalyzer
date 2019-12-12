#----------------------------// UmiVisualiser //--------------------------------
#
# A Shiny app for visualising data generated with the umierrorcorrect pipeline:
# https://github.com/tobbeost/umierrorcorrect
#
# This app uses and comes supplied with the umiAnalyzer package:
# https://github.com/ozimand1as/umiAnalyzer
#

# Quietly import required packages
library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(shinyWidgets, quietly = TRUE)
library(DT, quietly = TRUE)
library(shinydashboard, quietly = TRUE)
library(umiAnalyzer, quietly = TRUE)

# Maximum 5GB data upload
options(shiny.maxRequestSize=5000*1024^2)

# Define user interface
ui <- dashboardPage(
  dashboardHeader(
    title = 'umiVisualiser'
  ),
  # Define menu items on the sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = 'Dashboard',
        tabName = 'dashboard',
        icon = icon('dna')
      ),
      menuItem(
        text = 'Advanced',
        tabName = 'advanced',
        icon = icon('magic')
      ),
      menuItem(
        text = 'User Guide',
        tabName = 'vignette',
        icon = icon('book')
      ),
      menuItem(
        text = 'umierrorcorrect',
        icon = icon('git'),
        href = 'https://github.com/tobbeost/umierrorcorrect'
      )
    )
  ),

  # Define main windows
  dashboardBody(
    # List tab items ...
    tabItems(
      # ... each tab-item correponds to a menu-item in the sidebar
      tabItem(tabName = 'dashboard',
        fluidRow(
          # Box for primary data upload and selection
          box(
            title = "Input",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 420,
            # Tab box with two panels
            tabBox(
              width = 12,
              type = "tabs",
              # Panel 1: Data upload widgets
              tabPanel(
                title = "Upload data",
                icon = icon("upload"),
                # Encase i/o buttons in a fluid row environment
                fluidRow(
                  style = "margin-bottom: 10px;margin-left: 5px;margin-right: 5px;",
                  fileInput(
                    inputId = 'zipFile', width = "100%",
                    label = 'Choose a zip file (Max. 5 GB)',
                    multiple = FALSE,
                    accept = c('.zip')
                  ),
                  fileInput(
                    inputId = 'file', width = "100%",
                    label = 'Choose a file containing sample metadata',
                    multiple = FALSE,
                    accept = c('.txt','.csv','.tsv')
                  ),
                  dropdown(
                    label = "Load data",

                    shinyDirButton(
                      id = 'dir',
                      label = 'Choose directory',
                      title = 'Upload',
                      style = "margin-bottom: 10px;margin-left: 5px;",
                      icon = icon("folder-open")
                    ),
                    actionButton(
                      inputId = "importBam",
                      label = "Import .bam files",
                      style = "margin-bottom: 10px;margin-left: 5px;",
                      icon = icon("align-left")
                    ),
                    actionButton(
                      inputId = 'importTest',
                      label = 'Load test data',
                      icon = icon('file-import'),
                      style = "margin-bottom: 10px;margin-left: 5px;"
                    ),
                    downloadButton(
                      outputId = 'report',
                      label = 'Print Report',
                      style = "margin-bottom: 10px;margin-left: 5px;"
                    ),

                    icon = icon("gear"),
                    width = "100%"
                  )
                )
              ),
              # Panel 2: Data selection - the user chooses consensus depth for
              # filtering and which samples and assays to show.
              tabPanel(
                title = 'Data selection',
                icon = icon('edit'),

                selectInput(
                  inputId = 'consensus', width = "50%",
                  label = 'Consensus Depth:',
                  choices = c(1,2,3,4,5,10,20,30),
                  selected = 3
                ),
                selectInput(
                  inputId = 'samples', width = "100%",
                  label = 'Samples:',
                  choices = '',
                  multiple = TRUE
                ),
                selectInput(
                  inputId = 'assays', width = "100%",
                  label = 'Assays:',
                  choices = '',
                  multiple = TRUE
                )
              ),
              # Panel 3: Merging assays
              tabPanel(
                title = 'Merge assays',
                icon = icon('hubspot'),
                textInput(
                  inputId = 'new_name',
                  label = 'Merged assay name:'
                ),
                selectInput(
                  inputId = 'assay_list',
                  label = 'Assays to merge:',
                  choices = '',
                  multiple = TRUE
                ),
                actionButton(
                  inputId = 'mergeAssays',
                  label = 'Merge!',
                  icon = icon('cog')
                )
              )
            )
          ),
          # Box for interactive paramter selection
          box(
            title = "Parameters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 420,
            sliderInput(
              inputId = "minFreq",
              label = "Minimum Variant allele frequency:",
              min = 0, max = 1,
              value = 0, step = 0.01,
              post = "%", sep = ","
            ),
            sliderInput(
              inputId = "minCount",
              label = "Minimum Variant allele count:",
              min = 0, max = 10,
              value = 0, step = 1,
              post = " reads", sep = ","
            ),
            sliderInput(
              inputId = "famSize",
              label =  "Minimum and Maximum family size to show:",
              min = 0, max = 500,
              value = c(0,100), step = 1,
              post = " reads", sep = ","
            ),
            br(),
            materialSwitch(
              inputId = "abs_counts",
              label = "Absolute counts: ",
              status = "primary"
            )
          ),
          # View data tables in collapsable box
          box(
            title = "Data Viewer",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(width = 12,
              tabBox(width = 12,
                tabPanel(
                  title = "Data",
                  DT::dataTableOutput("dataTable"),
                  style = "font-size: 10px;"
                ),
                tabPanel(
                  title = "Sample info",
                  DT::dataTableOutput("metaDataTable")
                )
              )
            )
          ),
          # Show plots in collapsable box containing a tabBox with a tab for
          # each plot to be shown.
          box(
            title = 'Plots',
            status = 'primary',
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(
              width = 12,
              tabBox(
                type = 'tabs',
                width = 12,
                # Panel for the amplicon plots with download button
                tabPanel(
                  title = 'Amplicons',
                  icon = icon('chart-bar'),
                  style = 'margin-left: 10px;',
                  fluidRow(
                    # Option for plot customisation
                    dropdownButton(
                      tags$h3('Customise plot'),
                      selectInput(
                        inputId = 'colors',
                        label = 'Choose colour palette:',
                        choices = c('default','viridis','magma','plasma','inferno','cividis')
                      ),
                      selectInput(
                        inputId = 'direction',
                        label = 'Color palette direction:',
                        choices = c('default','reverse')
                      ),
                      selectInput(
                        inputId = 'theme',
                        label = 'Choose theme:',
                        choices = c('classic','gray','bw','minimal','light')
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '300px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotOutput(
                      outputId = 'amplicon_plot',
                      width = '100%'
                    ),
                    downloadButton(
                      outputId = 'download_plot',
                      label = 'Download figure'
                    )
                  )
                ),
                # Panel for quality control plot
                tabPanel(
                  title = 'QC Plot',
                  style = 'margin-left: 10px;',
                  fluidRow(
                    # Option for plot customisation
                    dropdownButton(
                      tags$h3('Customise plot'),
                      selectInput(
                        inputId = 'colors_qc',
                        label = 'Choose colour palette:',
                        choices = c('default','viridis','magma','plasma','inferno','cividis')
                      ),
                      selectInput(
                        inputId = 'direction_qc',
                        label = 'Color palette direction:',
                        choices = c('default','reverse')
                      ),
                      selectInput(
                        inputId = 'theme_qc',
                        label = 'Choose theme:',
                        choices = c('classic','gray','bw','minimal','light')
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '300px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotOutput('qcPlot'),
                    downloadButton(
                      outputId = 'download_qc_plot',
                      label = 'Download figure'
                    )
                  )
                ),
                # Panel for UMI count data
                tabPanel(
                  title = "UMI counts",
                  style = 'margin-left: 10px;',
                  fluidRow(
                    # Option for plot customisation
                    dropdownButton(
                      tags$h3('Customise plot'),
                      selectInput(
                        inputId = 'colors_umi',
                        label = 'Choose colour palette:',
                        choices = c('default','viridis','magma','plasma','inferno','cividis')
                      ),
                      selectInput(
                        inputId = 'direction_umi',
                        label = 'Color palette direction:',
                        choices = c('default','reverse')
                      ),
                      selectInput(
                        inputId = 'theme_umi',
                        label = 'Choose theme:',
                        choices = c('classic','gray','bw','minimal','light')
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '300px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotOutput("umiCounts"),
                    downloadButton(
                      outputId = 'download_umi_plot',
                      label = 'Download figure'
                    )
                  )
                ),
                # Panel for barcode family histogram
                tabPanel(
                  title = "Histogram",
                  plotOutput("histogram")
                )
              )
            )
          )
        )
      ),

      shinydashboard::tabItem(tabName = "advanced",
        shiny::fluidRow(
          shinydashboard::box(
            title = "Advanced data analysis",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            width = 6,
            shiny::actionButton(
              inputId = 'runVarCaller',
              label = 'Run variant caller'
            ),
            shiny::actionButton(
              inputId = 'mergeReplicates',
              label = "Merge Replicates"
            ),
            shiny::actionButton(
              inputId = 'timeSeries',
              label = "Analyse time series"
            ),
            shiny::selectInput(
              inputId = 'replicates',
              label = 'Replicates:',
              choices = '',
              multiple = FALSE
            ),
            shiny::selectInput(
              inputId = 'timeVar',
              label = 'Time variable:',
              choices = '',
              multiple = FALSE
            )
          ),

          shinydashboard::box(
            title = "Parameters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            width = 6,
            shiny::sliderInput(
              inputId = "minVarCount",
              label = "Minimum Variant allele count:",
              min = 0, max = 10,
              value = 0, step = 1,
              post = " reads", sep = ","
            ),
            shiny::sliderInput(
              inputId = "pVal",
              label =  "Minimum adjusted p-value:",
              min = 0, max = 1,
              value = 1, step = 0.05,
              sep = ","
            )
          ),

          shinydashboard::box(
            title = 'Plot Viewer',
            status = 'primary',
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(
              width = 12,
              tabBox(width = 8,
                tabPanel(
                  title = 'Normalization',
                  plotOutput('normPlot')
                ),
                tabPanel(
                  title = 'Stacked counts',
                  plotOutput('stackPlot')
                ),
                tabPanel(
                  title = 'Merged amplicons',
                  plotOutput('mergePlot')
                ),
                tabPanel(
                  title = 'Time series',
                  plotOutput('timeSeriesPlot')
                ),
                tabPanel(
                  title = 'Variant caller',
                  plotOutput('varPlot')
                )
              )
            )
          ),

          # View data tables in collapsable box
          shinydashboard::box(
            title = 'Data Viewer',
            status = 'primary',
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(width = 12,
              tabBox(
                width = 12,
                tabPanel(
                  title = 'Amplicons data'
                ),
                tabPanel(
                  title = 'Sample info'
                ),
                tabPanel(
                  title = 'Merged data',
                  DT::dataTableOutput('mergedDataTable'),
                  style = "font-size: 10px;"
                ),
                tabPanel(
                  title = 'Variant data',
                  DT::dataTableOutput('varDataTable'),
                  style = "font-size: 10px;"
                )
              )
            )
          )
        )
      ),

      tabItem(
        tabName = 'vignette',
        includeMarkdown("user_guide.Rmd")
        #includeHTML("user_guide.html")
      )
    )
  )
)

# Define server logic
server <- function(input, output, session, plotFun) {

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

  # Output pdf report upon button click
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

  # Upload zipped data
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

  metaData <- reactive({

    if(is.null(experiment())){
      return(NULL)
    }

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

  # Output the consensus data to screen, this will change depending on user input
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

  # Generate amplicon plots using umiAnalyzer
  output$amplicon_plot <- renderPlot({

    if(is.null(filteredData())){
      return(NULL)
    }

    shiny::withProgress(
      message = 'Rendering amplicon plot',
      value = 0.25, {

    umiAnalyzer::generateAmpliconPlots(
      object = filteredData(),
      do.plot = TRUE,
      amplicons = input$assays,
      samples = input$samples,
      abs.count = input$abs_counts,
      theme = input$theme,
      option = input$colors,
      direction = input$direction
    )

    shiny::incProgress(1, detail = paste("Rendering plot"))

    })

  })

  # Output the QC plot
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

    umiAnalyzer::plotUmiCounts(
      object = experiment(),
      do.plot = TRUE,
      amplicons = input$assays,
      samples = input$samples,
      theme = input$theme_umi,
      option = input$colors_umi,
      direction = direction
    )
  })

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
shiny::shinyApp(
  ui = ui,
  server = server
)
