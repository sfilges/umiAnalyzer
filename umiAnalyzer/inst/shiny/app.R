#----------------------------// UmiVisualiser //--------------------------------
#
# A Shiny app for visualising data generated with the umierrorcorrect pipeline:
# https://github.com/tobbeost/umierrorcorrect
#
# This app uses and comes supplied with the umiAnalyzer package:
# https://github.com/ozimand1as/umiAnalyzer
#
#
# Quietly import required packages
library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(shinyWidgets, quietly = TRUE)
library(DT, quietly = TRUE)
library(shinydashboard, quietly = TRUE)
library(umiAnalyzer, quietly = TRUE)

# Maximum 1GB data upload
options(shiny.maxRequestSize=1000*1024^2)

# Define user interface
ui <- dashboardPage(

  dashboardHeader(
    title = "umiVisualiser"
  ),

  # Define menu items on the sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = "Dashboard",
        tabName = "dashboard",
        icon = icon("dna")
      ),
      menuItem(
        text = "Advanced",
        tabName = "advanced",
        icon = icon("magic")
      ),
      menuItem(
        text = "User Guide",
        tabName = "vignette",
        icon = icon("book")
      )
    )
  ),

  # Define main windows
  dashboardBody(
    # List tab items ...
    tabItems(
      # ... each tab-item correponds to a menu-item in the sidebar
      tabItem(tabName = "dashboard",
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
              type = "tabs",
              width = 12,
              # Panel 1: Data upload widgets
              # TODO fix processing of zip files on server
              # TODO add options for cloud services?
              tabPanel(
                title = "Upload data",
                icon = icon("upload"),
                # Encase i/o buttons in a fluid row environment
                fluidRow(
                  fileInput(
                    inputId = 'zipFile',
                    label = 'Choose a zip file (Max. 1 GB)',
                    multiple = FALSE,
                    accept = c('.zip')
                  ),
                  fileInput(
                    inputId = 'file',
                    label = 'Choose a File containing sample info',
                    multiple = FALSE,
                    accept = c('.txt','.csv','.tsv')
                  ),
                  shinyDirButton(
                    id = 'dir',
                    label = 'Choose directory',
                    title = 'Upload',
                    icon = icon("folder-open")
                  ),
                  actionButton(
                    inputId = "importBam",
                    label = "Import .bam files",
                    icon = icon("align-left")
                  ),
                  actionButton(
                    inputId = "importTest",
                    label = "Load test data",
                    icon = icon("file-import")
                  ),
                  downloadButton(
                    outputId = "report",
                    label = "Print Report"
                  )
                )
              ),
              # Panel 2: Data selection - the user chooses consensus depth for
              # filtering and which samples and assays to show.
              tabPanel(
                title = "Data selection",
                icon = icon("edit"),

                selectInput(
                  inputId = 'consensus',
                  label = 'Consensus Depth:',
                  choices = c(1,2,3,4,5,10,20,30),
                  selected = 3
                ),
                selectInput(
                  inputId = 'samples',
                  label = 'Samples:',
                  choices = '',
                  multiple = TRUE
                ),
                selectInput(
                  inputId = 'assays',
                  label = 'Assays:',
                  choices = '',
                  multiple = TRUE
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
              label = "Absolute counts",
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
            mainPanel(
              tabBox(
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
            title = "Plots",
            status = "primary",
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(
              tabBox(
                type = "tabs",
                width = 12,
                # Panel for the amplicon plots with download button
                tabPanel(
                  title = "Amplicons",
                  icon = icon("chart-bar"),
                  fluidRow(
                    plotOutput(
                      outputId = "amplicon_plot",
                      width = "130%",
                      height = "600px"
                    ),
                    downloadButton(
                      outputId = "download_plot",
                      label = "Download figure"
                    )
                  )
                ),
                # Panel for quality control plot
                tabPanel(
                  title = "QC Plot",
                  plotOutput("qcPlot")
                ),
                # Panel for UMI count data
                tabPanel(
                  title = "UMI counts",
                  plotOutput("umiCounts")
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

      tabItem(tabName = "advanced",
        fluidRow(
          shinydashboard::box(
            title = "Advanced data analysis",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            width = 4,
            actionButton(
              inputId = 'runVarCaller',
              label = 'Run variant caller'
            ),
            actionButton(
              inputId = 'mergeReplicates',
              label = "Merge Replicates"
            ),
            actionButton(
              inputId = 'timeSeries',
              label = "Analyse time series"
            ),
            actionButton(
              inputId = "generateVCF",
              label = "Generate VCF file"
            ),
            selectInput(
              inputId = 'replicates',
              label = 'Replicates:',
              choices = '',
              multiple = FALSE
            ),
            selectInput(
              inputId = 'timeVar',
              label = 'Time variable:',
              choices = '',
              multiple = FALSE
            )
          ),

          shinydashboard::box(
            title = 'Plot Viewer',
            status = 'primary',
            solidHeader = TRUE,
            collapsible = TRUE,
            width = 12,
            mainPanel(
              tabBox(
                type = 'tabs',
                width = 12,
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
                  title = 'Time series'
                ),
                tabPanel(
                  title = 'Variant caller'
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
            mainPanel(
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
                )
              )
            )
          )
        )
      ),

      # TODO Including html vignette causes some issues as it seems the app size
      # becomes fixed to the size of vignette...
      tabItem(tabName = 'vignette',
              h4('User guide')
        #includeHTML(path = system.file('shiny', 'vignette.html',package = 'umiAnalyzer'))
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

  # TODO add function for processing fastq files
  fastq_main <- reactive({

    fastq_path <- input$fastqFiles$datapath

    if ( is.null(fastq_path) ) {
      return(NULL)
    } else {

      # unzip data to temporary directory
      temp_dir <- file.path(tempdir(), 'fastq')

      unzip(
        zipfile = fastq_path,
        list = FALSE,
        exdir = temp_dir,
        unzip = 'internal'
      )

      # execute run docker script in the folder containing the fastq files
      command <- sprintf('run_docker.sh %s', fastq_path)
      system(command = command)

      # return directory containing umierrorcorrect output
      return(temp_dir)
    }
  })

  # Output pdf report upon button click
  output$download_plot <- downloadHandler(
    filename = 'amplicon_plot.pdf',
    content = function(file) {

      if(is.null(filteredData())){
        return(NULL)
      }

      pdf(file)
        umiAnalyzer::generateAmpliconPlots(
          object = filteredData(),
          filter.name = 'user_filter',
          do.plot = TRUE,
          amplicons = input$assays,
          samples = input$samples)
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
        file = metaData
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
      data <- umiAnalyzer::createUmiExperiment(main)
      return(data)
    }
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

      # TODO keep only first two decimal points to keep output format reasonable
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

  # filteredData returns an updated version of the experimen() object containing
  # a single filter called "user_filter" which is used downstream
  filteredData <- reactive({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- umiAnalyzer::filterUmiObject(
      object = experiment(),
      minDepth = input$consensus,
      minCoverage = 100,
      minFreq = input$minFreq/100,
      minCount = input$minCount
    )

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
    orderClasses = T,
    pageLength = 5,
    lengthMenu = c(5, 10, 50, 100)
  ))


  output$metaDataTable <- DT::renderDataTable({

    if (is.null(metaData())){
      return(NULL)
    }

      metaData()

  }, options = list(
    orderClasses = T,
    pageLenght = 50,
    lengthMenu = c(10, 50, 100)
  ))


  # Generate amplicon plots using umiAnalyzer
  # TODO prettify output on generateAmpliconPlots
  output$amplicon_plot <- renderPlot({

    if(is.null(filteredData())){
      return(NULL)
    }

    umiAnalyzer::generateAmpliconPlots(
      object = filteredData(),
      do.plot = TRUE,
      amplicons = input$assays,
      samples = input$samples,
      abs.count = input$abs_counts
    )
  })

  # Output the QC plot
  # TODO QC plots need a lot more attention. Also consider adding raw read data if available
  output$qcPlot <- renderPlot({

    if(is.null(experiment())){
      return(NULL)
    }

    umiAnalyzer::generateQCplots(
      object = experiment(),
      do.plot = TRUE,
      group.by = 'assay',
      plotDepth = input$consensus,
      assays = input$assays,
      samples = input$samples
      )

  })

  #observeEvent(input$timeSeries, {
  #
  #  output$timeSeriesPlot <- renderPlot({
  #
  #    if(is.null(filteredData() || is.null(metaData()) )){
  #      return(NULL)
  #    } else {
  #      umiAnalyzer::analyzeTimeSeries(
  #        object = filteredData(),
  #        time.var = 'time',
  #        group.by = 'replicate'
  #      )
  #    }
  #  })
  #})

  output$umiCounts <- renderPlot({

    if(is.null(experiment())){
      return(NULL)
    }

    umiAnalyzer::plotUmiCounts(
      object = experiment(),
      do.plot = TRUE
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
        recursive = FALSE)
      # Parse each sample folder for .bam files containing consensus reads
      reads <- umiAnalyzer::parseBamFiles(
        mainDir = main,
        sampleNames = samples,
        consDepth = 0
      )
      # Output barcode family histogram
      output$histogram <- renderPlot({
        # Generate histogram plot using user defined parameters
        umiAnalyzer::plotFamilyHistogram(
          object = reads,
          xMin = input$famSize[1],
          xMax = input$famSize[2],
          samples = input$samples
        )
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

