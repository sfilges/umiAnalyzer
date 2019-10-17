
library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(shinyWidgets, quietly = TRUE)
library(DT, quietly = TRUE)
library(shinydashboard, quietly = TRUE)
library(umiAnalyzer, quietly = TRUE)

options(shiny.maxRequestSize=200*1024^2)


ui <- dashboardPage(
  dashboardHeader(title = "umiVisualiser"),

  # Define menu items on the sidebar
  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dna")),
      menuItem("Advanced", tabName = "advanced", icon = icon("magic")),
      menuItem("User Guide", tabName = "vignette", icon = icon("book"))
    )
  ),

  # Define main windows
  dashboardBody(
    # List tab items ...
    tabItems(
      # ... each tab-item correponds to a menu-item in the sidebar
      tabItem(tabName = "dashboard",
        fluidRow(
          box(title = "Input",status = "primary", solidHeader = TRUE, collapsible = FALSE,
              height = 420,

            shinyDirButton('dir',
                           'Choose directory',
                           'Upload'),

            actionButton("importBam", "Import BAM files"),

            actionButton("importTest", "Import test data"),

            downloadButton("report", label = "Print Report"),

            br(),

            selectInput('consensus',
                        label = 'Consensus Depth:',
                        choices = c(1,2,3,4,5,10,20,30),
                        selected = 3),

            selectInput('samples',
                        label = 'Samples:',
                        choices = '',
                        multiple = T),

            selectInput('assays',
                        label = 'Assays:',
                        choices = '',
                        multiple = T),

            fileInput('file',
                      'Choose a File containing sample info',
                      multiple = F,
                      accept = c('.txt','.csv','.tsv'))
          ),
          box(title = "Parameters",status = "primary", solidHeader = TRUE, collapsible = FALSE,
              height = 420,

            sliderInput("minFreq", "Minimum Variant allele frequency:",
                        min = 0, max = 1,
                        value = 0, step = 0.01,
                        post = "%", sep = ","),

            sliderInput("minCount", "Minimum Variant allele count:",
                        min = 0, max = 10,
                        value = 0, step = 1,
                        post = " reads", sep = ","),

            sliderInput("famSize", "Minimum and Maximum family size to show:",
                        min = 0, max = 500,
                        value = c(0,100), step = 1,
                        post = " reads", sep = ","),

            br(),

            materialSwitch(
              inputId = "id",
              label = "Absolute counts",
              right = TRUE,
              status = "primary")
          ),

          # View data tables in collapsable box
          box(title = "Data Viewer",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              width = 12,
              style = "font-size: 10px;",
            mainPanel(
              tabBox(
                type = "tabs", width = 12,
                tabPanel("Data", DT::dataTableOutput("dataTable")),
                tabPanel("Sample info", DT::dataTableOutput("metaDataTable"))
              )
            )
          ),

          # Show plots in collapsable box containing a tabBox with a tab for
          # each plot to be shown.
          box(title = "Plots",
              status = "primary",
              solidHeader = TRUE,
              collapsible = TRUE,
              width = 12,
            mainPanel(
              tabBox(type = "tabs",
                     width = 12,

                # Panel for the amplicon plots with download button
                tabPanel("Amplicons",
                  fluidRow(
                    plotOutput("amplicon_plot", width = "130%", height = "600px"),
                    downloadButton("download_plot", "Download figure")
                  )
                ),
                tabPanel("QC Plot", plotOutput("qcPlot")),
                tabPanel("UMI counts", plotOutput("umiCounts")),
                tabPanel("Histogram", plotOutput("histogram"))
              )
            )
          )
        )
      ),

      #
      tabItem(tabName = "advanced",
        h4("Advanced analysis"),

        fluidRow(
          mainPanel(

          )
        )
      ),

      # Including html vignette causes some issues as it seems the app size becomes
      # fixed to the size of vignette...
      tabItem(tabName = "vignette",
        h4("Include vignette")
        #includeHTML(path = system.file("shiny", "vignette.html", package = "umiAnalyzer"))
      )
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session, plotFun) {

  output$report <- downloadHandler(
    # For PDF output, change this to "report.pdf"
    filename = "report.html",
    content = function(file) {
      # Copy the report file to a temporary directory before processing it, in
      # case we don't have write permissions to the current working dir (which
      # can happen when deployed).
      tempReport <- file.path(tempdir(), "report.Rmd")
      file.copy("report.Rmd", tempReport, overwrite = TRUE)

      # Set up parameters to pass to Rmd document
      params <- list(data = filteredData(),
                     assays = input$assays,
                     samples = input$samples,
                     minDepth =  input$consensus)

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

  output$download_plot <- downloadHandler(
    filename = "amplicon_plot.pdf",
    content = function(file) {

      if(is.null(filteredData())){
        return(NULL)
      }

      pdf(file)
        umiAnalyzer::generateAmpliconPlots(
          object = filteredData(),
          filter.name = "user_filter",
          do.plot = TRUE,
          amplicons = input$assays,
          samples = input$samples)
      dev.off()
    }
  )

  # Define avalible volumes for shinyFiles
  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())

  shinyDirChoose(input, 'dir',
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"))

  shinyFileChoose(input, 'file', root=volumes, filetypes=c('.csv','.txt','.tsv'))

  metaData <- reactive({

    if(is.null(experiment())){
      return(NULL)
    }

    metaData <- input$file$datapath

    if (identical(metaData, character(0))) {
      return(NULL)
    } else {

      data <- umiAnalyzer::importDesign(object = experiment(), file = metaData)

      design <- umiAnalyzer::getMetaData(object = data, attributeName = "design")

      return(design)

    }
  })

  # Set up user_data_main
  user_data_main <- reactive({

    # Path selected by the user
    main <- parseDirPath(volumes, input$dir)

    # Create umiExperiment object
    if (identical(main, character(0))){
        return(NULL)
      } else {
        return(main)
      }
  })

  # Set up a test_data main

  test_data_main <- eventReactive(input$importTest,{
    main <- system.file("extdata", package = "umiAnalyzer")
    return(main)
  })

  # Create experiment

  experiment <- reactive({

    # select between main
    if(!is.null(user_data_main())){
      main <- user_data_main()
    } else {
      main <- test_data_main()
    }

    if (is.null(main)) {
      return(NULL)
    } else {

      samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

      data <- umiAnalyzer::createUmiExperiment(
        experimentName = "simsen",
        mainDir = main,
        sampleNames = samples)

      print(class(data))
      print(data@cons.data)

      return(data)
      }

  })

  # filteredData returns an updated version of the experimen() object containing
  # a single filter called "user_filter" which is used downstream
  filteredData <- reactive({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- umiAnalyzer::filterUmiobject(
      object = experiment(),
      name = "user_filter",
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

    updateSelectInput(session,
                      inputId = 'assays',
                      choices = unlist(strsplit(unique(data$Name), split = ',')),
                      selected = head(unlist(strsplit(unique(data$Name), split = ',')),1))

    updateSelectInput(session,
                      inputId = 'samples',
                      choices = unlist(strsplit(unique(data$`Sample Name`), split = ',')),
                      selected = head(unlist(strsplit(unique(data$`Sample Name`), split = ',')),1))

  })

  # Output the consensus data to screen, this will change depending on user input
  # TODO fix getFilterdData function to return a tibble and not a list!
  output$dataTable <- DT::renderDataTable({

    if (is.null(filteredData())){
      return(NULL)
    }

    filter <- umiAnalyzer::getFilterdData(object = filteredData(), name = "user_filter")
    filter <- filter['user_filter'][[1]]

    filter %>%
      dplyr::filter(.data$Name %in% input$assays) %>%
      dplyr::filter(.data$`Sample Name` %in% input$samples)

  }, options = list(
    orderClasses = T,
    pageLenght = 50,
    lengthMenu = c(10, 50, 100)
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

    umiAnalyzer::generateAmpliconPlots(object = filteredData(),
                          filter.name = "user_filter",
                          do.plot = TRUE,
                          amplicons = input$assays,
                          samples = input$samples)

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
      group.by = "assay",
      plotDepth = input$consensus,
      assays = input$assays,
      samples = input$samples
      )

  })

  output$umiCounts <- renderPlot({

    if(is.null(experiment())){
      return(NULL)
    }

    umiAnalyzer::plotUmiCounts(object = experiment(),do.plot = TRUE)

  })



  # Import consensus read bam file upon button click to generate histograms
  # of barcode distribution. It is possible to import directly into the umiExperiment object
  # by setting the importBam flag, but file and this might become very large, so this
  # function is outsourced to parseBamFiles
  observeEvent(input$importBam, {

    #main <- parseDirPath(volumes, input$dir)

    #test2 <- eventReactive(input$importTest,{
    #  main <- system.file("extdata", package = "umiAnalyzer")
    #})

    #if(!is.null(test2())){
    #  main <- test2()
    #}

    # select between main
    if(!is.null(user_data_main())){
      main <- user_data_main()
    } else {
      main <- test_data_main()
    }

    if (identical(main, character(0))) {
      return(NULL)
    } else {
      samples <- list.dirs(path = main,
                           full.names = FALSE,
                           recursive = FALSE)

      reads <- umiAnalyzer::parseBamFiles(mainDir = main,
                                          sampleNames = samples,
                                          consDepth = 0)

      output$histogram <- renderPlot({

        # TODO add min, max and samples arguments to be selected by user
        umiAnalyzer::plotFamilyHistogram(object = reads,
                                         xMin = input$famSize[1],
                                         xMax = input$famSize[2],
                                         samples = input$samples)
      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

