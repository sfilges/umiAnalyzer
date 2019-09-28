
if (!requireNamespace(c('tidyverse','shiny','shinyFiles','DT', 'shinydashboard'), quietly = TRUE)) {
  install.packages(c('tidyverse','shiny','shinyFiles','DT', 'shinydashboard'))
}

library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(DT, quietly = TRUE)
library(shinydashboard, quietly = TRUE)
library(umiAnalyzer, quietly = TRUE)

options(shiny.maxRequestSize=200*1024^2)


ui <- dashboardPage(
  dashboardHeader(title = "umiVisualiser"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Vignette", tabName = "vignette", icon = icon("book"))
    )
  ),

  dashboardBody(
    tabItems(
      tabItem(tabName = "dashboard",
        fluidRow(
          box(
            shinyDirButton('dir',
                           'Choose directory',
                           'Upload'),

            actionButton("importBam", "Import BAM files"),

            actionButton("importTest", "Import test data"),

            br(),

            selectInput('consensus',
                        label = 'Consensus Depth:',
                        choices = c(1,2,3,4,5,10,15,20,30),
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
          box(

            sliderInput("minFreq", "Minimum Variant allele frequency:",
                        min = 0, max = 1,
                        value = 0, step = 0.01,
                        post = "%", sep = ","),

            sliderInput("minCount", "Minimum Variant allele count:",
                        min = 0, max = 10,
                        value = 0, step = 1,
                        post = " reads", sep = ",")
          ),

          mainPanel(

            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(
              type = "tabs",
              tabPanel("Amplicons", plotOutput("amplicon_plot", width = "1024px", height = "800px")),
              tabPanel("Data", DT::dataTableOutput("dataTable")),
              tabPanel("Sample info", DT::dataTableOutput("metaDataTable")),
              tabPanel("QC Plot", plotOutput("qcPlot", width = "1024px", height = "800px")),
              tabPanel("Histogram", plotOutput("histogram", width = "1024px", height = "800px"))
            )
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
server <- function(input, output, session) {

  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())

  path <- shinyDirChoose(input, 'dir',
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

  # Set up a reactive umiExperiment object
  experiment <- reactive({

    # Path selected by the user
    main <- parseDirPath(volumes, input$dir)

    # If importTest is pressed the test data from the umianalyzer package
    # will be loaded.
    # TODO fix so this can't override user data
    test <- eventReactive(input$importTest,{

      main <- system.file("extdata", package = "umiAnalyzer")

    })

    if(!is.null(test())){
      main <- test()
    }

    # Create umiExperiment object
    if (identical(main, character(0))) {
      return(NULL)
    } else {
      samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

      simsen <- umiAnalyzer::createUmiExperiment(experimentName = "simsen",
                                                 mainDir = main,
                                                 sampleNames = samples)
      return(simsen)
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
                      selected = head(unlist(strsplit(unique(data$Name), split = ',')),2))

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

    umiAnalyzer::generateQCplots(object = experiment(),
                                 do.plot = TRUE,
                                 group.by = "assay",
                                 assays = input$assays,
                                 samples = input$samples)

  })

  # Import consensus read bam file upon button click to generate histograms
  # of barcode distribution. It is possible to import directly into the umiExperiment object
  # by setting the importBam flag, but file and this might become very large, so this
  # function is outsourced to parseBamFiles
  observeEvent(input$importBam, {

    main <- parseDirPath(volumes, input$dir)

    test2 <- eventReactive(input$importTest,{

      main <- system.file("extdata", package = "umiAnalyzer")

    })

    if(!is.null(test2())){
      main <- test2()
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

        umiAnalyzer::plotFamilyHistogram(object = reads)

      })
    }
  })
}

# Run the application
shinyApp(ui = ui, server = server)

