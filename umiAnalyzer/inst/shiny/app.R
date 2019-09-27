
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
  dashboardHeader(title = "Visualise UMI data"),

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
                        multiple = T)
          ),

          mainPanel(

            verbatimTextOutput('directorypath'),

            # Output: Tabset w/ plot, summary, and table ----
            tabsetPanel(
              type = "tabs",
              tabPanel("Amplicons", plotOutput("amplicon_plot", width = "1024px", height = "800px")),
              tabPanel("Data", DT::dataTableOutput("dataTable")),
              tabPanel("qcPlot", plotOutput("qcPlot", width = "1024px", height = "800px")),
              tabPanel("Histogram", plotOutput("histogram", width = "1024px", height = "800px"))
            )
          )
        )
      ),
      tabItem(tabName = "vignette",
              h4("Include vignette")
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


  experiment <- reactive({

    main <- parseDirPath(volumes, input$dir)

    test <- eventReactive(input$importTest,{

      main <- system.file("extdata", package = "umiAnalyzer")

    })

    if(!is.null(test())){
      main <- test()
    }

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

  filteredData <- reactive({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- umiAnalyzer::filterUmiobject(
      object = experiment(),
      name = "user_filter",
      minDepth = input$consensus,
      minCoverage = 100,
      minFreq = 0,
      minCount = 0
    )

    return(data)

    #filter <- umiAnalyzer::getFilterdData(object = data, name = "user_filter")
    #return(filter['user_filter'][[1]])

  })

  observe({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- saveConsData( experiment() )

    updateSelectInput(session,
                      inputId = 'assays',
                      choices = unlist(strsplit(unique(data$Name), split = ',')),
                      selected = head(unlist(strsplit(unique(data$Name), split = ',')),2))

    updateSelectInput(session,
                      inputId = 'samples',
                      choices = unlist(strsplit(unique(data$`Sample Name`), split = ',')),
                      selected = head(unlist(strsplit(unique(data$`Sample Name`), split = ',')),1))

  })

  #output$directorypath = renderPrint({  unique(saveConsData( experiment() )$`Sample Name`)  })

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

  output$amplicon_plot <- renderPlot({


    if(is.null(filteredData())){
      return(NULL)
    }

    generateAmpliconPlots(object = filteredData(),
                          filter.name = "user_filter",
                          do.plot = TRUE,
                          amplicons = input$assays,
                          samples = input$samples)

  })

  output$qcPlot <- renderPlot({


    if(is.null(experiment())){
      return(NULL)
    }

    generateQCplots(simsen, do.plot = TRUE, group.by = "assay")

  })

  observeEvent(input$importBam, {

    main <- parseDirPath(volumes, input$dir)

    if (identical(main, character(0))) {
      return(NULL)
    } else {
      samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

      reads <- parseBamFiles(mainDir = main, sampleNames = samples, consDepth = 0)

      output$histogram <- renderPlot({

        plotFamilyHistogram(object = reads)

      })
    }

  })

}

# Run the application
shinyApp(ui = ui, server = server)

