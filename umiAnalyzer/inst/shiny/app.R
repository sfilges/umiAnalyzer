
if (!requireNamespace(c('tidyverse','shiny','shinyFiles','DT'), quietly = TRUE)) {
  install.packages(c('tidyverse','shiny','shinyFiles','DT'))
}

library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(DT, quietly = TRUE)
library(shinydashboard)

options(shiny.maxRequestSize=200*1024^2)

# Define UI for application that draws a histogram
# Define UI for application that draws a histogram
ui <- dashboardPage(
  dashboardHeader(title = "Visualise UMI data"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Dashboard", tabName = "dashboard", icon = icon("dashboard")),
      menuItem("Settings", tabName = "settings", icon = icon("cog"))
    )
  ),

  dashboardBody(
    fluidRow(
      box(
        fileInput("file",
                  "Choose consensus file",
                  accept = c(".cons")
        ),

        shinyDirButton('dir',
                       'Choose directory',
                       'Upload'),

        br(),

        selectInput('consensus',
                    label = 'Consensus Depth:',
                    choices = 3,
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
          tabPanel("Amplicons", plotOutput("amplicon_plot")),
          tabPanel("Data", DT::dataTableOutput("dataTable")),
          tabPanel("qcPlot", plotOutput("qcPlot"))
        )
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

    if (identical(main, character(0))) {
      return(NULL)
    } else {
      samples <- list.dirs(path = main, full.names = FALSE, recursive = FALSE)

      simsen <- createUmiExperiment(experimentName = "simsen",
                                    mainDir = main,
                                    sampleNames = samples)
      return(simsen)
    }
  })

  observe({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- saveConsData( experiment() )

    updateSelectInput(session,
                      inputId = 'samples',
                      choices = unlist(strsplit(unique(data$`Sample Name`), split = ',')),
                      selected = head(unlist(strsplit(unique(data$`Sample Name`), split = ',')),1))

  })

  output$directorypath = renderPrint({  unique(saveConsData( experiment() )$`Sample Name`)  })

}

# Run the application
shinyApp(ui = ui, server = server)

