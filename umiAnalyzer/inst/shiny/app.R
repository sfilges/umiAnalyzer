#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

if (!requireNamespace(c('tidyverse','shiny','shinyFiles','DT'), quietly = TRUE)) {
  install.packages(c('tidyverse','shiny','shinyFiles','DT'))
}

library(tidyverse, quietly = TRUE)
library(shiny, quietly = TRUE)
library(shinyFiles, quietly = TRUE)
library(DT, quietly = TRUE)

options(shiny.maxRequestSize=200*1024^2)

check_me <- function(list, b){

  a = NULL

  for(i in seq_along(list))
  {
    assays_list <- as.list(strsplit(list[[i]], split = ",")[[1]])
    a[i] <- any(assays_list %in% b)
  }

  return(a)

}


# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("umiVizualisation"),

    # Sidebar with a slider input for number of bins
    sidebarLayout(
        sidebarPanel(

          shinyDirButton('dir',
                         'Choose directory',
                         'Upload'),

          fileInput('file',
                    'Choose Single File',
                    multiple = F,
                    accept = c('.cons')),

          fileInput('bed_file',
                    'Upload BED file',
                    multiple = F),

          #shinyFilesButton('file', 'Upload single file', 'upload', multiple = F),

          br(),

          checkboxGroupInput('samples',
                      label = 'Samples:',
                      choices = ''),

          selectInput('consensus',
                        label = 'Consensus Depth:',
                        choices = 3,
                        selected = 3),

          selectInput('assays',
                        label = 'Assays:',
                        choices = '',
                        multiple = T)

        ),

        # Show a plot of the generated distribution
        mainPanel(

          tabsetPanel(
            tabPanel("Plot",
                     plotOutput('freq_plot')),
            tabPanel("Datatable",
                     DT::dataTableOutput("cons_table")),
            tabPanel("QC plot", plotOutput('qc_plot'))
          )

        )
    )
)



server <- function(input, output, session) {

  volumes <- c(Home = fs::path_home(),
               "R Installation" = R.home(),
               getVolumes()())

  shinyDirChoose(input, "dir",
                 roots = volumes,
                 session = session,
                 restrictions = system.file(package = "base"))


  shinyFileChoose(input, 'file', root=volumes, filetypes=c('.cons'))

  # Updates selection if experiment is changed

  observe({

    if (is.null(experiment())){
      return(NULL)
    }

    data <- experiment()

    updateSelectInput(session,
                      inputId = 'consensus',
                      choices = unique(data$`Consensus group size`),
                      selected = 3)
    updateSelectInput(session,
                      inputId = 'assays',
                      choices = unlist(strsplit(unique(data$Name), split = ',')),
                      selected = head(unlist(strsplit(unique(data$Name), split = ',')),1))
    updateCheckboxGroupInput(session,
                      inputId = 'samples',
                      choices = unique(data$`Sample Name`),
                      selected = head(unique(data$`Sample Name`),1))
  })

  #updates assays if a subset of sample is choosen

  observeEvent(input$samples,{

    data <- cleaned_data() %>%
      filter(`Sample Name` %in% input$samples)

    updateSelectInput(session,
                      inputId = 'assays',
                      choices = unlist(strsplit(unique(data$Name), split = ',')),
                      selected = head(unlist(strsplit(unique(data$Name), split = ',')),1))

  })

  #See if file has been selected

  experiment <- reactive({
    main <- parseDirPath(volumes, input$dir)

    if (identical(main, character(0))) {
      if (is.null(input$file)) {
        return(NULL)
      }

      return(readr::read_tsv(
        file = input$file$datapath,
        col_types = cols(
          `Sample Name` = col_character(),
          Contig = col_character(),
          Position = col_double(),
          Name = col_character(),
          Reference = col_character(),
          A = col_double(),
          C = col_double(),
          G = col_double(),
          T = col_double(),
          I = col_double(),
          D = col_double(),
          N = col_double(),
          Coverage = col_double(),
          `Consensus group size` = col_double(),
          `Max Non-ref Allele Count` = col_double(),
          `Max Non-ref Allele Frequency` = col_double(),
          `Max Non-ref Allele` = col_character()
        )


      ))
    }

    sample.names <- list.dirs(path = main,
                              full.names = FALSE,
                              recursive = FALSE)

    exp1 <- tryCatch({
      #read_file(path = 'asdd')

      data <- umiAnalyzer::createUmiExperiment(experimentName = "exp1",
                                               mainDir = main,
                                               sampleNames = sample.names)

      data@cons.data

    },

    error = function(error_message) {
      return(NULL)
    },

    warning = function(error_message) {
      return(NULL)
    })

    exp1

  })



  #Filter data on consesus depth

  cleaned_data <- reactive({

    experiment() %>%
      filter(!is.na(Name)) %>%
      mutate(A = if_else(Reference == 'A', 0, A)) %>%
      mutate(C = if_else(Reference == 'C', 0, C)) %>%
      mutate(G = if_else(Reference == 'G', 0, G)) %>%
      mutate(`T` = if_else(Reference == 'T', 0, `T`)) %>%
      mutate(I = if_else(Reference == 'I', 0, I)) %>%
      mutate(D = if_else(Reference == 'D', 0, D)) %>%
      mutate(N = if_else(Reference == 'N', 0, N))


    if (is.null(experiment())) {
        return(NULL)
      }

    experiment() %>%
      filter(`Consensus group size` == input$consensus) %>%
      #filter(Coverage > 2) %>%
      filter(`Max Non-ref Allele Frequency` >= 0) %>%
      filter(`Max Non-ref Allele Count` >= 0)

  })

  #Filter based on user couses

  filtered_data <- reactive({


    if (is.null(experiment())) {
      return(NULL)
    }


    cleaned_data() %>%
      #filter(Name %in% input$assays) %>%
      filter(check_me(Name, input$assays)) %>%
      filter(`Sample Name` %in% input$samples)

  })


  output$freq_plot <- renderPlot({


    if(is.null(experiment())){
      return(NULL)}


    filtered_data() %>%
        ggplot() +
        geom_point(aes(x = as.factor(Position), y = `Max Non-ref Allele Frequency`, color = `Sample Name`)) +
        geom_col(aes(x = as.factor(Position), y = 5/Coverage, fill = `Sample Name`), position = position_identity(), alpha = 0.3) +
        #facet_wrap(~Name, scales = 'free', ncol = 1) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              legend.position = 'bottom')
    }

  )

  output$qc_plot <- renderPlot({

    if(is.null(experiment())){
       return(NULL)}

    experiment() %>%
      filter(`Sample Name` %in% input$samples) %>%
      filter(Name %in% input$assays) %>%
      filter(`Consensus group size` == input$consensus | `Consensus group size` == 0) %>%
      ggplot() +
      geom_point(aes(x = Name, y = Coverage, color = as.factor(`Consensus group size`))) +
      scale_y_log10() +
      facet_wrap(~`Sample Name`)

  })

  output$cons_table <- DT::renderDataTable({

    cons_data <- filtered_data() %>%
      select(Assay = Name,
           Chr = Contig,
           Position,
           Reference,
           Coverage,
           Alt = `Max Non-ref Allele`,
           alternativ = `Max Non-ref Allele Frequency`)

  }, options = list(
    orderClasses = T,
    pageLenght = 50,
    lengthMenu = c(10, 50, 100)
  ))

}

# Run the application
shinyApp(ui = ui, server = server)
