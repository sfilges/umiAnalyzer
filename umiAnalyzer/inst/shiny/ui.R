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


#----UI----

# Maximum 5GB data upload
options(shiny.maxRequestSize=5000*1024^2, shiny.reactlog=TRUE)

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
          # Box for interactive parameter selection
          box(
            title = "Parameters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 420,
            style = "margin-bottom: 10px;margin-left: 10px;margin-right: 10px;",
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
            title = 'Advanced data analysis',
            status = 'primary',
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

