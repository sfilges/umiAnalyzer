#----------------------------// UmiVisualizer //--------------------------------
#
# A Shiny app for visualizing data generated with the UMIErrorCorrect pipeline:
# https://github.com/tobbeost/umierrorcorrect
#
# This app uses and comes supplied with the umiAnalyzer package:
# https://github.com/sfilges/umiAnalyzer
#

# Quietly import required packages
suppressMessages(library(umiAnalyzer, quietly = TRUE))

suppressMessages(library(dplyr, quietly = TRUE))
suppressMessages(library(tidyr, quietly = TRUE))
suppressMessages(library(readr, quietly = TRUE))
suppressMessages(library(tibble, quietly = TRUE))
suppressMessages(library(plotly, quietly = TRUE))

suppressMessages(library(shiny, quietly = TRUE))
suppressMessages(library(shinyFiles, quietly = TRUE))
suppressMessages(library(shinyWidgets, quietly = TRUE))
suppressMessages(library(DT, quietly = TRUE))
suppressMessages(library(shinydashboard, quietly = TRUE))


#----UI----
# Maximum 5GB data upload
options(
  shiny.maxRequestSize=5000*1024^2,
  shiny.reactlog=TRUE
)

# Define user interface
ui <- dashboardPage(
  dashboardHeader(
    title = 'umiVisualizer'
  ),
  #--------- Define menu items on the sidebar -----------------
  dashboardSidebar(
    sidebarMenu(
      menuItem(
        text = 'Dashboard',
        tabName = 'dashboard',
        icon = icon('dna')
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

  #------------ Define main windows ----------
  dashboardBody(
    # List tab items ...
    tabItems(
      # ... each tab-item corresponds to a menu-item in the sidebar
      tabItem(tabName = 'dashboard',
        fluidRow(
          #------------- Box for data upload and selection ---------------
          box(
            title = "Input",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 460,
            width = 4,
                # Encase i/o buttons in a fluid row environment
                fluidRow(
                  style = "margin-bottom: 10px;margin-left: 5px;margin-right: 5px;",
                  fileInput(
                    inputId = 'zipFile', width = "80%",
                    label = 'Choose a zip file (Max. 5 GB)',
                    multiple = FALSE,
                    accept = c('.zip')
                  ),
                  fileInput(
                    inputId = 'file', width = "80%",
                    label = 'Choose a file containing sample metadata',
                    multiple = FALSE,
                    accept = c('.txt','.csv','.tsv')
                  ),
                  fileInput(
                    inputId = 'bed_file', width = "80%",
                    label = 'Choose a bed file with known mutations',
                    multiple = FALSE,
                    accept = c('.bed','.txt','.csv','.tsv')
                  ),
                  dropdown(
                    label = "Options",

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
                      outputId = 'report.html',
                      label = 'Print Report',
                      style = "margin-bottom: 10px;margin-left: 5px;"
                    ),
                    downloadButton(
                      outputId = 'template',
                      label = 'Get Template',
                      style = "margin-bottom: 10px;margin-left: 5px;"
                    ),

                    icon = icon("gear"),
                    width = "100%"
                  )
            )
          ),

          # Box 2: Data selection - the user chooses consensus depth for
          # filtering and which samples and assays to show.
          box(
            title = "Data filters",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 460,
            width = 4,
            style = "margin-bottom: 10px;margin-left: 10px;margin-right: 10px;",

            tabPanel(
              title = 'Data selection',
              icon = icon('edit'),

              selectInput(
                inputId = 'samples', width = "100%",
                label = 'Samples:',
                choices = '',
                size = 12,
                selectize=FALSE,
                multiple = TRUE
              ),
              selectInput(
                inputId = 'assays', width = "100%",
                label = 'Assays:',
                choices = '',
                multiple = TRUE
              )
            )
          ),

          #---------- Box for parameter selection ------------
          box(
            title = "Plotting options",
            status = "primary",
            solidHeader = TRUE,
            collapsible = FALSE,
            height = 460,
            width = 4,
            style = "margin-bottom: 10px;margin-left: 10px;margin-right: 10px;",
            # Tab box with two panels
            tabBox(
              width = 12,
              type = "tabs",
              # Panel 1: Basic plot types
              tabPanel(
                title = "Plot types",
                icon = icon("chart-line"),
                fluidRow(
                  column(12,
                       style = "margin-top: 10px;margin-left: 5px;margin-right: 5px;",
                       materialSwitch(
                         inputId = "abs_counts",
                         label = "Use absolute counts: ",
                         status = "primary"
                       ),
                       materialSwitch(
                         inputId = "stacked",
                         label = "Stacked plot: ",
                         status = "primary"
                       ),
                       materialSwitch(
                         inputId = "classic",
                         label = "Raw error plot: ",
                         status = "primary"
                       ),
                       materialSwitch(
                         inputId = "use_caller",
                         label = "Use variant caller: ",
                         value = TRUE,
                         status = "primary"
                       ),
                       materialSwitch(
                         inputId = "use_bed",
                         label = "Use bed mutations: ",
                         value = FALSE,
                         status = "primary"
                       ),
                       selectInput(
                         inputId = 'consensus', width = "50%",
                         label = 'Consensus Depth:',
                         choices = c(0,1,2,3,4,5,10,20,30),
                         selected = 3
                       )
                    )
                  )
              ),
              tabPanel(
                #---------- Box for parameter selection ------------
                  title = "Filters",
                  icon = icon("sort"),
                  fluidRow(
                    column(6,
                           sliderInput(
                             inputId = "minFreq",
                             label = "Minimum VAF (to plot):",
                             min = 0, max = 1,
                             value = 0, step = 0.01,
                             post = "%", sep = ","
                           ),
                           sliderInput(
                             inputId = "minCount",
                             label = "Min variant count (to plot):",
                             min = 0, max = 10,
                             value = 0, step = 1,
                             post = " reads", sep = ","
                           ),
                           sliderInput(
                             inputId = "manual_cutoff",
                             label =  "Manual cut-off:",
                             min = 0, max = 100,
                             value = 5, step = 1,
                             sep = ","
                           )
                    ),

                    column(6,
                           sliderInput(
                             inputId = "famSize",
                             label =  "Min and Max family size (histogram):",
                             min = 0, max = 500,
                             value = c(0,100), step = 1,
                             post = " reads", sep = ","
                           ),
                           sliderInput(
                             inputId = "fdr_cutoff",
                             label =  "FDR cut-off:",
                             min = 0, max = 0.2,
                             value = 0.05, step = 0.01,
                             sep = ","
                           ),
                           sliderInput(
                             inputId = "minCoverage",
                             label =  "Minimum coverage:",
                             min = 10, max = 1000,
                             value = 100, step = 10,
                             sep = ","
                           )
                    )
                  )

              ),
              # Panel 2: Merging assays
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
                #--------- Panel for the amplicon plots ------
                tabPanel(
                  title = 'Amplicons',
                  icon = icon('chart-bar'),
                  style = 'margin-left: 10px;',
                  fluidRow(
                    # Options for plot customisation:
                    #     - Select colour scheme
                    #     - Select orientation of colors
                    #     - Select ggplot plot theme
                    #     - Set y-axis range
                    dropdownButton(
                      tags$h3('Customise plot'),
                      column(6,
                             selectInput(
                               inputId = 'colors',
                               label = 'Choose colour palette:',
                               choices = c('default','viridis','magma','plasma','inferno','cividis',
                                           'Accent', 'Dark2', 'Paired', 'Pastel1', 'Pastel2',
                                           'Set1', 'Set2', 'Set3')
                             ),
                             selectInput(
                               inputId = 'direction',
                               label = 'Color palette direction:',
                               choices = c('default','reverse')
                             ),
                             selectInput(
                               inputId = 'theme',
                               label = 'Choose theme:',
                               choices = c('umiVisualiser', 'classic','gray','bw','minimal','light')
                             ),
                             sliderInput(
                               inputId = 'font_size_amplicons',
                               label = 'Font size',
                               value = 7, step = 1,
                               min = 1, max = 14
                             ),
                             sliderInput(
                               inputId = 'font_angle_amplicons',
                               label = 'Font angle',
                               value = 45, step = 45,
                               min = 0, max = 90
                             )
                      ),
                      column(6,
                             numericInput(
                               inputId = 'y_min',
                               label = 'y_min',
                               value = 0,
                               min = 0,
                               max = 100
                             ),
                             numericInput(
                               inputId = 'y_max',
                               label = 'y_max',
                               value = NULL,
                               min = 0,
                               max = 100
                             ),
                             numericInput(
                               inputId = 'amplicon_width',
                               label = 'Plot width:',
                               value = 12,
                               min = 1,
                               max = 20
                             ),
                             numericInput(
                               inputId = 'amplicon_height',
                               label = 'Plot height:',
                               value = 6,
                               min = 1,
                               max = 20
                             ),
                             selectInput(
                               inputId = 'amplicon_device',
                               label = 'File type:',
                               choices = c('pdf','png','svg','eps')
                             ),
                             shinyWidgets::materialSwitch(
                               inputId = 'plot_mutation',
                               label = 'Show mutant allele: ',
                               status = 'primary',
                               value = FALSE
                             ),
                             shinyWidgets::materialSwitch(
                               inputId = 'plot_reference',
                               label = 'Show reference base: ',
                               status = 'primary',
                               value = TRUE
                             )
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '500px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotly::plotlyOutput(
                      outputId = 'amplicon_plot',
                      width = '100%'
                    ),
                    downloadButton(
                      outputId = 'download_plot',
                      label = 'Download figure'
                    )
                  )
                ),
                #---------- Panel for quality control plot -------------
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
                        choices = c('viridis','default','magma','plasma','inferno','cividis')
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
                      sliderInput(
                        inputId = 'font_angle_qc',
                        label = 'Font angle',
                        value = 0, step = 45,
                        min = 0, max = 90
                      ),
                      shinyWidgets::materialSwitch(
                        inputId = "show_mean",
                        label = "Show mean depth: ",
                        status = "primary",
                        value = TRUE
                      ),
                      selectInput(
                        inputId = 'line_col_qc',
                        label = 'Choose line colour :',
                        choices = c('blue','red','green','black')
                      ),
                      selectInput(
                        inputId = 'centerpoint',
                        label = 'Choose center:',
                        choices = c('mean','median')
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '500px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotly::plotlyOutput('qcPlot'),
                    #plotOutput('qcPlot'),
                    downloadButton(
                      outputId = 'download_qc_plot',
                      label = 'Download figure'
                    )
                  )
                ),
                #-------------- Panel for mutation heatmap --------------
                tabPanel(
                  title = "Heatmap",
                  style = 'margin-left: 20px;',
                  fluidRow(
                    # Option for plot customization
                    dropdownButton(
                      tags$h3('Customise plot'),
                      selectInput(
                        inputId = 'heatmap_colors',
                        label = 'Choose colour palette:',
                        choices = c('Blues','Reds','Greens','YlOrRd','YlGnBu','RdPu',
                                    'Purples', 'PuBuGn', 'PuBu', 'OrRd', 'Greys',
                                    'GnBu', 'BuPu', 'BuGn', 'Spectral', 'RdYlBu')
                      ),
                      selectInput(
                        inputId = 'cluster_by',
                        label = 'Samples in:',
                        choices = c('columns','rows')
                      ),
                      numericInput(
                        inputId = 'font_size',
                        label = 'Font Size',
                        value = 10,
                        min = 1,
                        max = 30
                      ),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '500px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!')
                    ),
                    plotOutput("heatmap"),
                    downloadButton(
                      outputId = 'download_heatmap.pdf',
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
                    plotOutput('umiCounts'),
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
                ),
                #-------------- Panel for timeseries  --------------
                tabPanel(
                  title = "Time series",
                  style = 'margin-left: 20px;',
                  fluidRow(
                    # Option for plot customization
                    dropdownButton(
                      tags$h3('Customise plot'),
                      circle = FALSE,
                      status = 'default',
                      icon = icon('gear'),
                      width = '500px',
                      tooltip = tooltipOptions(title = 'Click to customise plot!'),
                      fluidRow(
                        column(6,
                               selectInput(
                                 inputId = 'columns', width = "100%",
                                 label = 'Columns',
                                 choices = '',
                                 multiple = FALSE
                               ),
                               selectInput(
                                 inputId = 'rows', width = "100%",
                                 label = 'Rows:',
                                 choices = '',
                                 multiple = FALSE
                               )
                        ),
                        column(6,
                               selectInput(
                                 inputId = 'color_var', width = "100%",
                                 label = 'Color variable:',
                                 choices = '',
                                 multiple = FALSE
                               ),
                               selectInput(
                                 inputId = 'time_var', width = "100%",
                                 label = 'Time variable:',
                                 choices = '',
                                 multiple = FALSE
                               )
                        )
                      )
                    ),
                    plotOutput("time_series"),
                    downloadButton(
                      outputId = 'download_time_series.pdf',
                      label = 'Download figure'
                    )
                  )
                )
              )
            )
          ),
          #-------------- View data tables in collapsable box ------------------
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
                  style = "font-size: 10px;height:500px; overflow-y: scroll;overflow-x: scroll;"
                ),
                tabPanel(
                  title = "Sample info",
                  DT::dataTableOutput("metaDataTable")
                )
              ),
              downloadButton("downloadData.csv", "Download")
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

