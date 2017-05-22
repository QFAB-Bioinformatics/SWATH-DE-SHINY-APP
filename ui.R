require(shiny)
require(shinythemes)
require(shinyjs)
shinyUI(fluidPage(theme= shinytheme("superhero"),
                  useShinyjs(),
                  titlePanel("SWATH differential expression"),
                  
                  div(id = "dselect",
                      sidebarPanel(width = 5,
                                   fileInput("dataFile", h4("Choose CSV File containing the data :"),
                                             accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                   "Note: The input file should be a .CSV file , containing a header and the first column being the protein names, 
                                   the other columns being the intensity for the differents samples.",
                                   
                                   br(),
                                   
                                   
                                   conditionalPanel(
                                     condition = "input.manual == false",
                                     fileInput("designFile", h5("Choose CSV File containing the experiment design matrix :"),
                                               accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv")),
                                     fileInput("contrastFile", h5("Choose CSV File containing the contrast matrix :"),
                                               accept=c("text/csv", "text/comma-separated-values,text/plain", ".csv"))
                                   ),
                                   
                                   selectInput("norm",  h5("normalization : "), choices = list('Mean centering and scaling' = "MEAN",
                                                                                               'quantile normalization' = "QUANTILE",
                                                                                               'median centering and scaling' = "MEDIAN",
                                                                                               'no normalization' = "NULL")),
                                   
                                   checkboxInput("manual","manual entry of the experiment design"),
                                   
                                   conditionalPanel(
                                     condition = "input.manual == true",
                                     fluidRow(
                                       column(3, numericInput('nbCond',h5('number of conditions : '),value = 2, min = 2)),
                                       column(7)
                                     ),
                                     br(),
                                     
                                     fluidRow(
                                       column(6, conditionalPanel(condition = "input.nbCond != 0", uiOutput("text"))),
                                       column(6, conditionalPanel(condition = "input.nbCond != 0",uiOutput("select")))
                                     )
                                   ),
                                   
                                   
                                   actionButton("submit","Submit")
                      )
                  ),
                  
                  hidden(div(id = "cselect",
                             sidebarPanel(
                               h4("choose condition to compare : "),
                               
                               conditionalPanel(
                                 condition= "input.manual == true",
                                 numericInput('cond1','condition 1',value = 1, min = 1),
                                 numericInput('cond2','condition 2',value = 2, min = 1)
                               ),
                               
                               conditionalPanel(
                                 condition= "input.manual == false",
                                 numericInput('comp','comparison',value = 1, min = 1)
                               ),
                               
                               
                               h4("choose threshold : "),
                               sliderInput("p", label = "P value", min = 0, max = 0.2, value = 0.05),
                               sliderInput("fc", label = "Fold change", min = 0, max = 10, value = 0.5, step = 0.1),
                               actionButton("back","back"),
                               
                               width = 2
                             )
                  )),
                  
                  mainPanel(
                    hidden(div(id = "plots",
                               tabsetPanel(
                                 tabPanel("Volcano plot", plotlyOutput("volcanoPlot", height = "720px", width = "1280px"),downloadButton('dlvp', 'Download Volcano plot (.svg)')),
                                 tabPanel("boxplot", plotOutput("unnormalizedPlot", height = "600px"), plotOutput("normalizedPlot", height = "600px"))
                               ),
                               br(), 
                               
                               tabsetPanel(
                                 tabPanel("all proteins", dataTableOutput("allProteins"),downloadButton('dld', 'Download the table with all proteins (.csv)')),
                                 tabPanel("significant proteins",dataTableOutput("significantProteins"),downloadButton('dlsd', 'Download the table with selected proteins (.csv)')),
                                 br()
                               )
                    ))
                  )
)) 
