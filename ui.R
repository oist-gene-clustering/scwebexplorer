################################################
#                                              #
#      SHINY APPLICATION: UI PART              #
#                                              #
################################################

##############
# Parameters #
##############

title = "Single Cell Data Web Explorer"

######
# UI #
######

tg <- tags$div()

appCSS <- "
#loading-content {
      position: absolute;
      background: #FFFFFF;
      opacity: 0.9;
      z-index: 100;
      left: 0;
      right: 0;
      height: 100%;
      text-align: center;
      color: #000000;
}
"

#___________________________#
#      Header               #
#___________________________#

header <- dashboardHeader(
  title = title,
  titleWidth = 400
)

#___________________________#
#      Sidebar              #
#___________________________#

sidebar <- dashboardSidebar(
  br(), br(), br(), br(), br(), br(),
  uiOutput("selectedOrganism"), 
  sidebarMenu(id = "sidebarMenu",
              fluidRow(
                column(width=1),
                column(width=1,
                       dropdownButton(
                         HTML("<font color='black'><p><h6>This application was developped
                              during an internship</h6></p> <p><h6>at the OIST 
                              by C. REDA in 2017. Please mail me at:</h6></p>
                              
                              <p><h6><center><b>creda at ens-paris-saclay dot fr</b></center></h6></p>
                              
                              <p><h6>for feedback and error reports.</h6></p></font>"),
                         circle = TRUE, status = "primary", icon = icon("comment"), width = "300px",
                         tooltip = tooltipOptions(title = "Help")
                         )
              )),
              br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(), br(),
              menuItem("Data visualization", tabName="daVisualization", icon=icon('file-image-o')),
              menuItem("Differential expression analysis", tabName="deAnalysis", icon=icon('bar-chart')),
              menuItem("GE pattern explorer", tabName="paExplorer", icon=icon('pie-chart'))
  )
)

#_______________________#
#      Body             #
#_______________________#

body <- dashboardBody(
  
  ## Cell selection ##
  box(title = "Choose the cells to compare", width = 13, solidHeader = TRUE, 
      status = "primary", height = "535px",
      fluidRow(column(width = 3), column(width = 9, HTML("<p><h4>Drag the cells into each bucket</h4></p>"))),
      fluidRow(
        column(width = 1),
        column(width = 9,
                      fluidRow(column(width = 3, h3("Cell group #1", align = "center")),
                               column(width = 3, uiOutput("selectedStage")),
                               column(width = 3, h3("Cell group #2", align = "center"))
                      ),
                      fluidRow(column(width = 2, br(), uiOutput("selectedCells1")),
                               column(width = 5, uiOutput("stageCells")),
                               column(width = 2, br(), uiOutput("selectedCells2")))
               ),
        column(width = 1, br(), br(), textOutput("selectedGenes"), br(), br(), 
               actionButton("deselectAllCells", "Empty both buckets", icon("trash"), class="act"))
      ),
      fluidRow(column(width=12,  
                      tags$hr(),
                      p(class = "text-muted", "Make sure to select at least one cell for each group."),
                      tags$hr(),
                      fluidRow(column(width=1),
                        column(width=1, uiOutput("nextbuttonDATA")),
                        column(width=1),
                        column(width=1, uiOutput("nextbuttonDE")),
                        column(width=1),
                        column(width=1, uiOutput("nextbuttonGE"))
                      )
                      ))
    ),
  
  ## Pop-up events ##
  uiOutput("selectionfromPopup"),
  
  tabItems(
  
  ## Data visualization tab ##
    tabItem("daVisualization",
            fluidRow(
              box(title = "PCA on all stage cells", 
                  width = 6, solidHeader = TRUE, 
                  status = "info", height = "500px",
                  br(),
                  plotOutput("plotVisu")
            ),
            box(title = "Expression value heatmap for selected-stage cells", 
                width = 6, solidHeader = TRUE, 
                status = "info", height = "500px",
                uiOutput("messagegeheatmap"),
                plotOutput("plotHeatmap")
            )
          )
    ),
  
  ## Gene expression pattern explorer tab ##
    tabItem("paExplorer",
              fluidRow(column(width = 3),
                       box(title = "Choice of gene expression level threshold", width = 12, solidHeader = TRUE, 
                           status = "primary", height = "150px",
                           fluidRow(
                             column(width = 1),
                             column(width = 3, 
                                    radioButtons(inputId = "selectThresholdType", label = "", 
                                                            choices = c("Binary", "Graded"), inline = TRUE),
                                    p(class = "text-muted", "The threshold is set in percentage of gene 
                      expression respect to the total gene expression in each cell.")
                                    ),
                             column(width = 8, uiOutput("thresholdIn")))
                       )),
              box(title = "Resulting gene expression patterns", width = 12, solidHeader = TRUE, 
                  status = "info", height = "485px",
                  uiOutput("messagegepattern"),
                  fluidRow(
                    column(width = 5, plotOutput("plotGenePattern")),
                    column(width = 1),
                    column(width = 5, plotOutput("plotGenePattern2"))
                  )
              )
    ),
    
  ## Differential analysis tab ##
    tabItem("deAnalysis",
            fluidRow(
              tabBox(width = 12,
                title = "",
                id = "tabset1", height = "950px",
                tabPanel("Gene selection", 
                         box(title = "Gene expression value table", width = 7, solidHeader = TRUE, 
                             status = "success", height = "850px",
                             p(class = "text-muted", "Make sure to select at least two genes."),
                             tags$hr(),
                             fluidRow(column(width = 3), column(width = 7, HTML("<p><h4>Select genes to study</h4></p>"))),
                             fluidRow(column(width = 1)),
                             fluidRow(
                               column(width = 3),
                               column(width = 1, actionButton("selectGenes", "Select range", icon("check-square"), class="act")),
                               column(width = 1),
                               column(width = 1, actionButton("deselectAllGenes", "Deselect all", icon("trash"), class="act")),
                               column(width = 1, uiOutput("selectionGenes")),
                               column(width = 1, uiOutput("runSelection"))
                             ),
                             br(), br(),
                             fluidRow(
                               column(width = 3, uiOutput("pvalueCutoffSelection")),
                               column(width = 3, br(), br(), uiOutput("pvalueSelection"))
                             ),
                             br(), br(),
                             DT::dataTableOutput("geneData", height = "50")
                         ),
                         box(title = "PCA biplot of cells according to listed genes", width = 5, solidHeader = TRUE, 
                             status = "warning", height = "500px",
                             uiOutput("biplotUI"),
                             plotOutput("plotDEbiplot"))
                ),
                tabPanel("DE heatmap",
                         column(width=3),
                         box(title = "Differential expression analysis heatmap", width = 6, solidHeader = TRUE, 
                             status = "warning", height = "700px",
                             uiOutput("messageheatmap"),
                             plotOutput("mastData")
                         ),
                         column(width=3)
                )
              )
            )
    )
  )
)

#_______________________#
#      Page             #
#_______________________#

fluidPage(
  useShinyjs(),
  inlineCSS(appCSS),
  div(id = "loading_page",
      withSpinner(HTML("<p><center><b><h1>Loading application...</h1></b></center></p>"))
  ),
  hidden(div(id = "main_content",
      dashboardPage(
        header,
        sidebar,
        body,
        title = title,
        skin = "blue"
      )
    )
  )
)