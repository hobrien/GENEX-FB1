#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinyBS)

# Application title
#titlePanel("Gene Expression in the Fetal Brain: Sex Biases"),

navbarPage("Gene Expression in the Fetal Brain: Sex Biases:",
           tabPanel("Plots",
                    sidebarLayout(
                      sidebarPanel(
                        textInput("geneID", "Gene symbol or Ensembl ID"),
                        radioButtons("ages", "Post-Conception Weeks", c('12', '13', '14', '15-16', '17-19', '12-19'), selected = '12-19', inline = FALSE,
                                     width = NULL)
                        #plotOutput("sampleSizeHist", height=200)
                      ),
                      
                      # Show a plot of the generated distribution
                      mainPanel(
                        plotOutput("distPlot"),
                        plotOutput("timeCourse")
                      )
                    )
           ),
           tabPanel("Sex Diffs",
                     sidebarLayout(
                       sidebarPanel(
                         radioButtons("Bias", "Bias Direction", c('Male Bias'='MaleUp', 'Female Bias'='FemaleUp', 'Both'), selected = 'Both', inline = FALSE,
                                      width = NULL),
                         checkboxGroupInput("ChrType", "Chromosome types", 
                                            choices = c('Autosomes'='autosomal', 'ChrX'='chrX', 'ChrY'='chrY'), selected = c('autosomal', 'chrX', 'chrY'),
                                            inline = FALSE, width = NULL),
                         radioButtons("p_type", "Maximum p-value", c('Uncorrected p-values' = 'pvalue', 'FDR corrected p-values (q-values)'= 'padj'), selected = 'padj', inline = FALSE,
                                      width = NULL),
                         sliderInput("pvalue", "p-value:", 
                                     min = 0, max = 1, value = 0.1, step= 0.01),
                         textInput("typedPval", "Type p-value", value=.1),
                         conditionalPanel(
                           'input.dataset === "12-19 PCW"',
                           downloadButton('download12_19', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist12_19", height=200)
                         ),
                         conditionalPanel(
                           'input.dataset === "12 PCW"',
                           downloadButton('download12', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist12", height=200)
                         ),
                         conditionalPanel(
                           'input.dataset === "13 PCW"',
                           downloadButton('download13', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist13", height=200)
                         ),
                         conditionalPanel(
                           'input.dataset === "14 PCW"',
                           downloadButton('download14', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist14", height=200)
                         ),
                         conditionalPanel(
                           'input.dataset === "15-16 PCW"',
                           downloadButton('download15_16', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist15_16", height=200)
                         ),
                         conditionalPanel(
                           'input.dataset === "17-19 PCW"',
                           downloadButton('download17_19', 'Download Table'),
                           HTML("<br><br>"),
                           plotOutput("sampleSizeHist17_19", height=200)
                         )
                       ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('12-19 PCW', DT::dataTableOutput('mytable1'),
                                    plotOutput("timepoint_12_19")
                                    ),
                           tabPanel('12 PCW', DT::dataTableOutput('mytable2'),
                                    plotOutput("timepoint_12")
                                    ),
                           tabPanel('13 PCW', DT::dataTableOutput('mytable3'),
                                    plotOutput("timepoint_13")
                                    ),
                           tabPanel('14 PCW', DT::dataTableOutput('mytable4'),
                                    plotOutput("timepoint_14")
                                    ),
                           tabPanel('15-16 PCW', DT::dataTableOutput('mytable5'),
                                    plotOutput("timepoint_15_16")
                                    ),
                           tabPanel('17-19 PCW', DT::dataTableOutput('mytable6'),
                                    plotOutput("timepoint_17_19")
                           )
                         )
                       )
                     )
            ),
           tabPanel("Expression Trajectory",
                    sidebarLayout(
                      sidebarPanel(
                        radioButtons("Direction", "Change with time", c('Upregulated'='MaleUp', 'Downregulated'='FemaleUp', 'Both'), selected = 'Both', inline = FALSE,
                                     width = NULL),
                        checkboxGroupInput("ChrTypePCW", "Chromosome types", 
                                           choices = c('Autosomes'='autosomal', 'ChrX'='chrX', 'ChrY'='chrY'), selected = c('autosomal', 'chrX', 'chrY'),
                                           inline = FALSE, width = NULL),
                        radioButtons("p_typePCW", "Maximum p-value", c('Uncorrected p-values' = 'pvalue', 'FDR corrected p-values (q-values)'= 'padj'), selected = 'padj', inline = FALSE,
                                     width = NULL),
                        sliderInput("pvaluePCW", "p-value:", 
                                    min = 0, max = 1, value = 0.1, step= 0.01),
                        textInput("typedPvalPCW", "Type p-value", value=.1),
                        HTML("<strong>Plot Expression Trajectory</strong><br>"),
                        actionButton("go", "Plot"),
                        HTML("<br><br><strong>Export Table</strong><br>"),
                        downloadButton('downloadPCW', 'Download')
       
                        
                      ),
                      mainPanel(
                        tabsetPanel(
                          id = 'dataset',
                          tabPanel('Expression shifts over development', DT::dataTableOutput('mytable7'),
                          #plotOutput("timeCourseRowNum"),
                          bsModal("timeCoursePlot", "Expression Trajectory", "go", size = "large",plotOutput("timeCourseRowNum"))
                          )
                        )
                      )
                    )
           )           
)
