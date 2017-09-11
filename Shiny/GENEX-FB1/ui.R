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

navbarPage("Fetal Brain Sequencing 1: Sex Differences",
           tabPanel("Sex Differences",
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
                           condition = 'input.tabs == "Gene-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("go", "Plot"),
                           HTML("<br><br><strong>Export Gene Table</strong><br>"),
                           downloadButton('downloadSEX', 'Download Table')
                         ),   
                         conditionalPanel(
                           condition = 'input.tabs == "Transcript-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("go2", "Plot individual transcript"), 
                           #HTML("Plot individual transcript:"), 
                           actionButton("go5", "Plot all transcripts from gene"),
                           #HTML("Plot all transcripts from gene"), 
                           HTML("<br><br><strong>Export Transcript Table</strong><br>"),
                           downloadButton('downloadSEX_tr', 'Download Table')
                         )   
                       ),
                       mainPanel(
                         tabsetPanel(
                           id = 'tabs',
                           tabPanel('Gene-level analysis', DT::dataTableOutput('mytable1'),
                                    bsModal("sexDiffsPlot", "Sex Differences", "go", size = "large",plotOutput("SEXdiffs"))
                                    ),
                           tabPanel('Transcript-level analysis', DT::dataTableOutput('mytable2'),
                                    bsModal("sexDiffsTrans", "Sex Differences", "go2", size = "large",plotOutput("SEXdiffs_trans")),
                                    bsModal("sexDiffsAllTrans", "Sex Differences", "go5", size = "large",plotOutput("SEXdiffs_all_trans"))
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
                        conditionalPanel(
                          condition = 'input.tabs2 == "Gene-level analysis"',
                          HTML("<strong>Select row to plot gene-level data</strong><br>"),
                          actionButton("go3", "Plot"),
                          HTML("<br><br><strong>Export Gene Table</strong><br>"),
                          downloadButton('downloadPCW', 'Download')
                        ),
                        conditionalPanel(
                          condition = 'input.tabs2 == "Transcript-level analysis"',
                          HTML("<strong>Select row to plot transcripts</strong><br>"),
                          actionButton("go4", "Plot"),
                          HTML("<br><br><strong>Export Transcript Table</strong><br>"),
                          downloadButton('downloadPCW_tr', 'Download')
                        )
                      ),
                      mainPanel(
                        tabsetPanel(
                          id = 'tabs2',
                          tabPanel('Gene-level analysis', DT::dataTableOutput('mytable3'),
                                    bsModal("timeCoursePlot", "Expression Trajectory", "go3", size = "large",plotOutput("timeCourseRowNum"))
                          ),
                          tabPanel('Transcript-level analysis', DT::dataTableOutput('mytable4'),
                                  bsModal("timeCoursePlot_tr", "Expression Trajectory", "go4", size = "large",plotOutput("timeCourseRowNum_tr"))
                          )
                        )
                      )
                    )
           )           
)
