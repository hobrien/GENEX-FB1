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
                         HTML("<strong>Plot Expression Trajectory</strong><br>"),
                         actionButton("go", "Plot"),
                         HTML("<br><br><strong>Export Table</strong><br>"),
                         downloadButton('download12_19', 'Download Table')
                        ),
                       mainPanel(
                         tabsetPanel(
                           id = 'dataset',
                           tabPanel('Genes exhibiting sex differences in fetal brain expression', DT::dataTableOutput('mytable1'),
                                    bsModal("sexDiffsPlot", "Sex Differences", "go", size = "large",plotOutput("timepoint_12_19"))
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
                        actionButton("go2", "Plot"),
                        HTML("<br><br><strong>Export Table</strong><br>"),
                        downloadButton('downloadPCW', 'Download')
       
                        
                      ),
                      mainPanel(
                        tabsetPanel(
                          id = 'dataset',
                          tabPanel('Genes exhibiting differences in fetal brain expression over development', DT::dataTableOutput('mytable7'),
                          #plotOutput("timeCourseRowNum"),
                          bsModal("timeCoursePlot", "Expression Trajectory", "go2", size = "large",plotOutput("timeCourseRowNum"))
                          )
                        )
                      )
                    )
           )           
)
