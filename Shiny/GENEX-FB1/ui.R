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

navbarPage("Fetal Brain Sequencing (FBSeq) 1: ",
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
                         sliderInput("pvalue", "p-value", 
                                     min = 0, max = 1, value = 1, step= 0.01),
                         textInput("typedPval", "Type p-value", value=1),
                         conditionalPanel(
                           condition = 'input.tabs == "Gene-level analysis"',
                           strong("Select row to plot data"),
                           br(),
                           actionButton("PlotSEXdiffs", "Plot"),
                           br(),
                           br(),
                           strong("Export Gene Table"),
                           br(),
                           downloadButton('downloadSEX', 'Download Table')
                         ),   
                         conditionalPanel(
                           condition = 'input.tabs == "Transcript-level analysis"',
                           HTML("<strong>Select row to plot data</strong><br>"),
                           actionButton("PlotSEXdiffsTr", "Plot individual transcript"), 
                           actionButton("PlotSEXdiffsAllTr", "Plot all transcripts from gene"),
                           HTML("<br><br><strong>Export Transcript Table</strong><br>"),
                           downloadButton('downloadSEX_tr', 'Download Table')
                         )   
                       ),
                       mainPanel(
                         tabsetPanel(
                           id = 'tabs',
                           tabPanel('Gene-level analysis', DT::dataTableOutput('SexDiffTable'),
                                    bsModal("sexDiffsPlot", "Sex Differences", "PlotSEXdiffs", size = "large",plotOutput("SEXdiffs"))
                                    ),
                           tabPanel('Transcript-level analysis', DT::dataTableOutput('SexDiffTrTable'),
                                    bsModal("sexDiffsTrans", "Sex Differences", "PlotSEXdiffsTr", size = "large",plotOutput("SEXdiffs_trans")),
                                    bsModal("sexDiffsAllTrans", "Sex Differences", "PlotSEXdiffsAllTr", size = "large",plotOutput("SEXdiffs_all_trans"))
                           )
                           
                         )
                       )
                     )
            )
          
)
