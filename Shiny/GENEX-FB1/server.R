#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
source("FormatGGplot.R")
library(DT)

# setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/Shiny/GENEX-FB1")

################################## Define functions ##################################
PlotExpressionRowNum<-function(row_num, counts, fittedPCW, target) {
  fit_params <- fittedPCW[row_num,]
  geneID <- fit_params$Id
  data <- counts %>% filter(Id == geneID) %>%  
    dplyr::select(-one_of('SYMBOL', 'Id', 'Chr', 'ChrType', 'GeneId')) %>%
    gather() %>%
    separate(key, into=c('norm', 'Sample'), sep='[.]') %>%
    dplyr::select(Sample, value) %>%
    left_join(target)
  mean_age<-mean(target$PCW)
  fit <- data.frame(PCW=seq(12,19)) %>% mutate(fit=fit_params$baseMean*2^(fit_params$log2FoldDiff*(PCW-mean_age)))
  title<-paste0(geneID, ' (', fit_params$SYMBOL, ')')#: log2 change/week = ', fit_params$log2FoldDiff, ', p=', fit_params$pvalue, ', q=', fit_params$padj)
  plot<-  ggplot(data, aes(x=PCW, y=value, colour=Sex)) + 
    geom_jitter(height = 0, width=.1, alpha=.75) + 
    geom_line(aes(y=fit), colour='black', data=fit) +
    scale_x_continuous(breaks=seq(12, 20)) +
    ylab("normalised counts") +
    xlab('post-conception weeks') +
    main_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(title) 
  plot
}
#PlotExpressionRowNum(9337, counts, fittedPCW, target) + scale_y_continuous(limits=c(0,10))

PlotTimepointRowNum<-function(row_num, counts, fitted, target) {
  selection <- fitted[row_num,]
  geneID <- selection$Id
  data <- counts %>% filter(Id == geneID) %>%  
    dplyr::select(-one_of('SYMBOL', 'Id', 'Chr', 'ChrType', 'GeneId')) %>%
    gather() %>%
    separate(key, into=c('norm', 'Sample'), sep='[.]') %>%
    dplyr::select(Sample, value) %>%
    left_join(target)
  mean <- selection %>% dplyr::select(Male, Female) %>%
    gather('Sex', 'mean')
  title<-paste0(selection$Id, ' (', selection$SYMBOL, ')')#: log2 change/week = ', selecton$log2FoldDiff, ', p=', selecton$pvalue, ', q=', selecton$padj)
  plot<-  ggplot(data, aes(x=Sex, colour=Sex)) + 
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    ggtitle(title) +
    scale_colour_brewer(type = "qual", palette = 6) 
  plot
}
#PlotTimepointRowNum(20, counts_tr, all_PCW_tr, target, '12-19')
  
PlotTranscriptsRowNum <- function(counts, selection, fitted, target) {
  #selection <- fitted[row_num,]
  geneID <- selection$GeneId
  symbol <- selection$SYMBOL
  mean <- filter(fitted, GeneId == geneID) %>% dplyr::select(Male, Female, Id, qval=padj) %>%
    gather('Sex', 'mean', -Id, -qval) %>%
    mutate(facet=paste0(Id, '\nFDR=', signif(qval, digits = 3)), Sex=ifelse(Sex=='Male', 'M', ifelse(Sex=='Female', 'F', NA)))

  data <- filter(counts, GeneId == geneID) %>% 
      dplyr::select(-SYMBOL, -GeneId, -Chr, -ChrType) %>%
      gather(key, value, -Id) %>%
      separate(key, into=c('norm', 'Sample'), sep='[.]') %>%
      dplyr::select(Sample, value, Id) %>%
      left_join(target) %>%
      mutate(Sex=ifelse(Sex=='Male', 'M', ifelse(Sex=='Female', 'F', NA))) %>%
      left_join(dplyr::select(mean, Id, qval) %>% group_by(Id) %>% dplyr::slice(1)) %>%
      mutate(facet=paste0(Id, '\nFDR=', signif(qval, digits = 3)))
  
  title<-paste0(selection$GeneId, ' (', selection$SYMBOL, ')')#: log2 change/week = ', selecton$log2FoldDiff, ', p=', selecton$pvalue, ', q=', selecton$padj)
  
  plot<-  ggplot(data, aes(x=Sex, colour=Sex)) + 
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    facet_wrap(~ facet, ncol=6) +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    ggtitle(title) +
    scale_colour_brewer(type = "qual", palette = 6) 
  plot
}
#PlotTranscriptsRowNum(20, counts_tr, fitted_tr, target, '12-19')
#PlotTranscriptsRowNum(counts_tr , filter(fitted_tr, Id=='ENST00000359939'), fitted_tr, target)

filter_table <- function(fitted, ChrTypeList, Bias, p_type, p_val) {
  fitted <- filter(fitted, !is.na(padj) & UQ(as.name(p_type)) <= p_val ) %>%
    arrange(UQ(as.name(p_type))) %>%
    filter(ChrType %in% ChrTypeList) %>% 
    dplyr::select(-ChrType)
  if (Bias == 'MaleUp') {
    fitted <- fitted %>% filter(log2FoldDiff > 0)
  } else if (Bias == 'FemaleUp') {
    fitted <- fitted %>% filter(log2FoldDiff < 0)
  }
  fitted <- mutate(fitted, pvalue = as.numeric(format(pvalue, digits=3)), padj = as.numeric(format(padj, digits=3)))
  fitted
}

add_links <-function(fitted) {
  mutate(fitted, SYMBOL=paste0("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=", SYMBOL, " target='_blank'>", SYMBOL, "</a>"),
         GTEx=paste0("<a href=https://gtexportal.org/home/gene/", Id, " target='_blank'>GTEx</a>"),
         Id=paste0("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", Id, " target='_blank'>", Id, "</a>")
  )         
}

add_links_tr <-function(fitted) {
  mutate(fitted, SYMBOL=paste0("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=", SYMBOL, " target='_blank'>", SYMBOL, "</a>"),
         GeneId=paste0("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", GeneId, " target='_blank'>", GeneId, "</a>")
  )         
}

################################## Load Data ##################################
target <- read_tsv("./Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 

counts <-  read_delim("./Data/counts12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-gene_type)

fitted <- read_delim("./Data/fitted.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange) %>% 
  dplyr::select(-gene_type) %>%
  arrange(padj)

fittedPCW <- read_delim("./Data/dropPCW.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange) %>% 
  dplyr::select(-gene_type) %>%
  arrange(padj)

counts_tr <-  read_delim("./Data/counts12_20_tr.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(-gene_type)
  
  
fitted_tr <- read_delim("./Data/fitted_tr.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange) %>%
  dplyr::select(-gene_type)
  
fittedPCW_tr <- read_delim("./Data/dropPCW_tr.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange) %>%
  dplyr::select(-gene_type)
  

################################## Run server ##################################
shinyServer(function(session, input, output) {
  observe({
    updateSliderInput(session, "pvalue", value = input$typedPval)
    updateSliderInput(session, "pvaluePCW", value = input$typedPvalPCW)
  })
  observe({
    updateSliderInput(session, "pvalue", value = input$typedPval)
    updateSliderInput(session, "pvaluePCW", value = input$typedPvalPCW)
  })

  ################################## Render plots ################################## 
  output$SEXdiffs<- renderPlot({
    validate(
      need(input$SexDiffTable_rows_selected != "", "Please select a row from the table")
    )
    PlotTimepointRowNum(input$SexDiffTable_rows_selected, counts , filter_table(fitted, input$ChrType, input$Bias, input$p_type, input$pvalue), target)
  })
  output$SEXdiffs_trans <- renderPlot({
    validate(
      need(input$SexDiffTrTable_rows_selected != "", "Please select a row from the table")
    )
    PlotTimepointRowNum(input$SexDiffTrTable_rows_selected, counts_tr, filter_table(fitted_tr, input$ChrType, input$Bias, input$p_type, input$pvalue), target)
  })
  output$SEXdiffs_all_trans <- renderPlot({
    validate(
      need(input$SexDiffTrTable_rows_selected != "", "Please select a row from the table")
    )
    PlotTranscriptsRowNum(counts_tr, filter_table(fitted_tr, input$ChrType, input$Bias, input$p_type, input$pvalue)[input$SexDiffTrTable_rows_selected,], fitted_tr, target)
  })
  output$distPlot <- renderPlot({
    req(input$geneID)
    PlotTimepoint(toupper(input$geneID), counts, fitted, target, input$ages)
  })
  output$timeCourseRowNum <- renderPlot({
    validate(
      need(input$PCWTable_rows_selected != "", "Please select a row from the table")
    )
    PlotExpressionRowNum(input$PCWTable_rows_selected, counts, filter_table(fittedPCW, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW), target)
  })
  output$timeCourseRowNum_tr <- renderPlot({
    validate(
      need(input$PCWTableTr_rows_selected != "", "Please select a row from the table")
    )
    PlotExpressionRowNum(input$PCWTableTr_rows_selected, counts_tr, filter_table(fittedPCW_tr, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW), target)
  })
  output$sampleSizeHist <- renderPlot({
    PlotSampleSize(target, input$ages)
  })

  ################################## Render tables ################################## 
  output$SexDiffTable <- DT::renderDataTable({
    DT::datatable(add_links(filter_table(fitted, input$ChrType, input$Bias, input$p_type, input$pvalue)), escape = FALSE, selection="single", caption = 'Genes exhibiting sex differences in fetal brain expression')
  })
  output$SexDiffTrTable <- DT::renderDataTable({
    DT::datatable(add_links_tr(filter_table(fitted_tr, input$ChrType, input$Bias, input$p_type, input$pvalue)), escape = FALSE, selection="single", caption = 'Genes exhibiting sex differences in fetal brain expression')
  })
  output$PCWTable <- DT::renderDataTable({
    DT::datatable(filter_table(fittedPCW, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW) %>% add_links(), escape = FALSE, selection="single", caption = 'Genes exhibiting differences in fetal brain expression over development')
  })
  output$PCWTableTr <- DT::renderDataTable({
    DT::datatable(filter_table(fittedPCW_tr, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW) %>% add_links_tr(), escape = FALSE, selection="single", caption = 'Genes exhibiting differences in fetal brain expression over development')
  })
  
  ################################## Download tables ################################## 
  output$downloadSEX <- downloadHandler(
    filename = function() { 'SEXdiffs.txt' },
    content = function(file) {
      filter_table(fitted, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }  
  )
  output$downloadSEX_tr <- downloadHandler(
    filename = function() { 'SEXdiffs_tr.txt' },
    content = function(file) {
      filter_table(fitted_tr, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }  
  )
  output$downloadPCW <- downloadHandler(
    filename = function() { 'PCWdiffs.txt' },
    content = function(file) {
      filter_table(fittedPCW, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW) %>%
        write_tsv(file)
    }
  )
  output$downloadPCW_tr <- downloadHandler(
    filename = function() { 'PCWdiffs_tr.txt' },
    content = function(file) {
      filter_table(fittedPCW_tr, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW) %>%
        write_tsv(file)
    }
  )
})
