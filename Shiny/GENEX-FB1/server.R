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

counts <-  read_delim("./Data/counts12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

fitted <- read_delim("./Data/fitted.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(pvalue = as.numeric(format(pvalue, digits=2)), padj = as.numeric(format(padj, digits=2))) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange)
target <- read_tsv("./Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 

fittedPCW <- read_delim("./Data/dropPCW.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(pvalue = as.numeric(format(pvalue, digits=2)), padj = as.numeric(format(padj, digits=2))) %>%
  dplyr::rename(log2FoldDiff = log2FoldChange)

target <- read_tsv("./Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 

PlotExpressionRowNum<-function(row_num, counts, fittedPCW, target) {
  fit_params <- fittedPCW[row_num,]
  geneID <- fit_params$SYMBOL
  data <- counts %>% filter(SYMBOL == geneID | Id == geneID) %>%  
    dplyr::select(-SYMBOL, -Id, -Chr, -ChrType) %>%
    gather() %>%
    separate(key, into=c('norm', 'Sample'), sep='[.]') %>%
    dplyr::select(Sample, value) %>%
    left_join(target)
  mean_age<-mean(target$PCW)
  fit <- data.frame(PCW=seq(12,19)) %>% mutate(fit=fit_params$baseMean*2^(fit_params$log2FoldDiff*(PCW-mean_age)))
  title<-paste0(geneID, ': log2 change/week = ', fit_params$log2FoldDiff, ', p=', fit_params$pvalue, ', q=', fit_params$padj)
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

PlotTimepointRowNum<-function(row_num, counts, fitted, target, ages) {
  selection <- fitted[row_num,]
  geneID <- selection$SYMBOL
  ageSplit <- strsplit(ages, '-')[[1]]
  min <- ageSplit[1]
  max <- ageSplit[length(ageSplit)]
  data <- counts %>% filter(SYMBOL == geneID | Id == geneID) %>%  
    dplyr::select(-SYMBOL, -Id, -Chr, -ChrType) %>%
    gather() %>%
    separate(key, into=c('norm', 'Sample'), sep='[.]') %>%
    dplyr::select(Sample, value) %>%
    left_join(target) %>%
    filter(PCW >= min & PCW <= max)
  mean <- selection %>% dplyr::select(Male, Female) %>%
    gather('Sex', 'mean')
  plot<-  ggplot(data, aes(x=Sex, colour=Sex)) + 
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    scale_colour_brewer(type = "qual", palette = 6) 
  plot
}

PlotSampleSize<-function(target, ages){
  ageSplit <- strsplit(ages, '-')[[1]]
  min <- ageSplit[1]
  max <- ageSplit[length(ageSplit)]
  target2 <- target %>% filter(PCW >= min & PCW <= max) %>% mutate(age_bin =ifelse(PCW == 16, 15, 
                                           ifelse(PCW > 16, 16, PCW)))
  
  plot <- ggplot(target2, aes(x=age_bin, fill=Sex)) +
    geom_bar() +
    facet_grid(Sex ~ .) +
    side_theme() +
    ggtitle("Sample Size") +
    scale_y_continuous(breaks=seq(0,100, 5)) +
    scale_x_continuous(limits=c(11, 17),
                       breaks=c(11,12,13,14,15,16),
                       labels=c('11','12','13','14','15-16','17-19')) +
    side_theme() +
    xlab("Post Conception Weeks") +
    ylab("Count") +
    scale_fill_brewer(type = "qual", palette = 6) 
  plot
}

shinyServer(function(session, input, output) {
  all_PCW = filter(fitted, ageBin=='12-19' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)

    observe({
    updateSliderInput(session, "pvalue", value = input$typedPval)
    updateSliderInput(session, "pvaluePCW", value = input$typedPvalPCW)
  })
  output$timepoint_12_19 <- renderPlot({
    req(input$mytable1_rows_selected)
    PlotTimepointRowNum(input$mytable1_rows_selected, counts , filter_table(all_PCW, input$ChrType, input$Bias, input$p_type, input$pvalue), target, '12-19')
  })
  output$distPlot <- renderPlot({
    req(input$geneID)
    PlotTimepoint(toupper(input$geneID), counts, fitted, target, input$ages)
  })
  output$timeCourseRowNum <- renderPlot({
    validate(
      need(input$mytable7_rows_selected != "", "Please select a row from the table")
    )
    PlotExpressionRowNum(input$mytable7_rows_selected, counts, filter_table(fittedPCW, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW), target)
  })
  output$sampleSizeHist <- renderPlot({
    PlotSampleSize(target, input$ages)
  })

  
  filter_table <- function(fitted, ChrTypeList, Bias, p_type, p_val) {
    fitted <- filter(fitted, UQ(as.name(p_type)) < p_val ) %>%
      arrange(UQ(as.name(p_type))) %>%
      filter(ChrType %in% ChrTypeList) %>% 
      dplyr::select(-ChrType)
    if (Bias == 'MaleUp') {
      fitted <- fitted %>% filter(log2FoldDiff > 0)
    } else if (Bias == 'FemaleUp') {
      fitted <- fitted %>% filter(log2FoldDiff < 0)
    }
    fitted
  }
  
  add_links <-function(fitted) {
    mutate(fitted, SYMBOL=paste0("<a href=http://www.genecards.org/cgi-bin/carddisp.pl?gene=", SYMBOL, " target='_blank'>", SYMBOL, "</a>"),
           GTEx=paste0("<a href=https://gtexportal.org/home/gene/", Id, " target='_blank'>GTEx</a>"),
           Id=paste0("<a href=http://www.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g=", Id, " target='_blank'>", Id, "</a>")
           )         
  }
  output$mytable1 <- DT::renderDataTable({
    DT::datatable(add_links(filter_table(all_PCW, input$ChrType, input$Bias, input$p_type, input$pvalue)), escape = FALSE, selection="single")
  })
  output$mytable7 <- DT::renderDataTable({
    DT::datatable(filter_table(fittedPCW, input$ChrTypePCW, input$Direction, input$p_typePCW, input$pvaluePCW) %>% add_links(), escape = FALSE, selection="single")
  })
  output$download12_19 <- downloadHandler(
    filename = function() { 'PCW12_19.txt' },
    content = function(file) {
      filter_table(all_PCW, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
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
  
})
