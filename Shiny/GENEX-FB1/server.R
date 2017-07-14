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

PlotExpression<-function(geneID, counts, target, ages) {
  geneID = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', geneID)
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
  plot<-  ggplot(data, aes(x=PCW, y=value, colour=Sex)) + 
    geom_jitter(height = 0, width=.1, alpha=.75) + 
    geom_smooth() +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    scale_colour_brewer(type = "qual", palette = 6) 
  plot
}

PlotTimepoint<-function(geneID, counts, fitted, target, ages) {
  geneID = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', geneID)
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
  mean <- filter(fitted, SYMBOL == geneID | Id == geneID) %>% filter(ageBin==ages)
  fc <- mean$log2FoldDiff[1] %>% format(digits=2)
  pval <- mean$pvalue[1]
  qval <- mean$padj[1]
  mean <- mean %>% dplyr::select(Male, Female) %>%
    gather('Sex', 'mean')
  title<-paste0(geneID, ': log2 difference = ', fc, ', p=', pval, ', q=', qval)
  plot<-  ggplot(data, aes(x=Sex, colour=Sex)) + 
    geom_jitter(aes(y=value), height = 0, width=.1, alpha=.75) + 
    geom_errorbar(aes(ymin=mean, ymax=mean), colour='black', size=1, width=.5, data=mean) +
    ylab("normalised counts") +
    xlab('') +
    main_theme() +
    scale_colour_brewer(type = "qual", palette = 6) +
    ggtitle(title) 
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

shinyServer(function(input, output) {
   
  output$distPlot <- renderPlot({
    PlotTimepoint(toupper(input$geneID), counts, fitted, target, input$ages)
  })
  output$timeCourse <- renderPlot({
    PlotExpression(toupper(input$geneID), counts, target, input$ages)
  })
  output$sampleSizeHist <- renderPlot({
    PlotSampleSize(target, input$ages)
  })
  output$summary <- renderPrint({
    summary(cars)
  })
  
  all_PCW = filter(fitted, ageBin=='12-19' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  PCW12 = filter(fitted, ageBin=='12' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  PCW13 = filter(fitted, ageBin=='13' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  PCW14 = filter(fitted, ageBin=='14' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  PCW15_16 = filter(fitted, ageBin=='15-16' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  PCW17_19 = filter(fitted, ageBin=='17-19' & !is.na(padj)) %>% dplyr::select(-ageBin) %>% arrange(padj)
  
  filter_table <- function(fitted, ChrTypeList, Bias, p_type, p_val) {
    fitted <- filter(fitted, UQ(as.name(p_type)) < p_val ) %>%
      arrange(UQ(as.name(p_type))) %>%
      filter(ChrType %in% ChrTypeList) %>% 
      select(-ChrType)
    if (Bias == 'MaleUp') {
      fitted <- fitted %>% filter(log2FoldDiff > 0)
    } else if (Bias == 'FemaleUp') {
      fitted <- fitted %>% filter(log2FoldDiff < 0)
    }  
    fitted
  }
  output$mytable1 <- DT::renderDataTable({
    DT::datatable(filter_table(all_PCW, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$mytable2 <- DT::renderDataTable({
    DT::datatable(filter_table(PCW12, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$mytable3 <- DT::renderDataTable({
    DT::datatable(filter_table(PCW13, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$mytable4 <- DT::renderDataTable({
    DT::datatable(filter_table(PCW14, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$mytable5 <- DT::renderDataTable({
    DT::datatable(filter_table(PCW15_16, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$mytable6 <- DT::renderDataTable({
    DT::datatable(filter_table(PCW17_19, input$ChrType, input$Bias, input$p_type, input$pvalue))
  })
  output$download12_19 <- downloadHandler(
    filename = function() { 'PCW12_19.txt' },
    content = function(file) {
      filter_table(all_PCW, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }  
  )
  output$download12 <- downloadHandler(
    filename = function() { 'PCW12.txt' },
    content = function(file) {
      filter_table(PCW12, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }
  )
  output$download13 <- downloadHandler(
    filename = function() { 'PCW13.txt' },
    content = function(file) {
      filter_table(PCW13, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }
  )
  output$download14 <- downloadHandler(
    filename = function() { 'PCW14.txt' },
    content = function(file) {
      filter_table(PCW14, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }
  )
  output$download15_16 <- downloadHandler(
    filename = function() { 'PCW15_16.txt' },
    content = function(file) {
      filter_table(PCW15_16, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }
  )
  output$download17_19 <- downloadHandler(
    filename = function() { 'PCW17_19.txt' },
    content = function(file) {
      filter_table(PCW17_19, input$ChrType, input$Bias, input$p_type, input$pvalue) %>%
        write_tsv(file)
    }
  )
  
})
