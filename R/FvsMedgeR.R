################################################################################
### R script based on SARTools by Hugo Varet
### Heath O'Brien
### 12 April, 2017
### designed to be executed with a modified version of SARTools 1.3.0 
### (https://github.com/hobrien/SARTools)
################################################################################
#setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/")

rm(list=ls())                                        # remove all the objects from the R session
library("optparse")

option_list <- list(
  make_option(c("-v", "--varInt"), type="character", default="Sex", 
              help="Variable of interest (either Sex or PCW)"),
  make_option(c("-m", "--min"), type="integer", default=12, 
              help="minimum age (PCW)", metavar="minimum age"),
  make_option(c("-x", "--max"), type="integer", default=13, 
              help="maximum age (PCW)", metavar="maximum age"),
  make_option(c("-r", "--rin"), type="numeric", default=NA, 
              help="minimum RIN", metavar="minimum RIN"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.1, 
              help="corrected pvalue cutoff", metavar="pvalue"),
  make_option(c("-c", "--cofactor"), type="character", default="", 
              help="Cofactors (either Sex or PCW)"),
  make_option(c("--male"), type="integer", default=NULL, 
              help="number of male samples"),
  make_option(c("--female"), type="integer", default=NULL, 
              help="number of female samples"),
  make_option(c("-b", "--batch"), type="character", default="RIN,ReadLength", 
              help="factors to be used as batch correction"),
  make_option(c("-t", "--tool"), type="character", default="DESeq", 
              help="Tool used for analysis (EdgeR, DESeq, DESeqLRT)"),
  make_option(c("-e", "--exclude"), type="character", default='', 
              help="Samples to exclude (comma separated list, no spaces)", metavar="excluded"),
  make_option(c("-s", "--sex_chromosomes"), action='store_true', type="logical", default=FALSE, 
              help="Exclude sex chromosomes", metavar="sex_chromosomes"),
  make_option(c("-i", "--interaction"), type="character", default="", 
              help="Cofactors to interact with varInt", metavar = 'interact')
  
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)
if (is.null(opt$min) | is.null(opt$max)){
  print_help(opt_parser)
  stop("Min and Max ages must be specified", call.=FALSE)
}

if ( opt$max-opt$min < 1 ) {
  print_help(opt_parser)
  stop("Min age must be less than max", call.=FALSE)
}

if (! opt$tool %in% c('EdgeR', 'DESeq', 'DESeqLRT')){
  print_help(opt_parser)
  stop("Tool must be one of EdgeR, DESeq or DESeqLRT", call.=FALSE)
}

if (! opt$interact == "" & ! opt$tool == 'DESeqLRT') {
  print_help(opt_parser)
  stop("Interactions can only be tested with LRTs", call.=FALSE)
}

PCW_cutoff <- c(opt$min, opt$max)
RIN_cutoff <- opt$rin
alpha <- opt$pvalue 
male <- opt$male
female <- opt$female
varInt <- opt$varInt  # factor of interest

ageBin <- ifelse(PCW_cutoff[2]-PCW_cutoff[1] > 1, 
                 paste(PCW_cutoff, 
                       collapse ='_'
                 ),
                 PCW_cutoff[1]
)
interact <- strsplit(opt$interact, ',')[[1]]
exclude <- strsplit(opt$exclude, ',')[[1]]
batch <- strsplit(opt$batch, ',')[[1]]
cofactor <- strsplit(opt$cofactor, ',')[[1]]

projectName <- paste0(varInt, 
       ifelse(length(cofactor)>0, paste0('_', cofactor, collapse = ''), ''),
       ifelse(length(interact)>0, paste0('_x_', interact, collapse = ''), ''),
       '_', ageBin,
       ifelse(!is.na(RIN_cutoff), paste0('_RIN', RIN_cutoff,'_'), '_'),
       ifelse(!is.null(male), 'Subsample_', ''),
       'FDR_', alpha,
       '_',
       opt$tool,
       ifelse(length(exclude > 0),
              paste(c('_excl', exclude), collapse='_', sep='_'), ''),
       ifelse(opt$sex_chromosomes, '_autosomes', '')
)                         # name of the project
print(paste("Saving output to", projectName))

author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- paste("Results", projectName, sep='/')      # working directory for the R session

rawDir <- "Counts"                                      # path to the directory containing raw counts files

targetFile <- "Data/SampleInfo.txt"

featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove


condRef <- "Female"                                      # reference biological condition

pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")


# EdgeR parameters
cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

# DESeq parameters
fitType <- "parametric"                              # mean-variance relationship: "parametric" (default) or "local"
#if numeric, features with maxCooks values above this number are removed 
cooksCutoff <-  0.75 #FALSE                          # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
# p-value adjustment method: "BH" (default) or "BY"
if (opt$tool == 'DESeqLRT') {
  testMethod <- 'LRT'
} else{
  testMethod <- 'Wald'
}
typeTrans <- "VST"                                   # transformation for PCA/clustering: "VST" or "rlog"
locfunc <- "median"                                  # "median" (default) or "shorth" to estimate the size factors

################################################################################
###                             running script                               ###
################################################################################

library(devtools)
load_all(pkg = "R/SARTools")
library(tidyverse)

# loading target file
LibraryInfo <- read_tsv(targetFile, 
                        col_types = cols(.default = col_character())
)
LibraryInfo <- mutate(LibraryInfo, Files=paste0(Sample, '.chr.counts.txt')) %>%
  select(Sample, Files, Sex, PCW, RIN, ReadLength = `Read Length`, Sequencer) %>%
  mutate(PCW=as.numeric(PCW), RIN=as.numeric(RIN)) %>%
  arrange(Sex)
if (!is.na(RIN_cutoff)) {
  LibraryInfo <- filter(LibraryInfo, RIN >= RIN_cutoff)
}
if (!is.null(PCW_cutoff)) {
  LibraryInfo <- filter(LibraryInfo, PCW >= PCW_cutoff[1] & PCW < PCW_cutoff[2])
}

if (length(exclude) > 0) {
  LibraryInfo <- dplyr::filter(LibraryInfo, !Sample == exclude)
}

LibraryInfo <- as.data.frame(LibraryInfo)
LibraryInfo$Sex=factor(LibraryInfo$Sex)

if (!is.null(male)) {
  MaleSamples <- filter(LibraryInfo, Sex=='Male')
  MaleSamples <- MaleSamples[sample(nrow(MaleSamples), male),]
  FemaleSamples <- filter(LibraryInfo, Sex=='Female')
  FemaleSamples <- FemaleSamples[sample(nrow(FemaleSamples), female),]
  LibraryInfo <- bind_rows(MaleSamples, FemaleSamples)
}

# loading counts
counts <- loadCountData(target=LibraryInfo, rawDir=rawDir, featuresToRemove=featuresToRemove)

if (opt$sex_chromosomes) {
  excludedFeatures <- read_tsv("Data/genes.txt") %>% 
    filter(seqid == 'chrX' | seqid == 'chrY') %>% 
    dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)
  
  counts <- counts[!rownames(counts) %in% excludedFeatures$Id, ]
}

dir.create(workDir)
setwd(workDir)

# description plots
majSequences <- descriptionPlots(counts=counts, group=LibraryInfo[,varInt], col=colors)

# edgeR analysis
if ( opt$tool == 'EdgeR' ) {
    # checking parameters
    checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

    out.edgeR <- run.edgeR(counts=counts, target=LibraryInfo, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

    # MDS + clustering
    exploreCounts(object=out.edgeR$dge, group=LibraryInfo[,varInt], gene.selection=gene.selection, col=colors)

    # summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
    summaryResults <- summarizeResults.edgeR(out.edgeR, group=LibraryInfo[,varInt], counts=counts, alpha=alpha, col=colors)
} else if ( opt$tool == 'DESeq' | opt$tool == 'DESeqLRT') {
# DEseq analysis
  checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                                              rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                                              condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                                              independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                                              typeTrans=typeTrans,locfunc=locfunc,colors=colors)

  if (testMethod=='Wald' ) {
    out.DESeq2 <- run.DESeq2(counts=counts, target=LibraryInfo, varInt=varInt, batch=batch, interact=interact,
                           locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                           cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
    if (is.numeric(cooksCutoff)) {
      cooksFiltered <- data.frame(out.DESeq2$results$Male_vs_Female)
      cooksFiltered$maxCooks <- apply(assays(out.DESeq2$dds)[["cooks"]], 1, max)
      cooksFiltered <- cooksFiltered %>% mutate(pvalue = ifelse(is.na(padj) | maxCooks > cooksCutoff, NA, pvalue))
      out.DESeq2$results$Male_vs_Female$pvalue <- cooksFiltered$pvalue
      out.DESeq2$results$Male_vs_Female$padj <- p.adjust(out.DESeq2$results$Male_vs_Female$pvalue, method="BH")
    }  
  } else if (testMethod=='LRT' ) {
    out.DESeq2 <- run.DESeq2.LRT(counts=counts, target=LibraryInfo, varInt=varInt, batch=batch, interact=interact,
                               locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod,
                               cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
  } else {
    stop("testMethod not recognised")
  }
  # PCA + clustering
  exploreCounts(object=out.DESeq2$dds, group=LibraryInfo[,varInt], typeTrans=typeTrans, col=colors)

  # summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
  summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=LibraryInfo[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, alpha=alpha)
} else {
  stop("Tool should be one of 'EdgeR', 'DESeq'")
}

################################################################################

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

gene_info <- read_tsv("../../Data/genes.txt") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid, gene_type)

if (opt$varInt == 'Sex') {
  upfile=paste0("tables/MaleUp", ageBin, ".txt")
  downfile=paste0("tables/FemaleUp", ageBin, ".txt")
  col_names <- c('Id', 'baseMean', 'Male', 'Female', 'FC', 'log2FoldChange', 'pvalue', 'padj')
} else {
  upfile=paste0("tables/Upregulated", ageBin, ".txt")
  downfile=paste0("tables/Downregulated", ageBin, ".txt")
  col_names <- c('Id', 'baseMean', 'FC', 'log2FoldChange', 'pvalue', 'padj')
}

Upregulated <- read.delim(paste('tables', list.files('tables', pattern = ".up.txt$") , sep= '/'), check.names=FALSE)  %>% 
    mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
    dplyr::select(one_of(col_names))
right_join(gene_info, Upregulated) %>% 
    write_tsv(upfile)
  
Downregulated <- read.delim(paste('tables', list.files('tables', pattern = ".down.txt$") , sep= '/'), check.names=FALSE)  %>% 
    mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
    dplyr::select(one_of(col_names))
right_join(gene_info, Downregulated) %>% 
    write_tsv(downfile)

# In the case of DESeq, some genes with high expression were filtered out on the basis of high Cooks distance
# This means that we can't use genes with padj == na to determine background
# For DESeq, I am now extracting the filtering threshold and using it to filter the background
Complete <- read.delim(paste('tables', list.files('tables', pattern = ".complete.txt$") , sep= '/'), check.names=FALSE) 
if (opt$tool == 'EdgeR') {
    Complete <- Complete %>% filter(! is.na(log2FoldChange))
} else{
    Complete <- Complete %>% filter(baseMean > metadata(out.DESeq2$results$Male_vs_Female)$filterThreshold)
}
Complete <- Complete %>% mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
  dplyr::select(one_of(col_names))

right_join(gene_info, Complete) %>% 
  write_tsv(paste0("tables/BG", ageBin, ".txt"))

# generating HTML report
if ( opt$tool == 'EdgeR' ) {
writeReport.edgeR(target=LibraryInfo, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=paste0("MvsF_", PCW_cutoff, collapse='_'), author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod)

} else {
  
  writeReport.DESeq2(target=LibraryInfo, counts=counts, out.DESeq2=out.DESeq2, summaryResults=summaryResults,
                   majSequences=majSequences, workDir=workDir, projectName=projectName, author=author,
                   targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                   condRef=condRef, batch=batch, fitType=fitType, cooksCutoff=cooksCutoff,
                   independentFiltering=independentFiltering, alpha=alpha, pAdjustMethod=pAdjustMethod,
                   typeTrans=typeTrans, locfunc=locfunc, colors=colors)
}


