################################################################################
### R script based on SARTools by Hugo Varet
### Heath O'Brien
### 12 April, 2017
### designed to be executed with a modified version of SARTools 1.3.0 
### (https://github.com/hobrien/SARTools)
################################################################################

rm(list=ls())                                        # remove all the objects from the R session
library("optparse")

option_list <- list (
#available: adghjlmnoquwyz
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
  make_option(c("--info"), type="character", default="Data/SampleInfo.txt", 
              help="Sample Info"),
  make_option(c("--male"), type="integer", default=NULL, 
              help="number of male samples"),
  make_option(c("--female"), type="integer", default=NULL, 
              help="number of female samples"),
  make_option(c("-b", "--batch"), type="character", default="RIN,ReadLength", 
              help="factors to be used as batch correction"),
  make_option(c("-t", "--tool"), type="character", default="DESeq", 
              help="Tool used for analysis (EdgeR, DESeq, DESeqLRT)"),
  make_option(c("-e", "--exclude"), type="character", default='none', 
              help="Samples to exclude (comma separated list, no spaces)", metavar="excluded"),
  make_option(c("-s", "--sex_chromosomes"), action='store_true', type="logical", default=FALSE, 
              help="Exclude sex chromosomes", metavar="sex_chromosomes"),
  make_option(c("-i", "--interaction"), type="character", default="", 
              help="Cofactors to interact with varInt", metavar = 'interact'),
  make_option(c("-f", "--feature"), type="character", default="genes", 
              help="Type of feature to Analyse (genes, junctions, transcripts)"),
  make_option(c("--sva"), type="integer", default=0, 
              help="Number of Surrogate Variables to estimate"),
  make_option(c("-k", "--kallisto"), action='store_true', type="logical", default=FALSE, 
              help="Use counts derived from Kallisto")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser, positional_arguments=TRUE)
if (is.null(opt$options$min) | is.null(opt$options$max)){
  print_help(opt_parser)
  stop("Min and Max ages must be specified", call.=FALSE)
}
if (! opt$options$feature %in% c('genes', 'junctions', 'transcripts')) {
  stop("feature type not recognised. Must be one of (genes, junctions, transcripts)")
}

if ( opt$options$max-opt$options$min < 1 ) {
  print_help(opt_parser)
  stop("Min age must be less than max", call.=FALSE)
}

if (! opt$options$tool %in% c('EdgeR', 'DESeq', 'DESeqLRT')){
  print_help(opt_parser)
  stop("Tool must be one of EdgeR, DESeq or DESeqLRT", call.=FALSE)
}

if (! opt$options$interact == "" & ! opt$options$tool == 'DESeqLRT') {
  print_help(opt_parser)
  stop("Interactions can only be tested with LRTs", call.=FALSE)
}

PCW_cutoff <- c(opt$options$min, opt$options$max)
RIN_cutoff <- opt$options$rin
alpha <- opt$options$pvalue 
male <- opt$options$male
female <- opt$options$female
varInt <- opt$options$varInt  # factor of interest

ageBin <- ifelse(PCW_cutoff[2]-PCW_cutoff[1] > 1, 
                 paste(PCW_cutoff, 
                       collapse ='_'
                 ),
                 PCW_cutoff[1]
)
interact <- strsplit(opt$options$interact, ',')[[1]]

exclude <- strsplit(opt$options$exclude, '_')[[1]]
if ( opt$options$batch == 'none') {
  batch <- strsplit(opt$options$cofactor, ',')[[1]]
} else{
    batch <- c(strsplit(opt$options$batch, ',')[[1]], strsplit(opt$options$cofactor, ',')[[1]])
}

projectName <- paste0(varInt, 
       ifelse(nchar(opt$options$cofactor)>0, paste0('_', opt$options$cofactor, collapse = ''), ''),
       ifelse(length(interact)>0, paste0('_x_', interact, collapse = ''), ''),
       '_', ageBin,
       ifelse(!is.na(RIN_cutoff), paste0('_RIN', RIN_cutoff,'_'), '_'),
       ifelse(!is.null(male), 'Subsample_', ''),
       'FDR_', alpha,
       '_',
       opt$options$tool,
       '_',
       opt$options$feature,
       ifelse(length(exclude > 0),
              paste(c('_excl', exclude), collapse='_', sep='_'), ''),
       ifelse(opt$options$sex_chromosomes, '_autosomes', ''),
       ifelse(opt$options$sva > 0, paste0("_sva", opt$options$sva), ""),
       ifelse(opt$options$kallisto, "_kallistoCounts", "")
)                         # name of the project
print(paste("Saving output to", projectName))

author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- paste("Results", projectName, sep='/')      # working directory for the R session

rawDir <- "Counts"                                      # path to the directory containing raw counts files

targetFile <- opt$options$info

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
cooksCutoff <-  FALSE #0.75                          # TRUE/FALSE to perform the outliers detection (default is TRUE)
independentFiltering <- TRUE                         # TRUE/FALSE to perform independent filtering (default is TRUE)
# p-value adjustment method: "BH" (default) or "BY"
if (opt$options$tool == 'DESeqLRT') {
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
library(stringr)
library(RColorBrewer)

# loading target file
LibraryInfo <- read_tsv(targetFile, 
                        col_types = cols(.default = col_character())
) %>%
  dplyr::select(one_of('Sample', varInt, batch))
  
if (opt$options$feature == 'genes') {
  LibraryInfo <- mutate(LibraryInfo, Files=paste0(Sample, '.chr.counts.txt')) 
} else if (opt$options$feature == 'junctions') {
  LibraryInfo <- mutate(LibraryInfo, Files=paste0(Sample, '.junctions.txt')) 
}

LibraryInfo <- LibraryInfo %>%
  mutate(PCW=as.numeric(PCW)) %>%
  arrange(Sex)
if ('RIN' %in% batch) {
  LibraryInfo <- LibraryInfo %>% mutate(RIN=as.numeric(RIN))
}
if (!is.na(RIN_cutoff)) {
  LibraryInfo <- filter(LibraryInfo, RIN >= RIN_cutoff)
}
if (!is.null(PCW_cutoff)) {
  LibraryInfo <- filter(LibraryInfo, PCW >= PCW_cutoff[1] & PCW < PCW_cutoff[2])
}

if (length(exclude) > 0 & exclude[1] != 'none') {
  LibraryInfo <- dplyr::filter(LibraryInfo, !Sample %in% exclude)
}

LibraryInfo$Sex=factor(LibraryInfo$Sex)

if (!is.null(male)) {
  MaleSamples <- filter(LibraryInfo, Sex=='Male')
  MaleSamples <- MaleSamples[sample(nrow(MaleSamples), male),]
  FemaleSamples <- filter(LibraryInfo, Sex=='Female')
  FemaleSamples <- FemaleSamples[sample(nrow(FemaleSamples), female),]
  LibraryInfo <- bind_rows(MaleSamples, FemaleSamples)
}
LibraryInfo <- as.data.frame(LibraryInfo)

# loading counts
if (opt$options$kallisto) {
    library(tximport)
    tx2gene <- read_tsv("Data/tx2gene.txt")
    files <- file.path("Kallisto", LibraryInfo$Sample, "abundance.tsv")
    names(files) <- str_split(files, '/')[[1]][2]
    if ( opt$options$feature == 'genes' ) {
      counts <- tximport(files, type = "kallisto", tx2gene = tx2gene, reader=read_tsv)
    } else if ( opt$options$feature == 'transcripts' ) {
      counts <- tximport(files, type = "kallisto", txOut = TRUE, reader=read_tsv)
    } else {
      stop("Kallisto output cannot be used to analyse junctions")
    }
} else {
    counts <- loadCountData(target=LibraryInfo, rawDir=rawDir, featuresToRemove=featuresToRemove)
}
if (opt$options$sex_chromosomes) {
  excludedFeatures <- read_tsv("Data/genes.txt") %>% 
    filter(seqid == 'chrX' | seqid == 'chrY') %>% 
    dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)
  
  if (opt$options$kallisto) {
    counts$counts <- counts$counts[!rownames(counts$counts) %in% excludedFeatures$Id, ]
    counts$length <- counts$length[!rownames(counts$length) %in% excludedFeatures$Id, ]
    counts$abundance <- counts$abundance[!rownames(counts$abundance) %in% excludedFeatures$Id, ]
    counts$counts <- counts$counts[!rownames(counts$counts) %in% excludedFeatures$Id, ]
  } else {  
    counts <- counts[!rownames(counts) %in% excludedFeatures$Id, ]
  }
}
if (varInt == 'PCW'){
  if (opt$options$feature == 'genes') {
    excludedFeatures <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/MalevsFemale.complete.txt",
                                   "\t", escape_double = FALSE, trim_ws = TRUE) %>% filter(is.na(padj))
  } else if (opt$options$feature == 'transcripts') {
    excludedFeatures <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/tables/MalevsFemale.complete.txt",
                                   "\t", escape_double = FALSE, trim_ws = TRUE) %>% filter(is.na(padj))
  }
  if (opt$options$kallisto) {
    counts$counts <- counts$counts[!rownames(counts$counts) %in% excludedFeatures$Id, ]
    counts$length <- counts$length[!rownames(counts$length) %in% excludedFeatures$Id, ]
    counts$abundance <- counts$abundance[!rownames(counts$abundance) %in% excludedFeatures$Id, ]
    counts$counts <- counts$counts[!rownames(counts$counts) %in% excludedFeatures$Id, ]
  } else {  
    counts <- counts[!rownames(counts) %in% excludedFeatures$Id, ]
  }
  independentFiltering <- FALSE
}


dir.create(workDir)
setwd(workDir)

# description plots
if (opt$options$kallisto) {
  majSequences <- descriptionPlots(counts=counts$counts, group=LibraryInfo[,varInt], col=colors)
} else {
  majSequences <- descriptionPlots(counts=counts, group=LibraryInfo[,varInt], col=colors)
}
# edgeR analysis
if ( opt$options$tool == 'EdgeR' ) {
    # checking parameters
    checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch, alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)

    out.edgeR <- run.edgeR(counts=counts, target=LibraryInfo, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

    # MDS + clustering
    exploreCounts(object=out.edgeR$dge, group=LibraryInfo[,varInt], gene.selection=gene.selection, col=colors)

    # summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
    summaryResults <- summarizeResults.edgeR(out.edgeR, group=LibraryInfo[,varInt], counts=counts, alpha=alpha, col=colors)
} else if ( opt$options$tool == 'DESeq' | opt$options$tool == 'DESeqLRT') {
# DEseq analysis
  checkParameters.DESeq2(projectName=projectName,author=author,targetFile=targetFile,
                                              rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                                              condRef=condRef,batch=batch,fitType=fitType,cooksCutoff=cooksCutoff,
                                              independentFiltering=independentFiltering,alpha=alpha,pAdjustMethod=pAdjustMethod,
                                              typeTrans=typeTrans,locfunc=locfunc,colors=colors)

  if (testMethod=='Wald' ) {
    out.DESeq2 <- run.DESeq2(counts=counts, target=LibraryInfo, varInt=varInt, batch=batch, interact=interact, num_sva=opt$options$sva, kallisto=opt$options$kallisto,
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
                               locfunc=locfunc, fitType=fitType, pAdjustMethod=pAdjustMethod, kallisto=opt$options$kallisto,
                               cooksCutoff=cooksCutoff, independentFiltering=independentFiltering, alpha=alpha)
  } else {
    stop("testMethod not recognised")
  }
  # PCA + clustering
  exploreCounts(object=out.DESeq2$dds, group=LibraryInfo[,varInt], typeTrans=typeTrans, col=colors)

  # summary of the analysis (boxplots, dispersions, diag size factors, export table, nDiffTotal, histograms, MA plot)
  summaryResults <- summarizeResults.DESeq2(out.DESeq2, group=LibraryInfo[,varInt], col=colors,
                                          independentFiltering=independentFiltering,
                                          cooksCutoff=cooksCutoff, kallisto=opt$options$kallisto, alpha=alpha)
} else {
  stop("Tool should be one of 'EdgeR', 'DESeq'")
}

################################################################################

write_tsv(as.data.frame(colData(out.DESeq2$dds)), "tables/col_data.txt")
# save image of the R session
save.image(file=paste0(projectName, ".RData"))

gene_info <- read_tsv("../../Data/genes.txt") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid, gene_type)

if ( opt$options$feature == 'transcripts' ) {
  gene_info <- dplyr::rename(tx2gene, Id = transcript_id) %>% right_join(gene_info, by=c("gene_id" = "Id"))
}
if (opt$options$varInt == 'Sex') {
  upfile=paste0("tables/MaleUp", ageBin, ".txt")
  downfile=paste0("tables/FemaleUp", ageBin, ".txt")
  col_names <- c('Id', 'baseMean', 'Male', 'Female', 'FC', 'log2FoldChange', 'pvalue', 'padj', 'maxCooks')
  
} else {
  upfile=paste0("tables/Upregulated", ageBin, ".txt")
  downfile=paste0("tables/Downregulated", ageBin, ".txt")
  col_names <- c('Id', 'baseMean', 'FC', 'log2FoldChange', 'pvalue', 'padj', 'maxCooks')
}

Upregulated <- read.delim(paste('tables', list.files('tables', pattern = ".up.txt$") , sep= '/'), check.names=FALSE)  %>% 
    dplyr::select(one_of(col_names))
if (opt$options$feature == 'junctions') {
  Upregulated <- Upregulated %>% mutate(originalId=Id) %>%
  separate_rows(Id, sep='\\+') %>% 
  mutate(Id= str_replace(Id, '\\..*', ''))
} else{
  Upregulated <- Upregulated %>% mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id))
}
right_join(gene_info, Upregulated) %>%
    write_tsv(upfile)
  
Downregulated <- read.delim(paste('tables', list.files('tables', pattern = ".down.txt$") , sep= '/'), check.names=FALSE)  %>% 
    dplyr::select(one_of(col_names))
if (opt$options$feature == 'junctions') {
  Downregulated <- Downregulated %>% mutate(originalId=Id) %>%
    separate_rows(Id, sep='\\+') %>% 
    mutate(Id= str_replace(Id, '\\..*', ''))
} else{
  Downregulated <- Downregulated %>% mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id))
}
right_join(gene_info, Downregulated) %>% 
    write_tsv(downfile)

# In the case of DESeq, some genes with high expression were filtered out on the basis of high Cooks distance
# This means that we can't use genes with padj == na to determine background
# For DESeq, I am now extracting the filtering threshold and using it to filter the background
Complete <- read.delim(paste('tables', list.files('tables', pattern = ".complete.txt$") , sep= '/'), check.names=FALSE) 
Complete <- Complete %>% filter(! is.na(padj))

Complete <- Complete %>% dplyr::select(one_of(col_names))
if (opt$options$feature == 'junctions') {
  Complete <- Complete %>% mutate(originalId=Id) %>%
    separate_rows(Id, sep='\\+') %>% 
    mutate(Id= str_replace(Id, '\\..*', ''))
} else{
  Complete <- Complete %>% mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id))
}
right_join(gene_info, Complete) %>% 
  write_tsv(paste0("tables/BG", ageBin, ".txt"))

# generating HTML report
if ( opt$options$tool == 'EdgeR' ) {
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


