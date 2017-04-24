################################################################################
### R script based on SARTools by Hugo Varet
### Heath O'Brien
### 12 April, 2017
### designed to be executed with a modified version of SARTools 1.3.0 
### (https://github.com/hobrien/SARTools)
################################################################################

getwd()
rm(list=ls())                                        # remove all the objects from the R session
library("optparse")
library(devtools)
load_all(pkg = "SARTools")
library(tidyverse)
source("../Shiny/GENEX-FB1/FormatGGplot.R")

option_list <- list(
  make_option(c("-m", "--min"), type="integer", default=NULL, 
              help="minimum age (PCW)", metavar="minimum age"),
  make_option(c("-x", "--max"), type="integer", default=NULL, 
              help="maximum age (PCW)", metavar="maximum age"),
  make_option(c("-r", "--rin"), type="numeric", default=NA, 
              help="minimum RIN", metavar="minimum RIN"),
  make_option(c("-p", "--pvalue"), type="numeric", default=0.1, 
              help="corrected pvalue cutoff", metavar="pvalue"),
  make_option(c("-a", "--age"), action='store_true', type="logical", default=FALSE, 
              help="Include age as a cofactor", metavar="age")
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

PCW_cutoff <- c(opt$min, opt$max)
RIN_cutoff <- opt$rin
alpha <- opt$pvalue 
BrainBank <- opt$brainbank
pcw <- opt$age
exclude_sex <- opt$sex_chromosomes


ageBin <- ifelse(PCW_cutoff[2]-PCW_cutoff[1] > 1, 
                 paste(PCW_cutoff, 
                       collapse ='_'
                 ),
                 PCW_cutoff[1]
)

projectName <- paste("MvsF", 
                     ageBin,
                     ifelse(pcw, "PCW_FDR", "FDR"),
                     alpha)                         # name of the project


author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- paste("../Results", projectName, sep='/')      # working directory for the R session

rawDir <- "../Counts/"                                      # path to the directory containing raw counts files

targetFile <- "../Shiny/GENEX-FB1/Data/target.txt"

featuresToRemove <- c("alignment_not_unique",        # names of the features to be removed
                      "ambiguous", "no_feature",     # (specific HTSeq-count information and rRNA for example)
                      "not_aligned", "too_low_aQual")# NULL if no feature to remove

if ( pcw ) {
  batch <- c("Centre", "RIN", "PCW")                # blocking factor: NULL (default) or "batch" for example
} else {
  batch <- c("Centre", "RIN")                # blocking factor: NULL (default) or "batch" for example
}

varInt <- "Sex"                                    # factor of interest
condRef <- "Female"                                      # reference biological condition

pAdjustMethod <- "BH"                                # p-value adjustment method: "BH" (default) or "BY"

cpmCutoff <- 1                                       # counts-per-million cut-off to filter low counts
gene.selection <- "pairwise"                         # selection of the features in MDSPlot
normalizationMethod <- "TMM"                         # normalization method: "TMM" (default), "RLE" (DESeq) or "upperquartile"

colors <- c("dodgerblue","firebrick1",               # vector of colors of each biological condition on the plots
            "MediumVioletRed","SpringGreen")


################################################################################
###                             running script                               ###
################################################################################

# loading target file
target <- read_delim("../Shiny/GENEX-FB1/Data/target.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(label=as.character(label))

dir.create(workDir)
setwd(workDir)

# checking parameters
checkParameters.edgeR(projectName=projectName,author=author,targetFile=targetFile,
                      rawDir=rawDir,featuresToRemove=featuresToRemove,varInt=varInt,
                      condRef=condRef,batch=batch,alpha=alpha,pAdjustMethod=pAdjustMethod,
                      cpmCutoff=cpmCutoff,gene.selection=gene.selection,
                      normalizationMethod=normalizationMethod,colors=colors)


# loading counts
counts <- loadCountData(target=target, rawDir=rawDir, featuresToRemove=featuresToRemove)

# description plots
majSequences <- descriptionPlots(counts=counts, group=target[,varInt], col=colors)

# edgeR analysis
out.edgeR <- run.edgeR(counts=counts, target=target, varInt=varInt, condRef=condRef,
                       batch=batch, cpmCutoff=cpmCutoff, normalizationMethod=normalizationMethod,
                       pAdjustMethod=pAdjustMethod)

# MDS + clustering
exploreCounts(object=out.edgeR$dge, group=target[,varInt], gene.selection=gene.selection, col=colors)

# summary of the analysis (boxplots, dispersions, export table, nDiffTotal, histograms, MA plot)
summaryResults <- summarizeResults.edgeR(out.edgeR, group=target[,varInt], counts=counts, alpha=alpha, col=colors)

# save image of the R session
save.image(file=paste0(projectName, ".RData"))

# generating HTML report
writeReport.edgeR(target=target, counts=counts, out.edgeR=out.edgeR, summaryResults=summaryResults,
                  majSequences=majSequences, workDir=workDir, projectName=paste0("MvsF_", PCW_cutoff, collapse='_'), author=author,
                  targetFile=targetFile, rawDir=rawDir, featuresToRemove=featuresToRemove, varInt=varInt,
                  condRef=condRef, batch=batch, alpha=alpha, pAdjustMethod=pAdjustMethod, colors=colors,
                  gene.selection=gene.selection, normalizationMethod=normalizationMethod)



gene_info <- read_tsv("~/BTSync/FetalRNAseq/GENEX-FB1/Data/genes.csv") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)

MalevsFemale.up <- read.delim("tables/MalevsFemale.up.txt", check.names=FALSE)  %>% 
  mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
  dplyr::select(Id, Female, Male, FC, log2FoldChange, pvalue, padj)
right_join(gene_info, MalevsFemale.up) %>% 
  write_tsv(paste0("tables/MaleUp", ageBin, ".txt"))

MalevsFemale.down <- read.delim("tables/MalevsFemale.down.txt", check.names=FALSE)  %>% 
  mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
  dplyr::select(Id, Female, Male, FC, log2FoldChange, pvalue, padj)
right_join(gene_info, MalevsFemale.down) %>% 
    write_tsv(paste0("tables/FemaleUp", ageBin, ".txt"))

MalevsFemale.down <- read.delim("tables/MalevsFemale.complete.txt", check.names=FALSE)  %>%
  filter(MalevsFemale.complete, ! is.na(padj)) %>%
  mutate(Id = sub("(ENSG[0-9]+)\\.[0-9]+", '\\1', Id)) %>%
  dplyr::select(Id, Female, Male, FC, log2FoldChange, pvalue, padj)
right_join(gene_info, MalevsFemale.down) %>% 
  write_tsv(paste0("tables/BG", ageBin, ".txt"))




