################################################################################
### R script to run WCGNA analyses, Inspired by SARTools (Hugo Varet)
### Heath O'Brien
### 10 Jan 2017
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
#setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/")
rm(list=ls())                                        # remove all the objects from the R session
library("optparse")

option_list <- list(
  make_option(c("-v", "--varInt"), type="character", default="Sex,PCW", 
              help="Variables of interest (Sex and PCW)"),
  make_option(c("-n", "--name"), type="character", default="Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts_RIN_excl17046_s50_p7_c30_m10", 
              help="Project name"),
  make_option(c("-o", "--outfolder"), type="character", default="/Volumes/FetalRNAseq/WGCNA", 
              help="Folder for output"),
  make_option(c("-i", "--input"), type="character", default="Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/counts_vst.txt", 
              help="Input counts"),
  make_option(c("-l", "--library"), type="character", default="Data/SampleInfo.txt", 
              help="Library info file (must contain batch effects)"),
  make_option(c("-b", "--batch"), type="character", default="ReadLength", 
              help="factors to be used as batch correction"),
  make_option(c("-c", "--covariate"), type="character", default="RIN", 
              help="numeric cofactors to be corrected for"),
  make_option(c("-e", "--exclude"), type="character", default='17046', 
              help="Samples to exclude (comma separated list, no spaces)", metavar="excluded"),
  make_option(c("-p", "--power"), type="integer", default=7, 
              help="power"),
  make_option(c("-m", "--min"), type="numeric", default=10, 
              help="minimum average transformed count (3.550576 = 0 untransformed)"),
  make_option(c("-c", "--cut"), type="numeric", default=.3, 
              help="mergeCutHeight"),
  make_option(c("-s", "--size"), type="numeric", default=50, 
              help="minimum cluster size"),
  make_option(c("-q", "--quantile"), action='store_true', type="logical", default=FALSE, 
              help="Use quantile normalisation")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)


projectName <- opt$name                        # name of the project
author <- "Heath O'Brien"                                # author of the statistical analysis/report

workDir <- opt$outfolder      # working directory for the R session

#inputDir <- "~/BTSync/FetalRNAseq/Counts/MvsF_12_20_noA_excl_15641_16491_18432_norm"

inputCounts <- opt$input #filter out genes with counts < 10 in more than 10% of samples

batch <- strsplit(opt$batch, ',')[[1]]
covar <- strsplit(opt$covariate, ',')[[1]] 
varInt <- strsplit(opt$varInt, ',')[[1]]
design=formula(paste("~", paste(varInt, collapse= ' + ')))
quantile_norm <- opt$quantile
exclude <- strsplit(opt$exclude, ',')[[1]]

################################################################################
###                             running script                               ###
################################################################################

library(WGCNA)
library(tidyverse)
#library(sva)
library(limma)
library(preprocessCore)
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

dir.create(workDir)
dir.create(paste(workDir, opt$name, sep='/'))
dir.create(paste(workDir, opt$name, 'figures', sep='/'))

# Prepare data for analysis
MalevsFemale <- read_delim(paste(inputCounts, sep='/'), "\t", escape_double = FALSE, trim_ws = TRUE)
if (length(exclude) > 0) {
  MalevsFemale <- select(MalevsFemale, -one_of(exclude))
}
#move Id to rownames
MalevsFemale <- as.data.frame(MalevsFemale)
rownames(MalevsFemale) <- MalevsFemale$Id
MalevsFemale <- MalevsFemale[,-1] 

#Add trait info
sample_info <- read_delim("Data/SampleInfo.txt",
                          "\t", escape_double = FALSE, col_types = cols(Sample = col_character()),
                          trim_ws = TRUE)
#sample_info <- select(sample_info, -one_of(exclude))
datTraits <- data.frame(Sample=as.character(colnames(MalevsFemale)))
datTraits <- left_join(datTraits, sample_info) %>% 
  as.data.frame()
rownames(datTraits) <- datTraits$Sample
datTraits <- datTraits[,-1]
collectGarbage()

setwd(paste(workDir, opt$name, sep='/'))

datTraits$Sex <- as.numeric(as.factor(datTraits$Sex))
datTraits$ReadLength <- as.numeric(as.factor(datTraits$ReadLength))
datTraits <- subset(datTraits, select=-c(Sequencer))

MalevsFemale <- MalevsFemale %>% filter( Reduce(`+`, .)/dim(MalevsFemale)[2] > opt$min)
#Correct for batch effects (see http://jtleek.com/genstats/inst/doc/02_13_batch-effects.html)
#if ('ReadLength' %in% batch) {
#    batchLevel <- datTraits$ReadLength
#    modcombat = model.matrix(~RIN, data=datTraits)
#    MalevsFemale = ComBat(dat=MalevsFemale, batch=batchLevel, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
#}

MalevsFemale <- removeBatchEffect(MalevsFemale, batch=datTraits[,batch], design=model.matrix(design, data=datTraits), covariates=datTraits[,covar]) 

if (quantile_norm) {
    MalevsFemale <- normalize.quantiles(as.matrix(MalevsFemale),copy=TRUE)
}

#transpose, remove metadata about genes
datExpr0 <- as.data.frame(t(MalevsFemale))
rownames(datExpr0) <- colnames(MalevsFemale)
colnames(datExpr0) <- rownames(MalevsFemale)


#test for gene and sample missingness
gsg = goodSamplesGenes(datExpr0, verbose = 3);
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)



traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
png(filename="figures/clusterSamples.png",width=1800,height=1800,res=300) 
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()
# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
# Plot the results:

# Scale-free topology fit index as a function of the soft-thresholding power
png(filename="figures/softThreshold.png",width=1800,height=1800,res=300) 
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

net = blockwiseModules(datExpr0, power =opt$power, maxBlockSize = 20000,
                       TOMType = "unsigned", minModuleSize = opt$size,
                       reassignThreshold = 0, mergeCutHeight = opt$cut,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "MalevsFemaleFilter_8TOM",
                       verbose = 3)

table(net$colors)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
png(filename="figures/modules.png",width=1800,height=1800,res=300) 
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)

png(filename="figures/heatmap.png",width=1800,height=1800,res=300) 
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

save.image(file=paste0(projectName, ".RData"))

