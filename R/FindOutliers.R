
HighCooks <- function(Id) {
  df <- as.data.frame(assays(out.DESeq2$dds)[['cooks']][Id,])
  colnames(df)<-'cooks_distance'
  df$SampleID<-rownames(df)
  df <- filter(df, cooks_distance>2) %>% arrange(SampleID)
  paste(df$SampleID, collapse='_')
}

load("/Users/heo3/BTSync/FetalRNAseq/Github/GENEX-FB1/Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts.RData")

x<-data.frame('Id'=c('ENSG00000270066', 'ENSG00000143546', 'ENSG00000102837', 
                     'ENSG00000263887', 'ENSG00000250182', 'ENSG00000167978', 
                     'ENSG00000234969', 'ENSG00000229807', 'ENSG00000234944',
                     'ENSG00000229851', 'ENSG00000215506', 'ENSG00000143552', 
                     'ENSG00000227289', 'ENSG00000250421', 'ENSG00000188399', 
                     'ENSG00000251838', 'ENSG00000182162', 'ENSG00000252766', 
                     'ENSG00000231535', 'ENSG00000250951', 'ENSG00000205936'))
excluded=c()
for (Id in x$Id) {
    excluded = c(excluded, HighCooks(Id))
}
x$excluded=excluded

for (duplicated in c('16024_16115', '12993')) {
  print(duplicated)
  x1<- filter(x, excluded==duplicated)
  x <- filter(x, excluded!=duplicated) %>% 
    bind_rows(data.frame(Id=paste(x1$Id, collapse='_'), excluded=duplicated))
}

#x <- mutate(x, excluded=paste0("'", excluded, "'"))
myList <- as.list(x$Id)
names(myList) <- as.character(x$excluded)

load("/Users/heo3/BTSync/FetalRNAseq/Github/GENEX-FB1/Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts.RData")
fitted_tr <- read_delim("~/BTSync/FetalRNAseq/Github/GENEX-FB1/Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/tables/MalevsFemale.complete.txt", 
                         "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(padj<.1) %>% 
  dplyr::select(Id, baseMean, Female, Male, FC, log2FoldChange, pvalue, padj, dispGeneEst, dispFit, dispMAP, dispersion, betaConv, maxCooks)

transcripts <- data.frame(Id=filter(fitted_tr, maxCooks<10 & maxCooks > 2)$Id)

excluded_tr=c()

# This scales VERY poorly. I'm not even going to attempt it for more than a few dozen genes
for (Id in transcripts$Id) {
  excluded_tr = c(excluded_tr, HighCooks(Id))
}
transcripts$excluded=excluded_tr

all_duplicated <- count(transcripts, excluded) %>% filter(n>1)
for (duplicated in all_duplicated$excluded) {
  print(duplicated)
  x1<- filter(transcripts, excluded==duplicated)
  transcripts <- filter(transcripts, excluded!=duplicated) %>% 
    bind_rows(data.frame(Id=paste(x1$Id, collapse='_'), excluded=duplicated))
}

#x <- mutate(x, excluded=paste0("'", excluded, "'"))
transcriptsList <- as.list(transcripts$Id)
names(transcriptsList) <- as.character(transcripts$excluded)
as.yaml(list(gene_level=myList, transcript_level=transcriptsList), indent=4) %>% write_file("config.yaml")


