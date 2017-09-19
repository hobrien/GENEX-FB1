library(tidyverse)
library(stringr)
library(biomaRt)
library(yaml)

# Oh yes, thought of a catchy (and hopefully memorable) name for the dataset: 
# GENEX-FB (for GENe EXpression in the Fetal Brain). 
# This dataset would be GENEX-FB1: Sex biases. 
# The larger dataset for eQTL / TWAS analysis (I reckon we'll get up to 150) would be GENEX-FB2: 
# Genotypic effects. (If we grew the sample, we'd call it FB3 etc). 
# We could also use GENEX for the adult brain samples (E.g. GENEX-AC (adult caudate)!
# setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/")

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "chromosome_name", "gene_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, Id = ensembl_transcript_id,
                     GeneId = ensembl_gene_id, SYMBOL = external_gene_name, Chr=chromosome_name, gene_type=gene_biotype)
t2g <- dplyr::mutate(t2g, Chr = paste0('chr', Chr), ChrType = ifelse(Chr == 'chrX' | Chr == 'chrY', Chr, 'autosomal')
  )

gene_info <- t2g %>% 
  dplyr::select(-Id) %>% 
  group_by(GeneId) %>% 
  dplyr::slice(1) %>%
  ungroup() %>%
  dplyr::rename(Id=GeneId)

dplyr::rename(gene_info, gene_id=Id, gene_name=SYMBOL, seqid=Chr) %>% write_tsv("Data/genes.txt")
# Results of gene level analyses
counts12_20 <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Id=str_extract(Id, '^[^.]+'))
counts12_20 <-right_join(gene_info, dplyr::select(counts12_20, Id, starts_with('norm')))

counts12_20<-as.data.frame(counts12_20)
rownames(counts12_20)<-counts12_20$Id

outliers <- yaml.load_file("config.yaml")
Ids <-c()
SampleIds<-c()
for (samples in names(outliers$gene_level)) {
  for (SampleId in str_split(samples, '_')[[1]]) {
    for (Id in str_split(outliers$gene_level[[samples]], '_')[[1]]) {
      counts12_20[Id, paste('norm', SampleId, sep='.')] <- NA
    }
  }
}

write_tsv(counts12_20, "Shiny/GENEX-FB1/Data/counts12_20.txt")

fittedBias <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/BG12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj, maxCooks) %>%
  mutate(pvalue= ifelse(maxCooks>10, 1, pvalue), padj= ifelse(maxCooks>10, 1, padj))
fittedBias <- right_join(gene_info, fittedBias)
write_tsv(fittedBias, "Shiny/GENEX-FB1/Data/fitted.txt")

fittedPCW <- read_tsv("Results/PCW_Sex_12_20_FDR_0.1_DESeqLRT_kallistoCounts/tables/BG12_20.txt") %>%
  dplyr::select(Id, baseMean, log2FoldChange, pvalue, padj, maxCooks) %>%
  mutate(pvalue= ifelse(maxCooks>10, 1, pvalue), padj= ifelse(maxCooks>10, 1, padj))
fittedPCW <- right_join(gene_info, fittedPCW)
write_tsv(fittedPCW, "Shiny/GENEX-FB1/Data/dropPCW.txt")

# Results of transcript level analyses

fittedBias_tr <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/tables/BG12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj, maxCooks) %>%
  mutate(pvalue= ifelse(maxCooks>10, 1, pvalue), padj= ifelse(maxCooks>10, 1, padj)) 
fittedBias_tr <- right_join(t2g, fittedBias_tr) 
write_tsv(fittedBias_tr, "Shiny/GENEX-FB1/Data/fitted_tr.txt")

fittedPCW_tr <- read_tsv("Results/PCW_Sex_12_20_FDR_0.1_DESeqLRT_transcripts_kallistoCounts/tables/BG12_20.txt") %>%
  dplyr::select(Id, baseMean, log2FoldChange, pvalue, padj, maxCooks) %>%
  mutate(pvalue= ifelse(maxCooks>10, 1, pvalue), padj= ifelse(maxCooks>10, 1, padj))
fittedPCW_tr <- right_join(t2g, fittedPCW_tr) 
write_tsv(fittedPCW_tr, "Shiny/GENEX-FB1/Data/dropPCW_tr.txt")

counts12_20_tr <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/tables/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Id=str_extract(Id, '^[^.]+'))
counts12_20_tr <- right_join(t2g, dplyr::select(counts12_20_tr, Id, starts_with('norm')))

counts12_20<-as.data.frame(counts12_20)
rownames(counts12_20)<-counts12_20$Id

for (samples in names(outliers$transcript_level)) {
  for (SampleId in str_split(samples, '_')[[1]]) {
    for (Id in str_split(outliers$transcript_level[[samples]], '_')[[1]]) {
      counts12_20_tr[Id, paste('norm', SampleId, sep='.')] <- NA
    }
  }
}

write_tsv(counts12_20_tr, "Shiny/GENEX-FB1/Data/counts12_20_tr.txt")


file.copy("Data/SampleInfo.txt", "Shiny/GENEX-FB1/Data/SampleInfo.txt", overwrite=TRUE)

#SampleInfo <- read_tsv("Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 
