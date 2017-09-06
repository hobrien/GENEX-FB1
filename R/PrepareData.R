library(tidyverse)
library(stringr)
library(biomaRt)

# Oh yes, thought of a catchy (and hopefully memorable) name for the dataset: 
# GENEX-FB (for GENe EXpression in the Fetal Brain). 
# This dataset would be GENEX-FB1: Sex biases. 
# The larger dataset for eQTL / TWAS analysis (I reckon we'll get up to 150) would be GENEX-FB2: 
# Genotypic effects. (If we grew the sample, we'd call it FB3 etc). 
# We could also use GENEX for the adult brain samples (E.g. GENEX-AC (adult caudate)!
# setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/")
gene_info <- read_tsv("Data/genes.txt", col_types = cols(
  seqid = col_character(),
  source = col_character(),
  feature = col_character(),
  start = col_integer(),
  end = col_integer(),
  score = col_character(),
  strand = col_integer(),
  frame = col_character(),
  gene_id = col_character(),
  gene_type = col_character(),
  gene_status = col_character(),
  gene_name = col_character()
)) %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id),
         ChrType = ifelse(seqid == 'chrX' | seqid == 'chrY', seqid, 'autosomal')
  ) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid, ChrType)


counts12_20 <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Id=str_extract(Id, '^[^.]+'))
right_join(gene_info, dplyr::select(counts12_20, Id, starts_with('norm'))) %>% 
  write_tsv("Shiny/GENEX-FB1/Data/counts12_20.txt")

fittedBias <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/tables/BG12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
  mutate(ageBin='12-19')

right_join(gene_info, fittedBias) %>% 
  write_tsv("Shiny/GENEX-FB1/Data/fitted.txt")

fittedPCW <- read_tsv("Results/PCW_Sex_12_20_FDR_0.1_DESeqLRT_kallistoCounts/tables/dropPCW.complete.txt") %>%
  dplyr::select(Id, baseMean, log2FoldChange, pvalue, padj) %>%
  mutate(Id = sub("\\.[0-9]+", "", Id))

right_join(gene_info, fittedPCW) %>% 
  write_tsv("Shiny/GENEX-FB1/Data/dropPCW.complete.txt")

file.copy("Data/SampleInfo.txt", "Shiny/GENEX-FB1/Data/SampleInfo.txt", overwrite=TRUE)

SampleInfo <- read_tsv("Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 


mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "chromosome_name", "gene_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     Id = ensembl_gene_id, SYMBOL = external_gene_name, Chr=chromosome_name, gene_type=gene_biotype)
Sleuth <- read_delim("Results/SleuthPCW_RIN.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
Sleuth <- right_join(t2g, mutate(Sleuth, target_id=str_replace(target_id, '\\.[0-9]+', '')))
write_tsv(Sleuth, "Results/SleuthPCW_RIN_annotated.txt")

AllKallisto <- read_delim("Counts/AllKallisto.txt", 
                     "\t", escape_double = FALSE, trim_ws = TRUE)
AllKallisto <- right_join(t2g, mutate(Sleuth, target_id=str_replace(target_id, '\\.[0-9]+', '')))
write_tsv(AllKallisto, "Counts/AllKallisto_annotated.txt")
