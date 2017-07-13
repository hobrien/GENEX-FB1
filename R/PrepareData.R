library(tidyverse)
library(stringr)

# Oh yes, thought of a catchy (and hopefully memorable) name for the dataset: 
# GENEX-FB (for GENe EXpression in the Fetal Brain). 
# This dataset would be GENEX-FB1: Sex biases. 
# The larger dataset for eQTL / TWAS analysis (I reckon we'll get up to 150) would be GENEX-FB2: 
# Genotypic effects. (If we grew the sample, we'd call it FB3 etc). 
# We could also use GENEX for the adult brain samples (E.g. GENEX-AC (adult caudate)!
setwd("~/BTSync/FetalRNAseq/Github/GENEX-FB1/")
gene_info <- read_tsv("Data/genes.txt") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)


counts12_20 <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq/tables/MalevsFemale.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(Id=str_extract(Id, '^[^.]+'))
right_join(gene_info, dplyr::select(counts12_20, Id, starts_with('norm'))) %>% 
  write_tsv("Shiny/GENEX-FB1/Data/counts12_20.txt")



fittedBias <- read_delim("Results/Sex_PCW_12_20_FDR_0.1_DESeq/tables/BG12_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
  mutate(ageBin='12-19') %>%
  bind_rows(
    read_delim("Results/Sex_12_FDR_0.1_DESeq/tables/BG12.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='12')
  ) %>%
  bind_rows(
    read_delim("Results/Sex_13_FDR_0.1_DESeq/tables/BG13.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='13')
  ) %>%
  bind_rows(
    read_delim("Results/Sex_14_FDR_0.1_DESeq/tables/BG14.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='14')
  ) %>%
  bind_rows(
    read_delim("Results/Sex_15_17_FDR_0.1_DESeq/tables/BG15_17.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='15-16')
  ) %>%
  bind_rows(
    read_delim("Results/Sex_17_20_FDR_0.1_DESeq/tables/BG17_20.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='17-19')
  )

right_join(gene_info, fittedBias) %>%
  write_tsv("Shiny/GENEX-FB1/Data/fitted.txt")

file.copy("Data/SampleInfo.txt", "Shiny/GENEX-FB1/Data/SampleInfo.txt", overwrite=TRUE)
