library(tidyverse)

# Oh yes, thought of a catchy (and hopefully memorable) name for the dataset: 
# GENEX-FB (for GENe EXpression in the Fetal Brain). 
# This dataset would be GENEX-FB1: Sex biases. 
# The larger dataset for eQTL / TWAS analysis (I reckon we'll get up to 150) would be GENEX-FB2: 
# Genotypic effects. (If we grew the sample, we'd call it FB3 etc). 
# We could also use GENEX for the adult brain samples (E.g. GENEX-AC (adult caudate)!
gene_info <- read_tsv("../Data/genes.csv") %>%
  mutate(gene_id = sub("\\.[0-9]+", "", gene_id)) %>%
  dplyr::select(Id = gene_id, SYMBOL=gene_name, Chr=seqid)

counts <- read_delim("../Results/MvsF12_19_FDR_0.1/BG12_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, starts_with('norm'))

right_join(gene_info, counts) %>% 
  write_tsv("../Shiny/GENEX-FB1/Data/counts.txt")

fittedBias <- read_delim("../Results/MvsF12_19_FDR_0.1/BG12_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
  dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
  mutate(ageBin='12-19') %>%
  bind_rows(
    read_delim("../Results/MvsF12_FDR_0.1/BG12.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='12')
  ) %>%
  bind_rows(
    read_delim("../Results/MvsF13_FDR_0.1/BG13.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='13')
  ) %>%
  bind_rows(
    read_delim("../Results/MvsF14_FDR_0.1/BG14.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='14')
  ) %>%
  bind_rows(
    read_delim("../Results/MvsF15_16_FDR_0.1/BG15_16.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='15-16')
  ) %>%
  bind_rows(
    read_delim("../Results/MvsF17_19_FDR_0.1/BG17_19.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
      dplyr::select(Id, Male, Female, log2FoldChange, pvalue, padj) %>% 
      mutate(ageBin='17-19')
  )

right_join(gene_info, fittedBias) %>% 
  write_tsv("../Shiny/GENEX-FB1/Data/fitted.txt")
