library(tidyverse)
library(stringr)


# Results of gene level analyses
counts <- read_delim("Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts/tables/sex2vssex1.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

fittedBias <- select(gene_info, Id=gene_id, SYMBOL=gene_name, Chr=seqid, gene_type, ChrType) %>% right_join(select(counts, Id, Female=sex1, Male=sex2, log2FoldDiff=log2FoldChange, pvalue, padj, maxCooks)) %>%
  filter(!is.na(padj) & !is.na(Chr) & maxCooks < 1) %>%
  dplyr::select(-gene_type, -maxCooks)

counts <- counts %>% select(Id, starts_with('norm.'))
colnames(counts) <- str_replace(colnames(counts), 'norm.', '')
 
write_tsv(counts, "Shiny/GENEX-FB1/Data/counts.txt")
write_tsv(fittedBias, "Shiny/GENEX-FB1/Data/fitted.txt")

# Results of transcript level analyses
counts_tr <- read_delim("Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts/tables/sex2vssex1.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

fittedBias_tr <- select(gene_info, gene_id, SYMBOL=gene_name, Chr=seqid, gene_type, ChrType) %>%
  left_join(tx2gene) %>%
  rename(Id=transcript_id, GeneId=gene_id) %>%
  right_join(select(counts_tr, Id, Female=sex1, Male=sex2, log2FoldDiff=log2FoldChange, pvalue, padj, maxCooks)) %>%
  filter(!is.na(padj) & !is.na(Chr) & maxCooks < 1) %>%
  dplyr::select(-gene_type, -maxCooks)

counts_tr <- counts_tr %>% select(Id, starts_with('norm.'))
colnames(counts_tr) <- str_replace(colnames(counts_tr), 'norm.', '')

write_tsv(fittedBias_tr, "Shiny/GENEX-FB1/Data/fitted_tr.txt")
write_tsv(counts_tr, "Shiny/GENEX-FB1/Data/counts_tr.txt")

SampleInfo <- read_tsv("Data/SampleInfo.txt", col_types = cols(Sample = col_character()))
SampleInfo <- SampleInfo[SampleInfo$Sample %in% colnames(counts),]

write_tsv(SampleInfo, "Shiny/GENEX-FB1/Data/SampleInfo.txt")

#SampleInfo <- read_tsv("Data/SampleInfo.txt", trim_ws = TRUE, col_names=TRUE, cols(Sample='c')) 
