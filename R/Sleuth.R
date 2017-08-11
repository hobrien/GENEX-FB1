library(tidyverse)
library(sleuth)

mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = 'ensembl.org')
t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                     "external_gene_name", "chromosome_name", "gene_biotype"), mart = mart)
t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                     Id = ensembl_gene_id, SYMBOL = external_gene_name, Chr=chromosome_name, gene_type=gene_biotype)


sample_id <- dir("Kallisto")
kal_dirs <- tibble(sample=sample_id, path=file.path("Kallisto", sample_id)) 
LibraryInfo <- read_tsv("Data/SampleInfo.txt", 
                        col_types=cols(Sample='c')
) %>%
  rename(sample=Sample)
LibraryInfo <- right_join(LibraryInfo, kal_dirs)

so <- sleuth_prep(LibraryInfo, num_cores=1, target_mapping = t2g, aggregation_column = 'Id')

so <- sleuth_fit(so, ~Sex + PCW + RIN, 'full')
so <- sleuth_fit(so, ~PCW + RIN, 'reduced')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)
sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) %>%
  write_tsv("Results/SleuthGenePCW_RIN.txt")

