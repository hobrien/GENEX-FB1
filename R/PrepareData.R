library(tidyverse)
library(stringr)


#setwd("../")
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
counts <- read_delim("Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts/tables/sex2vssex1.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

fittedBias <- gene_info %>% right_join(select(counts, Id, Female=sex1, Male=sex2, log2FoldDiff=log2FoldChange, pvalue, padj, maxCooks)) %>%
  filter(!is.na(padj) & !is.na(Chr) & maxCooks < 1) %>%
  dplyr::select(-gene_type, -maxCooks)

counts <- counts %>% select(Id, starts_with('norm.'))
colnames(counts) <- str_replace(colnames(counts), 'norm.', '')
 
write_tsv(counts, "Shiny/GENEX-FB1/Data/counts.txt")
write_tsv(fittedBias, "Shiny/GENEX-FB1/Data/fitted.txt")

# Results of transcript level analyses
counts_tr <- read_delim("Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts/tables/sex2vssex1.complete.txt", "\t", escape_double = FALSE, trim_ws = TRUE)

fittedBias_tr <- t2g %>%
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


#Prepare Biosample submission for SRA:
SampleInfo <- read_tsv("Data/SampleInfo.txt", col_types = cols(Sample='c'))
NumSamples <- nrow(SampleInfo)
Biosamples <- data.frame(sample_name=SampleInfo$Sample, 
                        sample_title=rep(NA, NumSamples), 
                        bioproject_accession=rep('PRJNA417945', NumSamples), 
                        organism=rep("Homo sapiens", NumSamples), 
                        isolate=SampleInfo$Sample, 
                        age=paste(SampleInfo$PCW, "weeks post-conception"), 
                        biomaterial_provider=rep('Human Developmental Biology Resource (HDBR), 30 Guilford St, London WC1N 1EH', NumSamples), 
                        sex=SampleInfo$Sex, 
                        tissue=rep('Fetal Brain', NumSamples), 
                        cell_line=rep(NA, NumSamples), 
                        cell_subtype=rep(NA, NumSamples), 
                        cell_type=rep(NA, NumSamples), 
                        culture_collection=rep(NA, NumSamples), 
                        dev_stage=rep("Fetal", NumSamples), 
                        disease=rep(NA, NumSamples), 
                        disease_stage=rep(NA, NumSamples), 
                        ethnicity=rep(NA, NumSamples), 
                        health_state=rep(NA, NumSamples), 
                        karyotype=rep(NA, NumSamples), 
                        phenotype=rep(NA, NumSamples), 
                        population=rep(NA, NumSamples), 
                        race=rep(NA, NumSamples), 
                        sample_type=rep(NA, NumSamples), 
                        treatment=rep(NA, NumSamples), 
                        description=rep(NA, NumSamples)
)

write_tsv(Biosamples, "Data/Human.1.0.tsv", col_names=FALSE, append=TRUE)
BioSampleObjects<- read_tsv("Data/BioSampleObjects.txt", col_types = cols(`Sample Name`='c', Isolate='c'))

sequences <- read_delim("Data/sequences.txt",
"\t", col_names=FALSE, escape_double = FALSE, trim_ws = TRUE)
sequences <- sequences %>% mutate(library_ID=str_replace(X2, '-\\d', '')) %>%
  dplyr::select(X1, library_ID) %>%
  group_by(library_ID) %>%
  mutate(read=paste0('filename', row_number())) %>% 
  ungroup() %>% 
  spread(read, X1) %>%
  dplyr::select(library_ID, filename=filename1, one_of(paste0('filename', seq(2, 12))))

SRA <- data.frame(bioproject_accession=rep('PRJNA417945', NumSamples),
                  title=rep('Human fetal brain RNA sequencing', NumSamples),
                  library_ID=SampleInfo$Sample,
                  design_description=rep('RNA was extracted from frozen human fetal brain homogenate followed by QIAGEN RNeasy MinElute RNA cleanup & ribosomal RNA depletion using Ribo-Zero and Illumina TruSeq Stranded Total RNA library preparation', NumSamples),
                  library_strategy=rep('RNA-Seq', NumSamples),
                  library_source=rep('TRANSCRIPTOMIC', NumSamples),
                  library_selection=rep('RANDOM', NumSamples),
                  library_layout=rep('PAIRED', NumSamples),
                  platform=rep('ILLUMINA', NumSamples),
                  filetype=rep('fastq', NumSamples)
)
SRA <- SampleInfo %>% mutate(instrument_model=ifelse(str_detect(ReadLength, '2x76bp'), "Illumina HiSeq 4000", "Illumina HiSeq 2500")) %>%
  dplyr::select(library_ID=Sample, instrument_model) %>%
  full_join(SRA)

SRA <- dplyr::select(BioSampleObjects, library_ID=`Sample Name`, biosample_accession=Accession) %>%
  full_join(SRA)

SRA <- full_join(SRA, sequences)

SRA$assembly<-NA

SRA <- dplyr::select(SRA, biosample_accession, bioproject_accession, title, library_ID, everything())

write_tsv(SRA, "Data/sra_metadata.txt", na='')
