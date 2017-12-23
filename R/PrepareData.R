library(tidyverse)
library(stringr)
library(biomaRt)
library(yaml)

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

fittedBias <- read_delim("Results/BGgenes.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
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

fittedBias_tr <- read_delim("Results/BGtranscripts.txt", "\t", escape_double = FALSE, trim_ws = TRUE) %>%
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

counts12_20_tr<-as.data.frame(counts12_20_tr)
rownames(counts12_20_tr)<-counts12_20_tr$Id

for (samples in names(outliers$transcript_level)) {
  for (SampleId in str_split(samples, '_')[[1]]) {
    for (Id in str_split(outliers$transcript_level[[samples]], '_')[[1]]) {
      print(c(Id, SampleId))
      counts12_20_tr[Id, paste('norm', SampleId, sep='.')] <- NA
    }
  }
}

write_tsv(counts12_20_tr, "Shiny/GENEX-FB1/Data/counts12_20_tr.txt")

file.copy("Data/SampleInfo.txt", "Shiny/GENEX-FB1/Data/SampleInfo.txt", overwrite=TRUE)

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
SRA <- dplyr::select(SampleInfo, library_ID=Sample, instrument_model=Sequencer) %>%
  mutate(instrument_model=ifelse(str_detect(instrument_model, '2500'), "Illumina HiSeq 2500", "Illumina HiSeq 4000")) %>%
  full_join(SRA)

SRA <- dplyr::select(BioSampleObjects, library_ID=`Sample Name`, biosample_accession=Accession) %>%
  full_join(SRA)

SRA <- full_join(SRA, sequences)

SRA$assembly<-NA

SRA <- dplyr::select(SRA, biosample_accession, bioproject_accession, title, library_ID, everything())

write_tsv(SRA, "Data/sra_metadata.txt", na='')
