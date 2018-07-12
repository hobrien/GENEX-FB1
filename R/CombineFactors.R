library(readr)
library(dplyr)
args <- commandArgs(trailingOnly = TRUE)


factors <- read_tsv(args[1], trim_ws = TRUE, col_types = cols(ID='c'))

SampleInfo <- read_tsv(args[2], trim_ws = TRUE, col_types = cols(Sample='c'))

SampleInfo <- full_join(SampleInfo, select(factors, Sample=ID, V5:V10))

write_tsv(SampleInfo, args[3])
