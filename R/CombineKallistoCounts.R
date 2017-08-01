library(readr)
library(dplyr)
library(stringr)
args = commandArgs(trailingOnly=TRUE)
print(args)
Counts <- data.frame('Feature'=c())
for (fileName in args[2:length(args)]) {
  input <- read_delim(fileName, 
             "\t", escape_double = FALSE, col_names = TRUE, 
             trim_ws = TRUE)
  input <- select(input, target_id, tpm)
  colnames(input) <- c('Id', str_replace(fileName, '.*/(.*).kallisto.txt', '\\1'))
  if (length(Counts) == 0) {
    Counts <- input
  } else {
    Counts <- full_join(Counts, input)
  }    
}
write_tsv(Counts, args[1] )
