library(tidyverse)
library(sleuth)

sample_id <- dir("Kallisto")
kal_dirs <- tibble(sample=sample_id, path=file.path("Kallisto", sample_id)) 
LibraryInfo <- read_tsv("Data/SampleInfo.txt", 
                        col_types=cols(Sample='c')
) %>%
  rename(sample=Sample)
LibraryInfo <- full_join(LibraryInfo, kal_dirs)

so <- sleuth_prep(LibraryInfo, extra_bootstrap_summary = TRUE)

so <- sleuth_fit(so, ~Sex, 'full')
so <- sleuth_fit(so, ~1, 'full')
so <- sleuth_lrt(so, 'reduced', 'full')
models(so)
sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE) %>%
  write_tsv("Results/Sleuth.txt")

