library(tidyr)
library(dplyr)
library(stringr)
library(readr)

folders <- list.dirs("../Mappings", recursive=FALSE)
folders<-folders[folders!="../Mappings/17193-1"]
folders<-folders[folders!="../Mappings/17193-2"]


# bam_stat.py
RSeQCstats <- data.frame()
for (folder in folders) {
  #print(folder)
  folder=paste(folder, str_extract(folder, '[^/]+$'), sep='/')
  if(file.exists(paste(folder, ".ex.stats.txt", sep=""))){
    temp <- read.delim(paste(folder, ".ex.stats.txt", sep=""), 
                       header=FALSE, stringsAsFactors=FALSE, skip=5, sep=':'
    )
  } else {
    temp <- read.delim(paste(folder, ".chr.stats.txt", sep=""), 
                       header=FALSE, stringsAsFactors=FALSE, skip=5, sep=':'
    )
  }
  temp[4,2] <-strsplit(temp[4,1], ' +')[[1]][4]
  temp[4,1] <- 'Non primary hits'
  temp <- temp[c(6,7,14),]
  temp2 <- read.delim(paste(folder, ".in.stats.txt", sep=""), 
                      header=FALSE, stringsAsFactors=FALSE, skip=5, sep=':'
  )
  temp2[4,2] <-strsplit(temp[4,1], ' +')[[1]][4]
  temp2[4,1] <- 'Non primary hits'
  temp <-rbind(temp, c("rDNA",sum(as.numeric(temp2[c(6,7),2]))))
  temp[,1] <- c("Multimapped", "Unique", "Paired", "rDNA")
  temp$sample <- tail(strsplit(folder, '/')[[1]], 1)
  RSeQCstats <- rbind(RSeQCstats, temp)
}
RSeQCstats <- mutate(RSeQCstats, V2=as.numeric(V2))
RSeQCstats <- spread(RSeQCstats, V1, V2)
ReadNumbers <- RSeQCstats[,c(1,5)]
RSeQCstats <- RSeQCstats[c(1,5,3,2,4)]
RSeQCstats <- RSeQCstats %>% mutate(sample = str_extract(sample, "^[^-]+")) %>% group_by(sample) %>% summarise_each("sum")
write_tsv(RSeQCstats, "../Figures/read_numbers.txt")

# read_distribution.py
RSeQCdistribution <- data.frame()
for (folder in folders) {
  #print(folder)
  folder=paste(folder, str_extract(folder, '[^/]+$'), sep='/')
  sample <- tail(strsplit(folder, '/')[[1]], 1)
  temp <- read.delim(paste(folder, ".dist.txt", sep=""), 
                     header=TRUE, stringsAsFactors=FALSE, skip=4, sep=''
  )
  temp <- temp[-11,]
  tag_total <- sum(temp[c(1,2,3,4,7,10),3])
  temp <- data.frame(V1 = c("Total Tags", "CDS", "UTR", "Intron", "Intergenic"), 
                     V2=c(
                       tag_total,
                       temp[1,3]/tag_total, 
                       sum(as.numeric(temp[c(2,3),3]))/tag_total,
                       temp[4,3]/tag_total,
                       sum(as.numeric(temp[c(7,10),3]))/tag_total
                     )
  )
  temp$sample <- sample
  RSeQCdistribution <- rbind(RSeQCdistribution, temp)
}
RSeQCdistribution$V2 <- round(RSeQCdistribution$V2, 3)
RSeQCdistribution <- spread(RSeQCdistribution, V1, V2)
RSeQCdistribution <- RSeQCdistribution[c(1,5,3,4,2,6)]
RSeQCdistribution <- RSeQCdistribution %>% mutate(sample = str_extract(sample, "^[^-]+")) %>% group_by(sample) %>% summarise_each("sum")
write_tsv(RSeQCdistribution, "../Figures/read_distribution.txt")

#infer_experiment.py
RSeQCexpt <- data.frame()
for (folder in folders) {
  folder=paste(folder, str_extract(folder, '[^/]+$'), sep='/')
  temp <- read.delim(paste(folder, ".expt.txt", sep=""), 
                     skip=3, header=FALSE, sep=':'
  )
  temp$sample <- tail(strsplit(folder, '/')[[1]], 1)
  temp$V1 <- c("Ambiguous", "First Strand", "Second Strand")
  RSeQCexpt <- rbind(RSeQCexpt, temp)
}
RSeQCexpt <- spread(RSeQCexpt, V1, V2)
# This assumes the same number of reads from each lane, which isn't true but is probably close enough
RSeQCexpt <- RSeQCexpt %>% mutate(sample = str_extract(sample, "^[^-]+")) %>% group_by(sample) %>% summarise_each("mean")
write_tsv(RSeQCexpt, "../Figures/read_strand.txt")

#inner_distance.py
RSeQCdistance <- data.frame()
for (folder in folders) {
  folder=paste(folder, str_extract(folder, '[^/]+$'), sep='/')
  #print(folder)
  temp <- read.delim(paste(folder, ".inner_distance_freq.txt", sep=""),
                     header=FALSE
  )
  sample <- tail(strsplit(folder, '/')[[1]], 1)
  temp$sample <- sample
  RSeQCdistance <- rbind(RSeQCdistance, temp)
}
RSeQCdistance <- RSeQCdistance %>% mutate(sample = str_extract(sample, "^[^-]+"), 
                                          size = (V1+V2)/2) %>% 
  group_by(sample, size) %>% summarise(count=sum(V3))
write_tsv(RSeQCdistance, "../Figures/read_distance.txt")

#junction_saturation.py
RSeQCsat <- data.frame()
for (folder in folders) {
  folder=paste(folder, str_extract(folder, '[^/]+$'), sep='/')
  eval(parse(file = paste(folder, ".junctionSaturation_plot.r", sep=""))[2:5])
  temp<-rbind(
    data.frame(percent_reads=x, junctions=z, Category='All'),
    data.frame(percent_reads=x, junctions=y, Category='Known'),
    data.frame(percent_reads=x, junctions=w, Category='Novel')
  )
  sample <- tail(strsplit(folder, '/')[[1]], 1)
  temp$sample <- sample
  RSeQCsat <- rbind(RSeQCsat, temp)
}
write_tsv(RSeQCsat, "../Figures/junction_sat.txt")
