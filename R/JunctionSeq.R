# set a default repository for packages                         
local({r <- getOption('repos')                                                                                                                                
      r['CRAN'] <- 'http://www.stats.bris.ac.uk/R/'                                                                                                          
      options(repos=r)                    
}) 

if(!require(JunctionSeq)) {
    source("https://bioconductor.org/biocLite.R")
    biocLite("JunctionSeq", ask=FALSE, suppressUpdates=FALSE)
    install.packages("JunctionSeq", dependencies=T)
}

library(JunctionSeq)

decoder <- read.table("Data/SampleInfo.txt",
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- "JunctionSeqF1/withNovel.forJunctionSeq.gff.gz"
countFiles<- paste0("JunctionSeqF1/",
                    decoder$Sample,
                    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$Sample,
                               condition=factor(decoder$Sex),
                               flat.gff.file = gff.file,
                               nCores = 5,
                               analysis.type = "junctionsAndExons"
)

buildAllPlots(jscs=jscs,
              outfile.prefix = "JunctionSeq/plotsJCT-2/",
              use.plotting.device = "png",
              FDR.threshold = 0.01)

buildAllPlots(jscs=jscs,
              outfile.prefix = "JunctionSeq/plotsJCT-3/",
              use.plotting.device = "png",
              FDR.threshold = 0.001)
save.image(file="JunctionSeq/junctionSeq.RData")
