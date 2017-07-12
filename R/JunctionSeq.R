library(JunctionSeq)

decoder <- read.table("JunctionSeqdecoder.bySample.txt",
                      header=TRUE,
                      stringsAsFactors=FALSE);
gff.file <- "JunctionSeq/withNovel.forJunctionSeq.gff.gz"
countFiles<- paste0("JunctionSeq/",
                    decoder$sample.ID,
                    "/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz")

jscs <- runJunctionSeqAnalyses(sample.files = countFiles,
                               sample.names = decoder$sample.ID,
                               condition=factor(decoder$group.ID),
                               flat.gff.file = gff.file,
                               nCores = 8,
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
