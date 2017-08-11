BASH and R scripts used for all steps of analyses, along with Shiny web app to view results.

To run mapping on all samples:
```
mkdir Mappings
mkdir FastQC
for sample in `cut -f2 Data/sequences.txt | sort | uniq`
do
    bash Bash/MappingPipeline.sh $sample
done
``` 

To analyse QC:
```
mkdir Tables
Rscript R/SummariseBamQC.R
```

To merge counts for samples sequenced on multiple lanes:
```
mkdir Counts
for BrainBankID in `cut -f 2 Data/sequences.txt | cut -d- -f 1 | sort |uniq`
do
    find Mappings/ -name $BrainBankID*.chr.counts.txt | xargs Rscript R/CombineCounts.R Counts/$BrainBankID.chr.counts.txt
done
```

To Run EdgeR:
```
mkdir Results
bash Bash/RunEdgeR.sh
```

To Run JunctionSeq:
```
mkdir JunctionSeq
find Mappings/ -name *.sort.bam | grep -v chr | sort | xargs -n 1 qsub Bash/QoRTs.sh
echo -e 'unique.ID\tsample.ID' > JunctionSeq/decoder.byUID.txt
ls Mappings | grep '-' | perl -pe 's/([^-]+)(.*)/$1$2\t$1/' >> JunctionSeq/decoder.byUID.txt 
qsub Bash/MergeQoRTs.sh JunctionSeq

echo -e 'sample.ID' > JunctionSeq/decoder.bySample.txt
cut -f 1 Data/SampleInfo.txt | tail -n +2 >> JunctionSeq/decoder.bySample.txt
qsub Bash/NovelSplices.sh JunctionSeq
```

To Run JunctionSeq on 16 PCW13 samples (using GTF filtered to >= 1 TPM in at least 56 samples)
```
cat /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf | python Python/FilterGTF.py 1 > Data/
genesFiltered1.gtf
mkdir JunctionSeqF1
for BrainBankID in `grep '\b13\b' Data/SampleInfo.txt |head -16 | cut -f 1`; do find Mappings/ -name $BrainBankID*.sort.bam | grep -v chr | sort | xargs -n 1 qsub Bash/QoRTs.sh; done 
```

To Run DESeq2 on Junctions
```
for BrainBankID in `cut -f 1 Data/SampleInfo.txt | tail -n +2`; do zcat JunctionSeq/$BrainBankID/QC.spliceJunctionAndExonCounts.withNovel.forJunctionSeq.txt.gz | grep -e ':N' -e ':J' > Counts/$BrainBankID.junctions.txt; done

```

To Run DEXSeq:
```
mkdir DEXSeq
find Mappings/ -name *.sort.bam | grep -v chr | sort | xargs -n 1 qsub Bash/dexseq-count.sh
```

To run Kallisto on all samples (combine all reads from each sample):

```
mkdir Kallisto
cut -f 1 Data/SampleInfo.txt | tail -n +2 | xargs -n 1 qsub Bash/Kallisto.sh 
``` 

To get transcript to gene mapping:
```
echo -e "transcript_id\tgene_id" > Data/tx2gene.txt 
cat /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf | grep transcript_id | perl -pe 's/.*(ENSGR?\d+).*(ENSTR?\d+).*/$2\t$1/' >> Data/tx2gene.txt 
```
