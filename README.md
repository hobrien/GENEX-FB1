BASH and R scripts used for all steps of analyses, along with Shiny web app to view results.

To run mapping on all samples:
```
mkdir Mappings
mkdir FastQC
for sample in `cut -f2 Data/sequences.txt | sort | uniq | grep 12116`
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
cut -f 2 Data/sequences.txt | grep '-' | sort | uniq | perl -pe 's/([^-]+)(.*)/$1$2\t$1/' >> JunctionSeq/decoder.byUID.txt 
qsub Bash/MergeQoRTs.sh JunctionSeq
```
