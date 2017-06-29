BASH and R scripts used for all steps of analyses, along with Shiny web app to view results.

To run mapping on all samples:
```
mkdir Mappings
mkdir FastQC
for sample in `cut -f2 /c8000xd3/rnaseq-heath/GENEX-FB1/Data/sequences.txt | sort | uniq | grep 12116`
do
    bash /c8000xd3/rnaseq-heath/GENEX-FB1/Bash/MappingPipeline.sh $sample
done
``` 

To merge counts for samples sequenced on multiple lanes:
```
mkdir Counts
for BrainBankID in `cut -f 2 Data/sequences.txt | cut -d- -f 1 | sort |uniq`
do
    find Mappings/ -name $BrainBankID*.chr.counts.txt | xargs Rscript R/CombineCounts.R Counts/$BrainBankID.chr.counts.txt
done
```

```
mkdir Tables
Rscript R/SummariseBamQC.R
```
