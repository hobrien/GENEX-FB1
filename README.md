BASH and R scripts used for all steps of analyses, along with Shiny web app to view results.

To run mapping on all samples:
```
for sample in `cut -f2 Data/sequences.txt`
do
    bash Bash/MappingPipeline.sh $sample
done
``` 

To merge counts for samples sequenced on multiple lanes:
```
for BrainBankID in `cut -f 2 Data/sequences.txt | grep '-' | cut -d- -f 1 | sort |uniq`
do
    find /c8000xd3/rnaseq-heath/Mappings/ -name $BrainBankID*.chr.counts.txt | xargs Rscript R/CombineCounts.R Counts/$SampleID.chr.counts.txt
done
```
