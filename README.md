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

To Run Differential Expression analysis:
```
mkdir Results
bash Bash/RunDEanalysis.sh
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

To run Differential Expression:
```
Bash Bash/RunDEanalyses.sh
```

To rerun Differential Expression analyses with outliers removed:
```
Rscript R/FindOutliers.R
source activate py35
snakemake --cluster "qsub -l h_vmem={params.maxvmem}" -j 50
```

