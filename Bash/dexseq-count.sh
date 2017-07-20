#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=8G
#

# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
filename=${1##*/}
filename=${filename%%.bam}
sampleID=${filename%%.*}

if [ ! -f DEXSeq/$sampleID/$filename.dex_counts.txt ]
then
    echo "Getting counts for $filename"
    mkdir DEXSeq/$sampleID/
    dexseq_count.py -p yes -r name -s reverse -f bam \
       /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gff \
      $1 \
      DEXSeq/$sampleID/$filename.dex_counts.txt
    if [ $? -eq 0 ]
    then
        echo "Finished getting counts for $filename"
    else
        echo "Could not get counts for $filename"
        exit 1
    fi    
fi



