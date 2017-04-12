#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#



echo "Getting counts for $1"
  
htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
  $1 \
  GRCh38/NCBI/GRCh38Decoy/Annotation/genes.gtf \
  > ${1%%.bam}_counts.txt
echo "Finished getting counts for $1"

