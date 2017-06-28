#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#



echo "Getting counts for $1"
  
htseq-count -f bam -s reverse -t exon -i gene_id -m intersection-strict \
  $1 \
   /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
  > ${1%%.bam}.counts.txt
echo "Finished getting counts for $1"

