#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#



tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
      --transcriptome-index /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.inx \
      --output-dir $3 \
      /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/Bowtie2Index/genome $1 $2
