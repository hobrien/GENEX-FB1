#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#


REFDIR=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy

tophat --keep-fasta-order --library-type fr-secondstrand --mate-inner-dist 500 --mate-std-dev 50 --num-threads 8 \
      --transcriptome-index $REFDIR/Annotation/Genes.gencode/genes.inx \
      --output-dir $3 \
      $REFDIR/Sequence/Bowtie2Index/genome $1 $2
