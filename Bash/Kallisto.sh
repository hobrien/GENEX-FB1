#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

REFDIR=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy

if [ ! -f $REFDIR//Annotation/Genes.gencode/genes.fa ]
then
    echo "Creating directory extracting transcript sequences"
    cut -f 1-12 \
      $REFDIR/Annotation/Genes.gencode/genes.bed > /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed12

    bedtools getfasta \
      -fi $REFDIR/Sequence/WholeGenomeFasta/genome.fa \
      -bed $REFDIR/Annotation/Genes.gencode/genes.bed12 \
      -split -name -s \
      | fold -w 60 \
      | perl -pe 's/\([+-]\)// \
      > $REFDIR//Annotation/Genes.gencode/genes.fa
    if [ $? -eq 0 ]
    then
        echo "Finished extracting transcript sequences"
    else
        echo "Could not extract transcript sequences"
        exit 1
    fi    
fi

kallisto index -i $REFDIR/Annotation/Genes.gencode/genes.fa