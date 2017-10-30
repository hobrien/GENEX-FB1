#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=6G
#

export PATH=$PATH:/share/apps/

BrainBankID=$1

echo "Starting Kallisto pipeline for $BrainBankID"

seq_folders=$(grep $BrainBankID Data/sequences.txt | cut -f 5 | sort | uniq)
sequences=$(for name in `grep $BrainBankID Data/sequences.txt | cut -f 1`; do find $seq_folders -name $name*f*q.gz -o -name $name*f*q; done)
set -- $sequences


REFDIR=/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy

if [ ! -f $REFDIR/Annotation/Genes.gencode/genes.fa ]
then
    echo "Extracting transcript sequences"
    cut -f 1-12 \
      $REFDIR/Annotation/Genes.gencode/genes.bed > /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed12

    bedtools getfasta \
      -fi $REFDIR/Sequence/WholeGenomeFasta/genome.fa \
      -bed $REFDIR/Annotation/Genes.gencode/genes.bed12 \
      -split -name -s \
      | fold -w 60 \
      | perl -pe 's/\([+-]\)//' \
      > $REFDIR/Annotation/Genes.gencode/genes.fa
    if [ $? -eq 0 ]
    then
        echo "Finished extracting transcript sequences"
    else
        echo "Could not extract transcript sequences"
        exit 1
    fi    
fi

if [ ! -f $REFDIR/Annotation/Genes.gencode/kallisto.inx ]
then
    echo "Creating kallisto index"
    kallisto index -i $REFDIR/Annotation/Genes.gencode/kallisto.inx $REFDIR/Annotation/Genes.gencode/genes.fa
    if [ $? -eq 0 ]
    then
        echo "Finished creating kallisto index"
    else
        echo "Could not create kallisto index"
        exit 1
    fi    
fi

if [ ! -f Kallisto/$SampleID/abundances.h5 ] | [ ! -f Kallisto/$SampleID/abundances.tsv ] | [ ! -f Kallisto/$SampleID/run_info.json ]
then
    echo "Running Kallisto on $sequences"
    kallisto quant -i $REFDIR/Annotation/Genes.gencode/kallisto.inx \
      -o Kallisto/$BrainBankID \
      --bias \
      -b 100 \
      --rf-stranded \
      $sequences
    if [ $? -eq 0 ]
    then
        echo "Finished running Kallisto on $sequences"
    else
        echo "Could not run Kallisto on $sequences"
        exit 1
    fi    
fi

echo "Finished Kallisto pipeline for $BrainBankID"

