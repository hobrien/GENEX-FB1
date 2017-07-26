#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=12G
#

filename=${1##*/}
filename=${filename%%.bam}
sampleID=${filename%%.*}

export PATH=$PATH:/share/apps/R-3.2.2/bin:/share/apps/

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
    echo "Running Kallisto on $SampleID"
    kallisto quant -i $REFDIR/Annotation/Genes.gencode/kallisto.inx \
      -o Kallisto/$SampleID \
      --bias \
      -b 100 \
      --rf-stranded \
      $@
    if [ $? -eq 0 ]
    then
        echo "Finished running Kallisto on $SampleID"
    else
        echo "Could not run Kallisto on $SampleID"
        exit 1
    fi    
fi
