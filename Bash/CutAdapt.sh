#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
outfile1=$(basename $1 .fastq.gz)_trimmed.fastq.gz
outfile2=$(basename $2 .fastq.gz)_trimmed.fastq.gz

cutadapt \
    -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
    -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
    -o $outfile1 -p $outfile2 $@

