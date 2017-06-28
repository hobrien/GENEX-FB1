#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#
echo "running fastQC on $@"
~/src/FastQC/fastqc --outdir=/c8000xd3/rnaseq-heath/GENEX-FB1/FastQC $@
echo "finished running fastQC on $@"
