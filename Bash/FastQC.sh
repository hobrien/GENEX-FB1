#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

BASEDIR=/c8000xd3/rnaseq-heath/GENEX-FB1

echo "running fastQC on $@"
~/src/FastQC/fastqc --outdir=$BASEDIR/FastQC $@
echo "finished running fastQC on $@"
