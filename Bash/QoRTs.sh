#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

INFILE=$1
OUTPUTDIR=$2
java -jar ~/src/QoRTs.jar QC \
--stranded \
--minMAPQ 50 \
$INFILE \
/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
$OUTPUTDIR
