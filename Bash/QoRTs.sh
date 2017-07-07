#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=12G
#

INFILE=$1
BASENAME=${INFILE##*/} 
OUTPUTDIR=JunctionSeq/${BASENAME%%.*}
_JAVA_OPTIONS='-Xms256M -Xmx8G -XX:ParallelGCThreads=1'

mkdir $OUTPUTDIR
java -jar ~/src/QoRTs.jar QC \
--stranded \
--minMAPQ 50 \
$INFILE \
/c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
$OUTPUTDIR
