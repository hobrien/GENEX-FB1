#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=12G
#

INFILE=$1
BASENAME=${INFILE##*/} 
OUTPUTDIR=JunctionSeqF1/${BASENAME%%.*}
_JAVA_OPTIONS='-Xms256M -Xmx8G -XX:ParallelGCThreads=1'

mkdir $OUTPUTDIR
java -jar ~/src/QoRTs.jar QC \
--stranded \
--minMAPQ 50 \
$INFILE \
Data/genesFiltered1.gtf \
$OUTPUTDIR
