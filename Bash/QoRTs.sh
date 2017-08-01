#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#

INFILE=$1
BASENAME=${INFILE##*/} 
OUTPUTDIR=JunctionSeqF1/${BASENAME%%.*}
_JAVA_OPTIONS='-Xms256M -Xmx16G -XX:ParallelGCThreads=1'

mkdir $OUTPUTDIR
java -jar ~/src/QoRTs.jar QC \
--stranded \
--minMAPQ 50 \
$INFILE \
Data/genesFiltered1.gtf \
$OUTPUTDIR
