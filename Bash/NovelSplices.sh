#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=20G
#


DIRECTORY=$1
_JAVA_OPTIONS='-Xms256M -Xmx16G -XX:ParallelGCThreads=1'

java -jar ~/src/QoRTs.jar mergeNovelSplices  \
                --minCount 6 \
                --stranded \
                $DIRECTORY/ \
                $DIRECTORY/decoder.bySample.txt \
                /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf \
                $DIRECTORY/
