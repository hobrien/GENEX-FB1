#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -pe smp 8
#$ -l h_vmem=8G
#

Rscript R/DEXSeq.R
