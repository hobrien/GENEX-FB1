#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

Rscript R/DEXSeq.R
