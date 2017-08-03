#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -l h_vmem=10G
#

Rscript R/Sleuth.R
