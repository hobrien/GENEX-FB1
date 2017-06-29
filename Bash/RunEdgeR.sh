#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash


#Rscript FvsMedgeR.r --min 12 --max 20 -a
Rscript R/FvsMedgeR.r --min 12 --max 13
Rscript R/FvsMedgeR.r --min 13 --max 14
Rscript R/FvsMedgeR.r --min 14 --max 15
Rscript R/FvsMedgeR.r --min 15 --max 17 -b Sequencer
Rscript R/FvsMedgeR.r --min 17 --max 20 -b Sequencer

# Rscript FvsMedgeR.r --min 12 --max 20 -a -r 5
# Rscript FvsMedgeR.r --min 12 --max 13 -r 5
# #Rscript FvsMedgeR.r --min 13 --max 14 -r 5
# Rscript FvsMedgeR.r --min 14 --max 15 -r 5
# Rscript FvsMedgeR.r --min 15 --max 17 -r 5
# Rscript FvsMedgeR.r --min 17 --max 20 -r 5

