#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

#cd `dirname $0`/../R # change to directory containing R scripts

#Rscript FvsMedgeR.r --min 12 --max 20 -a
Rscript FvsMedgeR.r --min 12 --max 13
Rscript FvsMedgeR.r --min 13 --max 14
Rscript FvsMedgeR.r --min 14 --max 15
Rscript FvsMedgeR.r --min 15 --max 17 -b Sequencer
Rscript FvsMedgeR.r --min 17 --max 20 -b Sequencer

# Rscript FvsMedgeR.r --min 12 --max 20 -a -r 5
# Rscript FvsMedgeR.r --min 12 --max 13 -r 5
# #Rscript FvsMedgeR.r --min 13 --max 14 -r 5
# Rscript FvsMedgeR.r --min 14 --max 15 -r 5
# Rscript FvsMedgeR.r --min 15 --max 17 -r 5
# Rscript FvsMedgeR.r --min 17 --max 20 -r 5

