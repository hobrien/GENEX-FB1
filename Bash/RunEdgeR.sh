#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash


Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW
Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --exclude 17046,18655,11971
Rscript R/FvsMedgeR.R --min 12 --max 13
Rscript R/FvsMedgeR.R --min 12 --max 13 --exclude 17046
Rscript R/FvsMedgeR.R --min 13 --max 14
Rscript R/FvsMedgeR.R --min 13 --max 14 --exclude 18655  
Rscript R/FvsMedgeR.R --min 14 --max 15
Rscript R/FvsMedgeR.R --min 15 --max 17
Rscript R/FvsMedgeR.R --min 17 --max 20
Rscript R/FvsMedgeR.R --min 17 --max 20 --exclude 11971

# Rscript R/FvsMedgeR.R --min 12 --max 20 -cofactor PCW -r 5
# Rscript R/FvsMedgeR.R --min 12 --max 13 -r 5
# Rscript R/FvsMedgeR.R --min 13 --max 14 -r 5
# Rscript R/FvsMedgeR.R --min 14 --max 15 -r 5
# Rscript R/FvsMedgeR.R --min 15 --max 17 -r 5
# Rscript R/FvsMedgeR.R --min 17 --max 20 -r 5

Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq
Rscript R/FvsMedgeR.R --min 12 --max 13 --tool DESeq
Rscript R/FvsMedgeR.R --min 13 --max 14 --tool DESeq
Rscript R/FvsMedgeR.R --min 14 --max 15 --tool DESeq
Rscript R/FvsMedgeR.R --min 15 --max 17 --tool DESeq
Rscript R/FvsMedgeR.R --min 17 --max 20 --tool DESeq

Rscript R/FvsMedgeR.R --min 12 --max 20 --interaction PCW --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --varInt PCW --interaction Sex --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --varInt PCW --cofactor Sex --tool DESeqLRT
