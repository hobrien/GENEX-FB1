#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash



Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq
Rscript R/FvsMedgeR.R --min 12 --max 13 --tool DESeq
Rscript R/FvsMedgeR.R --min 13 --max 14 --tool DESeq
Rscript R/FvsMedgeR.R --min 14 --max 15 --tool DESeq
Rscript R/FvsMedgeR.R --min 15 --max 17 --tool DESeq
Rscript R/FvsMedgeR.R --min 17 --max 20 --tool DESeq

# Use LRT rather than Wald test
Rscript R/FvsMedgeR.R --min 12 --max 20 --interaction PCW --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --varInt PCW --interaction Sex --tool DESeqLRT
Rscript R/FvsMedgeR.R --min 12 --max 20 --varInt PCW --cofactor Sex --tool DESeqLRT

# Exclude samples with anomalous depth 
Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --exclude 17046,18655,11971  --tool DESeq
Rscript R/FvsMedgeR.R --min 12 --max 13 --exclude 17046 --tool DESeq
Rscript R/FvsMedgeR.R --min 13 --max 14 --exclude 18655 --tool DESeq
Rscript R/FvsMedgeR.R --min 17 --max 20 --exclude 11971 --tool DESeq

# Exclude sex chromosomes
Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq --sex_chromosomes
Rscript R/FvsMedgeR.R --min 12 --max 13 --tool DESeq --sex_chromosomes
Rscript R/FvsMedgeR.R --min 13 --max 14 --tool DESeq --sex_chromosomes
Rscript R/FvsMedgeR.R --min 14 --max 15 --tool DESeq --sex_chromosomes
Rscript R/FvsMedgeR.R --min 15 --max 17 --tool DESeq --sex_chromosomes
Rscript R/FvsMedgeR.R --min 17 --max 20 --tool DESeq --sex_chromosomes


Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool EdgeR
Rscript R/FvsMedgeR.R --min 12 --max 13 --tool EdgeR
Rscript R/FvsMedgeR.R --min 13 --max 14 --tool EdgeR
Rscript R/FvsMedgeR.R --min 14 --max 15 --tool EdgeR
Rscript R/FvsMedgeR.R --min 15 --max 17 --tool EdgeR
Rscript R/FvsMedgeR.R --min 17 --max 20 --tool EdgeR

# Rscript R/FvsMedgeR.R --min 12 --max 20 -cofactor PCW -r 5
# Rscript R/FvsMedgeR.R --min 12 --max 13 -r 5
# Rscript R/FvsMedgeR.R --min 13 --max 14 -r 5
# Rscript R/FvsMedgeR.R --min 14 --max 15 -r 5
# Rscript R/FvsMedgeR.R --min 15 --max 17 -r 5
# Rscript R/FvsMedgeR.R --min 17 --max 20 -r 5
