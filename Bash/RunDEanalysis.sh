#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash

Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --feature transcripts
Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto
Rscript R/RunDE.R --min 12 --max 20 --varInt PCW --cofactor Sex --tool DESeqLRT --kallisto
Rscript R/RunDE.R --min 12 --max 20 --varInt PCW --cofactor Sex --tool DESeqLRT --kallisto --feature transcripts
