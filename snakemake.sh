#!/bin/bash

snakemake -p --use-conda --cluster-config cluster_config.yaml --cluster "qsub -pe smp {cluster.num_cores} -l h_vmem={cluster.maxvmem}" -j 50 $@ 2> snakemake.log
mail -s "snakemake finished" $EMAIL < snakemake.log
