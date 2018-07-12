# Workflow for Sex Differences Analyses in the Developing Human Brain

- This workflow consists of a [Snakemake](https://snakemake.readthedocs.io/en/stable) workflow
that can be executed using the supplied bash scripts:

    ```bash snakemake.sh```

- Software used:
    - [bedtools](http://bedtools.readthedocs.io/en/latest) v2.26.0
    - [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) v1.20.0
    - [kallisto](https://pachterlab.github.io/kallisto/) v0.43.0
    - [pandas](https://pandas.pydata.org) v0.21.0
    - [PEER](https://github.com/PMBio/peer/wiki) v1.3
    - [SARTools](https://github.com/PF2-pasteur-fr/SARTools) v1.3.2
    - [tidyverse](https://www.tidyverse.org) v1.1.1

- Genome References:
    - genome sequence ([GRCh38.p4](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.30) with [Decoy sequences](https://www.ncbi.nlm.nih.gov/assembly/GCA_000786075.2))
    - annotation (derived from [Gencode version 23 (GRCh38.p3)](https://www.gencodegenes.org/releases/23.html))
