configfile: "config.yaml"


rule all:
    input:
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts_report.html", excluded=config["high_cooks"])

rule high_cooks:
    input:
        sample_info="Data/SampleInfo.txt",
        vcf="Genotypes/{run}/hg19/chr{chr_num}.dose.vcf.gz"
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts_report.html"
    params:
        {excluded}
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --exclude {excluded}) 2> {log}"