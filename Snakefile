configfile: "config.yaml"

rule all:
    input:
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts_report.html", excluded=config["gene_level"]),
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts_report.html", excluded=config["transcript_level"])
rule gene_level:
    input:
        sample_info="Data/SampleInfo.txt"
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        out_name = "Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --exclude {params.exclude} --out {params.out_name}) 2> {log}"
        
rule transcript_level:
    input:
        sample_info="Data/SampleInfo.txt"
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        out_name = "Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --feature transcripts --exclude {params.exclude} --out {params.out_name}) 2> {log}"
