configfile: "config.yaml"


rule all:
    input:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts_report.html", 
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts_report.html"

rule format_bed:
    input:
        "$REFDIR/Annotation/Genes.gencode/genes.bed"
    output:
        "$REFDIR/Annotation/Genes.gencode/genes.bed12"
    params:
        maxvmem = "4G",
    shell:
        "cut -f 1-12 {input} > {output}"

rule transcript_seqs:
    input:
        bed = rules.format_bed.output,
        fasta = $REFDIR/Sequence/WholeGenomeFasta/genome.fa
    output:
        "$REFDIR/Annotation/Genes.gencode/genes.fa"
    params:
        maxvmem = "4G",
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed}" -split -name -s "
        "| fold -w 60 | perl -pe 's/\([+-]\)//' > {output}"

rule make_index:
    input:
        rules.transcript_seqs.output
    output:
        "$REFDIR/Annotation/Genes.gencode/kallisto.inx"
    params:
        maxvmem = "4G",
    log:
        "Logs/kallisto_index.txt"
    shell:
        "(kallisto index -i {output} {input}) 2> {log}"
        
rule run_kalliso:
    input:
        rules.make_index.output,
    output:
        "Kallisto/{sample}/abundances.h5",
        "Kallisto/{sample}/abundances.tsv",
        "Kallisto/{sample}/run_info.json"
    params:
        maxvmem = "4G",
        prefix = Kallisto/{sample},
        sample = {sample},
    log:
        "Logs/kallisto_quant_{sample}.txt"
    shell:
        "(kallisto quant -i {input.index} -o {params.prefix --bias -b 100 --rf-stranded $sequences) 2> {log}"

rule gene_level:
    input:
        expand(Kallisto/{sample}/abundance.tsv, sample=)
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts_report.html"
    params:
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts.txt"
    shell:
        "(Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto) 2> {log}"
        
rule transcript_level:
    input:
        expand(Kallisto/{sample}/abundance.tsv, sample=)
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts_report.html"
    params:
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/FvsMedgeR.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --feature transcripts) 2> {log}"
