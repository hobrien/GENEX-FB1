import GetSequences.py
configfile: "config.yaml"

files=get_sequences(config['seqfile'])

rule all:
    input:
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts_report.html", excluded=config["gene_level"]),
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts_report.html", excluded=config["transcript_level"])

rule format_bed:
    input:
        config['bedfile']
    output:
        config['bedfile'] + '12'
    shell:
        "cut -f 1-12 {input} > {output}"

rule transcript_seqs:
    input:
        bed = rules.format_bed.output,
        fasta = config['refseq']
    output:
        os.path.join(os.path.realpath(rules.format_bed.output, "genes.fa"))
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed}" -split -name -s "
        "| fold -w 60 | perl -pe 's/\([+-]\)//' > {output}"

rule make_index:
    input:
        rules.transcript_seqs.output
    output:
        os.path.join(os.path.realpath(rules.transcript_seqs.output, "kallisto.inx"))
    log:
        "Logs/kallisto_index.txt"
    shell:
        "(kallisto index -i {output} {input}) 2> {log}"
        
rule run_kalliso:
    input:
        index = rules.make_index.output,
        reads = lambda wildcards: files[wildcards.sample]
    output:
        "Kallisto/{sample}/abundances.h5",
        "Kallisto/{sample}/abundances.tsv",
        "Kallisto/{sample}/run_info.json"
    params:
        prefix = Kallisto/{sample}
    log:
        "Logs/kallisto_quant_{sample}.txt"
    shell:
        "(kallisto quant -i {input.index} -o {params.prefix} --bias -b 100 --rf-stranded {input.reads}) 2> {log}"
        
rule gene_level:
    input:
        sample_info="Data/SampleInfo.txt"
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --exclude {params.exclude}) 2> {log}"
        
rule transcript_level:
    input:
        sample_info="Data/SampleInfo.txt"
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --feature transcripts --exclude {params.exclude}) 2> {log}"
