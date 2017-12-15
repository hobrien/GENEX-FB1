import GetSequences.py
configfile: "config.yaml"

files=get_sequences(config['seqfile'])

rule all:
    input:
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts_report.html", excluded=config["gene_level"]),
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts_report.html", excluded=config["transcript_level"]),
        expand("Kallisto/{sample}/abundances.tsv", sample=set([string.split('-')[0] for string in files.keys()])), #remove run numbers from sample names
        expand("BAM/{seq_run}.sort.bam.bai", seq_run=files.keys())

# HiSat is run separately on each sequencing run from each sample
rule hisat:
    input:
        reads = lambda wildcards: files[wildcards.seq_run]
    output:
        temp("BAM/{seq_run}.bam")
    params:
        idx = config["index"],
        extra = '--known-splicesite-infile ' + config["splice_sites"]
    benchmark:
        "Benchmarks/{seq_run}.hisat.benchmark.txt"
    log:
        "Logs/{seq_run}_hisat_map.txt"
    threads: 8
    wrapper:
"0.17.4/bio/hisat2"

rule sort_bam:
    input:
        rules.hisat.output
    output:
        "BAM/{seq_run}.sort.bam"
    params:
        "-m 4G"
    threads: 8
    wrapper:
        "0.17.4/bio/samtools/sort"
        
rule samtools_index:
    input:
        rules.sort_bam.output
    output:
        "BAM/{seq_run}.sort.bam.bai"
    wrapper:
"0.17.4/bio/samtools/index"

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

# input is an ugly function to concatenate read files from all sequencing runs for the same sample 
rule run_kalliso:
    input:
        index = rules.make_index.output,
        reads = lambda wildcards: reduce(lambda x,y: x+y, [files[key] for key in filter(lambda key: key.startswith(wildcards.sample), files.keys())]) 
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
