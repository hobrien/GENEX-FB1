from GetSequences import get_sequences
configfile: "config.yaml"

files=get_sequences(config['seqfile'])

rule all:
    input:
        "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts_report.html",
        "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts_report.html",

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
        os.path.join(os.path.dirname(rules.format_bed.output[0]), "genes.fa")
    shell:
        "bedtools getfasta -fi {input.fasta} -bed {input.bed} -split -name -s "
        "| fold -w 60 | perl -pe 's/\([+-]\)//' > {output}"

rule make_index:
    input:
        rules.transcript_seqs.output
    output:
        os.path.join(os.path.dirname(rules.transcript_seqs.output[0]), "kallisto.inx")
    log:
        "Logs/kallisto_index.txt"
    shell:
        "(kallisto index -i {output} {input}) 2> {log}"
        
rule run_kalliso:
    input:
        index = rules.make_index.output,
        reads = lambda wildcards: files[wildcards.sample]
    output:
        "Kallisto/{sample}/abundance.h5",
        "Kallisto/{sample}/abundance.tsv",
        "Kallisto/{sample}/run_info.json"
    params:
        prefix = "Kallisto/{sample}"
    log:
        "Logs/kallisto_quant_{sample}.txt"
    shell:
        "(kallisto quant -i {input.index} -o {params.prefix} --bias -b 100 --rf-stranded {input.reads}) 2> {log}"

rule estimate_counts:
    input:
        expand("Kallisto/{sample}/abundance.tsv", sample=files.keys())
    output:
        "UncorrectedResults/uncorrected/tables/counts_vst.txt"
    params:
        sample_info=config['sample_info'],
        exclude = "none",
        out_name = "UncorrectedResults/uncorrected"
    log:
        "Logs/uncorrected_results.txt"
    shell:
        "(Rscript R/RunDE.R --info {params.sample_info} --min 12 --max 20 --cofactor PCW "
        "--tool DESeq --kallisto --exclude {params.exclude} "
        "--out {params.out_name} {input}) 2> {log}"

rule peer:
    input:
        counts = rules.estimate_counts.output
    output:
        "Peer/factors.txt"
    params:
        sample_info=config["sample_info"],
        residuals = "Peer/residuals.txt",
        alpha = "Peer/alpha.txt",
        num_peer = 10
    log:
        "Logs/PEER/peer.txt"
    conda:
        "env/peer.yaml"
    shell:
        "(Rscript R/PEER.R -n {params.num_peer} -c {input.counts} -b {params.sample_info} "
        "-f {output} -r {params.residuals} -a {params.alpha}) > {log}"

rule join_factors:
    input:
        factors = rules.peer.output,
        sample_info=config["sample_info"]
    output:
        "Data/SampleInfoPeer.txt"
    shell:
        "Rscript R/CombineFactors.R {input.factors} {input.sample_info} {output}"

rule format_cov:
    input:
        config['sample_info']
    output:
        "Data/covariates.txt"
    run:
        import pandas as pd
        cov=pd.read_csv(input[0], sep='\t')
        cov['V1'] = 'sex' + cov['V1'].astype(str) # change sex to string (categorical)
        cov['V4'] = 'batch' + cov['V4'].astype(str) # change batch to string (categorical)
        cov.columns = ['Sample', 'Sex', 'PCW', 'RIN', 'ReadLength', 'PC1', 'PC2', 'PC3'] + ['PEER' + str(i) for i in range(1, 11)]
        cov.to_csv(output[0], sep='\t', header=True, index=False)

rule gene_level:
    input:
        counts = expand("Kallisto/{sample}/abundance.tsv", sample=files.keys()),
        sample_info = rules.format_cov.output
    output:
        "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        batch = ','.join(['RIN', 'ReadLength', 'PC1', 'PC2', 'PC3'] + ['PEER' + str(i) for i in range(1, 11)]),
        out_name = "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --info {input.sample_info} --min 12 --max 20 --cofactor PCW "
        "--batch {params.batch} --tool DESeq --kallisto --exclude {params.exclude} "
        "--out {params.out_name} {input.counts}) 2> {log}"
        
rule transcript_level:
    input:
        counts = expand("Kallisto/{sample}/abundance.tsv", sample=files.keys()),
        sample_info = rules.format_cov.output
    output:
        "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        batch = ','.join(['RIN', 'ReadLength', 'PC1', 'PC2', 'PC3'] + ['PEER' + str(i) for i in range(1, 11)]),
        out_name = "Results/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_PEER_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --info {input.sample_info} --min 12 --max 20 --cofactor PCW "
        "--batch {params.batch} --tool DESeq --kallisto --exclude {params.exclude} "
        "--feature transcripts --out {params.out_name} {input.counts}) 2> {log}"