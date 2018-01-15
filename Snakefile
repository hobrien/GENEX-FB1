from GetSequences import get_sequences
configfile: "config.yaml"

files=get_sequences(config['seqfile'])

rule all:
    input:
        "Results/BGgenes.txt",
        "Results/BGtranscripts.txt",
        "Peer/factors.txt"

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

rule gene_level:
    input:
        expand("Kallisto/{sample}/abundance.tsv", sample=files.keys())
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts_report.html"
    params:
        sample_info=config['sample_info'],
        exclude = "{excluded}",
        out_name = "Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --info {params.sample_info} --min 12 --max 20 --cofactor PCW "
        "--tool DESeq --kallisto --exclude {params.exclude} "
        "--out {params.out_name} {input}) 2> {log}"
        
rule transcript_level:
    input:
        expand("Kallisto/{sample}/abundance.tsv", sample=files.keys())
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts_report.html"
    params:
        sample_info=config['sample_info'],
        exclude = "{excluded}",
        out_name = "Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --info {params.sample_info} --min 12 --max 20 --cofactor PCW "
        "--tool DESeq --kallisto --feature transcripts --exclude {params.exclude} "
        "--out {params.out_name} {input}) 2> {log}"

rule outlier_summary:
    input:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts_report.html",
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts_report.html",
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{excluded}_kallistoCounts_report.html", excluded=config["gene_level"]),
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{excluded}_kallistoCounts_report.html", excluded=config["transcript_level"]),
    output:
        "Results/BGgenes.txt",
        "Results/BGtranscripts.txt"
    params:
        transcripts = "Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_none_kallistoCounts",
        genes = "Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_none_kallistoCounts"
    run:
        import yaml
        import pandas as pd
        from subprocess import Popen, PIPE
        from os import path

        with open("config.yaml", 'r') as stream:
            try:
                high_cooks =(yaml.load(stream))
                updated_genes=pd.DataFrame(columns=('Id', 'SYMBOL', 'Chr', 'gene_type',
                                                    'baseMean', 'Male', 'Female', 'FC',
                                                    'log2FoldChange', 'pvalue', 'padj',
                                                    'maxCooks'))
                print("Updating genes with outliers removed:")
                for excluded in high_cooks['gene_level'].keys():
                    for Id in high_cooks['gene_level'][excluded].split('_'):
                        file = path.join("Results", 
                                         "Sex_PCW_12_20_FDR_0.1_DESeq_genes_excl_{}_kallistoCounts".format(excluded), 
                                         "tables", "BG12_20.txt")
                        p = Popen(["grep", Id, file], stdout=PIPE, stderr=PIPE)
                        result, err = p.communicate()
                        #print(result)
                        result=result.decode('ascii').strip().split('\t')
                        try:
                            updated_genes.loc[len(updated_genes)]=result
                        except ValueError:
                            print("Results could not be found for {}".format(excluded))
                updated_genes=updated_genes.set_index('Id')
                #print(updated_genes)
                FittedBias = pd.read_csv(path.join("Results", params['genes'],
                                                   "tables", "BG12_20.txt"), 
                                          sep='\t')
                FittedBias.set_index('Id', inplace=True)
                FittedBias.update(updated_genes)
                FittedBias.to_csv(path.join("Results", "BGgenes.txt"), sep='\t')
                updated_transcripts=pd.DataFrame(columns=('Id', 'gene_id', 'SYMBOL', 'Chr', 
                                                    'gene_type', 'baseMean', 'Male', 'Female', 
                                                    'FC', 'log2FoldChange', 'pvalue', 'padj',
                                                    'maxCooks'))
                print("Updating transcripts with outliers removed:")
                for excluded in high_cooks['transcript_level'].keys():
                    for Id in high_cooks['transcript_level'][excluded].split('_'):
                        file = path.join("Results", 
                                         "Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_excl_{}_kallistoCounts".format(excluded),
                                         "tables", "BG12_20.txt")
                        p = Popen(["grep", Id, file], stdout=PIPE, stderr=PIPE)
                        result, err = p.communicate()
                        result=result.decode('ascii').strip().split('\t')
                        try:
                            updated_transcripts.loc[len(updated_transcripts)]=result
                        except ValueError:
                            print("Results could not be found for {}".format(excluded))
                #print(updated_transcripts)
                updated_transcripts=updated_transcripts.set_index('Id')
                FittedBias_tr = pd.read_csv(path.join("Results", params['transcripts'],
                                                   "tables", "BG12_20.txt"), 
                                          sep='\t')
                FittedBias_tr.set_index('Id', inplace=True)
                FittedBias_tr.update(updated_transcripts)
                FittedBias_tr.to_csv(path.join("Results", "BGtranscripts.txt"), sep='\t')
            except yaml.YAMLError as exc:
                print(exc)
