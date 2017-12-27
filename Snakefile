from GetSequences import get_sequences
configfile: "config.yaml"

files=get_sequences(config['seqfile'])
rep_files=get_sequences(config['rep_seqfile'])

rule all:
    input:
        "Results/BGgenes.txt",
        "Results/BGtranscripts.txt",
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_none_genes_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_none_genes_kallistoCounts_report.html",
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_none_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_none_transcripts_kallistoCounts_report.html",
        expand("Kallisto/{sample}/abundance.tsv", sample=rep_files.keys())

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
        
rule gene_level:
    input:
        sample_info="Data/SampleInfo.txt",
        counts = expand("Kallisto/{sample}/abundance.tsv", sample=files.keys())
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_genes_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_genes_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --exclude {params.exclude}) 2> {log}"
        
rule transcript_level:
    input:
        sample_info="Data/SampleInfo.txt",
        counts = expand("Kallisto/{sample}/abundance.tsv", sample=files.keys())
    output:
        "Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts_report.html"
    params:
        exclude = "{excluded}",
        maxvmem = "20G"
    log:
        "Logs/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_kallistoCounts.txt"
    shell:
        "(Rscript R/RunDE.R --min 12 --max 20 --cofactor PCW --tool DESeq --kallisto --feature transcripts --exclude {params.exclude}) 2> {log}"

rule outlier_summary:
    input:
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_genes_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_genes_kallistoCounts_report.html", excluded=config["gene_level"]),
        expand("Results/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts/Sex_PCW_12_20_FDR_0.1_DESeq_excl_{excluded}_transcripts_kallistoCounts_report.html", excluded=config["transcript_level"]),
    output:
        "Results/BGgenes.txt",
        "Results/BGtranscripts.txt"
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
                                         "Sex_PCW_12_20_FDR_0.1_DESeq_excl_{}_kallistoCounts".format(excluded), 
                                         "tables", "BG12_20.txt")
                        p = Popen(["grep", Id, file], stdout=PIPE, stderr=PIPE)
                        result, err = p.communicate()
                        result=result.decode('ascii').strip().split('\t')
                        updated_genes.loc[len(updated_genes)]=result
                updated_genes=updated_genes.set_index('Id')
                print(updated_genes)
                FittedBias = pd.read_csv(path.join("Results",
                                           "Sex_PCW_12_20_FDR_0.1_DESeq_kallistoCounts",
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
                                         "Sex_PCW_12_20_FDR_0.1_DESeq_excl_{}_transcripts_kallistoCounts".format(excluded),
                                         "tables", "BG12_20.txt")
                        p = Popen(["grep", Id, file], stdout=PIPE, stderr=PIPE)
                        result, err = p.communicate()
                        result=result.decode('ascii').strip().split('\t')
                        try:
                            updated_transcripts.loc[len(updated_transcripts)]=result
                        except ValueError:
                            print("Results could not be found for {}".format(excluded))
                print(updated_transcripts)
                updated_transcripts=updated_transcripts.set_index('Id')
                FittedBias_tr = pd.read_csv(path.join("Results",
                                                   "Sex_PCW_12_20_FDR_0.1_DESeq_transcripts_kallistoCounts",
                                                   "tables", "BG12_20.txt"), 
                                          sep='\t')
                FittedBias_tr.set_index('Id', inplace=True)
                FittedBias_tr.update(updated_transcripts)
                FittedBias_tr.to_csv(path.join("Results", "BGtranscripts.txt"), sep='\t')
            except yaml.YAMLError as exc:
                print(exc)