configfile: "config.yaml"

rule all:
    input:
        "Results/BGgenes.txt",
        "Results/BGtranscripts.txt"

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