#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#

export PATH=/share/apps/R-3.2.2/bin:/share/apps/:$PATH

for dataset in $@
do
    folder_path=${dataset%/*}
    folder_path=${folder_path%/*}
    folder=${folder_path##*/}
    echo "Starting QC for $dataset"
# [BAMQC](https://github.com/s-andrews/BamQC)
#bamqc --outdir=$folder_path --gff /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.gtf $@

# [PicardTools](http://broadinstitute.github.io/picard/)
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CreateSequenceDictionary R=/home/heath/Ref/hg19.fa O=/home/heath/Ref/hg19.dict
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar ReorderSam INPUT=/home/heath/Mappings/15533_300/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300/accepted_hits_sorted.bam REFERENCE=/home/heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/WholeGenomeFasta/genome.fa
#/home/heath/bin/java -Xmx2g -jar /home/heath/src/picard-tools-2.1.1/picard.jar CollectRnaSeqMetrics REF_FLAT=/home/heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/refFlat.txt.gz STRAND_SPECIFICITY=SECOND_READ_TRANSCRIPTION_STRAND INPUT=/home/heath/Mappings/15533_300_secondstrand/accepted_hits.bam OUTPUT=/home/heath/Mappings/15533_300_secondstrand/RnaSeqMetrics.txt ASSUME_SORTED=false

#http://rseqc.sourceforge.net/
    if [ ! -f $folder_path/BAM/$folder.in.bam ] | [ ! -f $folder_path/BAM/$folder.ex.bam ]
    then
        echo "Splitting BAM for $folder"
        split_bam.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Sequence/AbundantSequences/humRibosomal.bed -i $dataset -o $folder_path/BAM/$folder   
        if [ $? -eq 0 ]
        then
            echo "Finished splitting BAM for $folder"
        else
            echo "Could not split BAM for $folder"
            exit 1
        fi 
    fi
    if [ ! -f $folder_path/$folder.in.stats.txt ]
    then
        echo "Running bam_stat on ribosomal reads for $folder"
        bam_stat.py -i $folder_path/BAM/$folder.in.bam > $folder_path/$folder.in.stats.txt
        if [ $? -eq 0 ]
        then
            echo "Finished running bam_stat on ribosomal reads for $folder"
        else
            echo "Could not run bam_stat on ribosomal reads for $folder"
            exit 1
        fi 
    fi
    if [ ! -f $folder_path/$folder.ex.stats.txt ]
    then
        echo "Running bam_stat on non-ribosomal reads for $folder"
        bam_stat.py -i $folder_path/BAM/$folder.ex.bam > $folder_path/$folder.ex.stats.txt
        if [ $? -eq 0 ]
        then
            echo "Finished running bam_stat on non-ribosomal reads for $folder"
        else
            echo "Could not run bam_stat on non-ribosomal reads for $folder"
            exit 1
        fi 
    fi
    if [ ! -f $folder_path/BAM/$folder.chr.bam ]
    then
        echo "Extracting chromosome reads for $folder"
        samtools index $folder_path/BAM/$folder.ex.bam
        samtools view -bh $folder_path/BAM/$folder.ex.bam \
          chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY \
          > $folder_path/BAM/$folder.chr.bam
        samtools index $folder_path/BAM/$folder.chr.bam

        if [ $? -eq 0 ]
        then
            echo "Finished extracting chromosome reads for $folder"
        else
            echo "Could not extract chromosome reads for $folder"
            exit 1
        fi 
    fi
        

    if [ ! -f $folder_path/$folder.chr.stats.txt ]
    then
        echo "calculating stats for $folder"
        bam_stat.py -i $folder_path/BAM/$folder.chr.bam > $folder_path/$folder.chr.stats.txt
        if [ $? -eq 0 ]
        then
            echo "Finished calculating stats for $folder"
        else
            echo "Could not calculate stats for $folder"
            exit 1
        fi 
    fi

    # determine the strand of experiment ("1++,1--,2+-,2-+" = first strand, "1+-,1-+,2++,2--" = second strand)
    if [ ! -f $folder_path/$folder.expt.txt ]
    then
        echo "determining strand for $folder"
        infer_experiment.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam > $folder_path/$folder.expt.txt
        if [ $? -eq 0 ]
        then
            echo "Finished determining strand for $folder"
        else
            echo "Could not determine strand for $folder"
            exit 1
        fi 
    fi

    # plot distribution of insert sizes (size - total read length)
    if [ ! -f $folder_path/$folder.inner_distance_freq.txt ]
    then
        echo "determining insert sizes for $folder"
        inner_distance.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder -u 1000 -s 10 >/dev/null
        if [ $? -eq 0 ]
        then
            echo "Finished determining insert sizes for $folder"
        else
            echo "Could not determine insert sizes for $folder"
            exit 1
        fi 
    fi

    # the necessary output from this is going to the log file, not to $folder.junction.txt
    if [ ! -f $folder_path/$folder.junction.txt ]
    then
        echo "annotating junctions for $folder"
        junction_annotation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder  >  $folder_path/$folder.junction.txt
        if [ $? -eq 0 ]
        then
            echo "Finished annotating junctions for $folder"
        else
            echo "Could not annotate junctions for $folder"
            exit 1
        fi 
    fi

    if [ ! -f $folder_path/$folder.junctionSaturation_plot.r ]
    then
        echo "Plotting junction saturation for $folder"
        junction_saturation.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder
        if [ $? -eq 0 ]
        then
            echo "Finished plotting junction saturation for $folder"
        else
            echo "Could not plot junction saturation for $folder"
            exit 1
        fi 
    fi

    if [ ! -f $folder_path/$folder.dist.txt ]
    then
        echo "Calculating read distribution for $folder"
        read_distribution.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam > $folder_path/$folder.dist.txt
        if [ $? -eq 0 ]
        then
            echo "Finished calculating read distribution for $folder"
        else
            echo "Could not calculate read distribution for $folder"
            exit 1
        fi 
    fi

    #read_duplication.py -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder

    #geneBody_coverage.py -r /c8000xd3/rnaseq-heath/Ref/Homo_sapiens/GRCh38/NCBI/GRCh38Decoy/Annotation/Genes.gencode/genes.bed -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder


    #deletion_profile.py -l 120 -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder
    #insertion_profile.py -s PE -i $folder_path/BAM?$folder.chr.bam -o $folder_path/$folder
    #mismatch_profile.py -l 120 -i $folder_path/BAM/$folder.chr.bam -o $folder_path/$folder

    echo "finished QC for $dataset"
done
exit $?
