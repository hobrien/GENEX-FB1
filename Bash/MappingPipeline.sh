#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#


# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution

BASEDIR=/c8000xd3/rnaseq-heath/GENEX-FB1
SampleID=$1 


echo "Starting mapping for $SampleID"

seq_folder=$(grep -P "\s$SampleID(\s|$)" $BASEDIR/Data/sequences.txt | cut -f 5 | head -1)
sequences=$(for name in `grep -P "\s$SampleID(\s|$)" $BASEDIR/Data/sequences.txt | cut -f 1`; do find $seq_folder -name $name*f*q.gz -o -name $name*f*q; done)
set -- $sequences

file1=${1##*/}
file1=${file1%.gz}
file1=${file1%.f*q}
file2=${2##*/}
file2=${file2%.gz}
file2=${file2%.f*q}

if [ ! -f $BASEDIR/FastQC/${file1}_fastqc.html ] | [ ! -f $BASEDIR/FastQC/${file1}_fastqc.zip ]
then
    echo "running FASTQC on $1"
    qsub $BASEDIR/Bash/FastQC.sh $1
fi
if [ ! -f $BASEDIR/FastQC/${file2}_fastqc.html ] | [ ! -f $BASEDIR/FastQC/${file2}_fastqc.zip ]
then
    echo "running FASTQC on $2"
    qsub $BASEDIR/Bash/FastQC.sh $2
fi

if [ ! -d $BASEDIR/Mappings/$SampleID ]
then
    echo "Creating directory $BASEDIR/Mappings/$SampleID"
    mkdir $BASEDIR/Mappings/$SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $BASEDIR/Mappings/$SampleID"
    else
        echo "Could not create directory for $BASEDIR/Mappings/$SampleID"
        exit 1
    fi    
fi

if [ ! -d $BASEDIR/Mappings/$SampleID/BAM ]
then
    echo "Creating directory $BASEDIR/Mappings/$SampleID/BAM"
    mkdir $BASEDIR/Mappings/$SampleID/BAM
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $BASEDIR/Mappings/$SampleID/BAM"
    else
        echo "Could not create directory for $BASEDIR/Mappings/$SampleID/BAM"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/Mappings/$SampleID/accepted_hits.bam ] || [ ! -f $BASEDIR/Mappings/$SampleID/unmapped.bam ]
then
    echo "$MAPPER mapping for $SampleID"
    if [ ! -z $sequences ]
    then
        echo "trying without gzip"
        sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  $BASEDIR/Data/sequences.txt | cut -f 1`; do find $seq_folder -name $name*f*q; done)
    fi
    if [ ! -z $sequences ]
    then
        echo "could not find sequence files for $SampleID"
        exit 1
    fi
    echo "Read files: $sequences"
    sequences="$sequences $BASEDIR/Mappings/$SampleID" # Add output folder to arguments
    qsub -N h${SampleID}_map $BASEDIR/Bash/Tophat2.sh $sequences
fi

if [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.sort.bam ]
then
    echo "Sorting $SampleID"
    qsub -N h${SampleID}_sort -hold_jid h${SampleID}_map \
      $BASEDIR/Bash/SamtoolsSort.sh \
      $BASEDIR/Mappings/$SampleID/accepted_hits.bam $BASEDIR/Mappings/$SampleID/BAM/$SampleID.sort.bam
fi   

if [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.sort.bam.bai ]
then
    echo "Indexing $SampleID"
    qsub -N h${SampleID}_index -hold_jid h${SampleID}_sort \
      $BASEDIR/Bash/SamtoolsIndex.sh \
      $BASEDIR/Mappings/$SampleID/BAM/${SampleID}.sort.bam     
fi

if [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.ex.bam ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.in.bam ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.in.stats.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.ex.stats.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.chr.bam ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.chr.stats.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.expt.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.inner_distance_freq.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.junction.txt ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.junctionSaturation_plot.r ] | \
   [ ! -f $BASEDIR/Mappings/$SampleID/$SampleID.dist.txt ]
then
    echo "Running RNAseqQC $SampleID"
    qsub -N h${SampleID}_stats -hold_jid h${SampleID}_index \
       $BASEDIR/Bash/RNAseqQC.sh $BASEDIR/Mappings/$SampleID/BAM/$SampleID.sort.bam   
fi

if [ ! -f $BASEDIR/Mappings/$SampleID/BAM/$SampleID.chr.counts.txt ]
then
    echo "Running htseq-count $SampleID"
    qsub -N h${SampleID}_count -hold_jid h${SampleID}_stats \
      $BASEDIR/Bash/htseq-count.sh $BASEDIR/Mappings/$SampleID/BAM/$SampleID.chr.bam   
fi

