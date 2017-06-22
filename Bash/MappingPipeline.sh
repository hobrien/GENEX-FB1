#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#


# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
echo "Starting mapping for $BASEDIR/$SampleID"

SampleID=$1 
BASEDIR=/c8000xd3/rnaseq-heath/Mappings
seq_folder=$(grep -P "\s$SampleID(\s|$)"  ~/GENEX-FB1/Data/sequences.txt | cut -f 5 | head -1)
sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/GENEX-FB1/Data/sequences.txt | cut -f 1`; do find $seq_folder -name $name*f*q.gz; done)
set -- $sequences

file1=${1##*/}
file2=${2##*/}
echo ${file1%%.*}_fastqc.html

if [ ! -f ${file1%%.*}_fastqc.html ]
then
    echo "running FASTQC on $1"
    qsub ~/GENEX-FB1/Bash/FastQC.sh $1
fi
if [ ! -f ${file2%%.*}_fastqc.html ]
then
    echo "running FASTQC on $2"
    qsub ~/GENEX-FB1/Bash/FastQC.sh $2
fi

if [ ! -d $BASEDIR/$SampleID ]
then
    echo "Creating directory $BASEDIR/$SampleID"
    mkdir $BASEDIR/$SampleID
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $BASEDIR/$SampleID"
    else
        echo "Could not create directory for $BASEDIR/$SampleID"
        exit 1
    fi    
fi

if [ ! -d $BASEDIR/$SampleID/BAM ]
then
    echo "Creating directory $BASEDIR/$SampleID/BAM"
    mkdir $BASEDIR/$SampleID/BAM
    if [ $? -eq 0 ]
    then
        echo "Finished creating directory for $BASEDIR/$SampleID/BAM"
    else
        echo "Could not create directory for $BASEDIR/$SampleID/BAM"
        exit 1
    fi    
fi

if [ ! -f $BASEDIR/$SampleID/accepted_hits.bam ] || [ ! -f $BASEDIR/$SampleID/unmapped.bam ]
then
    echo "$MAPPER mapping for $SampleID"
    if [ ! $sequences ]
    then
        echo "trying without gzip"
        sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/GENEX-FB1/Data/sequences.txt | cut -f 1`; do find $seq_folder -name $name*f*q; done)
    fi
    if [ ! $sequences ]
    then
        echo "could not find sequence files for $SampleID"
        exit 1
    fi
    echo "Read files: $sequences"
    sequences="$sequences $BASEDIR/$SampleID" # Add output folder to arguments
    qsub -N h${SampleID}_map ~/GENEX-FB1/Bash/Tophat2.sh $sequences
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.sort.bam ]
then
    echo "Sorting $SampleID"
    qsub -N h${SampleID}_sort -hold_jid h${SampleID}_map \
      ~/GENEX-FB1/Bash/SamtoolsSort.sh \
      $BASEDIR/$SampleID/accepted_hits.bam $BASEDIR/$SampleID/BAM/$SampleID.sort.bam
fi   

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.sort.bam.bai ]
then
    echo "Indexing $SampleID"
    qsub -N h${SampleID}_index -hold_jid h${SampleID}_sort \
      ~/GENEX-FB1/Bash/SamtoolsIndex.sh \
      $BASEDIR/$SampleID/BAM/${SampleID}.sort.bam     
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.ex.bam ] | \
   [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.in.bam ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.in.stats.txt ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.ex.stats.txt ] | \
   [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.chr.bam ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.chr.stats.txt ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.expt.txt ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.inner_distance_freq.txt ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.junction.txt ] | \
   [ ! -f $BASEDIR/$SampleID/$SampleID.junctionSaturation_plot.r ]
then
    echo "Running RNAseqQC $SampleID"
    qsub -N h${SampleID}_stats -hold_jid h${SampleID}_index \
       ~/GENEX-FB1/Bash/RNAseqQC.sh $BASEDIR/$SampleID/BAM/$SampleID.sort.bam   
fi

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.chr.counts.txt ]
then
    echo "Running htseq-count $SampleID"
    qsub -N h${SampleID}_count -hold_jid h${SampleID}_stats \
      ~/GENEX-FB1/Bash/htseq-count.sh $BASEDIR/$SampleID/BAM/$SampleID.chr.bam   
fi

