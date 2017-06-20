#!/bin/bash
#
#$ -cwd
#$ -j y
#$ -S /bin/bash
#


# see http://www.tldp.org/LDP/LG/issue18/bash.html for bash Parameter Substitution
SampleID=$1 
BASEDIR=/c8000xd3/rnaseq-heath/Mappings
seq_folder=$(grep -P "\s$SampleID(\s|$)"  ~/GENEX-FB1/Data/sequences.txt | cut -f 5 | head -1)
sequences=$(for name in `grep -P "\s$SampleID(\s|$)"  ~/GENEX-FB1/Data/sequences.txt | cut -f 1`; do find $seq_folder -name $name*f*q.gz; done)

echo "Starting mapping for $BASEDIR/$SampleID"
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

if [ ! -f $BASEDIR/$SampleID/BAM/accepted_hits.bam ] || [ ! -f $BASEDIR/$SampleID/BAM/unmapped.bam ]
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
      $BASEDIR/$SampleID/BAM/accepted_hits.bam $BASEDIR/$SampleID/BAM/$SampleID.sort.bam
fi   

if [ ! -f $BASEDIR/$SampleID/BAM/$SampleID.sort.bam.bai ]
then
    echo "Indexing $SampleID"
    qsub -N h${SampleID}_index -hold_jid h${SampleID}_sort \
      ~/GENEX-FB1/Bash/SamtoolsIndex.sh \
      $BASEDIR/$SampleID/${SampleID}_sort.bam     
fi
