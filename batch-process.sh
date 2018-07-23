#!/bin/bash

# This a general framework for running a pipeline in batch on bunch of files

# Read in tsv file of samples to be processed one at a time

tail -n +2 samples-list.tsv | while read line; do

#make sure necessary directories exist

mkdir $(pwd)/output
mkdir $(pwd)/logs

#set up variables

ID=`echo "$line" | cut -f 1`
SAMPLE=`echo "$line" | cut -f 2`
INPUT1=`echo "$line" | cut -f 3`
echo $INPUT1
INPUT2=`echo "$line" | cut -f 4`
echo $INPUT2
REFERENCE="/path/to/my/reference/reference.fasta"
OUTDIR="./output"
LOGDIR="./logs"

# Execute pipeline.  In this case it's mapping reads against a reference using BWA, followed by using samtools to convert to BAM and sorting

#### Step 1: Map with BWA (mem algorithm) #####
echo "mapping $ID"
bwa mem -t 8 -R "@RG\tID:$ID\tSM:$SAMPLE" $REFERENCE $INPUT1 $INPUT2 2> $LOGDIR/$ID.bwa.mem.log | samtools view -bS - > $OUTDIR/$ID.bam 2> $LOGDIR/$ID.bam.log

#### Step 2: Sort ####
echo "sorting $ID"
samtools sort --threads 4 -o $OUTDIR/$ID.sorted.bam $OUTDIR/$ID.bam 2> $LOGDIR/$ID.sorted.bam.log

done
