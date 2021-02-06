#!/bin/bash

# Pipeline for the GATK workflow for calling variants adapted for just fixing read groups and marking duplicates (Tuan Duong and Deon de Jager).
# Original source: https://github.com/Jverma/GATK-pipeline
# Created by Janu Verma on 10/24/14.
# jv367@cornell.edu

# This "fixRG_markdups.sh" script was copied to each folder containing the BAM file of that sample and the sample (directory) name changed as required under "BAM=..." (line 27).
# It is executed using script 7B

# Load modules
module add chpc/BIOMODULES
module add samtools/1.3.1

################################################################

f=$1
fn=${f%.bam}

################################################################

# Set up variables
GATKROOT=/mnt/lustre/users/djager/buffalo_refgenome
PICARDROOT=/home/djager/programs/gatk/GATK_Deon/JavaTools
PICARD="java -Xmx50g -Djava.io.tmpdir=$GATKROOT/tmp -jar $PICARDROOT/picard.jar"
SAMTOOLS="samtools"
REFERENCE="/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa"
BAM="/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/${f}" # Changed the sample (directory) name here

#################################################################

OUTDIR=./gatkout

if [ ! -d $OUTDIR ]; then
        mkdir $OUTDIR
fi


##################################################################
# Fix RG
FixedRGBAM=${BAM/.bam/_fixRG.bam}

function fixReadGroup(){

    $PICARD AddOrReplaceReadGroups \
            INPUT=$BAM \
            OUTPUT=$FixedRGBAM \
            RGID=${fn} \
            RGLB=${fn} \
            RGPL=illumina \
            RGPU=group1 \
            RGSM=${fn} \
            SORT_ORDER=coordinate
}

fixReadGroup


#############################################################
# Mark duplicate reads.
MarkedDupsBam="${BAM%.bam}_dups_marked.bam"

function markDups(){
echo "mark the duplicates in the bam file."

$PICARD MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=4000 MAX_RECORDS_IN_RAM=10000000 INPUT=$FixedRGBAM OUTPUT="${BAM%.bam}_dups_marked.bam" METRICS_FILE="${BAM%.bam}_dups_metrics.txt" REMOVE_DUPLICATES=false

echo "index the dup-marked bam file,"${BAM%.bam}_dups_marked.bam" "

$SAMTOOLS index "${BAM%.bam}_dups_marked.bam"

}

markDups


