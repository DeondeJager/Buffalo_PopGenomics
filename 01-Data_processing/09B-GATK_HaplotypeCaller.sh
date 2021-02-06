#!/bin/bash

# Pipeline for the GATK workflow for calling variants adapted for FABI by Tuan Duong and then by Deon de Jager for personal use.
# Original source: https://github.com/Jverma/GATK-pipeline
# Uses HaplotypeCaller.
# Created by Janu Verma on 10/24/14.
# jv367@cornell.edu

# This "GATK_HaplotypeCaller.sh" script was copied to each folder containing the BAM file of that sample and the sample (directory) name changed as required.
# It is executed using script 9A. 

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/

# Add required modules/software
module add chpc/BIOMODULES
module add samtools/1.3.1
module add picard/2.2.1
module add java/1.8.0_73

f=$1
fn=${f%.bam}

					
##################################################################
#   SOFTWARES
# Provide locations of the softwares to be used. 
OUTDIR=./gatkout

if [ ! -d $OUTDIR ]; then
	mkdir $OUTDIR 
fi

#################################################################
#  Set up environment variables

GATKROOT=/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14
JARROOT=$GATKROOT/JavaTools
DATAROOT=$GATKROOT/refs
BAM="$GATKROOT/${f}"
REFERENCE="$DATAROOT/buffalo.final.51scaffolds.fa"
PICARD="java -Xmx10g -Djava.io.tmpdir=$GATKROOT/tmp -jar $JARROOT/picard.jar"
SAMTOOLS="samtools"
GATK="java -Xmx40g -jar $JARROOT/GenomeAnalysisTK.jar"
FixedRGBAM=${BAM/.bam/_fixRG.bam} # From scripts 7B and 7C
MarkedDupsBam="${BAM%.bam}_dups_marked.bam" # From scripts 7B and 7C 


###########################################################################
# GATK Variant Calling -  HaplotypeCaller
# Set -nct, outmode, emit_thresh, call_threh,

outmode="EMIT_ALL_CONFIDENT_SITES"
#emit_thresh=20	# Threshold for tagging possible variants
call_thresh=30	# Threshold for tagging _good_ variants
hetrate=0.03	# Popgen heterozygosity rate (that is, for any two random chrom in pop, what is rate of mismatch). Human is ~0.01, so up buffalo to ~0.03
minBaseScore=20	# Minimum Phred base score to count a base (20 = 0.01 error, 30=0.001 error, etc)

echo "calling variants...."
$GATK \
    -T HaplotypeCaller \
    -R $REFERENCE \
    -I /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/A_243_14_aln-PE_sorted_dups_marked.bam \
    --emitRefConfidence GVCF \
    --variant_index_type LINEAR \
    --variant_index_parameter 128000 \
    -hets $hetrate \
    -mbq $minBaseScore \
    -stand_call_conf $call_thresh \
    -out_mode $outmode \
    -nct 24 \
    -ploidy 2 \
    -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/A_243_14_aln-PE_sorted_dups_marked_output.raw.snps.indels.g.vcf
