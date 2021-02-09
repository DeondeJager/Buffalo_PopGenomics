#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=15:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/subsample.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/subsample.err
#PBS -m abe
#PBS -N A_243_14_subsample
#PBS -M dejager4@gmail.com

# This script was copied and run for the three other samples for which this analysis was also conducted, namely
# B98_509, HC_32 and M_120_13, changing the sample and filenames as required.
# At after preparing the subsampled BAM files, the heterozygosity was calculated as before (script 12, using GL1) for each 
# subsampled BAM file of each sample.

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14

# Path to tools and files (set up environment variables)
REFERENCE='/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa'
INPUT=/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14
PICARD=/home/djager/programs/gatk/GATK_Deon/JavaTools

# Add required modules/software
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add samtools/1.6
module add java/1.8.0_73

# Subsample 50% of reads using samtools:
samtools view -s 1459.50 -@ 24 -b $INPUT/A_243_14_aln-PE_sorted_dups_marked.bam > $INPUT/A_243_14_sub50.bam

# Subsample 25% of reads using samtools:
samtools view -s 1459.25 -@ 24 -b $INPUT/A_243_14_aln-PE_sorted_dups_marked.bam > $INPUT/A_243_14_sub25.bam

# Index bam files (just in case it is necessary)
samtools index -@ 24 -b $INPUT/A_243_14_sub50.bam
samtools index -@ 24 -b $INPUT/A_243_14_sub25.bam


# Use Picard to calculate coverage statistics:
java -Xmx50g -jar  $PICARD/picard.jar CollectWgsMetrics \
          INPUT=$INPUT/A_243_14_sub50.bam \
          OUTPUT=$INPUT/A_243_14_sub50_collectwgsmetrics.txt \
          REFERENCE_SEQUENCE=$REFERENCE \
          MINIMUM_MAPPING_QUALITY=20 \
          MINIMUM_BASE_QUALITY=20

java -Xmx50g -jar  $PICARD/picard.jar CollectWgsMetrics \
          INPUT=$INPUT/A_243_14_sub25.bam \
          OUTPUT=$INPUT/A_243_14_sub25_collectwgsmetrics.txt \
          REFERENCE_SEQUENCE=$REFERENCE \
          MINIMUM_MAPPING_QUALITY=20 \
          MINIMUM_BASE_QUALITY=20

