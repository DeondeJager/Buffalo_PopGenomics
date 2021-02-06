#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=1:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/GATK_Deon/A_243_14/picard_metrics.log
#PBS -e /mnt/lustre/users/djager/buf_clean/GATK_Deon/A_243_14/picard_metrics.err
#PBS -m abe
#PBS -N A_243_14_picard_metrics
#PBS -M dejager4@gmail.com

# Collect alignment metrics using Picard

# This script was copied and run separately for each sample, changing the filenames as required.

# Add modules
module add chpc/BIOMODULES
module add gcc/6.2.0
module add R/3.3.1-gcc6.2.0

# Define file and folder paths
DATA=/mnt/lustre/users/djager/buf_clean/GATK_Deon

# Define software
PICARD="java -jar /home/djager/programs/picard/2.17.11/picard.jar"

# Navigate to relevant folder containing the data, or to parent folder
cd $DATA

# Have to create a sequence dictionary of the reference sequence, if not already done
#$PICARD CreateSequenceDictionary \
#  R=/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fas \
#  O=/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.dict

# Collect multiple metrics - Removed CollectGcBiasMetrics, as it was giving an error
$PICARD CollectMultipleMetrics \
  R=/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa \
  I=$DATA/A_243_14/A_243_14_aln-PE_sorted_dups_marked.bam \
  O=$DATA/A_243_14/A_243_14_aln-PE_sorted_dups_marked \
  PROGRAM=null \
  PROGRAM=CollectAlignmentSummaryMetrics \
  PROGRAM=CollectInsertSizeMetrics \
  PROGRAM=CollectBaseDistributionByCycle \
  PROGRAM=CollectQualityYieldMetrics

