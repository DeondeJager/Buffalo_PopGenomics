#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=12:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bwa_samtools_01.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bwa_samtools_01.err
#PBS -m abe
#PBS -N bwa_samtools_01
#PBS -M dejager4@gmail.com

# This script was copied and run separately for each sample, changing the read filenames as required.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/read_files/

# Add required modules/software
module add chpc/BIOMODULES
module add bwa/0.7.12
module add samtools/1.3.1

# Map reads using bwa and pass to samtools to sort by coordinates; save as bam file
bwa mem -M -a -t 24 /mnt/lustre/users/djager/buffalo_refgenome/index_bwa/buffalo.final.51scaffolds A_243_14_DSW37619_HCT32ALXX_L8_clean_1P.fq.gz A_243_14_DSW37619_HCT32ALXX_L8_clean_2P.fq.gz | samtools view -@ 24 -b - | samtools sort -@ 24 -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/A_243_14_aln-PE_sorted.bam -O bam /dev/stdin

# Move output file to a new directory, ready for the next analysis. This is not required; up to the user to decide
#mv /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/A_243_14_aln-PE_sorted.bam /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/


