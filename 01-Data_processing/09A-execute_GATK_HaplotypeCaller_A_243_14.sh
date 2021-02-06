#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=7:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/GATK_A_243_14_HC.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/GATK_A_243_14_HC.err
#PBS -m abe
#PBS -N GATK_A_243_14_HC
#PBS -M dejager4@gmail.com

# This script was copied and run separately for each sample, changing the filename as required. It was used to execute the "GATK_HaplotypeCaller.sh" script.
# The "GATK_HaplotypeCaller.sh" script was copied to each folder containing the BAM file of that sample.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/

# Execute GATK HaplotypeCaller script (9B) and provide the name of the input file
./GATK_HaplotypeCaller.sh A_243_14_aln-PE_sorted_dups_marked.bam
