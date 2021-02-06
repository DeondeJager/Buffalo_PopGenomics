#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=4:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/trimmomatic_out/trim_final_02.out
#PBS -e /mnt/lustre/users/djager/buf_clean/trimmomatic_out/trim_final_02.err
#PBS -m abe
#PBS -N trim_final_02
#PBS -M dejager4@gmail.com

# This script copied and run separately for each sequencing library changing the read file names as required.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buffalo_rawdata/all_raw_reads/

# Add required modules/software
module add chpc/BIOMODULES
module add trimmomatic/0.36
module add java/1.8.0_73

# Execute trimmomatic
java -jar /apps/chpc/bio/trimmomatic/0.36/bin/trimmomatic-0.36.jar PE -threads 24 -phred33 A_251_14_DSW37618_HCT32ALXX_L7_1.fq.gz A_251_14_DSW37618_HCT32ALXX_L7_2.fq.gz -baseout /mnt/lustre/users/djager/buf_clean/A_251_14_DSW37618_HCT32ALXX_L7_clean.fq.gz ILLUMINACLIP:/apps/chpc/bio/trimmomatic/0.36/bin/adapters/TruSeq3-PE-2.fa:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

