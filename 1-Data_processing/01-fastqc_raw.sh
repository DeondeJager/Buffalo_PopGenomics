#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=12:00:00
#PBS -o /mnt/lustre/users/djager/buffalo_rawdata/all_raw_reads/fastqc_raw.out
#PBS -e /mnt/lustre/users/djager/buffalo_rawdata/all_raw_reads/fastqc_raw.err
#PBS -m abe
#PBS -N fastqc_raw
#PBS -M dejager4@gmail.com

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buffalo_rawdata/all_raw_reads

# Add required modules/software
module add chpc/BIOMODULES
module add FastQC/0.11.5

# Execute fastqc on all fq.gz files in the directory
fastqc --noextract --format fastq -t 24 -o /mnt/lustre/users/djager/buffalo_rawdata/all_raw_reads/ *.fq.gz

