#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=2:00:00
#PBS -o /mnt/lustre/users/djager/buffalo_refgenome/bwa_index.out
#PBS -e /mnt/lustre/users/djager/buffalo_refgenome/bwa_index.err
#PBS -m abe
#PBS -N bwa_index
#PBS -M dejager4@gmail.com

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buffalo_refgenome/

# Add required modules/software
module add chpc/BIOMODULES
module add bwa/0.7.12

# Execute bwa to index the reference genome
bwa index -p buffalo.final.51scaffolds -a bwtsw buffalo.final.51scaffolds.fas
