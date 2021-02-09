#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=60:00:00
#PBS -o /mnt/lustre/users/djager/buffalo_refgenome/satsuma_synteny.out
#PBS -e /mnt/lustre/users/djager/buffalo_refgenome/satsuma_synteny.err
#PBS -m abe
#PBS -N satsuma_synteny
#PBS -M dejager4@gmail.com

# Identify X chromosome regions in buffalo genome, using syntenic alignment with Bos taurus X chromosome.
# This was done to mask the X chr regions in the buffalo genome (artificially concatenated) for the PSMC analysis.

# Add required modules/software
module add chpc/BIOMODULES
module add spines/1.1.15

# Set up environment variables
QUERYDIR=/mnt/lustre/users/djager/buffalo_refgenome
TARGET="/mnt/lustre/users/djager/bos_taurus_ref_genome/UMD3_1_chrX.fa"

# Navigate to directory with query fasta 
cd $QUERYDIR

# Run SatsumaSynteny with default parameters
SatsumaSynteny -q /mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fas \
  -t $TARGET \
  -o $QUERYDIR
