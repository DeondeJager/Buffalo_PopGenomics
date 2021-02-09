#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/ANGSD/angsd_VCF_to_beaglebinary.out
#PBS -e /mnt/lustre/users/djager/ANGSD/angsd_VCF_to_beaglebinary.err
#PBS -m abe
#PBS -N angsd_VCF_to_beaglebinary
#PBS -M dejager4@gmail.com

# Use ANGSD to convert VCF file to Beagle binary format for use in ngsF

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/ANGSD/

# Set up enviroment variables
ANGSDPATH=/home/djager/programs/angsd/angsd
OUTDIR=/mnt/lustre/users/djager/ANGSD
VCFPATH=/mnt/lustre/users/djager
REFERENCE='/mnt/lustre/users/djager/buffalo_refgenome/samtools_index/buffalo.final.51scaffolds.fa.fai'

# Add required modules/software
# Have to load bzip2 library for angsd to work
module add chpc/BIOMODULES
module add bzip2/1.0.6

# Convert vcf file into Beagle binary format for use in ngsF
$ANGSDPATH/angsd -vcf-gl $VCFPATH/buf_all_SNPs_allscaffolds_filtered_tabix.vcf.gz \
  -GL 0 \
  -out $OUTDIR/ngsF_beaglebinary \
  -fai $REFERENCE \
  -nThreads 24 \
  -doGlf 3 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -nind 40
