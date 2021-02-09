#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/ANGSD/ngsF/ngsF-HMM_reps.out
#PBS -e /mnt/lustre/users/djager/ANGSD/ngsF/ngsF-HMM_reps.err
#PBS -m abe
#PBS -N ngsF-HMM_reps
#PBS -M dejager4@gmail.com

# Identify inbreeding (identical by descent) tracts across individual genomes using ngsF-HMM
# See the following website for details: https://github.com/fgvieira/ngsF-HMM
# First need to convert the VCF file from GATK to Beagle format (script 15A) 

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/ANGSD/

# Set up environment variables
OUTDIR=/mnt/lustre/users/djager/ANGSD/ngsF/ngsF-HMM_reps
NGSFPATH=/home/djager/programs/ngsTools/ngsTools/ngsF-HMM
INPUTFILE=/mnt/lustre/users/djager/ANGSD/NGSadmix

# Add required modules/software
# NB: Have to load bzip2 library for angsd to work and gsl for ngsF-HMM
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add /apps/chpc/scripts/modules/bio/lib/gnu/gsl_2.4
module add samtools/1.3.1
module add R/3.4.1-gcc6.3.0

# Calculate the number of sites (SNPs) in the data (if you don't already know). Can also use this command to just double-check the file conversion (in 15A) was successful:
# In an interactive session:
#NSITES=zcat ngsF-HMM_beagle.mafs.gz | wc -l
#      = 3804536 - 1 (it has a header line)
#      = 3804535

# Then, estimate inbreeding tracts using the ngsF-HMM.sh executable provided with ngsF-HMM, which runs 10 repeats of the analysis to prevent convergence to local maxima:
$NGSFPATH/ngsF-HMM.sh --verbose 2 --n_ind 40 --n_sites 3804535 --geno $INPUTFILE/ngsF-HMM_beagle.beagle.gz --lkl --freq 0.1 --indF 0.1,0.1 --min_epsilon 1e-07 --out $OUTDIR/ngsF-HMM_vcf --n_threads 23 --seed 12345 --log 10
