#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/ANGSD/ngsF/ngsF_reps/ngsF_reps.out
#PBS -e /mnt/lustre/users/djager/ANGSD/ngsF/ngsF_reps/ngsF_reps.err
#PBS -m abe
#PBS -N ngsF_reps
#PBS -M dejager4@gmail.com

# Calculate individual inbreeding coefficients using ngsF
# Based on the tutorial at https://github.com/fgvieira/ngsF/tree/master/examples. 
# First need to convert VCF file from GATK to binary beagle format (script 14A)

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/ANGSD/ngsF

# Set up environment variables
OUTDIR=/mnt/lustre/users/djager/ANGSD/ngsF/ngsF_reps
NGSFPATH=/home/djager/programs/ngsTools/ngsTools/ngsF
INPUTFILE=/mnt/lustre/users/djager/ANGSD/ngsF

# Add required modules/software
# Have to load bzip2 library for angsd to work and gsl for ngsF
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add /apps/chpc/scripts/modules/bio/lib/gnu/gsl_2.4
module add samtools/1.3.1

# Calculate the number of sites (SNPs) in the data (if you don't already know). Can also use this command to just double-check the file conversion (in 14A) was successful:
# In an interactive session:
#N_SITES=$((`zcat ngsF_beaglebinary.mafs.gz | wc -l`-1))
#       = 3804536 -1 (it has a header line)
#       = 3804535

# Then, estimate inbreeding coefficients using the ngsF.sh executable provided with ngsF, which runs 10 repeats of the analysis to prevent convergence to local maxima:
# Initial (quick) run to get approximate values (as per the tutorial)
$NGSFPATH/ngsF.sh --glf $INPUTFILE/ngsF_beaglebinary.glf --n_ind 40 --n_sites 3804535 --min_epsilon 1e-05 --out $OUTDIR/testF_reps.approx_indF --approx_EM --seed 12345 --init_values r --n_threads 24

# Final (full) run using the initial values as priors
$NGSFPATH/ngsF.sh --glf $INPUTFILE/ngsF_beaglebinary.glf --n_ind 40 --n_sites 3804535 --min_epsilon 1e-07 --out $OUTDIR/inbreeding_values_ngsF_vcf_reps.indF --init_values $OUTDIR/testF_reps.approx_indF.pars --n_threads 24

