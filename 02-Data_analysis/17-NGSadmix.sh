#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/ANGSD/NGSadmix/NGSadmixK1.out
#PBS -e /mnt/lustre/users/djager/ANGSD/NGSadmix/NGSadmixK1.err
#PBS -m abe
#PBS -N NGSadmixK1
#PBS -M dejager4@gmail.com

# Estimate individual admixture proportions using NGSadmix.
# See the following website for details: http://www.popgen.dk/software/index.php/NgsAdmix#Installation
# And this one for the tutorial on which this script is based: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
# Uses the same input file as ngsF-HMM. See script 15A on how this file prepared.
# I just renamed the input file for consistency: mv ngsF-HMM_beagle.beagle.gz NGSadmix_beagle.beagle.gz

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/ANGSD/NGSadmix/

# Set up environment variables
OUTDIR=/mnt/lustre/users/djager/ANGSD/NGSadmix
NGSTOOLS=/home/djager/programs/ngsTools/ngsTools
NGSADMIX=/home/djager/programs/NGSadmix
INPUTFILE=/mnt/lustre/users/djager/ANGSD/NGSadmix

# Add required modules/software
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add R/3.4.1-gcc6.3.0

# Run NGSadmix for K = 1. I repeated the analysis for K 2-6 by just changing the value of K and filenames as required.
$NGSADMIX/NGSadmix -likes $INPUTFILE/NGSadmix_beagle.beagle.gz \
  -K 1 \
  -o $INPUTFILE/NGSadmix_K1 \
  -P 24 \
  -tolLike50 0.01 \
  -minMaf 0

# Plot using R
Rscript $NGSTOOLS/Scripts/plotAdmix.R -i $INPUTFILE/NGSadmix_K1.qopt -o $INPUTFILE/NGSadmix_K1.admix.pdf
