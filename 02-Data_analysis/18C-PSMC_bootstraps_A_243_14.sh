#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=12:00:00
#PBS -o /mnt/lustre/users/djager/PSMC/A_243_14/A_243_14_psmc_bootstrap.out
#PBS -e /mnt/lustre/users/djager/PSMC/A_243_14/A_243_14_psmc_bootstrap.err
#PBS -m abe
#PBS -N A_243_14_psmc_bootstrap
#PBS -M dejager4@gmail.com

# This script was copied and run for each sample, changing sample and filenames as required.
# This is the PSMC analysis with bootstraps, where the genome is divided into segments and PSMC is run on random segments, with replacement, to generate bootstrap replicates.
# Many thanks to Michael Westbury for helping me out with this script.

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/PSMC/A_243_14

# Add the required modules/software
module add chpc/BIOMODULES
module add samtools/1.3.1 # It seems loading the samtools module automatically loads the corresponding version of bcftools.
module add perl/5.16.3    # Add perl version that vcfutils.pl runs in. Don't know if this is necessary, but just in case.
module add psmc/0.6.5
module add gnu/plot-5.0.4

# Split long sequences in your psmcfa file (generated in script 18B) to be randomly sampled with replacement
splitfa A_243_14.psmcfa > A_243_14_split.psmcfa
# Bootstrap replicates (with multiple threads)
seq 100 | xargs -P 24 -i psmc -N25 -t15 -r5 -b -p "4+25*2+4+6" -o A_243_14_round-{}.psmc A_243_14_split.psmcfa | sh
# seq = number of bootstraps
# -P = number of threads
# -o = output file

#Combine the psmc output files (100 files)
cat A_243_14_PSMC_human_r5.psmc A_243_14_round-*.psmc > A_243_14_combinedboots.psmc

# Quick plot to check results
psmc_plot.pl -u 1.5e-08 -g 7.5 -p A_243_14_combinedboots A_243_14_combinedboots.psmc
# This was not how the final plot was generated- see the "03-Figures" folder for the R scripts to generate PSMC plots.
