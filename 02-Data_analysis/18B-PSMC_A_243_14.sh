#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/PSMC/A_243_14/psmc_A_243_14.out
#PBS -e /mnt/lustre/users/djager/PSMC/A_243_14/psmc_A_243_14.err
#PBS -m abe
#PBS -N A_243_14_psmc
#PBS -M dejager4@gmail.com

# This script was copied and run for each sample, changing sample and filenames as required.
# This is the "main" PSMC analysis of the entire genome (no bootstraps): 
# For plotting the PSMC results, we need a "main" result to plot as a dark line and then also the bootstraps (script 18C) to plot as thin lines around the dark "main" line.

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/PSMC/A_243_14

# Add the required modules/software
module add chpc/BIOMODULES
module add samtools/1.3.1 # It seems loading the samtools module automatically loads the corresponding version of bcftools.
module add perl/5.16.3    # Add perl version that vcfutils.pl runs in. Don't know if this is necessary, but just in case.
module add psmc/0.6.5
module add gnu/plot-5.0.4

# Set up environment variables
# "REGIONS" excludes Super_Scaffold41, which is the majority of the bovine chromosome X (as determined using satsuma_synteny)
REGIONS="/home/djager/scripts/PSMC/nonChrX.bed"
REFERENCE="/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa"
BAM="/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/A_243_14_aln-PE_sorted_dups_marked.bam"

# Generate the consensus sequence for one diploid individual
# Break up the pipe, so that I can troubleshoot vcfutils vcf2fq if need be
samtools mpileup -u -C 50 -l $REGIONS -f $REFERENCE $BAM | bcftools call -c -O v > A_243_14_cons.vcf

#Then run vcfutils separately
vcfutils.pl vcf2fq -Q 5 -d 2 -D 40 A_243_14_cons.vcf > A_243_14_cons.fq

# This is what the complete pipe would look like
#samtools mpileup -u -C 50 -l $REGIONS -f $REFERENCE $BAM | bcftools call -c | vcfutils.pl vcf2fq -Q 5 -d 2 -D 40 > A_243_14_cons.fq

# Convert consensus sequence to psmc format
fq2psmcfa -q20 A_243_14_cons.fq > A_243_14.psmcfa

# Run PSMC analysis
psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o A_243_14.psmc A_243_14.psmcfa

# Quick plot to check results
psmc_plot.pl -u 1.5e-8 -g 7.5 -p A_243_14_plot A_243_14.psmc
# This was not how the final plot was generated- see the "03-Figures" folder for the R scripts to generate PSMC plots.
