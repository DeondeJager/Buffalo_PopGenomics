#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=24:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/angsd_het.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14/angsd_het.err
#PBS -m abe
#PBS -N A_243_14_het
#PBS -M dejager4@gmail.com

# This script was copied and run for each sample, changing the sample and filenames as required.
# For the comparison of GL1 and GL2, the script was identical except for changing that parameter from GL1 to GL2.

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14

# Path to tools and files (set up environment variables)
REFERENCE='/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa'
INPUT=/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/A_243_14
BAM="$INPUT/A_243_14_aln-PE_sorted_dups_marked.bam"
ANGSD=/home/djager/programs/angsd

# Add required modules/software
# NB: Have to load bzip2 library for angsd to work!
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add R/3.4.1-gcc6.3.0

# Global estimate (SFS estimation for single samples) - see http://www.popgen.dk/angsd/index.php/Input#BAM.2FCRAM for explanation of parameters/options.
# Also this URL for the actual instructions: http://www.popgen.dk/angsd/index.php/Heterozygosity.
# -GL 2 means calculate genotype likelihood using the GATK algorithm (from 1st version of GATK) - see http://www.popgen.dk/angsd/index.php/Genotype_Likelihoods#GATK_genotype_likelihoods for more information
# But -GL 2 does no error modeling or correction, while -GL 1 (the SAMtools model) does. See Appendix G
# of Renaud et al. 2019 (https://www.genetics.org/content/genetics/212/3/587.full.pdf).

# Run the heterozygosity analysis:
# Create saf.idx file in ANGSD:
$ANGSD/angsd -i $BAM \
 -C 50 \
 -ref $REFERENCE \
 -anc $REFERENCE \
 -dosaf 1 \
 -fold 1 \
 -GL 1 \
 -nThreads 24 \
 -minmapq 30 \
 -minQ 20 \
 -uniqueOnly 1 \
 -remove_bads 1 \
 -only_proper_pairs 1

# This is followed by the actual estimation using the realSFS subprogram of ANGSD: 
$ANGSD/misc/realSFS angsdput.saf.idx -P 24 > $INPUT/est.ml

# Then, in an interactive session use R to calculate the proportion of heterozygous sites.
# By dividing the number of sites in the second entry of the sfs (the heterozygous sites) by
# the sum of the first (non-variant sites) and second entries (i.e. total number of sites)
R
a<-scan("est.ml")
He <- a[2]/sum(a)
write.table(He, file = "He_angsd.txt", sep = '\t', col.names = F, row.names = F)
quit()
