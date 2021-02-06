#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=2:00:00
#PBS -o /mnt/lustre/users/djager/buffalo_refgenome/scaffoldstitcher_scaffolds/scaffoldstitcher_scaffolds.out
#PBS -e /mnt/lustre/users/djager/buffalo_refgenome/scaffoldstitcher_scaffolds/scaffoldstitcher_scaffolds.err
#PBS -m abe
#PBS -N scaffoldstitcher_scaffolds
#PBS -M dejager4@gmail.com

# ScaffoldStitcher was used to concatenate the 442,402 scaffolds and contigs of the reference genome into artificial chromosomes, or "Super_Scaffolds",
# for use with GATK, which doesn't work well with >100-odd scaffolds in the reference genome.

# ScaffoldStitcher is available at: https://github.com/ameliahaj/ScaffoldStitcher

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buffalo_refgenome/

# Add required modules/software
module add chpc/BIOMODULES
module add python/3.5.2_gcc-6.1.0

# Execute ScaffoldStitcher
python /home/djager/programs/ScaffoldStitcher/ScaffoldStitcher.py \
	-fasta buffalo.final.fa \
	-identifier scaffold \
	-short 36 \
	-nlength 1000 \
	-maxlength 60000000 \
	> scaffoldstitcher_scaffolds/buffalo_scaffoldsjoined.fa

####The nlength was changed in this script and scaffoldstitcher_contigs until approximately 52 super scaffolds were obtained
####so that there would be ~52 super scaffolds, the number of chromosomes that buffalo have. 
