#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=3:00:00
#PBS -o /mnt/lustre/users/djager/genotypeGVCFs/filter_SNPs_tabix.out
#PBS -e /mnt/lustre/users/djager/genotypeGVCFs/filter_SNPs_tabix.err
#PBS -m abe
#PBS -N filter_SNPs
#PBS -M dejager4@gmail.com

# Navigate to the relevant directory
cd /mnt/lustre/users/djager/genotypeGVCFs

# Add required modules/software
module add chpc/BIOMODULES
module add bcftools/1.6.33

# Set up environment variables
INPUT="/mnt/lustre/users/djager/genotypeGVCFs/buf_all_raw_snps_indels_genotyping_allscaffolds.vcf.gz"

# Index the input VCF file (required for bcftools)
tabix $INPUT

# Filter for SNPs - see manuscript for details on filtering parameters
bcftools filter -e 'FORMAT/DP < 2 | FORMAT/DP > 40 | FORMAT/GQ < 30.0' --set-GTs . $INPUT --threads 24 -O u | bcftools view -U -i 'TYPE=="snp" & MAC >=2 & F_MISSING < 0.5' -M 2 --threads 24 -O u - | bcftools view -U -i 'AF < 0.9' --threads 24 -O u - | bcftools view -U -i 'INFO/DP < 500' --threads 24 -O z - > buf_all_SNPs_allscaffolds_filtered_tabix.vcf.gz

# Generate stats on retained sites after filtering
bcftools stats buf_all_SNPs_allscaffolds_filtered_tabix.vcf.gz > buf_all_SNPs_allscaffolds_filtered_tabix.vcf.gz_stats.txt
