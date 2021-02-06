#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=24:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/catvariants.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/catvariants.err
#PBS -m abe
#PBS -N catvariants
#PBS -M dejager4@gmail.com

# This script was used to join all the Super_Scaffold VCF files output from GATK's GenotypeGVCFs (script 10A) into one VCF file

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/GenotypeGVCFs/

# Add required modules/software
module add chpc/BIOMODULES
module add java/1.8.0_73

# Join all VCF files from GenotypeGVCFs output. Use catvariants tool, 
# which works a bit differently to other GATK tools and thus has to be called 
# in a different way (for more details see: 
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_CatVariants.php)
# Execute catvariants
java -cp ../JavaTools/GenomeAnalysisTK.jar org.broadinstitute.gatk.tools.CatVariants \
 -R ../refs/buffalo.final.51scaffolds.fa \
 -V scaffolds_vcf.list \
 -out buf_all_raw_snps_indels_genotyping_allscaffolds.vcf \
 -assumeSorted 

# Compress output file VCF file with bgzip or gzip.