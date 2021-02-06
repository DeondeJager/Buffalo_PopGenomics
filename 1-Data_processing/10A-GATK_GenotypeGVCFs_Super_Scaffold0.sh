#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=6:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/genotypeGVCFs_SS0.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/genotypeGVCFs_SS0.err
#PBS -m abe
#PBS -N genotypeGVCFs_SS0
#PBS -M dejager4@gmail.com

# This script was copied and run separately for each Super_Scaffold, changing the filename as required.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/

# Add required modules/software
module add chpc/BIOMODULES
module add java/1.8.0_73

# Generate the list of the gvcf (sample) files that you would like to do genotyping (those created by HaplotypeCaller; scripts 9A and 9B):
#ls *.g.vcf > gvcfs.list ## I just did this outside of the script, hence commenting it out.

# Execute GATK GenotypeGVCFs to call genotypes across all samples in the gvcfs.list, for Super_Scaffold0:
java -jar JavaTools/GenomeAnalysisTK.jar \
 -T GenotypeGVCFs \
 -R refs/buffalo.final.51scaffolds.fa \
 -L Super_Scaffold0 \
 --variant gvcfs.list \
 --includeNonVariantSites \
 -o buf_all_raw_snps_indels_genotyping_SS0.vcf
