#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=96:00:00
#PBS -o /mnt/lustre/users/djager/ANGSD/PCA/angsd_PCA.out
#PBS -e /mnt/lustre/users/djager/ANGSD/PCA/angsd_PCA.err
#PBS -m abe
#PBS -N angsd_PCA
#PBS -M dejager4@gmail.com

# Perform principal component analysis using ngsCovar, of the ngsTools package
# Following tutorial from this website: https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md
# Except I used the VCF file with SNPs already filtered from the GATK-BCFtools pipeline.

# Navigate to the folder with the data files
cd /mnt/lustre/users/djager/ANGSD/

# Set up environment variables
ANGSDPATH=/home/djager/programs/angsd/angsd
OUTDIR=/mnt/lustre/users/djager/ANGSD/PCA
VCFPATH=/mnt/lustre/users/djager
REFERENCE='/mnt/lustre/users/djager/buffalo_refgenome/samtools_index/buffalo.final.51scaffolds.fa.fai'
NGSCOVAR=/home/djager/programs/ngsTools/ngsTools/ngsPopGen
NGSTOOLS=/home/djager/programs/ngsTools/ngsTools

# Add required modules/software
# NB: Have to load bzip2 library for angsd to work!
module add chpc/BIOMODULES
module add bzip2/1.0.6
module add R/3.4.1-gcc6.3.0
module add gcc/6.3.0

# Create genotype probability file in binary format using ANGSD (-doGeno 32) from VCF file:
$ANGSDPATH/angsd -vcf-gl $VCFPATH/buf_all_SNPs_allscaffolds_filtered_tabix.vcf \
  -out $OUTDIR/PCA \
  -nind 40 \
  -fai $REFERENCE \
  -nThreads 24 \
  -doGeno 32 \
  -doMajorMinor 1 \
  -doMaf 1 \
  -doPost 1

# Double-check the number of sites to make sure the conversion was successful (if you want)
# In an interactive session:
#NSITES=zcat PCA.mafs.gz | wc -l
#      = 3804536 - 1 (it contains a header line
#      = 3804535

# Unzip the the geno output file from angsd:
gunzip $OUTDIR/PCA.geno.gz

# Then estimate the covariance matrix between individuals based on genotype probabilities using ngsCovar:
$NGSCOVAR/ngsCovar -probfile $OUTDIR/PCA.geno \
  -outfile $OUTDIR/PCA.covar \
  -nind 40 \
  -nsites 3804535 \
  -call 0 \
  -norm 0

# Finally, we perform an eigenvector decomposition and plot the resulting map:
# First we need to create a plink cluster-like file defining the labelling (population) for each sample
Rscript -e 'write.table(cbind(seq(1,40),rep(1,40),c(rep("AENP",5),rep("KNP",15),rep("HiP",15),rep("MNP",5))), row.names=F, sep="\t", col.names=c("FID","IID","CLUSTER"), file="/mnt/lustre/users/djager/ANGSD/PCA/PCA.clst", quote=F)'

#Run and plot:
Rscript $NGSTOOLS/Scripts/plotPCA.R -i $OUTDIR/PCA.covar -c 1-2 -a $OUTDIR/PCA.clst -o $OUTDIR/PCA4.pdf
