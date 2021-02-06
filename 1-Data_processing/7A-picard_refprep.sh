#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=0:30:00
#PBS -o /mnt/lustre/users/djager/buffalo_refgenome/picardrefprep.out
#PBS -e /mnt/lustre/users/djager/buffalo_refgenome/picardrefprep.err
#PBS -m abe
#PBS -N refprep_picard
#PBS -M dejager4@gmail.com

# Prepare the reference genome for Picard

# Note: This script takes about 30 seconds to run, so you can really just do it in an interactive session.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buffalo_refgenome

# Add required modules/software
module add chpc/BIOMODULES
module add samtools/1.3.1

# Set up variables
GATKROOT=/mnt/lustre/users/djager/buffalo_refgenome
PICARDROOT=/home/djager/programs/gatk/GATK_Deon/JavaTools
PICARD="java -Xmx50g -Djava.io.tmpdir=$GATKROOT/tmp -jar $PICARDROOT/picard.jar"
SAMTOOLS="samtools"
REFERENCE="/mnt/lustre/users/djager/buffalo_refgenome/buffalo.final.51scaffolds.fa"

# Run commands to create sequence dictionary of reference genome:
$PICARD CreateSequenceDictionary \
            REFERENCE=$REFERENCE \
            OUTPUT="${REFERENCE%.fa}.dict"

$SAMTOOLS faidx $REFERENCE

