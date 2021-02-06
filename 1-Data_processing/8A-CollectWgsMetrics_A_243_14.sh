#!/bin/bash
#PBS -l select=1:ncpus=24:nodetype=haswell_reg
#PBS -P CBBIXXXX
#PBS -q smp
#PBS -l walltime=2:00:00
#PBS -o /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/collectwgsmetrics_dupsMarked_A_243_14.out
#PBS -e /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/collectwgsmetrics_dupsMarked_A_243_14.err
#PBS -m abe
#PBS -N A_243_14_collectwgsmetrics_dupsMarked
#PBS -M dejager4@gmail.com

# Collect genome coverage statistics using Picard

# This script was copied and run separately for each sample, changing the filename as required.

# Navigate to relevant directory
cd /mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon/


######################################
# Set up files and folders
DATAROOT=/mnt/lustre/users/djager/buf_clean/alignment_files/bam_sorted/gatk/GATK_Deon
REFROOT=$DATAROOT/refs
PICARDROOT=$DATAROOT/JavaTools
BAM="$DATAROOT/A_243_14/A_243_14_aln-PE_sorted_dups_marked.bam"
REFERENCE="$REFROOT/buffalo.final.51scaffolds.fa"
PICARD="java -Xmx50g -Djava.io.tmpdir=$DATAROOT/tmp -jar $PICARDROOT/picard.jar"

#####################################
# Add required modules/software
module add chpc/BIOMODULES
module add java/1.8.0_73

# Define and set up new function called "picardstats"
function picardstats(){

  echo "collecting WGS metrics... Started at $(date)"
  $PICARD CollectWgsMetrics \
          INPUT=$BAM \
          OUTPUT=$DATAROOT/A_243_14/A_243_14_aln-PE_sorted_dups_marked_stats_picard.txt \
          REFERENCE_SEQUENCE=$REFERENCE \
          MINIMUM_MAPPING_QUALITY=20 \
          MINIMUM_BASE_QUALITY=20
}

# Execute function "picardstats"
picardstats

# Print completion time to screen
echo "Picard completed at $(date)"


