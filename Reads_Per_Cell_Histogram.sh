#!/bin/bash

# Damon Polioudakis
# 2016-06-22
# Code to submit all samples to Drop-seq BAMTagHistogram

# To submit this script to Geschwind Short Que:
# qsub -cwd -o logs/Reads_Per_Cell_Histogram_$(date +%Y%m%d).log -e logs/Reads_Per_Cell_Histogram_$(date +%Y%m%d).error -S /bin/bash -V -N DRH_QSUB -q geschwind_short.q Reads_Per_Cell_Histogram.sh

# Reminder: make /logs directory in code directory

################################################################################
echo ""
echo "Starting Reads_Per_Cell_Histogram.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# Path to directory that contains star_gene_exon_tagged.bams
workParentDir=../data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC

################################################################################

# Loop through sequencing lanes in parent directory

# Alter Sxa* accordingly to match sequencing lane directories
for workDir in ${workParentDir}/Sxa*/N*; do

  echo "Submitting Drop-seq BAMTagHistogram script for:"
  echo ${workDir}

  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar BAMTagHistogram I=${workDir}/star_gene_exon_tagged.bam O=${workDir}/out_cell_readcounts.txt.gz TAG=XC

done
################################################################################

echo ""
echo "End of Reads_Per_Cell_Histogram.sh... "$(date)
