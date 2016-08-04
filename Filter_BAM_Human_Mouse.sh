#!/bin/bash

# Damon Polioudakis
# 2016-06-22
# Code to submit all samples to Drop-seq FilterBAM to subset to human or mouse

# To submit this script to Geschwind Short Que:
# qsub -cwd -o logs/Filter_BAM_Human_Mouse_$(date +%Y%m%d).log -e logs/Filter_BAM_Human_Mouse_$(date +%Y%m%d).error -S /bin/bash -V -N DsFltBAM -q geschwind_short.q Filter_BAM_Human_Mouse.sh

# Reminder: make /logs directory in code directory

################################################################################
echo ""
echo "Starting Filter_BAM_Human_Mouse.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# Path to directory that contains star_gene_exon_tagged.bams
workParentDir=../data/bam/human_mouse

################################################################################

# Loop through sequencing lanes in parent directory

# Alter Sxa* accordingly to match sequencing lane directories
for workDir in ${workParentDir}/Sxa*/N*; do

  echo "Submitting Drop-seq BAMTagHistogram script for:"
  echo ${workDir}

  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar FilterBAM I=${workDir}/star_gene_exon_tagged.bam O=${workDir}/star_gene_exon_tagged_human.bam REF_SOFT_MATCHED_RETAINED=HUMAN

  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar FilterBAM I=${workDir}/star_gene_exon_tagged.bam O=${workDir}/star_gene_exon_tagged_mouse.bam REF_SOFT_MATCHED_RETAINED=MOUSE

done
################################################################################

echo ""
echo "End of Filter_BAM_Human_Mouse.sh... "$(date)
