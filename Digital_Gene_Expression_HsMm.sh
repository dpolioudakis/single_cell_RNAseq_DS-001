#!/bin/bash

# Damon Polioudakis
# 2016-06-20
# Code to submit all aligned Drop-seq samples to Digital Gene Expression

# To submit this script:
# qsub -cwd -o logs/Digital_Gene_Expression_HsMm_$(date +%Y%m%d).log -e logs/Digital_Gene_Expression_HsMm_$(date +%Y%m%d).error -S /bin/bash -V -N DgtGE -q geschwind.q -l h_data=8G,h_rt=12:00:00 Digital_Gene_Expression_HsMm.sh

# To submit this script to Geschwind Short Que:
# qsub -cwd -o logs/Digital_Gene_Expression_HsMm_$(date +%Y%m%d).log -e logs/Digital_Gene_Expression_HsMm_$(date +%Y%m%d).error -S /bin/bash -V -N DgtGE -q geschwind_short.q Digital_Gene_Expression_HsMm.sh

# Reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting Digital_Gene_Expression_HsMm.sh... "$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# Path to directory that contains fastqs
inParentDir=/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/bam/human_mouse
# Path for parent directory to output alignment, will create subdirectories
# below for each sample ID
outParentDir=/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/digital_gene_expression/human_mouse
# List of cell barcodes derived from Reads_Per_Cell_Histogram_HsMm.R
# Made with:
# for dir in *; do gzip -cd ${dir}/out_cell_readcounts.txt.gz | head -400 | cut -f2 > ${dir}/out_cell_subset_barcodes_list.txt; done
cellBarcodeFile=out_cell_subset_barcodes_list.txt
# Reads mapped to human output from Filter_BAM_Human_Mouse.sh
inHsBamName=star_gene_exon_tagged_human.bam
# Reads mapped to mouse output from Filter_BAM_Human_Mouse.sh
inMmBamName=star_gene_exon_tagged_mouse.bam

mkdir -p ${outParentDir}
################################################################################

# Loop through sequencing lanes and samples in parent directory

# Alter Sxa* accordingly to match sequencing lane directories
for inDir in ${inParentDir}/Sxa*/N*; do

  # Extract lane from directory path
  seqLane=$(basename $(dirname ${inDir}))

  # Extract sample ID from directory path
  sampleID=$(basename ${inDir})
  sampleID=${sampleID%%[._]R*}
  echo ""
  echo "Sample ID:"
  echo ${sampleID}

  # Designate an output directory for the alignment specific to that sample ID
  outDir=${outParentDir}/${seqLane}/${sampleID}
  echo "Alignment output directory:"
  echo ${outDir}

  # if [ ! -d ${outDir} ]; then

  mkdir -p ${outDir}

  echo "Running Drop-seq DigitalExpression for:"
  echo ${inDir}/star_gene_exon_tagged.bam

  ## Run Drop-seq McCarroll lab processing and alignment tool

  # Human
  # Collapse UMIs
  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar DigitalExpression INPUT=${inDir}/${inHsBamName} OUTPUT=${outDir}/out_gene_exon_tagged_human.dge.txt.gz SUMMARY=${outDir}/out_gene_exon_tagged_human.dge.summary.txt CELL_BC_FILE=${inDir}/${cellBarcodeFile}
  # Output number of reads (not collapsed by UMI)
  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar DigitalExpression INPUT=${inDir}/${inHsBamName} OUTPUT=${outDir}/out_gene_exon_tagged_human.counts.txt.gz SUMMARY=${outDir}/out_gene_exon_tagged_human.counts.summary.txt CELL_BC_FILE=${inDir}/${cellBarcodeFile} OUTPUT_READS_INSTEAD=true

  # Mouse
  # Collapse UMIs
  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar DigitalExpression INPUT=${inDir}/${inMmBamName} OUTPUT=${outDir}/out_gene_exon_tagged_mouse.dge.txt.gz SUMMARY=${outDir}/out_gene_exon_tagged_mouse.dge.summary.txt CELL_BC_FILE=${inDir}/${cellBarcodeFile}
  # Output number of reads (not collapsed by UMI)
  /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar DigitalExpression INPUT=${inDir}/${inMmBamName} OUTPUT=${outDir}/out_gene_exon_tagged_mouse.counts.txt.gz SUMMARY=${outDir}/out_gene_exon_tagged_mouse.counts.summary.txt CELL_BC_FILE=${inDir}/${cellBarcodeFile} OUTPUT_READS_INSTEAD=true

done
################################################################################

echo ""
echo "End of Digital_Gene_Expression_HsMm.sh... "$(date)
