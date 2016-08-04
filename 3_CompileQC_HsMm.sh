#!/bin/bash

# Damon Polioudakis
# 2016-06-29
# Script to compile QC information output from QC_and_Index.sh

# Suggested script call:
# sh 3_CompileQC.sh 2>&1 | tee logs/3_CompileQC_$(date +%Y%m%d).log

# To submit via QSUB:
# qsub -cwd -V -N QC_comp -S /bin/bash -q geschwind.q -l h_data=16G,h_rt=3:00:00 3_CompileQC.sh

# NOTE:
#   The headers from Picards statistics do not match the output columns
#   Code for compiling alignment stats needs to be changed depending on if
#   RNAseq is paired end or single end

# Adopted from Jill Haney's 3_CompileQC.sh script, most QC derived from
# Neel Parikshak's script, CompileQCOutput.sh
################################################################################

# Set variables: Location for QC files and desired output directory
inParentDir=../data/bam/human_mouse/SxaQSEQsXap108L1
outQCdir=../data/QC/human_mouse/SxaQSEQsXap108L1
################################################################################

mkdir -p ${outQCdir}
## Start the RNA-seq summary file

for inDir in ${inParentDir}/*/Picard_Stats; do

  echo ""
  echo "Input directory:"
  echo ${inDir}
  sampleID=${inDir##${inParentDir}/}
  echo "Sample ID:"
  echo ${sampleID}

  if [ -f ${inDir}/gcbias_stats.txt ]; then

    echo "Getting rnaseq stats"
    if [ -f ${outQCdir}/rnaseq_stats.txt ]; then
    	var=`sed -n 8p ${inDir}/rnaseq_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/rnaseq_stats.txt
    else
      var=`sed -n 7p ${inDir}/rnaseq_stats.txt`
    	echo "SAMPLE" ${var} > ${outQCdir}/rnaseq_stats.txt
    	var=`sed -n 8p ${inDir}/rnaseq_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/rnaseq_stats.txt
    fi
    var=`sed -n 11,112p ${inDir}/rnaseq_stats.txt`
    echo ${sampleID} ${var} >> ${outQCdir}/rnaseq_stats_Transcript_Coverage.txt

  	echo "Getting gcbias summary"
    if [ -f ${outQCdir}/gcbias_summary.txt ]; then
    	var=`sed -n 8p ${inDir}/gcbias_summary.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/gcbias_summary.txt
    else
      var=`sed -n 7p ${inDir}/gcbias_summary.txt`
    	echo "SAMPLE" ${var} > ${outQCdir}/gcbias_summary.txt
      var=`sed -n 8p ${inDir}/gcbias_summary.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/gcbias_summary.txt
    fi

    echo "Getting gcbias stats"
  	var=`sed -n 8,1088p ${inDir}/gcbias_stats.txt`
  	echo ${sampleID} ${var} >> ${outQCdir}/gcbias_stats.txt

  	echo "Getting alignment summary"
    if [ -f ${outQCdir}/alignment_stats.txt ]; then
      # Change depending on if RNAseq is PE or single end
    	# var=`sed -n 10p ${inDir}/alignment_stats.txt`
      var=`sed -n 8p ${inDir}/alignment_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/alignment_stats.txt
    else
      var=`sed -n 7p ${inDir}/alignment_stats.txt`
    	echo "SAMPLE" ${var} > ${outQCdir}/alignment_stats.txt
      # Change depending on if RNAseq is PE or single end
      # var=`sed -n 10p ${inDir}/alignment_stats.txt`
      var=`sed -n 8p ${inDir}/alignment_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/alignment_stats.txt
    fi

  	echo "Getting duplication summary"
    if [ -f ${outQCdir}/duplication_stats.txt ]; then
    	var=`sed -n 9p ${inDir}/duplication_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/duplication_stats.txt
    else
      var=`sed -n 7p ${inDir}/duplication_stats.txt`
    	echo "SAMPLE" ${var} > ${outQCdir}/duplication_stats.txt
      var=`sed -n 9p ${inDir}/duplication_stats.txt`
    	echo ${sampleID} ${var} >> ${outQCdir}/duplication_stats.txt
    fi

  else
  	echo "No file for" ${inDir}
  fi

done
