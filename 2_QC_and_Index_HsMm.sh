#!/bin/bash

# Damon Polioudakis
# 2016-06-29
# Script to merge single end Drop-seq samples STAR output from multiple lanes (obtained in Stage 1 of the pipeline), collect QC information, and mark duplicates.

# This script should be fed into 2_QC_and_Index_HsMm_QSUB.sh
# Most QC is derived from Jill Haney's 4_QC.sh script and Neel Parikshak's script, RunPicardScripts.sh
################################################################################
echo "Starting 2_QC_and_Index_HsMm.sh... "$(date)
################################################################################

# Define Input Variables and Functions

inParentDir=$1
outDir=$2
bnSampleDir=$3
# outDir=../test/bam/S37-D4

JAV=/usr/java/jdk1.7.0_51/bin
PIC=/geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar
# Fasta is from Macosko Dropseq
inFasta=/geschwindlabshares/RNAseq_singlecellfetal/source/Macosko_Dropseq/hg19_mm10_transgenes.fasta
refFlat=/geschwindlabshares/RNAseq_singlecellfetal/source/Macosko_Dropseq/hg19_mm10_transgenes.refFlat

inDirs=${inParentDir}/Sxa*/${bnSampleDir}
echo "In dirs:"
echo ${inDirs}
echo "Out dir:"
echo ${outDir}
mkdir -p ${outDir}
################################################################################

### Merge bams

nLanes=$(echo ${inDirs} | wc -w)
echo "Number of lanes:"
echo ${nLanes}

if [ ! -f ${outDir}/star_gene_exon_tagged.bam  ] && [ ${nLanes} -gt 1 ]; then

	echo "Merging:"
	echo ${inDirs}/star_gene_exon_tagged.bam
	samtools merge -f ${outDir}/star_gene_exon_tagged.bam ${inDirs}/star_gene_exon_tagged.bam
	echo "Done merging..."

fi
################################################################################

### QC on Bamfiles using Picard Tools and RSeqQC

## Execute scripts if either QC output file is not present
# if [ ! -f ${outDir}/rnaseq_stats.txt ] || [ ! -f ${outDir}/gcbias_summary.txt ] || [ ! -f ${outDir}/PEmatched_markdup_sorted.bai  ]; then

echo "Running QC stats scripts from PicardTools..."
      mkdir -p ${outDir}/tmp
  	mkdir -p ${outDir}/tmp2

## Reorder the BAM file according to the reference at hand
# For multiple lanes after merging:
if [ ! -f ${outDir}/sorted_reads.bai  ] && [ ${nLanes} -gt 1 ]; then
	echo "Sorting merged BAM..."
	# Sort the BAM
  ${JAV}/java -Xmx4g -Djava.io.tmpdir=${outDir}/tmp -jar ${PIC} SortSam INPUT=${outDir}/star_gene_exon_tagged.bam OUTPUT=${outDir}/sorted_reads.bam SORT_ORDER=coordinate TMP_DIR=${outDir}/tmp
	echo "sorted_reads.bam is now created"
	# Index the sorted reads file
	${JAV}/java -Xmx2g -jar ${PIC} BuildBamIndex INPUT=${outDir}/sorted_reads.bam
	echo "sorted reads file is now indexed"
# For one lane:
elif [ ! -f ${outDir}/sorted_reads.bai  ]; then
	echo "Sorting BAM..."
	# Sort the BAM
	inDir=$(echo ${inDirs})
  ${JAV}/java -Xmx4g -Djava.io.tmpdir=${outDir}/tmp -jar ${PIC} SortSam INPUT=${inDir}/star_gene_exon_tagged.bam OUTPUT=${outDir}/sorted_reads.bam SORT_ORDER=coordinate TMP_DIR=${outDir}/tmp
	echo "sorted_reads.bam is now created"
	# Index the sorted reads file
	${JAV}/java -Xmx2g -jar ${PIC} BuildBamIndex INPUT=${outDir}/sorted_reads.bam
	echo "sorted reads file is now indexed"
else
    echo ".bam file already sorted"
fi

## Reorder the .bam file according to the reference at hand
if [ ! -f ${outDir}/reordered_reads.bam  ]; then
	${JAV}/java -Xmx4g -Djava.io.tmpdir=${outDir}/tmp -jar ${PIC} ReorderSam INPUT=${outDir}/sorted_reads.bam OUTPUT=${outDir}/reordered_reads.bam REFERENCE=${inFasta} TMP_DIR=${outDir}/tmp
	# Reorder the .bam file according to the reference at hand
	echo "reordered_reads.bam is now created"
else
	echo ".bam file already reordered"
fi

## Collect alignment metrics if the file is not present
if [ ! -f ${outDir}/alignment_stats.txt ]; then
	${JAV}/java -Xmx2g -jar ${PIC} CollectAlignmentSummaryMetrics REFERENCE_SEQUENCE=${inFasta} INPUT=${outDir}/reordered_reads.bam OUTPUT=${outDir}/alignment_stats.txt ASSUME_SORTED=false ADAPTER_SEQUENCE=null
	# Collect alignment metrics if the file is not present
	echo "alignment_stats.txt is now created"
else
	echo ".bam file already analyzed for alignment metrics"
fi

## Collect sequencing metrics if the file is not present
if [ ! -f ${outDir}/rnaseq_stats.txt ]; then
	${JAV}/java -Xmx2g -jar ${PIC} CollectRnaSeqMetrics REFERENCE_SEQUENCE=${inFasta} INPUT=${outDir}/reordered_reads.bam OUTPUT=${outDir}/rnaseq_stats.txt STRAND_SPECIFICITY=NONE REF_FLAT=${refFlat} ASSUME_SORTED=false
	# Collect sequencing metrics if the file is not present
	echo "rnaseq_stats.txt is now created"
else
	echo ".bam file already analyzed for RNA seq metrics"
fi

## Collect gc bias metrics if the file is not present
if [ ! -f ${outDir}/gcbias_summary.txt ]; then
	${JAV}/java -Xmx2g -jar ${PIC} CollectGcBiasMetrics REFERENCE_SEQUENCE=${inFasta} INPUT=${outDir}/reordered_reads.bam OUTPUT=${outDir}/gcbias_stats.txt ASSUME_SORTED=false CHART_OUTPUT=${outDir}/gcbias_chart.pdf SUMMARY_OUTPUT=${outDir}/gcbias_summary.txt
	# Collect gc bias metrics if the file is not present
	echo "gcbias_summary stats are now created"
else
	echo ".bam file already analyzed for GC bias"
fi

## Collect read duplication metrics if the file is not present
if [ ! -f ${outDir}/markdup_sorted.bai ]; then
	# Collect read duplication metrics if the file is not present, output the marked duplicates file AND keep duplicates for future expression analysis
	${JAV}/java -Xmx4g -Djava.io.tmpdir=${outDir}/tmp -jar ${PIC} MarkDuplicates INPUT=${outDir}/reordered_reads.bam METRICS_FILE=${outDir}/duplication_stats.txt ASSUME_SORTED=false OUTPUT=${outDir}/reordered_duplication_marked_reads.bam REMOVE_DUPLICATES=FALSE TMP_DIR=${outDir}/tmp
	echo "duplicates are now marked"

	# # Sort the marked-duplicates file
	# ${JAV}/java -Xmx4g -Djava.io.tmpdir=${outDir}/tmp2 -jar ${PIC} SortSam INPUT=${outDir}/reordered_duplication_marked_reads.bam OUTPUT=${outDir}/markdup_sorted.bam SORT_ORDER=coordinate TMP_DIR=${outDir}/tmp2
	# echo "marked-duplicates file is now sorted"
	#
	# # Index the marked-duplicates file
	# ${JAV}/java -Xmx2g -jar ${PIC} BuildBamIndex INPUT=${outDir}/markdup_sorted.bam
	# echo "marked-duplicates file is now indexed"
else
	echo ".bam file already analyzed for duplicates and processed for deduplication"
fi
# else
#   echo "RNA seq QC metric files already present"
# fi

# Save space... if QC outputs are present, delete the .bam files
if [ -f ${outDir}/reordered_duplication_marked_reads.bam ] && [ -f ${outDir}/gcbias_summary.txt ]; then
	echo "cleaning up the extra .bam files"
	rm ${outDir}/reordered_reads.bam #
	rm ${outDir}/reordered_duplication_marked_reads.bam
	rm ${outDir}/sorted_reads.bam
	rm ${outDir}/sorted_reads.bai
        # rm -rf ${PARENTDIR}/${DIRNAME}_*_L00*
	# Remove merged intermediate bams if multiple lanes were merged
	if [ ${nLanes} -gt 1 ]; then
		# rm ${outDir}/star_gene_exon_tagged.bam
		rm ${outDir}/Aligned.sortedByCoord.out.bai
	fi
fi
################################################################################

echo "End of 2_QC_and_Index_HsMm.sh... "$(date)
