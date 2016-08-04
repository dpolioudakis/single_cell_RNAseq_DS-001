#!/bin/bash

# Damon Polioudakis
# 2016-06-17
# Code to submit all samples for Drop-seq STAR alignment

# To submit this script:
#   qsub -cwd -o logs/1_Process_And_Align_$(date +%Y%m%d).log -e logs/1_Process_And_Align_$(date +%Y%m%d).error -S /bin/bash -V -N DAHM_QSUB -q geschwind.q -l h_data=8G,h_rt=12:00:00 1_Process_And_Align.sh

# To submit this script to Geschwind Short Que:
#   qsub -cwd -o logs/1_Process_And_Align_$(date +%Y%m%d).log -e logs/1_Process_And_Align_$(date +%Y%m%d).error -S /bin/bash -V -N DAHM_QSUB -q geschwind_short.q 1_Process_And_Align.sh

# Reminder: make /logs directory in code directory
################################################################################
echo ""
echo "Starting 1_Process_And_Align.sh..."$(date)
echo ""
################################################################################

# Define Input Variables and Functions

# STAR Index is mouse index I built for CNTNAP2 project
idxGenomeDir=/geschwindlabshares/RNAseq_singlecellfetal/source/IndexedGenome_GRCh37.75_NoERCC
# Fasta is from Macosko Dropseq
inFasta=/geschwindlabshares/RNAseq_singlecellfetal/source/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC.fasta
rnaStar=/share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR
dropseq_root=/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/code/Drop-seq_tools-1.12

# Path to directory that contains fastqs
inParentDir=/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/fastq
# Path for parent directory to output alignment, will create subdirectories
# below for each sample ID
outParentDir=/geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC

mkdir -p ${outParentDir}
################################################################################

# Loop through sequencing lanes in parent directory

# Alter Sxa* accordingly to match sequencing lane directories
for inDir in ${inParentDir}/Sxa*; do

  # Extract sequencing lane from file path
  echo ""
  echo "#######################################################################"
  seqLane=$(basename ${inDir})
  echo "Sequencing Lane:"
  echo ${seqLane}

  # alter 'R3*_R1_001.fastq.gz' accordingly to match FASTQ files, but make sure that *_R1 is part of regExard term
  regEx=*R1.fastq*

  # Number of files for each lane
  nFiles=$(ls ${inDir}/${regEx} | wc -l)
  echo "Total files: ${nFiles}"

  for inFastq1 in ${inDir}/${regEx}; do

    echo ""
    echo "fastq Read 1 path:"
    echo ${inFastq1}
    # fastq Read 2 paths
    inFastq2=${inFastq1//R1/R2}
    echo "fastq Read 2 path:"
    echo ${inFastq2}

    # Extract sample ID from file path
    sampleID=$(basename ${inFastq1})
    sampleID=${sampleID%%[._]R*}
    echo "Sample ID:"
    echo ${sampleID}

    # Designate an output directory for the alignment specific to that sample ID
    outDir=${outParentDir}/${seqLane}/${sampleID}
    echo "Alignment output directory:"
    echo ${outDir}
    # Tmp directory for the alignment specific to that sample ID
    # RNA STAR stats output here
    tmpDir=${outDir}/tmp

    # if [ ! -d ${outDir} ]; then

      mkdir -p ${outDir}
      mkdir -p ${tmpDir}

      echo "Converting SAM to BAM for:"
      echo ${inFastq1}
      echo ${inFastq2}
      # Convert FASTQ to BAM
      /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=${tmpDir} -jar /geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar FastqToSam \
        FASTQ=${inFastq1} \
        FASTQ2=${inFastq2} \
        OUTPUT=${outDir}/unaligned_data.bam \
        QUALITY_FORMAT=Illumina \
        SAMPLE_NAME=${sampleID} #Sample name to insert into the read group header  Required.

      echo "Submitting Drop-seq_alignment.sh script for:"
      echo ${outDir}/unaligned_data.bam

      # Run Drop-seq McCarroll lab processing and alignment tool
      # sh Drop-seq_tools-1.12/Drop-seq_alignment.sh -p -g ${idxGenomeDir} -r ${inFasta} -o ${outDir} -t ${tmpDir} -s ${rnaStar} -d ${dropseq_root} ${outDir}/unaligned_data.bam

      # File names for QSUB stdout and error output
      outSampLog=logs/1_Process_And_Align_$(date +%Y%m%d)_${seqLane}_${sampleID}.log
      outSampErr=logs/1_Process_And_Align_$(date +%Y%m%d)_${seqLane}_${sampleID}.error

      # Note name of job (-N) will only display 1st 10 characters
      qsub -cwd -o ${outSampLog} -e ${outSampErr} -V -S /bin/bash -N DA_${sampleID} -q geschwind.q -l h_data=64G,h_rt=12:00:00 -pe shared 8 Drop-seq_tools-1.12/Drop-seq_alignment.sh -p -g ${idxGenomeDir} -r ${inFasta} -o ${outDir} -t ${tmpDir} -s ${rnaStar} -d ${dropseq_root} ${outDir}/unaligned_data.bam

    # else
    #   echo "${outDir} already exists..."
    # fi

  done

done
################################################################################

echo ""
echo "End of 1_Process_And_Align.sh... "$(date)







# mkdir ../data/bam
# mkdir ../tmp
#
# # Convert FASTQ to BAM
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar /geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar FastqToSam \
#   FASTQ=../data/test/N701.R1.fastq \
#   FASTQ2=../data/test/N701.R2.fastq \
#   OUTPUT=../data/test/my_unaligned_data.bam \
#   QUALITY_FORMAT=Illumina \
#   SAMPLE_NAME=N701 #Sample name to insert into the read group header  Required.
#
# # Run Drop-seq McCarroll lab processing and alignment tool
# sh Drop-seq_tools-1.12/Drop-seq_alignment.sh -g /geschwindlabshares/RNAseq_singlecellfetal/source/IndexedGenome_GRCh37.75_NoERCC -r /geschwindlabshares/RNAseq_singlecellfetal/source/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC.fa -o ../data/bam -t ../tmp -s /share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR -p ../data/test/my_unaligned_data.bam
#
#
#
#
#
#
#
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar TagBamWithReadSequenceExtended \
#   INPUT=../data/test/my_unaligned_data.bam \
#   OUTPUT=../data/test/unaligned_tagged_Cell.bam \
#   SUMMARY=../data/test/unaligned_tagged_Cellular.bam_summary.txt \
#   BASE_RANGE=1­-12 \
#   BASE_QUALITY=10 \
#   BARCODED_READ=1 \
#   DISCARD_READ=False \
#   TAG_NAME=XC \
#   NUM_BASES_BELOW_QUALITY=1
#
#
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar TagBamWithReadSequenceExtended \
#   INPUT=../data/test/unaligned_tagged_Cell.bam \
#   OUTPUT=../data/test/unaligned_tagged_CellMolecular.bam \
#   SUMMARY=../data/test/unaligned_tagged_Molecular.bam_summary.txt \
#   BASE_RANGE=13­-20 \
#   BASE_QUALITY=10 \
#   BARCODED_READ=1 \
#   DISCARD_READ=True \
#   TAG_NAME=XM \
#   NUM_BASES_BELOW_QUALITY=1
#
#
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar FilterBAM \
#   TAG_REJECT=XQ \
#   INPUT=../data/test/unaligned_tagged_CellMolecular.bam \
#   OUTPUT=../data/test/unaligned_tagged_filtered.bam
#
# # This Drop­seq program is one of two sequence cleanup programs designed to trim away any extra sequence that might have snuck it’s way into the reads. In this case, we trim the SMART Adapter that can occur 5’ of the read. In our standard run, we look for at least 5 contiguous bases (NUM_BASES) of the SMART adapter (SEQUENCE) at the 5’ end of the read with no errors (MISMATCHES) , and hard clip those bases off the read.
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar TrimStartingSequence \
#   INPUT=../data/test/unaligned_tagged_filtered.bam \
#   OUTPUT=../data/test/unaligned_tagged_trimmed_smart.bam \
#   OUTPUT_SUMMARY=../data/test/adapter_trimming_report.txt \
#   SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
#   MISMATCHES=0 \
#   NUM_BASES=5
#
# # This Drop­seq program is the second sequence cleanup program designed to trim away trailing polyA tails from reads. It searches for at least 6 (NUM_BASES) contiguous A’s in the read with 0 mismatches (MISMATCHES), and hard clips the read to remove these bases and all bases 3’ of the polyA run.
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar PolyATrimmer \
#   INPUT=../data/test/unaligned_tagged_trimmed_smart.bam \
#   OUTPUT=../data/test/unaligned_mc_tagged_polyA_filtered.bam \
#   OUTPUT_SUMMARY=../data/test/polyA_trimming_report.txt \
#   MISMATCHES=0 \
#   NUM_BASES=6
#
#
# #Now that your data has had the cell and molecular barcodes extracted, the reads have been cleaned of SMARTSeq primer and polyA tails, and the data is now unpaired reads, it’s time to align. To do this, w e e x t r a c t t h e F A S T Q f i l e s u s i n g P i c a r d ’ s S​a m T o F a s t q ​p r o g r a m .
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar /geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar SamToFastq \
#   INPUT=../data/test/unaligned_mc_tagged_polyA_filtered.bam \
#   FASTQ=../data/test/unaligned_mc_tagged_polyA_filtered.fastq
#
# # Alignment STAR
# /share/apps/STAR_2.4.0j/bin/Linux_x86_64/STAR --genomeDir /geschwindlabshares/RNAseq_singlecellfetal/source/IndexedGenome_GRCh37.75_NoERCC --readFilesIn /geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/test/unaligned_mc_tagged_polyA_filtered.fastq --outFileNamePrefix /geschwindlabshares/RNAseq_singlecellfetal/DS-001_DP/data/test/star
#
# # T h i s p​i c a r d p r o g r a m ​i s i n v o k e d a f t e r a l i g n m e n t , t o g u a r a n t e e t h a t t h e o u t p u t f r o m a l i g n m e n t i s s o r t e d in queryname order. As a side bonus, the output file is a BAM (compressed) instead of SAM (uncompressed.)
# /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar /geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar SortSam \
#   I=../data/test/starAligned.out.sam \
#   O=../data/test/aligned.sorted.bam \
#   SORT_ORDER=queryname
#
# # This Picard program merges the sorted alignment output from STAR (ALIGNED_BAM) with the unaligned BAM that had been previously tagged with molecular/cell barcodes (UNMAPPED_BAM). This recovers the BAM tags that were “lost” during alignment. The REFERENCE_SEQUENCE argument refers to the fasta metadata file. We ignore secondary alignments, as we want only the best alignment from STAR (or another aligner), instead of assigning a single sequencing read to multiple locations on the genome.
# /usr/java/jdk1.7.0_51/bin/java -Xmx4g -Djava.io.tmpdir=../data/test/tmp -jar /geschwindlabshares/CrossDisorder_transcriptome_comparison/1_RNAseq/bin/picard-master/dist/picard.jar MergeBamAlignment \
# REFERENCE_SEQUENCE=../../source/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC.fa \
# UNMAPPED_BAM=../data/test/unaligned_mc_tagged_polyA_filtered.bam \
# ALIGNED_BAM=../data/test/aligned.sorted.bam \
# OUTPUT=../data/test/merged.bam \
# INCLUDE_SECONDARY_ALIGNMENTS=false \
# PAIRED_RUN=false
#
# # This is a Drop­seq program that adds a BAM tag “GE” onto reads when the read overlaps the exon of a gene. This tag contains the name of the gene, as reported in the annotations file. You can use either a GTF or a RefFlat annotation file with this program, depending on what annotation data source you find most useful. This is used later when we extract digital gene expression (DGE) from the BAM.
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar TagReadWithGeneExon \
#   I=../data/test/merged.bam \
#   O=../data/test/star_gene_exon_tagged.bam \
#   ANNOTATIONS_FILE=../../source/Hs.v19.refFlat \
#   TAG=GE
#
# # Detecting and repairing barcode synthesis errors
# /usr/java/jdk1.7.0_51/bin/java -Xmx2g -Djava.io.tmpdir=../data/test/tmp -jar Drop-seq_tools-1.12/jar/dropseq.jar DetectBeadSynthesisErrors \
#   I=../data/test/my.bam \
#   O=../data/test/my_clean.bam \
#   OUTPUT_STATS=../data/test/my.synthesis_stats.txt \
#   SUMMARY=../data/test/my.synthesis_stats.summary.txt \
#   NUM_BARCODES=2400 \
#   PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC
