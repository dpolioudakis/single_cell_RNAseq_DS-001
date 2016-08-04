#!/bin/bash

# Damon Polioudakis
# 2016-06-28
# Code to count number of lines in star_gene_exon_tagged.bam with "GE" tag added by drop-seq TagReadWithGeneExon

echo -e "SAMPLE\tNUMBER_OF_TAGGED_EXONS" > ../data/QC/human_mouse/SxaQSEQsXap108L1/Tagged_Exon_Number.txt

for inDir in ../data/bam/*/*/N70*; do
  echo -ne ${inDir}"\t" >> ../data/QC/Tagged_Exon_Number.txt
  samtools view ${inDir}/star_gene_exon_tagged.bam \
  | grep -c '.*\s*GE:.*' >> ../data/QC/Tagged_Exon_Number.txt
done
