#!/usr/bin/perl
#$ -l time=10:00:00

##Convert Read 1 (first pair of the paired end) into a FASTQ file

$filename = $ARGV[0];

#use strict;
#use warnings;

open(INPUT_FILE, $filename);
        #       or die "Couldn't open $filename for reading you GSPs!";
        
while (<INPUT_FILE>) {
  chomp;
  my ($seqIndex, $seq1, $seq2, $qualIndex, $qual1, $qual2, $PF, $header)=split /\t/, $_;
  my ($instrumentname, $runID, $lanenumber, $tilenumber, $xcoord, $ycoord)=split /_/, $header;

  # According to the instruction manual from BSCRC sequencing core
  # the PHRED scores are +64 and this is default for subsequent programs so not changing this
  # Convert Illumina's quality values to true Phred scale
  # $qual1 =~ tr/\x40-\xff\x00-\x3f/\x21-\xe0\x21/;

  print "\@$instrumentname:$runID:$runID:$lanenumber:$tilenumber:$xcoord:$ycoord 1:N:0:$seqIndex\n$seq1\n+\n$qual1\n"
}
