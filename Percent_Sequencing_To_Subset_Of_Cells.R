# Damon Polioudakis
# 2016-07-18
# Number of exonic reads for subset of cell barcodes versus all exonic reads
################################################################################

rm(list=ls())
sessionInfo()

require(reshape2)
require(ggplot2)

## Load data and assign variables

# dirs <- c("N702", "N703")
dirs <- list.files("../data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/")

inParent <- "../data/digital_gene_expression"

dir.create("../analysis/graphs", recursive = TRUE)
dir.create("../analysis/tables/", recursive = TRUE)

tReadsDF <- read.table("../data/QC/Reads_Per_Sample.txt", fill = TRUE
                       , header = TRUE)

# Set ggplot2 theme
theme_set(theme_bw())
theme_set(theme_get() + theme(text=element_text(size=18)))
theme_update(plot.title = element_text(size = 16))
################################################################################

# Convert FASTQ line count to read count (divide by 4)
tReadsDF$TOTAL_FRAGMENTS <- tReadsDF$LINES_IN_FASTQ / 4

## Calculate percentages for each sample

# Loop through sample directories
pctExonicL <- list()
pctTotalL <- list()
mnGeneDetL <- list()
nCellsL <- list()
for (dir in dirs) {
  
  # Load human counts pre UMI collapse from human mouse mapping to filter
  hsDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
                             , "/out_gene_exon_tagged_human.counts.summary.txt"), header = TRUE)
  
  print("")
  print(dir)
  str(hsDat)
  # List of barcodes with > X human exonic reads
  print(dim(hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]))
  hsBarcodes <- hsDat[hsDat$NUM_TRANSCRIPTS > 10000, ]$CELL_BARCODE
  
  # Load mouse counts pre UMI collapse from human mouse mapping to filter
  mmDat <- read.table(paste0(inParent, "/human_mouse/SxaQSEQsXap108L1/", dir
                             , "/out_gene_exon_tagged_mouse.counts.summary.txt"), header = TRUE)
  str(mmDat)
  # List of barcodes with < X mouse exonic reads
  print(dim(mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]))
  mmBarcodes <- mmDat[mmDat$NUM_TRANSCRIPTS < 1000, ]$CELL_BARCODE
  
  hsCountsStats <- read.table(paste0(inParent
      , "/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/"
      , dir, "/out_gene_exon_tagged_human.counts.summary.txt")
      , header = TRUE)
  
  # Number of exonic reads after filtering to list of real cell barcodes
  bcExonic <- sum(hsCountsStats$NUM_TRANSCRIPTS)
  
  # Number of exonic reads after filtering to list of cell barcodes with > X
  # human and < X mouse exonic reads
  idxHs <- hsCountsStats$CELL_BARCODE %in% mmBarcodes & hsCountsStats$CELL_BARCODE %in% hsBarcodes
  print(table(idxHs))
  bcFtExonic <- sum(hsCountsStats$NUM_TRANSCRIPTS[idxHs])
  
  # Number of exonic reads in sample
  exonicDF <- read.table("../data/QC/Tagged_Exon_Number.txt", header = TRUE)
  tExonic <- exonicDF[exonicDF$SAMPLE == "../data/bam/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/N701", ]$NUMBER_OF_TAGGED_EXONS
  tReads <- tReadsDF$TOTAL_FRAGMENTS[tReadsDF$SAMPLE == paste0("../data/fastq/SxaQSEQsXap108L1/", dir, ".R1.fastq.gz")]
  
  # Percent exonic reads remaining after filtering to list of cell barcodes
  pctExonicL[[dir]] <- c(bcExonic / tExonic * 100, bcFtExonic / tExonic * 100)
  pctTotalL[[dir]] <- c(bcExonic / tReads * 100, bcFtExonic / tReads * 100)
  
  ## Number of genes detected
  dgeDF <- read.table(paste0(inParent
      , "/Homo_sapiens.GRCh37.75.dna.primary_assembly_NoERCC/SxaQSEQsXap108L1/"
      , dir, "/out_gene_exon_tagged_human.dge.txt.gz")
                      , header = TRUE, stringsAsFactors = FALSE, row.names = 1)
  
  # Genes detected: Filtered list of cell barcodes
  dgeHsFtDF <- dgeDF[ ,idxHs]
  dgeHsFtDF <- dgeHsFtDF >= 1
  # Convert to numeric
  dgeHsFtDF <- dgeHsFtDF + 0
  nGenesFt <- colSums(dgeHsFtDF)
  mnFt <- round(mean(nGenesFt, na.rm = TRUE), 2)
  
  # Genes detected: Real cell barcodes
  dgeHsDF <- dgeDF >= 1
  # Convert to numeric
  dgeHsDF <- dgeHsDF + 0
  nGenes <- colSums(dgeHsDF)
  mn <- round(mean(nGenes, na.rm = TRUE), 2)
  
  mnGeneDetL[[dir]] <- c(mn, mnFt)
  
  ## Number of cells
  nCellsL[[dir]] <- c(dim(hsCountsStats)[1], sum(idxHs))
}
# Format percent
pctDF <- data.frame(Sample = melt(data.frame(pctTotalL))[ ,1]
      , BarcodeFilter = rep(c("Read Depth", "Read Depth + Exonic Read Depth"), 3)
      , Total = melt(data.frame(pctTotalL))[ ,2]
      , Exonic = melt(data.frame(pctExonicL))[ ,2]
      , nGenes = melt(data.frame(mnGeneDetL))[ ,2]
      , nCells = melt(data.frame(nCellsL))[ ,2])

## Write table
write.table(pctDF, "../analysis/tables/Percent_Sequencing_To_Subset_Of_Cells.txt"
            , quote = FALSE, sep = "\t")